""" Import Apr 2022 murine in vitro data """
function importMurineInVitro(fn = "CHO-mFcgR-apr2022.csv")
    df = CSV.File(joinpath(dataDir, fn), comment = "#") |> DataFrame
    df = df[df."Experiment" .!= 5, :]   # throw away Exp. 5
    df[!, ["IgG1", "IgG2b", "IgG2c"]] .-= df[!, "TNP-BSA"]  # subtract TNP-BSA control from meas
    df = df[df."Receptor" .!= "CHO", Not("TNP-BSA")]
    df = dropmissing(stack(df, ["IgG1", "IgG2b", "IgG2c"], variable_name = "Subclass", value_name = "Value"))
    baseline = combine(groupby(df, "Experiment"), "Value" => geomean => "Baseline")
    df = innerjoin(df, baseline, on = "Experiment")
    df[!, "Value"] ./= df[!, "Baseline"]    # normalize fluorescence by daily geomean
    df = df[!, Not(["Experiment", "Baseline"])]
    return sort!(df, ["Receptor", "Subclass"])
end

""" Priors for murine affinities. Not including IgG3"""
@memoize function murineKavDist_nested(; regularKav = false, retdf = true)
    Kav = importKav(; murine = true, retdf = true)
    Kav = Kav[Kav."IgG" .!= "IgG3", :]
    function retDist(x; regularKav = regularKav)
        x = maximum([1e4, x])
        if regularKav
            return x
        end
        return inferLogNormal(x, x * 3)
    end
    Kav[!, Not("IgG")] = retDist.(Kav[!, Not("IgG")], regularKav = regularKav)
    if !retdf
        return Matrix(Kav[!, Not("IgG")])
    end
    return Kav
end

murineKavDist(; kwargs...) = deepcopy(murineKavDist_nested(; kwargs...))

const InVitroMurineRcpExp = Dict("FcgRI" => 1e6, "FcgRIIB" => 1e6, "FcgRIII" => 1e6, "FcgRIV" => 1e6)

function predictMurine(dfr::DataFrameRow; KxStar = KxConst, recepExp = InVitroMurineRcpExp)
    if (dfr."Subclass" == "TNP-BSA") || (dfr."Affinity" <= 0.0)
        return 0.0
    end
    if dfr."Affinity" >= Inf
        return Inf
    end
    return polyfc(1e-9, KxStar, dfr."Valency", [recepExp[dfr."Receptor"] * dfr."Expression" / 100], [1.0], reshape([dfr."Affinity"], 1, 1)).Lbound
end

function predictMurine(df::AbstractDataFrame; Kav = murineKavDist(; regularKav = true), kwargs...)
    # Add affinity information to df
    Kav = deepcopy(Kav)
    Kav[Kav."IgG" .== "IgG2a", "IgG"] .= "IgG2c"
    dft = innerjoin(df, Kav, on = "Subclass" => "IgG")
    dft = stack(dft, murineFcgR)
    rename!(dft, "value" => "Affinity")
    dft = dft[dft."Receptor" .== dft."variable", Not("variable")]
    sort!(dft, ["Receptor", "Subclass"])
    @assert df."Value" == dft."Value"   # must ensure the right order
    preds = Vector(undef, size(dft)[1])
    Threads.@threads for i = 1:size(dft)[1]
        preds[i] = predictMurine(dft[i, :]; kwargs...)
    end
    dft."Predict" = preds
    #@assert all(isfinite(dft[!, "Predict"]))
    dft[dft."Predict" .<= 0.0, "Predict"] .= 1e-12
    dft[!, "Predict"] ./= geomean(dft."Predict") / geomean(dft."Value")
    
    return dft
end

@model function murineFit(df, values; KxStar = KxConst)
    # Rtot sample
    Rtot = Vector(undef, length(murineFcgR))
    for ii in eachindex(Rtot)
        ref = InVitroMurineRcpExp[murineFcgR[ii]]
        Rtot[ii] ~ inferLogNormal(ref, ref * 1e2)
        # Treat receptor amount as unknown with a wide prior
    end
    Rtotd = Dict([murineFcgR[ii] => Rtot[ii] for ii = 1:length(Rtot)])

    # Kav sample
    Kav_dist = Matrix(murineKavDist()[:, Not("IgG")])
    Kavd = murineKavDist(; regularKav = true)
    Kavd = Kavd[Kavd."IgG" .!= "IgG3", :]
    Kav = Matrix(undef, size(Kav_dist)...)
    for ii in eachindex(Kav)
        Kav[ii] ~ Kav_dist[ii]
    end
    Kavd[!, Not("IgG")] = Kav

    # fit predictions
    if all(Rtot .>= 0.0) && all(0.0 .<= Kav .< Inf)
        df = predictMurine(deepcopy(df); recepExp = Rtotd, Kav = Kavd, KxStar = KxStar)
    else
        df = deepcopy(df)
        df."Predict" .= Inf
    end

    stdv = std(log.(df."Predict") - log.(values))
    values ~ MvLogNormal(log.(df."Predict"), stdv * I)
    nothing
end

function runMurineMCMC(fname = "murine_ADVI_0423.dat"; KxStar = KxConst, mcmc_iter = 1_000)
    if isfile(fname)
        return deserialize(fname)
    end
    df = importMurineInVitro()
    m = murineFit(df, df."Value"; KxStar = KxStar)
    # use MAP estimate as starting point
    opts = Optim.Options(iterations = 500, show_every = 10, show_trace = true)
    opt = optimize(m, MAP(), LBFGS(; m = 20), opts)
    c = sample(m, NUTS(), mcmc_iter, init_params = opt.values.array)
    #q = vi(m, ADVI(10, 1000))
    f = serialize(fname, c)
    return c
end

function plot_murineMCMC_dists(c::Union{Chains, MultivariateDistribution} = runMurineMCMC(); bincount = 20)
    setGadflyTheme()

    # Plot Kav's
    Kav = murineKavDist()
    ligg, lfcr = size(Kav)
    lfcr -= 1
    Kav_pls = Matrix{Plot}(undef, ligg, lfcr)
    if c isa MultivariateDistribution
        cc = rand(c, 100_000)
        @assert size(c)[1] == 18
    end
    for ii in eachindex(Kav_pls)
        IgGname = Kav."IgG"[(ii - 1) % ligg + 1]
        FcRname = names(Kav)[2:end][(ii - 1) ÷ ligg + 1]
        name = IgGname * " to " * FcRname
        dat = (c isa Chains) ? c["Kav[$ii]"].data : cc[4+ii, :]  # 5-16 are Kav's
        Kav_pls[ii] = plotHistPriorDist(dat, Kav[(ii-1)%ligg+1, FcRname], name; bincount = bincount)
    end
    Kav_plot = plotGrid((ligg, lfcr), permutedims(Kav_pls, (2, 1)); sublabels = false)
    draw(SVG("MCMCmurine_Kav.svg", 16inch, 13inch), Kav_plot)
    draw(PDF("MCMCmurine_Kav.pdf", 16inch, 13inch), Kav_plot)

    # Plot Rtot's
    Rtot_pls = Vector{Plot}(undef, lfcr)
    Rtot_dist = [inferLogNormal(InVitroMurineRcpExp[fcr], InVitroMurineRcpExp[fcr] * 1e2) for fcr in names(Kav)[2:end]]
    for ii in eachindex(Rtot_pls)
        FcRname = names(Kav)[2:end][ii]
        dat = (c isa Chains) ? c["Rtot[$ii]"].data : cc[ii, :]   # 1-4 are Rtot's
        Rtot_pls[ii] = plotHistPriorDist(dat, Rtot_dist[ii], FcRname; bincount = bincount)
    end
    Rtot_plot = plotGrid((1, lfcr), Rtot_pls; sublabels = false)
    draw(SVG("MCMCmurine_Rtot.svg", 16inch, 4inch), Rtot_plot)
    draw(PDF("MCMCmurine_Rtot.pdf", 16inch, 4inch), Rtot_plot)
end

function plot_murineMCMC_predict(c::Union{Chains, MultivariateDistribution} = runMurineMCMC(), 
        df = importMurineInVitro(); title = nothing, KxStar = KxConst, kwargs...)
    
    if c isa Chains
        # Extract MCMC results
        Rtot = [median(c["Rtot[$i]"].data) for i = 1:length(murineFcgR)]
        Rtotd = Dict([murineFcgR[ii] => Rtot[ii] for ii = 1:length(murineFcgR)])

        Kavd = murineKavDist(; regularKav = true)
        Kav = [median(c["Kav[$i]"].data) for i = 1:length(murineKavDist(; regularKav = true, retdf = false))]
        Kavd[!, Not("IgG")] = typeof(Kav[1, 1]).(reshape(Kav, size(Kavd)[1], :))
        title = (title === nothing) ? "Murine fitting with NUTS (median)" : title
    else
        # Extract ADVI results
        x = median(rand(c, 100000), dims = 2)
        Rtot = x[1:4]
        Rtotd = Dict([FcRegression.murineFcgR[ii] => Rtot[ii] for ii = 1:length(FcRegression.murineFcgR)])
        Kav = x[5:16]
        Kavd = FcRegression.murineKavDist(; regularKav = true)
        Kavd[!, Not("IgG")] = typeof(Kav[1, 1]).(reshape(Kav, size(Kavd)[1], :))
        title = (title === nothing) ? "Murine fitting with ADVI (sampled median)" : title
    end

    ndf = predictMurine(df; Kav = Kavd, KxStar = KxStar, recepExp = Rtotd)
    plotPredvsMeasured(ndf; xx = "Value", yy = "Predict", color = "Receptor", shape = "Subclass", 
        clip2one = false, R2pos = (0, -1.2), title = title, kwargs...)
end

function MAPmurineLikelihood(df = importMurineInVitro(); KxStar = KxConst)
    m = murineFit(df, df."Value"; KxStar = KxStar)
    opts = Optim.Options(iterations = 500, show_every = 10, show_trace = true)
    opt = optimize(m, MAP(), LBFGS(; m = 20), opts)
    x = opt.values

    Rtot = Dict([rcp => x[Symbol("Rtot[$i]")] for (i, rcp) in enumerate(murineFcgR)])
    Kav = murineKavDist(; regularKav = true, retdf = true)
    Kav[!, Not("IgG")] = reshape([x[Symbol("Kav[$i]")] for i = 1:12], 3, 4)

    ndf = predictMurine(df; Kav = Kav, recepExp = Rtot, KxStar = KxStar)
    pl = plotPredvsMeasured(
        ndf;
        xx = "Value",
        yy = "Predict",
        color = "Receptor",
        shape = "Subclass",
        clip2one = false,
        R2pos = (-0.5, -1.2),
        title = "Murine MAP fitting results",
    )
    return pl, [Rtot, Kav]
end




function importMurineLeukocyte(fn = "leukocyte-apr2022.csv")
    df = CSV.File(joinpath(dataDir, fn), comment = "#") |> DataFrame
    #df = df[df."Cell" .!= "Tcell", :]   # throw away T cells
    df = df[df."Experiment" .!= 7, :]   # throw away Exp. 7
    df[!, ["IgG1", "IgG2b", "IgG2c"]] .-= df[!, "TNP-BSA"]  # subtract TNP-BSA control from meas
    df = df[!, Not("TNP-BSA")]
    df = dropmissing(stack(df, ["IgG1", "IgG2b", "IgG2c"], variable_name = "Subclass", value_name = "Value"))
    df[df."Value" .< 1.0, "Value"] .= 1.0   # clip values to 1.0
    baseline = combine(groupby(df, "Experiment"), "Value" => geomean => "Baseline")
    df = innerjoin(df, baseline, on = "Experiment")
    df[!, "Value"] ./= df[!, "Baseline"]    # normalize fluorescence by daily geomean
    df = df[!, Not(["Experiment", "Baseline"])]
    return sort!(df, ["Cell", "Subclass", "Valency"])
end

function predictLeukocyte(dfr::DataFrameRow; KxStar = KxConst, Kav = importKav(; murine = true, retdf = true))
    igg = dfr."Subclass"
    igg = (igg == "IgG2c") ? "IgG2a" : igg
    Kav_vs = Matrix(Kav[Kav."IgG" .== igg, Not("IgG")])
    Rtot_vs = Vector(dfr[murineFcgR])
    return polyfc(1e-9, KxStar, dfr."Valency", Rtot_vs, [1.0], Kav_vs).Lbound
end

function predictLeukocyte(df::AbstractDataFrame = importMurineLeukocyte(); average = false, title = "Murine Leukocyte Predictions", kwargs...)
    if average
        df = combine(
            groupby(df, Not("Value")),
            "Value" => geomean => "Value",
            "Value" => (xs -> quantile(xs, 0.25)) => "xmin",
            "Value" => (xs -> quantile(xs, 0.75)) => "xmax",
        )
    end

    Rtot = importRtot(; murine = true, retdf = true, cellTypes = unique(df."Cell"))
    Rtot = stack(Rtot, Not("Receptor"), variable_name = "Cell", value_name = "Abundance")
    Rtot = dropmissing(unstack(Rtot, "Cell", "Receptor", "Abundance"))
    rdf = innerjoin(df, Rtot, on = "Cell")
    sort!(rdf, ["Cell", "Subclass", "Valency"])
    @assert df."Value" == rdf."Value"
    preds = Vector(undef, size(rdf)[1])
    Threads.@threads for i = 1:size(rdf)[1]
        preds[i] = predictLeukocyte(rdf[i, :]; kwargs...)
    end
    rdf."Predict" = preds
    @assert all(isfinite(rdf[!, "Predict"]))
    rdf[rdf."Predict" .<= 0.0, "Predict"] .= 1e-8
    rdf = rdf[!, Not(murineFcgR)]

    rdf[rdf."Predict" .<= 0.0, "Predict"] .= 1e-12
    rdf."Predict" ./= geomean(rdf."Predict") / geomean(rdf."Value")
    pldf = deepcopy(rdf)
    pldf."Cell" = replace.(pldf."Cell", cellTypeFullName...)
    pl = plotPredvsMeasured(pldf; xx = "Value", yy = "Predict", 
        color = "Cell", shape = "Valency", clip2one = false, 
        R2pos = (-0.5, -2), title = title)
    return rdf, pl
end
