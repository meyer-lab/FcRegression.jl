""" Import Apr 2022 murine in vitro data """
function importMurineInVitro(fn = "CHO-mFcgR-apr2022.csv")
    df = CSV.File(joinpath(dataDir, fn), comment = "#") |> DataFrame
    df = df[df."Experiment" .!= 5, :]   # throw away Exp. 5
    df[!, ["IgG1", "IgG2b", "IgG2c"]] .-= df[!, "TNP-BSA"]  # subtract TNP-BSA control from meas
    df = df[df."Receptor" .!= "CHO", Not("TNP-BSA")]
    df = dropmissing(stack(df, ["IgG1", "IgG2b", "IgG2c"]))
    rename!(df, "variable" => "Subclass")
    rename!(df, "value" => "Value")
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
        return inferLogNormal(x, x * 10)
    end
    Kav[!, Not("IgG")] = retDist.(Kav[!, Not("IgG")], regularKav = regularKav)
    if !retdf
        return Matrix(Kav[!, Not("IgG")])
    end
    return Kav
end

murineKavDist(; kwargs...) = deepcopy(murineKavDist_nested(; kwargs...))

const InVitroMurineRcpExp = Dict("FcgRI" => 1e5, "FcgRIIB" => 1e7, "FcgRIII" => 1e6, "FcgRIV" => 1e5)

function predictMurine(dfr::DataFrameRow; KxStar = KxConst, recepExp = InVitroMurineRcpExp)
    if dfr."Subclass" == "TNP-BSA" || dfr."Affinity" <= 0.0
        return 0.0
    end
    return polyfc(1e-9, KxStar, dfr."Valency", [recepExp[dfr."Receptor"] * dfr."Expression" / 100], [1.0], reshape([dfr."Affinity"], 1, 1)).Lbound
end

function predictMurine(df::AbstractDataFrame; Kav = murineKavDist(; regularKav = true), conv = 20561, kwargs...)
    # Add affinity information to df
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
    @assert all(isfinite(dft[!, "Predict"]))
    dft[!, "Predict"] ./= conv
    dft[dft."Predict" .<= 0.0, "Predict"] .= 1e-6
    return dft
end

@model function murineFit(df, values)
    # Rtot sample
    Rtot = Vector(undef, length(murineFcgR))
    for ii in eachindex(Rtot)
        ref = InVitroMurineRcpExp[murineFcgR[ii]]
        Rtot[ii] ~ truncated(inferLogNormal(ref, ref * 1e2), 1e3, 1e9)
        # Treat receptor amount as unknown with a wide prior
    end
    Rtotd = Dict([murineFcgR[ii] => Rtot[ii] for ii = 1:length(Rtot)])

    # Kav sample
    Kav_dist = Matrix(murineKavDist()[:, Not("IgG")])
    Kavd = murineKavDist(; regularKav = true)
    Kavd = Kavd[Kavd."IgG" .!= "IgG3", :]
    Kav = Matrix(undef, size(Kav_dist)...)
    for ii in eachindex(Kav)
        Kav[ii] ~ truncated(Kav_dist[ii], 1e2, 1E10)
    end
    Kavd[!, Not("IgG")] = Kav

    KxStar ~ truncated(KxStarDist, 1E-18, 1E-9)
    # conversion factor: 20561
    conv ~ truncated(inferLogNormal(20561, 20561), 1e-6, 1e8)

    # fit predictions
    if all(Rtot .> 0.0) && all(Kav .> 0.0)
        df = predictMurine(deepcopy(df); recepExp = Rtotd, Kav = Kavd, KxStar = KxStar, conv = conv)
    else
        df = deepcopy(df)
        df."Predict" .= -1000.0
    end

    stdv = std(log.(df."Predict") - log.(values))
    values ~ MvLogNormal(log.(df."Predict"), stdv * I)
    nothing
end

function runMurineMCMC(fname = "murine_ADVI_0423.dat")
    if isfile(fname)
        return deserialize(fname)
    end
    df = importMurineInVitro()
    m = murineFit(df, df."Value")
    q = vi(m, ADVI(10, 1000))
    f = serialize(fname, q)
    return q
end

function plot_murineMCMC_dists(c::Union{Chains, MultivariateDistribution} = runMurineMCMC())
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
        Kav_pls[ii] = plotHistPriorDist(dat, Kav[(ii-1)%ligg+1, FcRname], name)
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
        Rtot_pls[ii] = plotHistPriorDist(dat, Rtot_dist[ii], FcRname)
    end
    Rtot_plot = plotGrid((1, lfcr), Rtot_pls; sublabels = false)
    draw(SVG("MCMCmurine_Rtot.svg", 16inch, 4inch), Rtot_plot)
    draw(PDF("MCMCmurine_Rtot.pdf", 16inch, 4inch), Rtot_plot)

    # Plot conv, KxStar
    other_pls = Vector{Plot}(undef, 2)
    KxStar_dat = (c isa Chains) ? c["KxStar"].data : cc[17, :]
    other_pls[1] = plotHistPriorDist(KxStar_dat, KxStarDist, "K<sub>x</sub><sup>*</sup>")
    conv_dat = (c isa Chains) ? c["conv"].data : cc[18, :]
    other_pls[2] = plotHistPriorDist(conv_dat, inferLogNormal(20561, 20561), "Conversion factor")
    other_plot = plotGrid((1, 2), other_pls; sublabels = false)
    draw(SVG("MCMCmurine_others.svg", 12inch, 4inch), other_plot)
    draw(PDF("MCMCmurine_others.pdf", 12inch, 4inch), other_plot)
end

function murineMCMC_params_predict_plot(c = runMurineMCMC(), df = importMurineInVitro(); kwargs...)
    # Extract MCMC results
    Rtot = [median(c["Rtot[$i]"].data) for i = 1:length(murineFcgR)]
    Rtotd = Dict([murineFcgR[ii] => Rtot[ii] for ii = 1:length(murineFcgR)])

    Kavd = murineKavDist(; regularKav = true)
    Kav = [median(c["Kav[$i]"].data) for i = 1:length(murineKavDist(; regularKav = true, retdf = false))]
    Kavd[!, Not("IgG")] = typeof(Kav[1, 1]).(reshape(Kav, size(Kavd)[1], :))

    KxStar = median(c["KxStar"].data)
    conv = median(c["conv"].data)

    ndf = predictMurine(df; Kav = Kavd, conv = conv, KxStar = KxStar, recepExp = Rtotd)
    plotPredvsMeasured(ndf; xx = "Value", yy = "Predict", color = "Receptor", shape = "Subclass", kwargs...)
end

function MAPmurineLikelihood()
    df = importMurineInVitro()
    m = murineFit(df, df."Value")
    opts = Optim.Options(iterations = 1000, show_every = 10, show_trace = true)
    opt = optimize(m, MAP(), LBFGS(; m = 20), opts)
    x = opt.values

    Rtot = Dict([rcp => x[Symbol("Rtot[$i]")] for (i, rcp) in enumerate(murineFcgR)])
    Kav = murineKavDist(; regularKav = true, retdf = true)
    Kav[!, Not("IgG")] = reshape([x[Symbol("Kav[$i]")] for i = 1:12], 3, 4)
    KxStar, conv = x[:KxStar], x[:conv]

    ndf = predictMurine(df; Kav = Kav, recepExp = Rtot, KxStar = KxStar, conv = conv)
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
    return pl, [Rtot, Kav, KxStar, conv]
end

function murineADVI_predict_plot(q, df = importMurineInVitro())
    x = median(rand(q, 100000), dims = 2)
    Rtot = x[1:4]
    Rtotd = Dict([FcRegression.murineFcgR[ii] => Rtot[ii] for ii = 1:length(FcRegression.murineFcgR)])
    Kav = x[5:16]
    Kavd = FcRegression.murineKavDist(; regularKav = true)
    Kavd[!, Not("IgG")] = typeof(Kav[1, 1]).(reshape(Kav, size(Kavd)[1], :))
    KxStar, conv = x[17], x[18]

    ndf = predictMurine(df; Kav=Kavd, recepExp = Rtotd, KxStar = KxStar, conv = conv)
    return plotPredvsMeasured(ndf; xx = "Value", yy = "Predict", 
        color = "Receptor", shape = "Subclass", clip2one = false, 
        R2pos = (0, -1.2), title = "Murine fitting with ADVI (sampled median)")
end
