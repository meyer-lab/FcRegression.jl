function importMurineLeukocyte(fn = "leukocyte-apr2022.csv"; average = true)
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
    if average
        df = combine(
            groupby(df, Not("Value")),
            "Value" => geomean => "Value",
            "Value" => (xs -> quantile(xs, 0.25)) => "xmin",
            "Value" => (xs -> quantile(xs, 0.75)) => "xmax",
        )
    end
    return sort!(df, ["Cell", "Subclass", "Valency"])
end

function predictLeukocyte(dfr::DataFrameRow; KxStar = KxConst, Kav = murineKavDist(; regularKav = true, retdf = true), f = [4, 33])
    igg = dfr."Subclass"
    igg = (igg == "IgG2a") ? "IgG2c" : igg
    Kav_vs = Matrix(Kav[Kav."IgG" .== igg, Not("IgG")])
    Rtot_vs = Vector(dfr[murineFcgR])
    val = (dfr."Valency" > 12) ? f[2] : f[1]
    return polyfc(1e-9, KxStar, val, Rtot_vs, [1.0], Kav_vs).Lbound
end

function predictLeukocyte(df::AbstractDataFrame = importMurineLeukocyte(; average = true); Rtot = nothing, kwargs...)
    # transpose Rtot
    cellTypes = unique(df."Cell")
    if Rtot === nothing
        Rtot = importRtot(; murine = true, retdf = true, cellTypes = cellTypes)
    else
        Rtot = Rtot[!, ["Receptor"; names(Rtot)[in(cellTypes).(names(Rtot))]]]
    end
    Rtot = stack(Rtot, Not("Receptor"), variable_name = "Cell", value_name = "Abundance")
    Rtot = dropmissing(unstack(Rtot, "Cell", "Receptor", "Abundance"))
    rdf = innerjoin(df, Rtot, on = "Cell")
    sort!(rdf, ["Cell", "Subclass", "Valency"])
    @assert df."Value" == rdf."Value"
    preds = Vector(undef, size(rdf)[1])
    Threads.@threads for i = 1:size(rdf)[1]
        preds[i] = predictLeukocyte(rdf[i, :]; kwargs...)
    end
    rdf."Predict" = typeof(preds[1]).(preds)
    rdf = rdf[!, Not(murineFcgR)]

    # one conversion factor per valency
    rdf[rdf."Predict" .< 0.0, "Predict"] .= 0.0
    for val in unique(rdf."Valency")
        if any(rdf[rdf."Valency" .== val, "Predict"] .> 0.0)
            rows = (rdf."Valency" .== val) .& (rdf."Predict" .> 0.0)
            rdf[(rdf."Valency" .== val), "Predict"] ./= geomean(rdf[rows, "Predict"]) / geomean(rdf[rows, "Value"])
        end
    end
    return rdf
end

function extractLeukocyteKav(c::Chains; retdf = false)
    Kavd = murineKavDist(; regularKav = true, retdf = true)
    Kav = retdf ? [c["Kav[$i]"].data for i = 1:12] : [median(c["Kav[$i]"].data) for i = 1:12]
    Kavd[!, Not("IgG")] = typeof(Kav[1, 1]).(reshape(Kav, size(Kavd)[1], :))
    return Kavd
end

function predictLeukocyte(c::Chains, df::AbstractDataFrame; Kavd = nothing)
    @assert (Symbol("Kav[1]") in c.name_map[1]) ⊻ (Kavd isa AbstractDataFrame)    # either providing Kavd, or from chain
    Rtotd = importCellRtotDist()
    Rtotd = Rtotd[!, ["Receptor"; names(Rtotd)[in(unique(df."Cell")).(names(Rtotd))]]]
    Rtot = [median(c["Rtot[$i]"].data) for i = 1:24]
    Rtotd[!, Not("Receptor")] = typeof(Rtot[1, 1]).(reshape(Rtot, size(Rtotd)[1], :))

    f4 = median(c["f4"])
    f33 = median(c["f33"])
    KxStar = median(c["KxStar"])
    return predictLeukocyte(df; Rtot = Rtotd, Kav = ((Kavd !== nothing) ? Kavd : extractLeukocyteKav(c)), 
        KxStar = KxStar, f = [f4, f33])
end


function plot_murine_leukocyte(ndf; kwargs...)
    @assert "Predict" in names(ndf)
    pldf = deepcopy(ndf)
    pldf."Cell" = replace.(pldf."Cell", cellTypeFullName...)
    return plotPredvsMeasured(pldf; xx = "Value", yy = "Predict", color = "Cell", shape = "Valency", R2pos = (0, -1.5), kwargs...)
end

@model function mLeukocyteModel(df, values; Kavd::Union{Nothing, AbstractDataFrame} = nothing)
    # sample Rtot
    Rtotd = importCellRtotDist()
    Rtotd = Rtotd[!, ["Receptor"; names(Rtotd)[in(unique(df."Cell")).(names(Rtotd))]]]
    Rtot_dist = Matrix(Rtotd[:, Not("Receptor")])
    Rtot = Matrix(undef, size(Rtot_dist)...)
    for ii in eachindex(Rtot)
        Rtot[ii] ~ Rtot_dist[ii]
    end
    Rtotd[!, Not("Receptor")] = typeof(Rtot[1, 1]).(Rtot)

    if Kavd === nothing
        # sample Kav
        Kav_dist = Matrix(murineKavDist()[:, Not("IgG")])
        Kavd = murineKavDist(; regularKav = true)
        Kavd = Kavd[Kavd."IgG" .!= "IgG3", :]
        Kav = Matrix(undef, size(Kav_dist)...)
        for ii in eachindex(Kav)
            Kav[ii] ~ Kav_dist[ii]
        end
        Kavd[!, Not("IgG")] = typeof(Kav[1, 1]).(Kav)
    end

    # sample valency
    f4 ~ f4Dist
    f33 ~ f33Dist
    KxStar ~ KxStarDist

    # fit predictions
    if all(0.0 .<= Rtot .< Inf) && all(0.0 .<= Matrix(Kavd[!, Not("IgG")]) .< Inf) && all(0.0 .< [f4, f33, KxStar] .< Inf)
        df = predictLeukocyte(deepcopy(df); Rtot = Rtotd, Kav = Kavd, KxStar = KxStar, f = [f4, f33])
    else
        df = deepcopy(df)
        df."Predict" .= Inf
    end

    stdv = std(log.(df."Predict") - log.(values))
    values ~ MvLogNormal(log.(df."Predict"), stdv * I)
    nothing
end

function fitLeukocyteMCMC(fname = "leukNUTSfit_0509.dat"; mcmc_iter = 1_000, Kavd = nothing)
    if (fname !== nothing) && isfile(fname)
        return deserialize(fname)
    end
    df = importMurineLeukocyte(; average = false)
    m = mLeukocyteModel(df, df."Value"; Kavd = Kavd)
    # use MAP estimate as starting point
    opts = Optim.Options(iterations = 200, show_every = 10, show_trace = true)
    opt = optimize(m, MAP(), LBFGS(; m = 20), opts)
    c = sample(m, NUTS(), mcmc_iter, init_params = opt.values.array)

    if fname !== nothing    # don't save if file name provided is nothing
        f = serialize(fname, c)
    end
    return c
end


function plot_MCMCLeuk_dists(c = fitLeukocyteMCMC())
    setGadflyTheme()

    # Plot Kav's
    mIgGs = murineKavDist(; retdf = true)."IgG"
    ligg, lfcr = length(mIgGs), length(murineFcgR)
    Kav_dist = murineKavDist(; retdf = false)
    Kav_pls = Matrix{Plot}(undef, ligg, lfcr)
    for ii in eachindex(Kav_pls)
        IgGname = mIgGs[(ii - 1) % ligg + 1]
        FcRname = murineFcgR[(ii - 1) ÷ ligg + 1]
        name = "m$IgGname to m$FcRname"
        Kav_pls[ii] = plotHistPriorDist(c["Kav[$ii]"].data, Kav_dist[ii], name)
    end
    Kav_plot = plotGrid((ligg, lfcr), permutedims(Kav_pls, (2, 1)); sublabels = false)
    draw(PDF("MCMCLeuk_Kav.pdf", 11inch, 8inch), Kav_plot)

    # Plot Rtot's
    cellTypes = unique(importMurineLeukocyte(; average = true)."Cell")
    lcell, lfcr = length(cellTypes), length(murineFcgR)
    Rtotd = importCellRtotDist(; retdf = true)
    Rtotd = Matrix(Rtotd[!, names(Rtotd)[in(cellTypes).(names(Rtotd))]])
    Rtot_pls = Matrix{Plot}(undef, lcell, lfcr)
    for ii in eachindex(Rtot_pls)
        cellname = cellTypes[(ii - 1) % lcell + 1]
        FcRname = murineFcgR[(ii - 1) ÷ lcell + 1]
        name = "m$FcRname on $cellname"
        Rtot_pls[ii] = plotHistPriorDist(c["Rtot[$ii]"].data, Rtotd[ii], name)
    end
    Rtot_plot = plotGrid((lcell, lfcr), permutedims(Rtot_pls, (2, 1)); sublabels = false)
    draw(PDF("MCMCLeuk_Rtot.pdf", 11inch, 14inch), Rtot_plot)

    # Plot f4, f33, KxStar
    other_pls = Vector{Plot}(undef, 3)
    other_pls[1] = plotHistPriorDist(c["f4"].data, f4Dist, "f = 4 effective valency")
    other_pls[2] = plotHistPriorDist(c["f33"].data, f33Dist, "f = 33 effective valency")
    other_pls[3] = plotHistPriorDist(c["KxStar"].data, f33Dist, "KxStar")
    other_plot = plotGrid((1, 3), other_pls; sublabels = false)
    draw(PDF("MCMCLeuk_others.pdf", 9inch, 3inch), other_plot)
end
