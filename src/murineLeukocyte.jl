function importMurineLeukocyte(fn = "leukocyte-apr2022.csv"; average = true)
    df = CSV.File(joinpath(dataDir, fn), comment = "#") |> DataFrame
    #df = df[df."ImCell" .!= "Tcell", :]   # throw away T cells
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
    return sort!(df, ["ImCell", "Subclass", "Valency"])
end

function predictLeukocyte(dfr::DataFrameRow; KxStar = KxConst, Kav = importKavDist(; murine = true, regularKav = true, retdf = true), f = [4, 33])
    igg = dfr."Subclass"
    igg = (igg == "IgG2a") ? "IgG2c" : igg
    Kav_vs = Matrix(Kav[Kav."IgG" .== igg, Not("IgG")])
    Rtot_vs = Vector(dfr[murineFcgR])
    val = (dfr."Valency" > 12) ? f[2] : f[1]
    return polyfc(1e-9, KxStar, val, Rtot_vs, [1.0], Kav_vs).Lbound
end

function predictLeukocyte(df::AbstractDataFrame = importMurineLeukocyte(; average = true); Rtot = nothing, kwargs...)
    # transpose Rtot
    cellTypes = unique(df."ImCell")
    if Rtot === nothing
        Rtot = importRtot(; murine = true, retdf = true, cellTypes = cellTypes)
    else
        Rtot = Rtot[!, ["Receptor"; names(Rtot)[in(cellTypes).(names(Rtot))]]]
    end
    Rtot = stack(Rtot, Not("Receptor"), variable_name = "ImCell", value_name = "Abundance")
    Rtot = dropmissing(unstack(Rtot, "ImCell", "Receptor", "Abundance"))
    rdf = innerjoin(df, Rtot, on = "ImCell")
    sort!(rdf, ["ImCell", "Subclass", "Valency"])
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

function predictLeukocyte(c::Chains, df::AbstractDataFrame; Kavd = nothing)
    p = extractMCMC(c; murine = true)
    @assert (p["Kav"] !== nothing) ⊻ (Kavd isa AbstractDataFrame)    # either providing Kavd, or from chain
    if Kavd === nothing
        Kavd = p["Kav"]
    end

    return predictLeukocyte(df; Rtot = p["Rtot"], Kav = Kavd, KxStar = p["KxStar"], f = [p["f4"], p["f33"]])
end


function plot_murine_leukocyte(ndf; kwargs...)
    @assert "Predict" in names(ndf)
    pldf = deepcopy(ndf)
    pldf."ImCell" = replace.(pldf."ImCell", cellTypeFullName...)
    return plotPredvsMeasured(pldf; xx = "Value", yy = "Predict", color = "ImCell", shape = "Valency", R2pos = (0, -1.5), kwargs...)
end

@model function mLeukocyteModel(df, values; Kavd::Union{Nothing, AbstractDataFrame} = nothing)
    # sample Rtot
    Rtotd = importCellRtotDist()
    Rtotd = Rtotd[!, ["Receptor"; names(Rtotd)[in(unique(df."ImCell")).(names(Rtotd))]]]
    Rtot_dist = Matrix(Rtotd[:, Not("Receptor")])
    Rtot = Matrix(undef, size(Rtot_dist)...)
    for ii in eachindex(Rtot)
        Rtot[ii] ~ Rtot_dist[ii]
    end
    Rtotd[!, Not("Receptor")] = typeof(Rtot[1, 1]).(Rtot)

    if Kavd === nothing
        # sample Kav
        Kav_dist = Matrix(importKavDist(; murine = true)[:, Not("IgG")])
        Kavd = importKavDist(; murine = true, regularKav = true)
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
