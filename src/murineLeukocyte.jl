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

function predictLeukocyte(dfr::DataFrameRow; KxStar = KxConst, 
        Kav = murineKavDist(; regularKav = true, retdf = true), f = [4, 33])
    igg = dfr."Subclass"
    igg = (igg == "IgG2c") ? "IgG2a" : igg
    Kav_vs = Matrix(Kav[Kav."IgG" .== igg, Not("IgG")])
    Rtot_vs = Vector(dfr[murineFcgR])
    val = (dfr."Valency" > 12) ? f[2] : f[1]
    return polyfc(1e-9, KxStar, val, Rtot_vs, [1.0], Kav_vs).Lbound
end

function predictLeukocyte(df::AbstractDataFrame = importMurineLeukocyte(; average = true); kwargs...)
    # transpose Rtot
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

function predictLeukocyte(c::Chains, df::AbstractDataFrame; Kav = nothing, kwargs...)   
    f4 = median(c["f4"])
    f33 = median(c["f33"])

    Kavd = murineKavDist(; regularKav = true, retdf = true)
    Kav = [median(c["Kav[$i]"].data) for i = 1:12]
    Kavd[!, Not("IgG")] = typeof(Kav[1, 1]).(reshape(Kav, size(Kavd)[1], :))
    
    return predictLeukocyte(df; Kav = Kavd, f = [f4, f33], kwargs...)
end


function plot_murine_leukocyte(ndf; kwargs...)
    @assert "Predict" in names(ndf)
    pldf = deepcopy(ndf)
    pldf."Cell" = replace.(pldf."Cell", cellTypeFullName...)
    return plotPredvsMeasured(pldf; xx = "Value", yy = "Predict", 
        color = "Cell", shape = "Valency", R2pos = (0, -1.5), kwargs...)
end

@model function mLeukocyteModel(df, values; KxStar = KxConst)
    # fit Kav
    Kav_dist = Matrix(murineKavDist()[:, Not("IgG")])
    Kavd = murineKavDist(; regularKav = true)
    Kavd = Kavd[Kavd."IgG" .!= "IgG3", :]
    Kav = Matrix(undef, size(Kav_dist)...)
    for ii in eachindex(Kav)
        Kav[ii] ~ Kav_dist[ii]
    end
    Kavd[!, Not("IgG")] = typeof(Kav[1, 1]).(Kav)

    # fit valency
    f4 ~ f4Dist
    f33 ~ f33Dist
    
    # fit predictions
    if all(0.0 .< [f4, f33] .< Inf)
        df = predictLeukocyte(deepcopy(df); Kav = Kavd, KxStar = KxStar, f = [f4, f33])
    else
        df = deepcopy(df)
        df."Predict" .= Inf
    end

    stdv = std(log.(df."Predict") - log.(values))
    values ~ MvLogNormal(log.(df."Predict"), stdv * I)
    nothing
end

function fitLeukocyteMCMC(fname = "leukNUTSfit_0505.dat"; mcmc_iter = 1_000, KxStar = KxConst)
    if isfile(fname)
        return deserialize(fname)
    end
    df = importMurineLeukocyte(; average = false)
    m = mLeukocyteModel(df, df."Value"; KxStar = KxStar)
    # use MAP estimate as starting point
    opts = Optim.Options(iterations = 200, show_every = 10, show_trace = true)
    opt = optimize(m, MAP(), LBFGS(; m = 20), opts)
    c = sample(m, NUTS(), mcmc_iter, init_params = opt.values.array)

    f = serialize(fname, c)
    return c
end

#=
function validateLeukocyte(c = runMurineMCMC(); KxStar = KxConst)
    df = importMurineLeukocyte(; average = true)

    Kav1 = importKav(; murine = true, retdf = true)
    ndf1 = predictLeukocyte(deepcopy(df); KxStar = KxStar, 
        f = fitLeukocyteMCMC(; KxStar = KxStar))
    pl1 = plot_murine_leukocyte(ndf1; title = "Murine Leukocyte Raw Predictions")

    Kavd = murineKavDist(; regularKav = true)
    Kav = [median(c["Kav[$i]"].data) for i = 1:12]
    Kavd[!, Not("IgG")] = typeof(Kav[1, 1]).(reshape(Kav, size(Kavd)[1], :))

    ndf2 = predictLeukocyte(deepcopy(df); Kav = Kavd, KxStar = KxStar, 
        f = fitLeukocyteMCMC(; KxStar = KxStar, Kav = Kavd))
    pl2 = plot_murine_leukocyte(ndf2; title = "Murine Leukocyte Updated Predictions")
    return pl1, pl2
end
=#
