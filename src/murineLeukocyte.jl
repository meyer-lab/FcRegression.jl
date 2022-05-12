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


@model function mLeukocyteModel(df, values; Kavd::Union{Nothing, AbstractDataFrame} = nothing)
    # sample Rtot
    Rtotd = importRtotDist(:mLeuk; retdf = true)
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
        df = predMix(deepcopy(df); Rtot = Rtotd, Kav = Kavd, KxStar = KxStar, fs = [f4, f33])
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
