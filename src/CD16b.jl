function importCD16b()
    df = CSV.File(joinpath(dataDir, "dec2022-FcgR3b-binding.csv"), delim = ",", comment = "#") |> DataFrame
    df = stack(df, Not(["Valency", "Cell", "Subclass"]), variable_name = "Experiment", value_name = "Value")
    df = dropmissing(df)
    df[!, "Value"] = convert.(Float64, df[!, "Value"])
    df[(df[!, "Value"]) .< 1.0, "Value"] .= 1.0

    expression = DataFrame(Dict("Experiment" => unique(df."Experiment"),
                            "CHO-FcgRIIIB-NA1" => [77.1, 90.8, 94.1, 90.8, 95.8, 98.0], 
                            "CHO-FcgRIIIB-NA2" => [53.3, 81.7, 91.8, 86.6, 92.4, 92.9]))
    # Expression data for "07.10.19 AM" was not available; used geomean of others to impute
    expression = stack(expression, Not("Experiment"), variable_name = "Cell", value_name = "Expression")
    df = innerjoin(df, expression, on = ["Cell", "Experiment"])

    rename!(df, "Cell" => "Receptor")
    replace!(df."Receptor", "CHO-FcgRIIIB-NA1" => "FcgRIIIB-NA1")
    replace!(df."Receptor", "CHO-FcgRIIIB-NA2" => "FcgRIIIB-NA2")
    return df
end

@model function gmodelCD16b(df, values)
    CD16bRtot = Dict("FcgRIIIB-NA1" => 100_000.0, "FcgRIIIB-NA2" => 100_000.0)

    # sample Kav
    Kavd = importKavDist(; murine = false, regularKav = true, retdf = true, CD16b = true)[:, [1,8,9]]
    Kav_dist = importKavDist(; murine = false, regularKav = false, retdf=false, CD16b = true)[:, [7,8]]
    Kav = Matrix(undef, size(Kav_dist)...)
    for ii in eachindex(Kav)
        Kav[ii] ~ Kav_dist[ii]
    end
    Kavd[!, Not("IgG")] = typeof(Kav[1, 1]).(Kav)

    # sample f4, f33, KxStar
    f4 ~ f4Dist
    f33 ~ f33Dist
    KxStar ~ KxStarDist

    # fit predictions
    if all(0.0 .<= Matrix(Kavd[!, Not("IgG")]) .< Inf) && all(0.0 .< [f4, f33, KxStar] .< Inf)
        df = predMix(deepcopy(df); Rtot = CD16bRtot, Kav = Kavd, KxStar = KxStar, fs = [f4, f33])
    else
        df = deepcopy(df)
        df."Predict" .= Inf
    end

    stdv = std(log.(df."Predict") - log.(values))
    values ~ MvLogNormal(log.(df."Predict"), stdv * I)
    nothing
end

function inferCD16b(; mcmc_iter = 1000)
    df = importCD16b()
    m = gmodelCD16b(df, df."Value")
    opts = Optim.Options(iterations = 500, show_every = 10, show_trace = true)
    opt = optimize(m, MAP(), LBFGS(; m = 20), opts)
    c = sample(m, NUTS(), mcmc_iter, init_params = opt.values.array)
    return c

    Kavd = importKavDist(; murine = false, regularKav = true, retdf = true, CD16b = true)[:, [1,8,9]]
    Kav = [median(c["Kav[$i]"].data) for i = 1:8]
    Kavd[!, Not("IgG")] = typeof(Kav[1, 1]).(reshape(Kav, size(Kavd)[1], :))
end

"""
 Row │ IgG      FcgRIIIB-NA1   FcgRIIIB-NA2        
─────┼───────────────────────────────────────
   1 │ IgG1         3.47807e5      3.75184e5
   2 │ IgG2     18524.8        25543.6
   3 │ IgG3         9.80651e5      8.82411e5
   4 │ IgG4      9481.83       21827.7
"""


function plotCD16bCHOaff(c)
    # extract
    Kav_priors = importKavDist(; murine = false, retdf = true, CD16b = true)[:, [1,8,9]]
    len = length(Matrix(Kav_priors[!, Not("IgG")]))
    Kav_posts = deepcopy(Kav_priors)
    Kav = [c["Kav[$i]"].data for i = 1:len]
    Kav_posts[!, Not("IgG")] = typeof(Kav[1, 1]).(reshape(Kav, size(Kav_posts)[1], :))

    pls = Vector{Union{Gadfly.Plot, Context}}(undef, size(Kav_priors)[2]-1)
    for (ii, fcr) in enumerate(names(Kav_priors)[2:end])
        priors = reshape(Matrix(Kav_priors[!, [fcr]]), :)
        posts = DataFrame(
            hcat([reshape(Kav_posts[i, fcr], :) for i = 1:(size(Kav_posts)[1])]...),
            Kav_posts."IgG",
        )
        fcr = replace(fcr, "FcgRIIIB" => "FcγRIIIB")
        pls[ii] = plot_distribution_violins(
            posts, 
            priors; 
            y_range = (4, 7),
            title = "$fcr Affinity Distributions",
            legend = (ii == length(names(Kav_priors))-1),
        )
    end
    return pls
end


@model function gmodelCD16bNeu(df, values; 
        Kav = extractNewHumanKav(), 
        Rtot = importRtot(; murine = false, retdf = true),
    )

    CD16bdist = importKavDist(; murine = false, regularKav = false, retdf = true, CD16b = true)[!, "FcgRIIIB-NA1"]
    CD16bKav = Vector(undef, length(CD16bdist))
    for ii in eachindex(CD16bKav)
        CD16bKav[ii] ~ CD16bdist[ii]
    end
    Kav[!, "FcgRIIIB"] .= CD16bKav

    f4 ~ f4Dist
    f33 ~ f33Dist
    KxStar ~ KxStarDist

    # fit predictions
    if all(0.0 .<= CD16bKav .< Inf) && all(0.0 .< [f4, f33, KxStar] .< Inf) && (f4 < f33)
        pred = predictLbound(Kav, Rtot; specificRcp = false, L0 = 1e-9, fs = [f4, f33], KxStar = KxStar,)
        dft = innerjoin(df, pred, on = ["Valency" => "Valency", "Subclass" => "IgG", "Cell" => "Cell"], makeunique=true)
    else
        dft = deepcopy(df)
        dft."Lbound" .= Inf
    end 

    stdv = std(log.(dft."Lbound") - log.(values))
    values ~ MvLogNormal(log.(dft."Lbound"), stdv * I)
    nothing
end


function CD16bHumanBlood(; mcmc_iter = 1000)
    # Only fit FcgRIIIB affinity using Neutrophils, with no allotype (NA1/NA2) distinction.
    df = importEffectorBind(; avg = false)
    df = df[df."Cell" .== "Neu", :]
    df[df."Value" .<= 1.0, "Value"] .= 1.0

    m = gmodelCD16bNeu(df, df."Value")
    opts = Optim.Options(iterations = 500, show_every = 10, show_trace = true)
    opt = optimize(m, MAP(), LBFGS(; m = 20), opts)
    c = sample(m, NUTS(), mcmc_iter, init_params = opt.values.array)
    return c

    
    affs = [median(c["CD16bKav[$i]"].data) for i = 1:4]
    Kav = extractNewHumanKav()
    Kav."FcgRIIIB" = [median(c["CD16bKav[$i]"].data) for i = 1:4]


    Kav1 = extractNewHumanKav()
    Rtot = importRtot(; murine = false, retdf = true)[!, ["Receptor", "Neu"]]
    df1 = predictLbound(Kav1, Rtot)

    f4 = median(c["f4"].data)
    f33 = median(c["f33"].data)
    KxStar = median(c["KxStar"].data)


    predictLbound(Kav, Rtot; specificRcp = false, L0 = 1e-9, fs = [f4, f33], KxStar = KxStar,)


    gmodelCD16bNeu(df, values)
end