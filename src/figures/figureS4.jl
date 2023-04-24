""" Supplemental figure to related to CD16b """

function importCD16b()
    df = CSV.File(joinpath(dataDir, "CHO_FcgR3b_binding.csv"), delim = ",", comment = "#") |> DataFrame
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
end


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


function figureS4(; ssize = (6inch, 3inch), widths = [3 4])
    setGadflyTheme()

    # FcgRIIIB: IgG2 ~= 60000, IgG4 ~= 80000

    c = inferCD16b(; mcmc_iter = 5000)
    af1, af2 = plotCD16bCHOaff(c[4000:end])

    draw(
        PDF("output/figureS4.pdf", ssize...),
        plotGrid(
            (1, 2),
            [af1, af2];
            widths = widths,
        )
    )
end