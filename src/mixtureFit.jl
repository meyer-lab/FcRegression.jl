""" Fitting mixture measurements """

using Optim

function mixSqLoss(df; logscale = true)
    # square root differences of model prediction and adjusted measurements
    if "Adjusted" in names(df)
        adj = logscale ? log.(df."Adjusted") : df."Adjusted"
    else
        adj = logscale ? log.(df."Value") : df."Value"
    end
    pred = logscale ? log.(df."Predict") : df."Predict"
    return sum((adj .- pred) .^ 2)
end

function fitExperiment(df; recepExp = measuredRecepExp, KxStar = KxConst, Kav::DataFrame = importKav(; murine = false, retdf = true))
    df = predictMix(df; recepExp = recepExp, KxStar = KxStar, Kav = Kav)
    if !("Experiment" in names(df))
        return df
    end
    exps = sort(unique(df."Experiment"))
    df[!, "Value"] = convert.(typeof(KxStar), df[!, "Value"])

    if !("Adjusted" in names(df))
        df[!, "Adjusted"] .= df[!, "Value"]
    end

    for exp in exps
        factor = ols(df[df."Experiment" .== exp, "Adjusted"], df[df."Experiment" .== exp, "Predict"])
        df[df."Experiment" .== exp, "Adjusted"] .*= factor
    end
    return df
end

function fitMixFunc(
    x::Vector, 
    df; 
    fitRs = true, 
    fitKav = false, 
    Kav::DataFrame = importKav(; murine = false, retdf = true),
    xres = nothing,
)
    ## assemble numbers to a vector for optimization
    # order: log(Rtot), log(valency), log(Kx*), log(Kav)
    # xres is the results from fitting just Rtot, valency, and Kx*

    cells = sort(unique(df."Cell"))
    @assert fitRs | fitKav
    RsLen = (fitRs ? length(cells) + 3 : 0)
    @assert length(x) == RsLen + (fitKav ? size(Kav)[1] * (size(Kav)[2] - 1) : 0)
    x = exp.(x)

    if fitRs
        xres = x[1:RsLen]
    else
        @assert length(xres) == length(cells) + 3
    end

    recepExp = Dict([cell => xres[ii] for (ii, cell) in enumerate(cells)])
    # adjust the valency, given there are only 4 and 33 originally
    df[!, "NewValency"] .= xres[(length(cells) + 1)]
    df[df."Valency" .== 33, "NewValency"] .= xres[(length(cells) + 2)]
    KxStar = xres[(length(cells) + 3)]

    if fitKav
        Kav[!, Not("IgG")] .= 0.0
        Kav[!, Not("IgG")] = reshape(x[(RsLen + 1):end], (size(Kav)[1], :))
    end

    return fitExperiment(df; recepExp = recepExp, KxStar = KxStar, Kav = Kav)
end

function fitMixMaster(
    df = loadMixData();
    fitRs = true,
    fitKav = false,
    recepExp = measuredRecepExp,
    Kav::DataFrame = importKav(; murine = false, retdf = true),
    xres = nothing,
    show_trace = false,
)
    # order: log(Rtot), log(valency), log(Kx*), log(Kav)
    # x0 for Rtot, valency, Kx*
    @assert fitRs | fitKav
    x0, x_ub, x_lb = Vector{Float64}(), Vector{Float64}(), Vector{Float64}()

    # fitting Rtot, valency, and Kx*
    if fitRs
        cells = sort(unique(df."Cell"))
        append!(x0, [log(recepExp[cell]) for cell in cells])
        append!(x0, log.([4, 33]))
        append!(x0, log(KxConst))
        # bounds
        x_ub = x0 .+ 3.0
        x_lb = x0 .- 3.0
        x_ub[(length(cells) + 1):(length(cells) + 2)] .= x0[(length(cells) + 1):(length(cells) + 2)] .+ 0.001   # valency ub is the original valency
        x_lb[(length(cells) + 1):(length(cells) + 2)] .= x0[(length(cells) + 1):(length(cells) + 2)] .- 1.0    # valency lb is exp(1) below 
    end

    # fitting Kav
    if fitKav
        x0aff = reshape(Matrix(Kav[!, Not("IgG")]), (:))
        x0aff[x0aff .<= 0.0] .= 1.0
        append!(x0, log.(x0aff))
        # bounds
        append!(x_ub, log.(x0aff) .+ 6.0)
        append!(x_lb, log.(x0aff) .- 6.0)
    end

    println(x0, x_ub, x_lb)
    f = x -> mixSqLoss(fitMixFunc(x, df; fitRs = fitRs, fitKav = fitKav, Kav = Kav, xres = xres))

    dfc = TwiceDifferentiableConstraints(x_lb, x_ub)
    res = optimize(f, dfc, x0, IPNewton(), Optim.Options(iterations = 100, show_trace = show_trace); autodiff = :forward)
    ndf = fitMixFunc(Optim.minimizer(res), averageMixData(df); fitRs = fitRs, fitKav = fitKav, Kav = Kav, xres = xres)
    if fitKav
        RsLen = (fitRs ? length(cells) + 3 : 0)
        Kav[!, Not("IgG")] .= 0.0
        Kav[!, Not("IgG")] = reshape(exp.(res.minimizer[(RsLen + 1):end]), (size(Kav)[1], :))
        nKav = unstack(stack(Kav, Not("IgG")), "variable", "IgG", "value")
        CSV.write(joinpath(dataDir, "fitted_human_new_affinity.csv"), nKav)
        return res, ndf, Kav
    end
    return res, ndf
end

function loadFittedKav(; retdf = true)
    df = CSV.File(joinpath(dataDir, "fitted_human_new_affinity.csv"), comment = "#") |> DataFrame
    df = stack(df; variable_name = "IgG", value_name = "Kav")
    df = unstack(df, "variable", "Kav")
    df = df[in(humanIgG).(df.IgG), :]
    if retdf
        return deepcopy(df[!, ["IgG"; humanFcgR]])
    else
        return deepcopy(Matrix{Float64}(df[!, humanFcgR]))
    end
end
