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

function fitMixFunc(x::Vector, df; fitKav = false)
    ## assemble numbers to a vector for optimization
    # order: log(Rtot), log(valency), log(Kx*), log(Kav)

    Kav = importKav(; murine = false, retdf = true)
    cells = sort(unique(df."Cell"))
    @assert length(x) == length(cells) + 3 + (fitKav ? size(Kav)[1] * (size(Kav)[2] - 1) : 0)
    x = exp.(x)

    recepExp = Dict([cell => x[ii] for (ii, cell) in enumerate(cells)])
    # adjust the valency, given there are only 4 and 33 originally
    df[!, "NewValency"] .= x[(length(cells) + 1)]
    df[df."Valency" .== 33, "NewValency"] .= x[(length(cells) + 2)]
    KxStar = x[(length(cells) + 3)]
    if fitKav
        Kav[!, Not("IgG")] .= 0.0
        Kav[!, Not("IgG")] = reshape(x[(length(cells) + 4):end], (size(Kav)[1], :))
    end

    return fitExperiment(df; recepExp = recepExp, KxStar = KxStar, Kav = Kav)
end

function fitMixMaster(df = loadMixData(); fitKav = false, recepExp = measuredRecepExp)
    # order: log(Rtot), log(valency), log(Kx*), log(Kav)
    # x0 for Rtot, valency, Kx*
    cells = sort(unique(df."Cell"))
    x0 = [log(recepExp[cell]) for cell in cells]
    append!(x0, log.([4, 33]))
    append!(x0, log(KxConst))

    # x0 for Kav
    if fitKav
        Kav = importKav(; murine = false, retdf = true)
        x0aff = reshape(Matrix(Kav[!, Not("IgG")]), (:))
        x0aff[x0aff .<= 0.0] .= 1.0
        append!(x0, log.(x0aff))
    end

    # set bounds for optimization
    x_ub = x0 .+ 3.0
    x_lb = x0 .- 3.0
    x_ub[(length(cells) + 1):(length(cells) + 2)] .= x0[(length(cells) + 1):(length(cells) + 2)] .+ 0.01   # valency ub is the original valency
    x_lb[(length(cells) + 1):(length(cells) + 2)] .= x0[(length(cells) + 1):(length(cells) + 2)] .- 1.0    # valency lb is exp(1) below
    println(x0, x_ub, x_lb)
    f = x -> mixSqLoss(fitMixFunc(x, df; fitKav = fitKav))

    dfc = TwiceDifferentiableConstraints(x_lb, x_ub)
    res = optimize(f, dfc, x0, IPNewton(), Optim.Options(iterations = 100, show_trace = true); autodiff = :forward)
    ndf = fitMixFunc(Optim.minimizer(res), averageMixData(df); fitKav = fitKav)
    if fitKav
        Kav[!, Not("IgG")] .= 0.0
        Kav[!, Not("IgG")] = reshape(exp.(res.minimizer[(length(cells) + 4):end]), (size(Kav)[1], :))
        return res, ndf, Kav
    end
    return res, ndf
end
