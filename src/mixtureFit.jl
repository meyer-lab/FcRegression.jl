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

function fitExperiment(df; recepExp = measuredRecepExp, KxStar = KxConst)
    df = predictMix(df; recepExp = recepExp, KxStar = KxStar)
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

function fitMixFunc(x, df)
    ## assemble numbers to a vector for optimization
    # order: log(Rtot), log(valency), log(Kx*)

    cells = sort(unique(df."Cell"))
    recepExp = Dict([cell => exp(x[ii]) for (ii, cell) in enumerate(cells)])
    # adjust the valency, given there are only 4 and 33 originally
    df[!, "NewValency"] .= exp(x[(length(cells) + 1)])
    df[df."Valency" .== 33, "NewValency"] .= exp(x[(length(cells) + 2)])
    KxStar = exp(x[(length(cells) + 3)])
    return fitExperiment(df; recepExp = recepExp, KxStar = KxStar)
end

function fitMixMaster(df = loadMixData())
    # order: log(Rtot), log(valency), log(Kx*)

    cells = sort(unique(df."Cell"))
    x0 = [log(measuredRecepExp[cell]) for cell in cells]
    append!(x0, log.([4, 33]))
    append!(x0, log(KxConst))
    x_ub = x0 .+ 3.0
    x_lb = x0 .- 3.0
    x_ub[(length(cells)+1):(length(cells)+2)] .= x0[(length(cells)+1):(length(cells)+2)] .+ 0.01   # valency ub is the original valency
    x_lb[(length(cells)+1):(length(cells)+2)] .= x0[(length(cells)+1):(length(cells)+2)] .- 1.0    # valency lb is exp(1) below
    println(x0, x_ub, x_lb)
    f = x -> mixSqLoss(fitMixFunc(x, df))

    dfc = TwiceDifferentiableConstraints(x_lb, x_ub)
    res = optimize(f, dfc, x0, IPNewton(), Optim.Options(iterations = 100, show_trace = true); autodiff = :forward)
    ndf = fitMixFunc(Optim.minimizer(res), averageMixData(df))
    return res, ndf
end
