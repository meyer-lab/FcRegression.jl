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
    exps = sort(unique(df."Experiment"))
    df = predictMix(df; recepExp = recepExp, KxStar = KxStar)

    df[!, "Value"] = convert.(typeof(KxStar), df[!, "Value"])

    if !("Adjusted" in names(df))
        df[!, "Adjusted"] .= df[!, "Value"]
    end

    for (ii, exp) in enumerate(exps)
        factor = ols(df[df."Experiment" .== exp, "Adjusted"], df[df."Experiment" .== exp, "Predict"])
        df[df."Experiment" .== exp, "Adjusted"] .*= factor
    end
    return df
end

function fitMixFunc2(x, df)
    ## assemble numbers to a vector for optimization
    # order: log(Rtot), log(Kx*)
    # no valency yet

    cells = sort(unique(df."Cell"))
    recepExp = Dict([cell => exp(x[ii]) for (ii, cell) in enumerate(cells)])
    KxStar = exp(x[(length(cells) + 1)])
    return mixSqLoss(fitExperiment(df; recepExp = recepExp, KxStar = KxStar))
end

function fitMixMaster()
    # order: Rtot, Kx*
    # no valency yet

    df = loadMixData()
    cells = sort(unique(df."Cell"))
    x0 = [log(measuredRecepExp[cell]) for cell in cells]
    append!(x0, log(KxConst))
    x_ub = x0 .+ 3.0
    x_lb = x0 .- 3.0
    f = x -> fitMixFunc2(x, df)
    res = optimize(f, x_lb, x_ub, x0, Fminbox(GradientDescent()), Optim.Options(iterations = 5, show_trace=true); autodiff = :forward)
    return res
end
