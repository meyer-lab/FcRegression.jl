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
    factors = zeros(length(exps))
    df = predictMix(df; recepExp = recepExp, KxStar = KxStar)

    if !("Adjusted" in names(df))
        df[!, "Adjusted"] .= df[!, "Value"]
    end

    for (ii, exp) in enumerate(exps)
        factors[ii] = ols(df[df."Experiment" .== exp, "Adjusted"], df[df."Experiment" .== exp, "Predict"])[1]
        df[df."Experiment" .== exp, "Adjusted"] .*= factors[ii]
    end
    return factors, df
end

function fitMixFunc2(x, df)
    ## assemble numbers to a vector for optimization
    # order: log(Rtot), log(Kx*)
    # no valency yet

    cells = sort(unique(df."Cell"))
    recepExp = Dict([cell => exp(x[ii]) for (ii, cell) in enumerate(cells)])
    KxStar = exp(x[(length(cells)+1)])
    return mixSqLoss(fitExperiment(df; recepExp = recepExp, KxStar = KxStar)[2])
end

function fitMixMaster()
    # order: Rtot, Kx*
    # no valency yet

    df = loadMixData()
    cells = sort(unique(df."Cell"))
    x0 = [log(measuredRecepExp[cell]) for cell in cells]
    append!(x0, log(KxConst))
    x_ub = x0 .+ 2.0
    x_lb = x0 .- 2.0
    println(x0)
    println(x_ub)
    println(x_lb)
    f = x -> fitMixFunc2(x, df)
    od = OnceDifferentiable(f, x0; autodiff = :forward)
    #res = optimize(f, x0, GradientDescent(), Optim.Options(show_trace=true, iterations = 100))
    res = optimize(f, x_lb, x_ub, x0; iterations = 100)#, Fminbox(GradientDescent()); show_trace=true, iterations = 100)
    #res = optimize(f, x0, x_lb, x_ub; show_trace=true, iterations = 100)
    return res
end
