
function MAPLikelihood(df; robinett = false)
    model = sfit(df, df."Value"; robinett = robinett)
    opts = Optim.Options(iterations = 1000, show_every = 10, show_trace = true)

    opt = optimize(model, MAP(), LBFGS(; m = 20), opts)
    x = opt.values.array

    Rtot = Dict([humanFcgRiv[ii] => x[ii] for ii = 1:length(humanFcgRiv)])
    
    Kav = deepcopy(importKavDist(; murine = false, regularKav = true, retdf = true))
    Kav[!, Not("IgG")] = reshape(x[10:33], length(humanIgG), length(humanFcgRiv))

    return predMix(df; Rtot = Rtot, Kav = Kav, KxStar = x[9], fs = [x[7], x[8]])
end


