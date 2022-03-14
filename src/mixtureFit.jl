""" Fitting the mixture measurements. """

using Optim
import Turing: optimize, MAP, MLE


function mixturePredictions(
    df = loadMixData();
    Rtot = measuredRecepExp,
    Kav = importKav(; murine = false, invitro = true, retdf = true),
    KxStar = KxConst,
    vals = [4.0, 33.0],
)
    @assert all(isfinite(df."Value"))
    df[!, "NewValency"] .= vals[1]
    df[df."Valency" .== 33, "NewValency"] .= vals[2]

    ndf = predictMix(df; recepExp = Rtot, KxStar = KxStar, Kav = Kav)
    @assert all(isfinite(ndf."Predict"))

    # Least squares with one var and no intercept
    scale = sum(ndf."Value" .* ndf."Predict") / sum(ndf."Value" .* ndf."Value")
    ndf."Predict" *= scale

    return ndf
end


function MAPLikelihood(df; robinett=false)
    model = sfit(df, df."Value"; robinett=robinett)
    opts = Optim.Options(iterations=1000, show_every=10, show_trace=true)

    opt = optimize(model, MAP(), LBFGS(; m = 20), opts)
    x = opt.values.array

    Rtot = Dict([humanFcgRiv[ii] => x[ii] for ii = 1:length(humanFcgRiv)])
    Kav = deepcopy(importKav(; murine = false, invitro = true, retdf = true))
    Kav[!, Not("IgG")] = reshape(x[10:33], length(humanIgG), length(humanFcgRiv))
    
    return mixturePredictions(df; Rtot = Rtot, Kav = Kav, KxStar = x[9], vals = [x[7], x[8]])
end
