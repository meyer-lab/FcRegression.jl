""" Fitting mixture measurements """

using Optim
import Turing: optimize, MAP, MLE

##### 
# Below are for the MLE approach
##### 

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


function dismantle_x0(x::Vector)
    # order: Rtot, vals, KxStar, Kav
    @assert all(x .>= 0.0) "In dismantle(): Received a negative number."
    x = deepcopy(x)
    Rtot = Dict([humanFcgRiv[ii] => popfirst!(x) for ii = 1:length(humanFcgRiv)])
    vals = [popfirst!(x), popfirst!(x)]
    KxStar = popfirst!(x)
    Kav = deepcopy(importKav(; murine = false, invitro = true, retdf = true))
    Kav[!, Not("IgG")] = convert.(eltype(x), Kav[:, Not("IgG")])
    Kav[!, Not("IgG")] = reshape(x, length(humanIgG), length(humanFcgRiv))
    return Rtot, vals, KxStar, Kav
end


function MAPLikelihood(df; robinett=false)
    model = sfit(df, df."Value"; robinett=robinett)
    opts = Optim.Options(iterations=1000, show_every=10, show_trace=true)

    opt = optimize(model, MAP(), LBFGS(; m = 20), opts)

    Rtot, vals, KxStar, Kav = FcRegression.dismantle_x0(opt.values.array[1:33])
    return FcRegression.mixturePredictions(df; Rtot = Rtot, Kav = Kav, KxStar = KxStar, vals = vals)
end
