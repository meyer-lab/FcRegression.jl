""" Fitting the mixture measurements. """

using Optim
import Turing: optimize, MAP, MLE


function mixturePredictions(
    df = loadMixData();
    Rtot = measuredRecepExp,
    Kav = importKav(; murine = false, invitro = true, retdf = true),
    KxStar = KxConst,
    vals = [4.0, 33.0],
    convs = [2.27, 3.26],
)
    df[!, "NewValency"] .= vals[1]
    df[df."Valency" .== 33, "NewValency"] .= vals[2]

    ndf = predictMix(df; recepExp = Rtot, KxStar = KxStar, Kav = Kav)

    ndf[ndf."Valency" .== 4, "Predict"] ./= convs[1]
    ndf[ndf."Valency" .== 33, "Predict"] ./= convs[2]
    return ndf
end

const f4conv_dist = LogNormal(log(2.27), 0.5)   # std ~= 1.37
const f33conv_dist = LogNormal(log(3.26), 0.5)  #std ~= 1.95

function MAPLikelihood(df; robinett = false)
    model = sfit(df, df."Value"; robinett = robinett)
    opts = Optim.Options(iterations = 1000, show_every = 10, show_trace = true)

    opt = optimize(model, MAP(), LBFGS(; m = 20), opts)
    x = opt.values.array

    Rtot = Dict([humanFcgRiv[ii] => x[ii] for ii = 1:length(humanFcgRiv)])
    Kav = deepcopy(importKav(; murine = false, invitro = true, retdf = true))
    Kav[!, Not("IgG")] = reshape(x[10:33], length(humanIgG), length(humanFcgRiv))

    return mixturePredictions(df; Rtot = Rtot, Kav = Kav, KxStar = x[9], vals = [x[7], x[8]])
end
