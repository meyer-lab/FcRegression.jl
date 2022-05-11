import Serialization: serialize, deserialize
using LinearAlgebra


const f4Dist = LogNormal(log(4), 0.2)   # std ~= 0.82
const f33Dist = LogNormal(log(33), 0.2)   # std ~= 6.80
const KxStarDist = LogNormal(log(KxConst), 2.0)   # std ~= 4.37 in Robinett

@model function sfit(df, values; robinett = false, Kavd::AbstractDataFrame = importKavDist(; murine = false, regularKav = true, retdf = true))
    Rtot_dist = importInVitroRtotDist(robinett)
    Rtot = Vector(undef, length(Rtot_dist))

    # Order of distribution definitions here matches MAPLikelihood()
    for ii in eachindex(Rtot)
        Rtot[ii] ~ Rtot_dist[ii]
    end
    Rtotd = Dict([humanFcgRiv[ii] => Rtot[ii] for ii = 1:length(humanFcgRiv)])

    if !robinett    # Don't fit affinity for Robinett data
        Kav_dist = importKavDist(; murine = false)
        Kav_dist = Matrix(Kav_dist[:, Not("IgG")])
        Kav = Matrix(undef, size(Kav_dist)...)

        for ii in eachindex(Kav)
            Kav[ii] ~ Kav_dist[ii]
        end
        Kavd[!, Not("IgG")] = typeof(Kav[1, 1]).(Kav)
    end

    f4 ~ f4Dist
    f33 ~ f33Dist
    KxStar ~ KxStarDist

    # fit predictions
    if all(Rtot .>= 0.0) && all(0.0 .<= Matrix(Kavd[!, Not("IgG")]) .< Inf) && all(0.0 .< [f4, f33, KxStar] .< Inf)
        df = predictMix(deepcopy(df); recepExp = Rtotd, Kav = Kavd, KxStar = KxStar, vals = [f4, f33])
    else
        df = deepcopy(df)
        df."Predict" .= Inf
    end

    stdv = std(log.(df."Predict") - log.(values))
    values ~ MvLogNormal(log.(df."Predict"), stdv * I)
    nothing
end

function runMCMC(fname = "humanNUTSfit_0505.dat"; mcmc_iter = 1_000)
    if isfile(fname)
        return deserialize(fname)
    end
    df = loadMixData()

    # only use single IgG
    df = df[(df."%_1" .== 1.0) .| (df."%_2" .== 1.0), :]
    m = sfit(df, df."Value")

    opts = Optim.Options(iterations = 500, show_every = 10, show_trace = true)
    opt = optimize(m, MAP(), LBFGS(; m = 20), opts)
    c = sample(m, NUTS(), mcmc_iter, init_params = opt.values.array)

    f = serialize(fname, c)
    return c
end

function MAPLikelihood(df; robinett = false)
    model = sfit(df, df."Value"; robinett = robinett)
    opts = Optim.Options(iterations = 1000, show_every = 10, show_trace = true)

    opt = optimize(model, MAP(), LBFGS(; m = 20), opts)
    x = opt.values.array

    Rtot = Dict([humanFcgRiv[ii] => x[ii] for ii = 1:length(humanFcgRiv)])
    
    Kav = deepcopy(importKavDist(; murine = false, regularKav = true, retdf = true))
    Kav[!, Not("IgG")] = reshape(x[10:33], length(humanIgG), length(humanFcgRiv))

    return predictMix(df; recepExp = Rtot, Kav = Kav, KxStar = x[9], vals = [x[7], x[8]])
end


