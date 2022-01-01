using Turing
using DataFrames


@model function sfit(measurements = FcRegression.loadMixData(; discard_small = true)."Value")
    Rtot_dist = FcRegression.importInVitroRtotDist()
    Kav_dist = FcRegression.importKavDist(; inflation = 0.1)
    Kav_dist = Matrix(Kav_dist[:, Not("IgG")])

    lRtot = Vector(undef, length(Rtot_dist))
    lKav = Matrix(undef, size(Kav_dist)...)
    
    
    for ii in eachindex(lRtot)
        lRtot[ii] ~ Rtot_dist[ii]
    end
    for ii in eachindex(lKav)
        lKav[ii] ~ Kav_dist[ii]
    end

    lf4 ~ Normal(log(4), 0.2*log(4))
    lf33 ~ Normal(log(33), 0.2*log(33))
    lKxStar ~ Normal(log(FcRegression.KxConst), 4.37)
    
    Rtotd, vals, KxStar, Kav = FcRegression.dismantle_x0(exp.(vcat(lRtot, [lf4, lf33, lKxStar], reshape(lKav, :))))
    sigma = 0.1
    ps = FcRegression.mixturePredictions(;
        Rtot = Rtotd,
        Kav = Kav,
        KxStar = KxStar,
        vals = vals,
    )."Predict"
    for ii in eachindex(measurements)
        measurements[ii] ~ Normal(ps[ii], ps[ii] * sigma)
    end
end



m = sfit()

c = sample(m, MH(), 1000, discard_initial = 500);

c2 = sample(m, Gibbs(MH()), 1000, discard_initial = 500);

c = sample(m, NUTS(), MCMCThreads(), 1000, 4)

using Serialization
f = serialize("MCMC_run20211231.dat", c)

c′ = deserialize("MCMC_run20211231.dat")