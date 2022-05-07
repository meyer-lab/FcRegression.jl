using ForwardDiff
using Random
import Turing: sample, MH, NUTS

@testset "Building the MCMC model can work" begin
    rng = MersenneTwister(1234)
    df = FcRegression.loadMixData()
    model = FcRegression.sfit(df, df."Value")
    model(rng)

    df = FcRegression.MAPLikelihood(df)
end

@testset "Check murine fitting model can work" begin
    df = FcRegression.importMurineInVitro()
    model = FcRegression.murineFit(df, df."Value")
    df = FcRegression.MAPmurineLikelihood(df)
end
