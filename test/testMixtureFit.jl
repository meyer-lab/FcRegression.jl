using ForwardDiff
using Random
import Turing: sample, MH, NUTS

@testset "Building the MCMC model can work" begin
    rng = MersenneTwister(1234);
    df = FcRegression.loadMixData()
    model = FcRegression.sfit(df, df."Value")
    model(rng)
end
