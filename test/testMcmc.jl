import Turing: sample, Priors

@testset "mcmc.jl tests" begin
    @testset "Test the MCMC model can make sane Priors" begin
        df = FcRegression.loadMixData()
        m = FcRegression.sfit(df, df."Value")
        c = sample(m, Prior(), 2)
    end
end
