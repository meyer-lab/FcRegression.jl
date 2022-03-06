using ForwardDiff
import Turing: sample, MH, NUTS

@testset "Building the MCMC model can work" begin
    model = FcRegression.sfit()
end
