using ForwardDiff
import Turing: sample, MH, NUTS

@testset "MLELikelihood() can take ForwardDiff" begin
    x0 = log.(FcRegression.assemble_x0())
    f = lx -> -FcRegression.totalLikelihood(exp.(lx))
    @test f(x0 .+ 0.1) isa Float64
    @test ForwardDiff.gradient(f, x0) isa Vector{Float64}
end

@testset "Building the MCMC model can work" begin
    model = FcRegression.sfit()
end
