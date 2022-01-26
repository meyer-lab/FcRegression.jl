using ForwardDiff

@testset "MLELikelihood() can take ForwardDiff" begin
    x0 = log.(FcRegression.assemble_x0())
    f = lx -> -FcRegression.totalLikelihood(exp.(lx))
    @test f(x0 .+ 0.1) isa Float64
    @test ForwardDiff.gradient(f, x0) isa Vector{Float64}
end

@testset "dismantle_x0() and assemble_x0() work" begin
    x = FcRegression.assemble_x0(vals = [1.0, 2.0])
    FcRegression.dismantle_x0(x)
end