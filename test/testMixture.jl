@testset "Test TwoDFit()" begin
    X = rand(5, 10)
    p = rand(5)
    q = rand(10)
    Y = X .* reshape(p, :, 1) .* reshape(q, 1, :)
    p0, q0 = FcRegression.TwoDFit(X, Y)
    @test all(p ./ p[1] .≈ p0)
    @test all(q ./ q[1] .≈ q0)
end
