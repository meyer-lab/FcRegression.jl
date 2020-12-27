"""function TwoDFit(X::Matrix, Y::Matrix)
    """
    Fit X_ij * p_i * q_j ≈ Y_ij
    p[1] == 1 (fixed)
    length(p) == size(X, 1) - 1 == m
    length(q) == size(X, 2) == n
    v == vcat(p, q)
    """
    @assert size(X) == size(Y)
    m, n = size(X)
    f(p::Vector, q::Vector) = sum((reshape([1.0; p], :, 1) .* X .* reshape(q, 1, :) .- Y) .^ 2)
    f(v::Vector) = f(v[1:(m - 1)], v[m:end])
    init_v = ones(m + n - 1)
    od = OnceDifferentiable(f, init_v; autodiff = :forward)
    res = optimize(od, init_v, BFGS()).minimizer
    return [1.0; res[1:(m - 1)]], res[m:end]
end

@testset "Test TwoDFit()" begin
    X = rand(5, 10)
    p = rand(5)
    q = rand(10)
    Y = X .* reshape(p, :, 1) .* reshape(q, 1, :)
    p0, q0 = TwoDFit(X, Y)
    @test all(p ./ p[1] .≈ p0)
    @test all(q .* p[1] .≈ q0)
end"""
