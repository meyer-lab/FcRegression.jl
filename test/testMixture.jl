function TwoDFit(X::Matrix, Y::Matrix)
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
    X = rand(5, 10) .* 10 .+ 1
    p = rand(5) .* 2 .+ 1
    q = rand(10) .* 2 .+ 1
    Y = X .* reshape(p, :, 1) .* reshape(q, 1, :)
    df = convert(DataFrame, X)

    df[!, "Valency"] .= 1:5
    df = stack(df, 1:10)
    rename!(df, "variable" => "Experiment")
    rename!(df, "value" => "Value")
    df[!, "Predict"] .= vec(Y)
    xdf = FcRegression.MixtureFitLoss(df, p, q)[2]
    @test all(xdf[!, "Predict"] .≈ xdf[!, "Adjusted"]) 

    res = FcRegression.MixtureFit(df)
    @test norm(p ./ p[1] .- [res["ValConv"][i] for i in 1:5], 2) < 1e-12
    @test norm(q .* p[1] .- [res["ExpConv"]["x" * string(i)] for i in 1:10], 2) < 1e-12
end
