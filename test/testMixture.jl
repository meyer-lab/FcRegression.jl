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
    @test norm(p ./ p[1] .- [res["ValConv"][i] for i = 1:5], 2) < 1e-12
    @test norm(q .* p[1] .- [res["ExpConv"]["x" * string(i)] for i = 1:10], 2) < 1e-12
end
