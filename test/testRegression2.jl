@testset "Test genRegPred() can take ForwardDiff." begin
    for data in ("ITP", "melanoma", "blood", "bone")
        df = FcRegression.importDepletion(data)
        X, _ = FcRegression.genModelPred(df; L0=1e-9, f=4, murine=true)
        f = x -> sum(FcRegression.genRegPred(X, [1.0, 2, 3, 4, 5], x))
        @test isa(ForwardDiff.gradient(f, [1.0, 1, -1, 1]), Vector{Float64})
        f = x -> sum(FcRegression.genRegPred(X, x, [1.0, 1, -1, 1]))
        @test isa(ForwardDiff.gradient(f, [1.0, 2, 3, 4, 5]), Vector{Float64})
    end
end
