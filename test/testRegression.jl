@testset "" begin

end



#=@testset "Test genRegPred() can take ForwardDiff." begin
    for data in ("ITP", "melanoma", "blood", "bone")
        df = FcRegression.importDepletion(data)
        Xo, Yo = old_regGenData(df; L0 = 1e-10, f = 6, murine = true)

        Xfc, Xdf, Yex = FcRegression.modelPred(df; L0 = 1e-10, f = 6, murine = true)
        extra = Xdf[!, in(["C1q", "Neutralization"]).(names(Xdf))]
        weights = rand(size(Xfc, 1) + size(extra, 2))
        Xmat, _ = FcRegression.regressionPred(Xfc, Xdf, weights, FcRegression.murineActI; showXmat = true)
        @test all(Matrix(Xmat) .≈ Xo)
        @test all(Yex .≈ Yo)

        f = x -> sum(FcRegression.regressionPred(Xfc, Xdf, weights, x))
        @test isa(ForwardDiff.gradient(f, [1.0, 1, -1, 1]), Vector{Float64})
        f = x -> sum(FcRegression.regressionPred(Xfc, Xdf, x, [1.0, 1, -1, 1]))
        @test isa(ForwardDiff.gradient(f, weights), Vector{Float64})
    end
end

@testset "Test regressionResult()" begin
    for data in ("ITP", "blood", "bone", "melanoma", "HIV", "Bcell")
        res, odf, effects, ActI_df = FcRegression.regressionResult(data; L0 = 1e-9, f = 6, murine = true)
        @test "Fitted" in names(odf)
        @test "LOOPredict" in names(odf)
        @test all(effects[!, "Q10"] .<= effects[!, "Median"])
        @test all(effects[!, "Q90"] .>= effects[!, "Median"])
    end
    for data in ("blood", "spleen", "bone")
        res, odf, effects, ActI_df = FcRegression.regressionResult(data; L0 = 1e-9, f = 6, murine = false)
        @test "Fitted" in names(odf)
        @test "LOOPredict" in names(odf)
        @test all(effects[!, "Q10"] .<= effects[!, "Median"])
        @test all(effects[!, "Q90"] .>= effects[!, "Median"])
    end
end=#
