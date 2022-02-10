@testset "Test each depletion data can run through regression successfully" begin

    for dataType in ["ITP", "blood", "bone", "melanoma", "HIV"]
        df = FcRegression.importDepletion(dataType)
        @test sum(Matrix(ismissing.(df))) == 0
        Xdf = FcRegression.modelPred(df; L0 = 1e-9, f = 4, murine = true)
        res = FcRegression.fitRegNNLS(Xdf; murine = true)
    end
    for dataType in ["blood", "spleen", "bone"]
        df = FcRegression.importHumanized(dataType)
        @test sum(Matrix(ismissing.(df))) == 0
        Xdf = FcRegression.modelPred(df; L0 = 1e-9, f = 4, murine = false)
        res = FcRegression.fitRegNNLS(Xdf; murine = false)
    end

end
