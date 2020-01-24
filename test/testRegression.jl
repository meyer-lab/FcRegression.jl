""" Generate simulated Y given X and parameters. """

@testset "Check that all combinations of the actual regression converge." begin
    for data in ("ITP", "melanoma", "blood", "bone")
        fit = FcgR.fitRegression(FcgR.importDepletion(data), FcgR.proportion_loss)
        @test Optim.converged(fit)
    end
end

@testset "Test cross validation and bootstrapping." begin
    for data in ("ITP", "melanoma", "blood", "bone")
        for loss in (FcgR.quadratic_loss, FcgR.proportion_loss)
            @test all(Optim.converged.(FcgR.LOOCrossVal(data, loss)))
            @test all(Optim.converged.(FcgR.bootstrap(data, loss)))
        end
    end
end
