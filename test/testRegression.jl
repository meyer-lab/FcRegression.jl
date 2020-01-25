""" Generate simulated Y given X and parameters. """

@testset "Check that all combinations of the actual regression converge." begin
    for data in ("ITP", "melanoma", "blood", "bone")
        for loss in (FcgR.quadratic_loss, FcgR.proportion_loss)
            fit = FcgR.fitRegression(FcgR.importDepletion(data), loss)
            @test Optim.converged(fit)
        end
    end
end

@testset "Test cross validation and bootstrapping." begin
    @test all(Optim.converged.(FcgR.LOOCrossVal("ITP", FcgR.proportion_loss)))
    @test all(Optim.converged.(FcgR.bootstrap("ITP", FcgR.proportion_loss)))
end
