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

@testset "Test for similar results with row shuffling." begin
	dataA = FcgR.importDepletion("ITP")
	dataB = dataA[shuffle(1:size(dataA, 1)), :]

	fitA = FcgR.fitRegression(dataA, FcgR.proportion_loss)
	fitB = FcgR.fitRegression(dataB, FcgR.proportion_loss)

	@test Optim.converged(fitA)
	@test Optim.converged(fitB)
	@test Optim.minimum(fitA) ≈ Optim.minimum(fitB)
end
