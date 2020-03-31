""" Generate simulated Y given X and parameters. """

# @testset "Regression: Check that all combinations of the actual regression converge." begin
#     for data in ("ITP", "melanoma", "blood", "bone")
#         for loss in (FcgR.quadratic_loss, FcgR.proportion_loss)
#             fit = FcgR.fitRegression(FcgR.importDepletion(data), loss)
#             @test Optim.converged(fit)
#         end
#     end
# end

@testset "Regression: Test for similar results with row shuffling." begin
    for data in ("ITP", "melanoma", "blood", "bone")
        dataA = FcgR.importDepletion(data)
        dataB = dataA[shuffle(1:size(dataA, 1)), :]

        diffFit = FcgR.fitRegression(dataA, L0 = 1.0e-9, f = 4, murine = true)
        diffFit -= FcgR.fitRegression(dataB, L0 = 1.0e-9, f = 4, murine = true)

        @test norm(diffFit) < 0.001
    end
end
