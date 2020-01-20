""" Generate simulated Y given X and parameters. """

@testset "Check that all combinations of the actual regression converge." begin
    for data in ("ITP", "melanoma", "blood", "bone")
        for method in (FcgR.exponential, FcgR.gompertz)
            fit = FcgR.fitRegression(data, method)
            @test Optim.converged(fit)
        end
    end
end
