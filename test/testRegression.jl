""" Generate simulated Y given X and parameters. """

@testset "Check that all combinations of the actual regression converge." begin
    for data in ("ITP", "melanoma", "blood", "bone")
        for method in (FcgR.exponential, FcgR.gompertz)
            for L0f in (false, true)
                fit = FcgR.fitRegression(data, method; wL0f = L0f)
                @test fit.converged
            end
        end
    end
end
