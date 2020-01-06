""" Generate simulated Y given X and parameters. """
function regSimY(X, regMethod::Function, p)
    N = size(X)[1]
    Y = regMethod(X, p) .+ 0.01 * randn(N)

    # Clip to be within a feasible range
    Y[Y .>= 0.999] .= 0.999
    Y[Y .<= 0.001] .= 0.001
    return Y
end


@testset "Test regression parts with synthetic data." begin
    IgGCs = [1 1; 1 0; 0 1; 1 1; 1 0; 0 1; 1 1; 1 0; 0 1; 1 1; 1 0; 0 1]
    Rcpon = [1 1 1; 1 1 1; 1 1 1; 1 1 0; 1 1 0; 1 1 0; 1 0 1; 1 0 1; 1 0 1; 0 1 1; 0 1 1; 0 1 1]
    X = FcgR.regGenX(IgGCs, Rcpon)

    Y_expo = regSimY(X, FcgR.exponential, [80, 150, 230, 340])
    Y_gomp = regSimY(X, FcgR.gompertz, [1.5, 80, 150, 130, 140])

    ## fitting without L0 and f
    fit1 = curve_fit(FcgR.exponential, X, Y_expo, [zeros(4);], lower = [zeros(4);]; autodiff = :forwarddiff)
    @test fit1.converged
    fit2 = curve_fit(FcgR.gompertz, X, Y_gomp, [1; zeros(4)], lower = [zeros(5);]; autodiff = :forwarddiff)
    @test fit2.converged

    ## use LsqFit for parameters and L0, f
    fit3 = curve_fit(FcgR.reg_wL0f_expo, hcat(IgGCs, Rcpon), Y_expo, [-10; 6; zeros(4)], lower = [-12; 1; zeros(4)]; autodiff = :forwarddiff)
    @test fit3.converged
    fit4 = curve_fit(FcgR.reg_wL0f_gomp, hcat(IgGCs, Rcpon), Y_gomp, [-10; 6; 1; zeros(4)], lower = [-12; 1; 0; zeros(4)]; autodiff = :forwarddiff)
    @test fit4.converged
end


@testset "Check that all combinations of the actual regression converge." begin
    for data in ("ITP", "melanoma")
        for method in (FcgR.exponential, FcgR.gompertz)
            for L0f in (false, true)
                fit = FcgR.fitRegression(data, method; wL0f = L0f)
                @test fit.converged
            end
        end
    end
end
