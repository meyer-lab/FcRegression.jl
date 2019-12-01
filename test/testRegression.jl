
function regSimY(X, regMethod::Function, p; randr=0.05)
    """ Generate simulated Y given X and parameters """
    N = size(X)[1]
    Y = regMethod(X, p) .+ randr * randn(N)
    if sum(Y.>=0.99) + sum(Y.<=0.01) > 0.2 * N
        @warn "Too many y's out of range (0 to 1)!"
    end
    Y[Y.>=0.99] .= 0.99
    Y[Y.<=0.01] .= 0.01
    return Y
end


@testset "Example of running regression." begin
	IgGCs = [1 1; 1 0; 0 1; 1 1; 1 0; 0 1; 1 1; 1 0; 0 1; 1 1; 1 0; 0 1]
	Rcpon = [1 1 1; 1 1 1; 1 1 1; 1 1 0; 1 1 0; 1 1 0; 1 0 1; 1 0 1; 1 0 1; 0 1 1; 0 1 1; 0 1 1]
	X = regGenX(IgGCs, Rcpon)

	Y_expo = regSimY(X, exponential, [80, 150, 230, 340])
	Y_gomp = regSimY(X, gompertz, [1.5, 80, 150, 130, 140])

	## fitting without L0 and f
	fit1 = curve_fit(FcgR.exponential, FcgR.X, FcgR.Y_expo, [zeros(4);], lower=[zeros(4);]; autodiff=:forwarddiff)
	fit2 = curve_fit(FcgR.gompertz, FcgR.X, FcgR.Y_gomp, [1; zeros(4);], lower=[zeros(5);]; autodiff=:forwarddiff)

    ## use LsqFit for parameters and L0, f
	fit3 = curve_fit(FcgR.reg_wL0f_expo, hcat(FcgR.IgGCs, FcgR.Rcpon), FcgR.Y_expo, [-10; 6; zeros(4);], lower=[-12; 1; zeros(4);]; autodiff=:finiteforward)
	fit4 = curve_fit(FcgR.reg_wL0f_gomp, hcat(FcgR.IgGCs, FcgR.Rcpon), FcgR.Y_gomp, [-10; 6; 1; zeros(4);], lower=[-12; 1; 0; zeros(4);]; autodiff=:finiteforward)
end
