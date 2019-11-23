using Test
using Profile
using LsqFit
using FcgR

@testset "Example of running regression." begin
	## fitting without L0 and f
	fit1 = curve_fit(FcgR.exponential, FcgR.X, FcgR.Y_expo, [zeros(4);], lower=[zeros(4);]; autodiff=:forwarddiff)
	fit2 = curve_fit(FcgR.gompertz, FcgR.X, FcgR.Y_gomp, [1; zeros(4);], lower=[zeros(5);]; autodiff=:forwarddiff)

    ## use LsqFit for parameters and L0, f
	fit3 = curve_fit(FcgR.reg_wL0f_expo, hcat(FcgR.IgGCs, FcgR.Rcpon), FcgR.Y_expo, [-10; 6; zeros(4);], lower=[-12; 1; zeros(4);]; autodiff=:finiteforward)
	fit4 = curve_fit(FcgR.reg_wL0f_gomp, hcat(FcgR.IgGCs, FcgR.Rcpon), FcgR.Y_gomp, [-10; 6; 1; zeros(4);], lower=[-12; 1; 0; zeros(4);]; autodiff=:finiteforward)
end

include("testSynergy.jl")
