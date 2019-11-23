using LsqFit
using fcBindingModel

exponential(t) = -expm1.(-t)
gompertz(t::Real, shape) = -expm1.( -shape .* expm1.(t) )
exponential(X, p) = -expm1.( -X * p )
gompertz(X::Array, p) = -expm1.( -p[1] .* expm1.(X * p[2:end]) )

function regGenX(IgGCs, Rcpon;
    L0 = 1e-9,
    f = 4,
    KxStar = 10^-12.2,
    Rpho = [100 300 600; 400 300 0; 400 500 100; 1000 10 900],
    Kav = rand(2,3) .* 1e6 .* [1, 1.5, 0.8]',
    ActI = [ -1, 1, 1])
    """
    Input:
    IgGCs is a N * ni size matrix. Each row is a vec of IgG composition
    Rcpon is a N * nr size matrix. Each row is a binary vec for existing receps
    Default kw_input:
    Rpho is a nct * nr size matrix, the abundance matrix
    Output:
    A N * nct size matrix. Each row is a condition, and entries are ActV
    """

    N = size(IgGCs)[1]
    @assert N == size(Rcpon)[1]
    @assert size(Kav) == (size(IgGCs)[2], size(Rcpon)[2])
    nct = size(Rpho)[1]

    X = Array{Float64,2}(undef, N, nct)
    for xi in 1:N
        Kav_n = Kav .* Rcpon[xi, :]'
        IgGC = IgGCs[xi, :]
        vals = []
        for ict in 1:nct
            Rtot = Rpho[ict, :]
            X[xi, ict] = fcBindingModel.polyfc(L0, KxStar, f, Rtot, IgGC, Kav_n, ActI)["ActV"]
        end
    end
    return X
end

IgGCs = [1 1; 1 0; 0 1; 1 1; 1 0; 0 1; 1 1; 1 0; 0 1; 1 1; 1 0; 0 1]
Rcpon = [1 1 1; 1 1 1; 1 1 1; 1 1 0; 1 1 0; 1 1 0; 1 0 1; 1 0 1; 1 0 1; 0 1 1; 0 1 1; 0 1 1]
X = regGenX(IgGCs, Rcpon)

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

Y_expo = regSimY(X, exponential, [80, 150, 230, 340])
Y_gomp = regSimY(X, gompertz, [1.5, 80, 150, 130, 140])

## fitting without L0 and f
fit1 = curve_fit(exponential, X, Y_expo, [zeros(4);], lower=[zeros(4);]; autodiff=:forwarddiff)
fit2 = curve_fit(gompertz, X, Y_gomp, [1; zeros(4);], lower=[zeros(5);]; autodiff=:forwarddiff)

function reg_wL0f(Xcond, ps; regMethod::Function)
    L0 = 10^ps[1]
    f = ps[2]
    p = ps[3:end]
    X = regGenX(Xcond[:, 1:2], Xcond[:, 3:5]; L0 = L0, f = f)
    return regMethod(X, p)
end

reg_wL0f_expo = (Xcond, ps) -> reg_wL0f(Xcond, ps; regMethod = exponential)
reg_wL0f_gomp = (Xcond, ps) -> reg_wL0f(Xcond, ps; regMethod = gompertz)

## use LsqFit for parameters and L0, f
fit3 = curve_fit(reg_wL0f_expo, hcat(IgGCs, Rcpon), Y_expo, [-10; 6; zeros(4);], lower=[-12; 1; zeros(4);]; autodiff=:finiteforward)
fit4 = curve_fit(reg_wL0f_gomp, hcat(IgGCs, Rcpon), Y_gomp, [-10; 6; 1; zeros(4);], lower=[-12; 1; 0; zeros(4);]; autodiff=:finiteforward)
