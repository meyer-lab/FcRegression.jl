using LsqFit

exponential(t) = -expm1.(-t)
gompertz(t::Real, shape) = -expm1.( -shape .* expm1.(t) )
exponential(X, p) = -expm1.( -X * p )
gompertz(X::Array, p) = -expm1.( -p[1] .* expm1.(X * p[2:end]) )

function regGenData(dataType;
    L0 = 1e-9,
    f = 4,
    KxStar = KxConst,
    Rtot = importRtot(),
    ActI = murineActI)

    df = importDepletion(dataType)
    ndpt = size(df,1)
    if :Concentration in names(df)
        df[!, :Concentration] .*= L0
    else
        insertcols!(df, 3, :Concentration => L0)
    end

    resX = Matrix{Float64}(undef, ndpt, size(Rtot,2))
    for i in 1:ndpt
        row = df[i, :]
        Kav = convert(Vector{Float64}, row[murineFcgR])
        Kav = reshape(Kav, 1, :)
        subActV = polyfc_ActV(row[:Concentration], KxStar, f, Rtot, [1.], Kav, ActI)
        resX[i, :] = subActV
    end

    resX[ df[:, :Background].=="NeuKO", cellTypes .== :Neu ] .= 0.0
    resX[ df[:, :Background].=="ncMOKO", cellTypes .== :ncMO ] .= 0.0

    @assert all(isfinite.(resX))
    return (resX, df[!, :Target])
end


function reg_wL0f(Xcond, ps, regMethod::Function, dataType)
    (X, Y) = regGenData(dataType; L0 = 10.0^ps[1], f = ps[2])
    return regMethod(X, ps[3:end])
end

function fitRegression(dataType, regMethod::Function; wL0f=false)
    (X, Y) = regGenData(dataType)
    if regMethod == exponential
        p_init = [ones(Float64, size(X,2));]
        p_lower = [zeros(size(X,2));]
        p_upper = [ones(Float64, size(X,2)) .* 1e5;]
    elseif regMethod == gompertz
        p_init = [ones(Float64, size(X,2)+1);]
        p_lower = [zeros(size(X,2)+1);]
        p_upper = [100; ones(Float64, size(X,2)) .* 1e5;]
    end

    # to fit L0 and f
    if wL0f
        fitMethod = (Xcond, ps) -> reg_wL0f(Xcond, ps, regMethod=regMethod)
        p_init = vcat(-9, 4, p_init)
        p_lower = vcat(-16, 2, p_lower)
        p_upper = vcat(-6, 12, p_upper)
    else
        fitMethod = regMethod
    end

    fit = curve_fit(fitMethod, X, Y, p_init; lower=p_lower, upper=p_upper, autodiff=:forwarddiff)
    if !fit.converged
        @warn "Fitting did not converge"
    end
    return fit
end


""" Build the input data matrix based on the binding parameters. """
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

    Xtype = promote_type(eltype(Rpho), typeof(L0), typeof(KxStar), typeof(f))
    X = Array{Xtype,2}(undef, N, nct)

    for xi in 1:N
        Kav_n = Kav .* Rcpon[xi, :]'
        for ict in 1:nct
            Rtot = Rpho[ict, :]
            X[xi, ict] = polyfc(L0, KxStar, f, Rtot, IgGCs[xi, :], Kav_n, ActI)["ActV"]
        end
    end
    return X
end


function reg_wL0f(Xcond, ps; regMethod::Function)
    X = regGenX(Xcond[:, 1:2], Xcond[:, 3:5]; L0=10.0^ps[1], f=ps[2])
    @assert all(isfinite.(X))
    return regMethod(X, ps[3:end])
end

reg_wL0f_expo = (Xcond, ps) -> reg_wL0f(Xcond, ps; regMethod=exponential)
reg_wL0f_gomp = (Xcond, ps) -> reg_wL0f(Xcond, ps; regMethod=gompertz)
