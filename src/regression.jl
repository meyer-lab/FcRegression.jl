using LsqFit

exponential(X, p) = -expm1.(-X * p)
gompertz(X::Array, p) = -expm1.(-p[1] .* expm1.(X * p[2:end]))

function regGenData(dataType; L0 = 1e-9, f = 4, KxStar = KxConst, Rtot = importRtot(), ActI = murineActI)
    df = importDepletion(dataType)

    if :Concentration in names(df)
        df[!, :Concentration] .*= L0
    else
        insertcols!(df, 3, :Concentration => L0)
    end

    resX = Matrix(undef, size(df, 1), size(Rtot, 2))
    for i = 1:size(df, 1)
        Kav = convert(Vector{Float64}, df[i, murineFcgR])
        Kav = reshape(Kav, 1, :)
        resX[i, :] = polyfc_ActV(df[i, :Concentration], KxStar, f, Rtot, [1.0], Kav, ActI)
    end

    resX[df[:, :Background] .== "NeuKO", cellTypes .== :Neu] .= 0.0
    resX[df[:, :Background] .== "ncMOKO", cellTypes .== :ncMO] .= 0.0

    @assert all(isfinite.(resX))
    return (resX, df[!, :Target])
end


function reg_wL0f(Xcond, ps, regMethod::Function, dataType)
    (X, Y) = regGenData(dataType; L0 = 10.0^ps[1], f = ps[2])
    return regMethod(X, ps[3:end])
end


function fitRegression(dataType, regMethod::Function; wL0f = false)
    (X, Y) = regGenData(dataType)
    if regMethod == exponential
        p_init = [ones(Float64, size(X, 2));]
        p_lower = [zeros(size(X, 2));]
        p_upper = [ones(Float64, size(X, 2)) .* 1e5;]
        autod = :forwarddiff
    elseif regMethod == gompertz
        p_init = [ones(Float64, size(X, 2) + 1);]
        p_lower = [zeros(size(X, 2) + 1);]
        p_upper = [100; ones(Float64, size(X, 2)) .* 1e5]
        autod = :finiteforward
    end

    # to fit L0 and f
    if wL0f
        fitMethod = (Xcond, ps) -> reg_wL0f(Xcond, ps, regMethod, dataType)
        p_init = vcat(-9, 4, p_init)
        p_lower = vcat(-16, 1, p_lower)
        p_upper = vcat(-7, 6, p_upper)
    else
        fitMethod = regMethod
    end

    fit = curve_fit(fitMethod, X, Y, p_init; lower = p_lower, upper = p_upper, autodiff = autod)
    if !fit.converged
        @warn "Fitting did not converge"
    end
    return fit
end
