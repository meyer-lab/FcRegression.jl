
exponential(X, p) = -expm1.(-X * p)
gompertz(X::Array, p) = -expm1.(-p[1] .* expm1.(X * p[2:end]))

function regGenData(dataType; L0, f, KxStar = KxConst, Rtot = importRtot(), ActI = murineActI)
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
    Y = df[!, :Target]

    @assert all(isfinite.(resX))
    @assert all(isfinite.(Y))
    return (resX, Y)
end


function reg_wL0f(ps, regMethod::Function, dataType)
    (X, Y) = regGenData(dataType; L0 = 10.0^ps[1], f = ps[2])
    return norm(regMethod(X, ps[3:end]) - Y)
end


function fitRegression(dataType, regMethod::Function)
    Np = size(X, 2)

    p_init = 0.1*ones(Float64, Np)
    p_lower = zeros(Float64, Np)
    p_upper = ones(Float64, Np)

    if regMethod == gompertz
        p_init = [1.0; p_init]
        p_lower = [0.1; p_lower]
        p_upper = [10; p_upper]
    end

    fitMethod = (ps) -> reg_wL0f(ps, regMethod, dataType)
    g! = (G, ps) -> ForwardDiff.gradient!(G, fitMethod, ps)
    p_init = vcat(-9, 4, p_init)
    p_lower = vcat(-16, 1, p_lower)
    p_upper = vcat(-7, 6, p_upper)

    fit = optimize(fitMethod, g!, p_lower, p_upper, p_init, Fminbox(GradientDescent()))
    if !fit.converged
        @warn "Fitting did not converge"
    end
    return fit
end
