
exponential(X::Matrix, p::Vector) = Distributions.cdf.(Distributions.Exponential(), X * p)
weibull(X::Matrix, p::Vector) = Distributions.cdf.(Distributions.Weibull(p[1]), X * p[2:end])

function regGenData(df; L0, f, KxStar = KxConst, murine = true)
    df = copy(df)
    Rtot = importRtot(murine = murine)

    if murine
        ActI = murineActI
    else
        ActI = humanActI
    end

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


function reg_wL0f(ps::Vector{T}, regMethod::Function, df)::T  where {T<:Real}
    (X, Y) = regGenData(df; L0 = 10.0^ps[1], f = ps[2])
    return Distances.sqeuclidean(regMethod(X, ps[3:end]), Y)
end


function fitRegression(dataType, regMethod::Function)
    df = importDepletion(dataType)
    (X, Y) = regGenData(df; L0 = 1.0e-9, f = 4)
    Np = size(X, 2)

    fitMethod = (ps) -> reg_wL0f(ps, regMethod, df)
    g! = (G, ps) -> ForwardDiff.gradient!(G, fitMethod, ps)

    p_init = 0.1*ones(Float64, Np)
    p_lower = zeros(Float64, Np)
    p_upper = ones(Float64, Np)

    if regMethod == weibull
        p_init = [1.0; p_init]
        p_lower = [0.1; p_lower]
        p_upper = [10; p_upper]
    end
    
    p_init = vcat(-9, 4, p_init)
    p_lower = vcat(-16, 1, p_lower)
    p_upper = vcat(-7, 6, p_upper)

    fit = optimize(fitMethod, g!, p_lower, p_upper, p_init, Fminbox())
    if !Optim.converged(fit)
        @warn "Fitting did not converge"
    end
    return fit
end
