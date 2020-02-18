import MLBase.LOOCV
import StatsBase.sample
import Statistics.std
import Distributions: cdf, Exponential

exponential(X::Matrix, p::Vector) = cdf.(Exponential(), X * p)

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

    resX = Matrix{Number}(undef, size(df, 1), size(Rtot, 2))
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


function quadratic_loss(Y0::Vector, Y::Vector)
    return Distances.sqeuclidean(Y0, Y)
end


function proportion_loss(p::Vector, Y::Vector)
    res = sum( (Y .- p) .^ 2 ./ (p .* (1 .- p) .+ 0.25) )
    @assert all(isfinite.(res))
    return res
end


function loss_wL0f(df, ps::Vector{T}, lossFunction::Function)::T where {T <: Real}
    (X, Y) = regGenData(df; L0 = 10.0^ps[1], f = ps[2])
    Y0 = exponential(X, ps[3:end])
    return lossFunction(Y0, Y)
end


function fitRegression(df, lossFunction::Function; wL0f = false)
    ## this method only supports expoential distribution due to param choice

    (X, Y) = regGenData(df; L0 = 1.0e-9, f = 4)
    Np = size(X, 2)
    if wL0f
        fitMethod = (ps) -> loss_wL0f(df, ps, lossFunction)
    else
        fitMethod = (ps) -> lossFunction(exponential(X, ps), Y)
    end
    g! = (G, ps) -> ForwardDiff.gradient!(G, fitMethod, ps)

    p_init = 1.0 * ones(Float64, Np)
    p_lower = zeros(Float64, Np)
    p_upper = 10.0 * ones(Float64, Np)

    if wL0f
        p_init = vcat(-9, 4, p_init)
        p_lower = vcat(-16, 1, p_lower)
        p_upper = vcat(-7, 6, p_upper)
    end

    fit = optimize(fitMethod, g!, p_lower, p_upper, p_init, Fminbox())
    if !Optim.converged(fit)
        @warn "Fitting did not converge"
    end
    return fit
end


function LOOCrossVal(dataType, lossFunction::Function; wL0f = false)
    df = importDepletion(dataType)
    n = size(df, 1)
    fitResults = Vector(undef, n)
    LOOindex = LOOCV(n)
    for (i, idx) in enumerate(LOOindex)
        fitResults[i] = fitRegression(df[idx, :], lossFunction, wL0f = wL0f)
    end
    return fitResults
end


function bootstrap(dataType, lossFunction::Function; nsample = 100, wL0f = false)
    df = importDepletion(dataType)
    n = size(df, 1)
    fitResults = Vector(undef, nsample)
    for i = 1:nsample
        fitResults[i] = fitRegression(df[sample(1:n, n, replace = true), :], lossFunction, wL0f = wL0f)
    end
    return fitResults
end


function CVResults(dataType, lossFunction::Function = proportion_loss)
    df = importDepletion(dataType)
    fit_out = fitRegression(df, lossFunction)
    loo_out = LOOCrossVal(dataType, lossFunction)
    btp_out = bootstrap(dataType, lossFunction)

    (X, Y) = regGenData(df; L0 = 1.0e-9, f = 4)
    fit_w = fit_out.:minimizer

    odf = df[!, [:Condition, :Background]]
    odf[!, :Y] = Y
    odf[!, :Fitted] = exponential(X, fit_w)
    odf[!, :LOOPredict] = vcat([exponential(X[[i], :], loo_out[i].:minimizer) for i in 1:length(loo_out)]...)

    effects = X .* fit_w'
    btp_ws = cat([X .* (a.:minimizer)' for a in btp_out]..., dims=(3))
    btp_std = std(btp_ws; dims=3)

    return (odf, effects, btp_std)
end
