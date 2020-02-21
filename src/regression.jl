import MLBase.LOOCV
import StatsBase.sample
import Statistics.std
import Distributions: cdf, Exponential

exponential(X::Matrix, p::Vector) = cdf.(Exponential(), X * p)

function regGenData(df; L0, f, KxStar = KxConst, murine = true, c1q = false)
    df = copy(df)
    Rtot = importRtot(murine = murine)
    ActI = murine ? murineActI : humanActI

    if :Concentration in names(df)
        df[!, :Concentration] .*= L0
    else
        insertcols!(df, 3, :Concentration => L0)
    end

    nct = size(Rtot, 2)
    resX = Matrix{Number}(undef, size(df, 1), (c1q ? nct+1 : nct))
    for i = 1:size(df, 1)
        Kav = convert(Vector{Float64}, df[i, murineFcgR])
        Kav = reshape(Kav, 1, :)
        resX[i, 1:nct] = polyfc_ActV(df[i, :Concentration], KxStar, f, Rtot, [1.0], Kav, ActI)
        if c1q
            resX[i, nct+1] = polyfc(df[i, :Concentration], KxStar, 1, [C1qConc], [1.0], hcat(df[i, :C1q]))["Lbound"]
        end
    end

    components = c1q ? [cellTypes; :C1q] : cellTypes
    resX[df[:, :Background] .== "NeuKO", components .== :Neu] .= 0.0
    resX[df[:, :Background] .== "ncMOKO", components .== :ncMO] .= 0.0
    Y = df[!, :Target]

    @assert all(isfinite.(resX))
    @assert all(isfinite.(Y))
    return (resX, Y)
end


function quadratic_loss(Y0::Vector, Y::Vector)
    return Distances.sqeuclidean(Y0, Y)
end


function proportion_loss(p::Vector, Y::Vector)
    res = sum((Y .- p) .^ 2 ./ (p .* (1 .- p) .+ 0.25))
    @assert all(isfinite.(res))
    return res
end


function loss_wL0f(df, ps::Vector{T}, lossFunction::Function; c1q = false)::T where {T <: Real}
    (X, Y) = regGenData(df; L0 = 10.0^ps[1], f = ps[2], c1q = c1q)
    Y0 = exponential(X, ps[3:end])
    return lossFunction(Y0, Y)
end


function fitRegression(df, lossFunction::Function; wL0f = false, c1q = false)
    ## this method only supports expoential distribution due to param choice

    (X, Y) = regGenData(df; L0 = 1.0e-9, f = 4, c1q = c1q)
    Np = size(X, 2)
    if wL0f
        fitMethod = (ps) -> loss_wL0f(df, ps, lossFunction; c1q = c1q)
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


function LOOCrossVal(dataType, lossFunction::Function; wL0f = false, c1q = false)
    df = importDepletion(dataType; c1q = c1q)
    n = size(df, 1)
    fitResults = Vector(undef, n)
    LOOindex = LOOCV(n)
    for (i, idx) in enumerate(LOOindex)
        fitResults[i] = fitRegression(df[idx, :], lossFunction; wL0f = wL0f, c1q = c1q)
    end
    return fitResults
end


function bootstrap(dataType, lossFunction::Function; nsample = 100, wL0f = false, c1q = false)
    df = importDepletion(dataType; c1q = c1q)
    n = size(df, 1)
    fitResults = Vector(undef, nsample)
    for i = 1:nsample
        for j = 1:5
            fit = try
                fitRegression(df[sample(1:n, n, replace = true), :], lossFunction; wL0f = wL0f, c1q = c1q)
            catch e
                @warn "This bootstrapping set failed at fitRegression"
                nothing
            end
            if fit != nothing
                fitResults[i] = fit
                break
            end
        end
    end
    @assert all(isa.(fitResults, Optim.MultivariateOptimizationResults))
    return fitResults
end


function CVResults(dataType, lossFunction::Function = proportion_loss; c1q = false)
    df = importDepletion(dataType; c1q = c1q)
    fit_out = fitRegression(df, lossFunction; c1q = c1q)
    loo_out = LOOCrossVal(dataType, lossFunction; c1q = c1q)
    btp_out = bootstrap(dataType, lossFunction; c1q = c1q)

    (X, Y) = regGenData(df; L0 = 1.0e-9, f = 4, c1q = c1q)
    fit_w = fit_out.:minimizer

    odf = df[!, [:Condition, :Background]]
    odf[!, :Y] = Y
    odf[!, :Fitted] = exponential(X, fit_w)
    odf[!, :LOOPredict] = vcat([exponential(X[[i], :], loo_out[i].:minimizer) for i = 1:length(loo_out)]...)

    effects = X .* fit_w'
    btp_ws = cat([X .* (a.:minimizer)' for a in btp_out]..., dims = (3))
    btp_std = dropdims(std(btp_ws; dims = 3), dims = 3)
    @assert size(effects) == size(btp_std)

    return (fit_w, odf, effects, btp_std)
end
