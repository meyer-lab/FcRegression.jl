import MLBase.LOOCV
import Statistics: mean, std
import Distributions: cdf, Exponential
import StatsBase: sample, mode
using NonNegLeastSquares

exponential(X::Matrix, p::Vector) = cdf.(Exponential(), X * p)
inv_exponential(y::Real) = -log(1 - y)

function regGenData(df; L0, f, murine::Bool, retdf = false)
    df = copy(df)
    FcRecep = murine ? murineFcgR : humanFcgR
    ActI = murine ? murineActI : humanActI

    if :Concentration in names(df)
        df[!, :Concentration] .*= L0
    else
        insertcols!(df, 3, :Concentration => L0)
    end

    X = DataFrame(repeat([Float64], length(cellTypes)), cellTypes)
    for i = 1:size(df, 1)
        Kav = convert(Vector{Float64}, df[i, FcRecep])
        Kav = reshape(Kav, 1, :)
        Rtot = importRtot(; murine = murine, genotype = murine ? "NA" : df[i, :Genotype])
        push!(X, polyfc_ActV(df[i, :Concentration], KxConst, f, Rtot, [1.0], Kav, ActI))
    end

    if :C1q in names(df)
        X[!, :C1q] = df[!, :C1q] .* df[!, :Concentration]
    end
    if :Neutralization in names(df)
        X[!, :Neutralization] = df[!, :Neutralization]
    end

    if :Background in names(df)
        X[df[:, :Background] .== "NeuKO", :Neu] .= 0.0
        X[df[:, :Background] .== "ncMOKO", :ncMO] .= 0.0
    end
    Y = df[!, :Target]

    @assert all(isfinite.(Matrix(X)))
    @assert all(isfinite.(Y))
    if retdf
        return (X, Y)
    else
        return (Matrix{Float64}(X), Y)
    end
end


function quadratic_loss(Y0::Vector, Y::Vector)
    return Distances.sqeuclidean(Y0, Y)
end

function proportion_loss(p::Vector, Y::Vector)
    res = sum((Y .- p) .^ 2 ./ (p .* (1 .- p) .+ 0.25))
    @assert all(isfinite.(res))
    return res
end

function loss_wL0f(df, ps::Vector{T}, lossFunc::Function; murine::Bool)::T where {T <: Real}
    (X, Y) = regGenData(df; L0 = 10.0^ps[1], f = ps[2], murine = murine)
    Y0 = exponential(X, ps[3:end])
    return lossFunc(Y0, Y)
end

function fitRegression(df, lossFunc::Function = proportion_loss; L0, f, murine::Bool, wL0f = false)
    ## this method only supports expoential distribution due to param choice

    (X, Y) = regGenData(df; L0 = L0, f = f, murine = murine)
    Np = size(X, 2)
    if wL0f
        fitMethod = (ps) -> loss_wL0f(df, ps, lossFunc)
    else
        fitMethod = (ps) -> lossFunc(exponential(X, ps), Y)
    end
    g! = (G, ps) -> ForwardDiff.gradient!(G, fitMethod, ps)

    p_init = vec(nonneg_lsq(X, inv_exponential.(Y)))
    p_lower = zeros(Float64, Np)
    p_upper = maximum(p_init) .* ones(Float64, Np)

    if wL0f
        p_init = vcat(-9, 4, p_init)
        p_lower = vcat(-16, 1, p_lower)
        p_upper = vcat(-7, 6, p_upper)
    end

    return p_init
end

function LOOCrossVal(df, lossFunc::Function; L0, f, murine, wL0f = false)
    n = size(df, 1)
    fitResults = Vector(undef, n)
    LOOindex = LOOCV(n)
    for (i, idx) in enumerate(LOOindex)
        fitResults[i] = fitRegression(df[idx, :], lossFunc; L0 = L0, f = f, murine = murine, wL0f = wL0f)
    end
    return fitResults
end

function bootstrap(df, lossFunc::Function; nsample = 100, L0, f, murine, wL0f = false)
    n = size(df, 1)
    fitResults = Vector(undef, nsample)
    for i = 1:nsample
        for j = 1:5
            fit = try
                fitRegression(df[sample(1:n, n, replace = true), :], lossFunc; L0 = L0, f = f, murine = murine, wL0f = wL0f)
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
    @assert all(isa.(fitResults, Array{<:Real}))
    return fitResults
end


function CVResults(df, lossFunc::Function = proportion_loss; L0, f, murine::Bool)
    fit_w = fitRegression(df, lossFunc; L0 = L0, f = f, murine = murine)
    loo_out = LOOCrossVal(df, lossFunc; L0 = L0, f = f, murine = murine)
    btp_out = bootstrap(df, lossFunc; L0 = L0, f = f, murine = murine)

    (X, Y) = regGenData(df; L0 = L0, f = f, murine = murine, retdf = true)
    components = names(X)
    @assert length(fit_w) == length(components)
    
    # Need labels for plotting the HIV case
    if :Neutralization in names(df)
        odf = df[!, [:Condition, :Background, :Label]]
    else
        odf = df[!, [:Condition, :Background]]
    end

    odf[!, :Concentration] .= (:Concentration in names(df)) ? (df[!, :Concentration] .* L0) : L0
    odf[!, :Y] = Y
    odf[!, :Fitted] = exponential(Matrix(X), fit_w)
    odf[!, :LOOPredict] = vcat([exponential(Matrix(X[[i], :]), loo_out[i]) for i = 1:length(loo_out)]...)

    selected = (odf[!, :Background] .== "wt") .& (odf[!, :Concentration] .== mode(odf[!, :Concentration]))
    X = Matrix(X[selected, :])

    effects = X .* fit_w'
    btp_ws = cat([X .* a' for a in btp_out]..., dims = (3))
    btp_std = dropdims(std(btp_ws; dims = 3), dims = 3)
    @assert size(effects) == size(btp_std)

    wdf = DataFrame(Weight = vec(effects), BtpStdev = vec(btp_std))
    wdf[!, :Condition] = vec(repeat(odf[selected, :Condition], 1, length(fit_w)))
    wdf[!, :Component] = vec(repeat(reshape(components, 1, :), size(effects, 1), 1))
    @assert size(wdf, 1) == size(X, 1) * length(fit_w)

    return fit_w, odf, wdf
end
