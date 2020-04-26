import MLBase.LOOCV
import Statistics: mean, quantile
import Distributions: cdf, Exponential
import StatsBase: sample, mode
using NonNegLeastSquares

mutable struct fitResult{T}
    x::Array{T}     # the param for best fit
    intercept::T
    r::T            # best fit residue
end
exponential(X::Matrix, p::Vector) = cdf.(Exponential(), X * p)
exponential(X::Matrix, p::fitResult) = cdf.(Exponential(), X * p.x .+ p.intercept)
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

function fitRegression(df, intercept = false, preset_W::Union{Vector,Nothing} = nothing;  L0, f, murine::Bool)
    (X, Y) = regGenData(df; L0 = L0, f = f, murine = murine)
    Xo = copy(X)
    if preset_W != nothing
        @assert size(X, 2) == length(preset_W)
        X = X * reshape(preset_W, :, 1)
    end
    if intercept
        X = hcat(ones(size(X, 1), 1), X)
    end

    cY = inv_exponential.(Y)
    W = vec(nonneg_lsq(X, cY))
    res = fitResult{promote_type(typeof(L0), typeof(f))}([NaN], 0, NaN)
    if intercept
        res.intercept = W[1]
        W = W[2:end]
    end
    res.x = preset_W != nothing ? W[1] .* preset_W : W
    res.r = norm(Xo * res.x .+ res.intercept - cY, 2) / length(Y)
    return res
end

function LOOCrossVal(df, intercept, preset_W; L0, f, murine)
    n = size(df, 1)
    fitResults = Vector(undef, n)
    LOOindex = LOOCV(n)
    for (i, idx) in enumerate(LOOindex)
        fitResults[i] = fitRegression(df[idx, :], intercept, preset_W; L0 = L0, f = f, murine = murine)
    end
    return fitResults
end

function bootstrap(df, intercept, preset_W; nsample = 100, L0, f, murine)
    n = size(df, 1)
    fitResults = Vector(undef, nsample)
    for i = 1:nsample
        for j = 1:5
            fit = try
                fitRegression(df[sample(1:n, n, replace = true), :], intercept, preset_W; L0 = L0, f = f, murine = murine).x
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


function CVResults(df, intercept = false, preset_W = nothing; L0, f, murine::Bool)
    fit_res = fitRegression(df, intercept, preset_W; L0 = L0, f = f, murine = murine)
    loo_out = LOOCrossVal(df, intercept, preset_W; L0 = L0, f = f, murine = murine)
    btp_out = bootstrap(df, intercept, preset_W; L0 = L0, f = f, murine = murine)

    (X, Y) = regGenData(df; L0 = L0, f = f, murine = murine, retdf = true)
    @assert length(fit_res.x) == length(names(X))

    odf = df[!, in([:Condition, :Background, :Genotype, :Label]).(names(df))]
    odf[!, :Concentration] .= (:Concentration in names(df)) ? (df[!, :Concentration] .* L0) : L0
    odf[!, :Y] = Y
    odf[!, :Fitted] = exponential(Matrix(X), fit_res)
    odf[!, :LOOPredict] = vcat([exponential(Matrix(X[[i], :]), loo_out[i]) for i = 1:length(loo_out)]...)

    wildtype = copy(importKav(; murine = murine, c1q = (:C1q in names(df)), IgG2bFucose = (:IgG2bFucose in df.Condition), retdf = true))
    wildtype[!, :Background] .= "wt"
    wildtype[!, :Target] .= 0.0
    if !murine
        wildtype[!, :Genotype] .= "ZZZ"
    end
    rename!(wildtype, :IgG => :Condition)
    wtX, _ = regGenData(wildtype; L0 = L0, f = f, murine = murine, retdf = true)

    comp = in(names(wtX)).(names(X))
    fit_res.x = fit_res.x[comp]    # remove neutralization from HIV
    effects = wtX .* fit_res.x'
    effects[!, :Condition] .= wildtype[!, :Condition]
    btp_ws = cat([Matrix(wtX) .* a[comp]' for a in btp_out]..., dims = (3))
    btp_qtl = mapslices(x -> quantile(x, [0.1, 0.5, 0.9]), btp_ws, dims = [3])

    wdf = stack(effects, Not(:Condition))
    rename!(wdf, :variable => :Component)
    rename!(wdf, :value => :Weight)
    wdf[!, :Q10] .= vec(btp_qtl[:, :, 1])
    wdf[!, :Median] .= vec(btp_qtl[:, :, 2])
    wdf[!, :Q90] .= vec(btp_qtl[:, :, 3])

    return fit_res, odf, wdf
end
