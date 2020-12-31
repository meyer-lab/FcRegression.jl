import MLBase.LOOCV
import Statistics: mean, quantile
import Distributions: cdf, Exponential
import StatsBase: sample, mode
using NonNegLeastSquares

exponential(x::Real) = cdf(Exponential(), x)
exponential(X::Array) = cdf.(Exponential(), X)
exponential(X::Matrix, p::Vector) = cdf.(Exponential(), X * p)
inv_exponential(y::Real) = -log(1 - y)

mutable struct optResult{T}
    cellWs::Array{T}
    ActI::Array{T}
    residual::T
end

function modelPred(df; L0, f)
    """
    Return model predictions on depletion data
    Xfc = 3D array, cellTypes * FcRecep * items
    Y = Vector, items
    """
    df = copy(df)

    # Xfc = cellTypes * FcRecep * items
    Xfc = Array{Float64}(undef, length(cellTypes), length(murineFcgR), size(df, 1))
    for k = 1:size(Xfc, 3)    # dataframe row
        Kav = convert(Vector{Float64}, df[k, murineFcgR])
        Kav = reshape(Kav, 1, :)
        Rtot = copy(importRtot(; murine = true))
        if "Background" in names(df)
            if df[k, "Background"] == "NeuKO"
                Rtot[:, "Neu" .== cellTypes] .= 0.0
            elseif df[k, "Background"] == "ncMOKO"
                Rtot[:, "ncMO" .== cellTypes] .= 0.0
            end
        end
        for i = 1:size(Xfc, 1)    # cell type
            Xfc[i, :, k] = polyfc(L0, KxConst, f, Rtot[:, i], [1.0], Kav).Rmulti_n
        end
    end

    Y = df[!, "Target"]
    return Xfc, Y
end

function regressionPred(Xfc, cellWeights, recepActI; showXmat = false)
    """ Derive Ypred from all Xs, not exponentially transformed. Used only when having weights ready """
    ansType = promote_type(eltype(Xfc), eltype(cellWeights), eltype(recepActI))

    @assert length(recepActI) == size(Xfc, 2)
    Xmat = DataFrame(repeat([ansType], length(cellWeights)), cellTypes)
    Ypred = Array{ansType}(undef, size(Xfc, 3))
    for k = 1:size(Xfc, 3)
        recepEff = Xfc[:, :, k] * recepActI
        recepEff[recepEff .< 0.0] .= 0.0
        push!(Xmat, recepEff)
        Ypred[k] = sum(cellWeights .* recepEff)      # exponential?
    end
    if showXmat
        return Xmat, Ypred
    else
        return Ypred
    end
end

regressionPred(Xfc, fit::optResult; showXmat = false) = regressionPred(Xfc, fit.cellWs, fit.ActI; showXmat = showXmat)

function nnls_fit(Xfc, Y, ActI)
    cY = inv_exponential.(Y)
    ansType = promote_type(eltype(Xfc), eltype(Y), eltype(ActI))
    Xmat = Matrix{ansType}(undef, size(Xfc, 3), size(Xfc, 1))
    for i = 1:size(Xfc, 1)
        for j = 1:size(Xfc, 3)
            Xmat[j, i] = Xfc[i, :, j]' * ActI
        end
    end
    Xmat[Xmat .< 0.0] .= 0.0
    w = vec(nonneg_lsq(Xmat, cY))
    Yr = Xmat * w
    residual = norm(cY - Yr, 2) / length(Y)

    return w, residual
end

function fitRegression(Xfc, Y; ActI::Union{Nothing, Vector} = nothing)
    upper = ones(length(murineActI)) .* 4.0
    lower = ones(length(murineActI)) .* -4.0
    init = Float64.(murineActI)

    if ActI == nothing
        func = x -> nnls_fit(Xfc, Y, x)[2]
        opt = optimize(func, lower, upper, init)
        ActI = opt.minimizer
    end

    cellWs, residual = nnls_fit(Xfc, Y, ActI)
    return optResult(cellWs, ActI, residual)
end


function LOOCrossVal(Xfc, Y; ActI::Union{Nothing, Vector} = nothing)
    n = size(Xfc, 3)
    fitResults = Vector{optResult}(undef, n)
    LOOindex = LOOCV(n)
    for (i, idx) in enumerate(LOOindex)
        fitResults[i] = fitRegression(Xfc[:, :, idx], Y[idx], ActI = ActI)
    end
    return fitResults
end


function regressionResult(dataType; L0, f)
    df = importDepletion(dataType)

    Xfc, Y = modelPred(df; L0 = L0, f = f)
    res = fitRegression(Xfc, Y)
    loo_res = LOOCrossVal(Xfc, Y; ActI = res.ActI)

    odf = df[!, in(["Condition", "Background"]).(names(df))]
    odf[!, "Concentration"] .= ("Concentration" in names(df)) ? (df[!, "Concentration"] .* L0) : L0
    odf[!, "Y"] = Y
    odf[!, "Fitted"] = exponential(regressionPred(Xfc, res))
    odf[!, "LOOPredict"] = exponential([regressionPred(Xfc[:, :, ii], loo_res[ii])[1] for ii = 1:length(loo_res)])

    # Prepare for cell type weights in wildtype
    wildtype = copy(importKav(; murine = true, IgG2bFucose = (:IgG2bFucose in df.Condition), retdf = true))
    wildtype[!, "Background"] .= "wt"
    wildtype[!, "Target"] .= 0.0
    rename!(wildtype, "IgG" => "Condition")

    wtXfc, _ = modelPred(wildtype; L0 = L0, f = f)
    Xmat, _ = regressionPred(wtXfc, res; showXmat = true)
    Xmat = Xmat .* res.cellWs'
    Xmat[!, "Condition"] = wildtype[!, "Condition"]
    Cell_df = stack(Xmat, Not("Condition"))
    rename!(Cell_df, "value" => "Weight")
    rename!(Cell_df, "variable" => "Component")

    # ActI interval
    ActI_df = DataFrame(Receptor = murineFcgR, Activity = res.ActI)

    return res, odf, Cell_df, ActI_df
end 
