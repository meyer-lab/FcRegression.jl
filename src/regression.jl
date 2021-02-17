import MLBase.LOOCV
import Statistics: mean, quantile, std
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

function modelPred(df; L0, f, murine::Bool = true)
    df = copy(df)
    FcRecep = murine ? murineFcgR : humanFcgR
    cellTypes = murine ? murineCellTypes : humanCellTypes

    if "Concentration" in names(df)
        df[!, "Concentration"] .*= L0
    else
        insertcols!(df, 3, "Concentration" => L0)
    end
    # Xfc = cellTypes * FcRecep * items
    Xfc = Array{Float64}(undef, length(cellTypes), length(FcRecep), size(df, 1))
    for k = 1:size(Xfc, 3)    # dataframe row
        Kav = convert(Vector{Float64}, df[k, FcRecep])
        Kav = reshape(Kav, 1, :)
        Rtot = copy(importRtot(; murine = murine, genotype = murine ? "NA" : df[k, :Genotype]))
        if "Background" in names(df)
            if df[k, "Background"] == "NeuKO"
                Rtot[:, "Neu" .== cellTypes] .= 0.0
            elseif df[k, "Background"] == "ncMOKO"
                Rtot[:, "ncMO" .== cellTypes] .= 0.0
            end
        end
        for i = 1:size(Xfc, 1)    # cell type
            Rtotc = Rtot[:, i]
            Xfc[i, :, k] = polyfc(df[k, :Concentration], KxConst, f, Rtotc, [1.0], Kav).Rmulti_n
        end
    end
    Xdf = df[!, in(["Condition", "Concentration", "Genotype", "C1q", "Neutralization"]).(names(df))]
    if "C1q" in names(df)
        Xdf[!, "C1q"] = Xdf[!, "C1q"] .* Xdf[!, "Concentration"]
    end

    Y = df[!, "Target"]
    return Xfc, Xdf, Y
end


function regressionPred(Xfc, Xdf::Union{DataFrame, Nothing}, cellWeights, recepActI; showXmat = false, murine = true)
    ansType = promote_type(eltype(Xfc), eltype(cellWeights), eltype(recepActI))
    noextra = true
    if Xdf == nothing
        extra = nothing
        noextra = true
    else
        extra = Xdf[!, in(["C1q", "Neutralization"]).(names(Xdf))]
        noextra = size(extra, 2) == 0
    end

    @assert length(cellWeights) == size(Xfc, 1) + (noextra ? 0 : size(extra, 2))
    @assert length(recepActI) == size(Xfc, 2)
    cellTypes = murine ? murineCellTypes : humanCellTypes
    if !noextra
        @assert size(Xfc, 3) == size(extra, 1)
        Xmat = DataFrame(repeat([ansType], length(cellWeights)), vcat(cellTypes, names(extra)))
    else
        Xmat = DataFrame(repeat([ansType], length(cellWeights)), cellTypes)
    end

    Ypred = Array{ansType}(undef, size(Xfc, 3))
    for k = 1:size(Xfc, 3)
        recepEff = Xfc[:, :, k] * recepActI
        recepEff[recepEff .< 0.0] .= 0.0
        if !noextra
            append!(recepEff, extra[k, :])
        end
        push!(Xmat, recepEff)
        Ypred[k] = sum(cellWeights .* recepEff)      # exponential?
    end
    if showXmat
        return Xmat, Ypred
    else
        return Ypred
    end
end

regressionPred(Xfc, Xdf, fit::optResult; showXmat = false, murine = true) =
    regressionPred(Xfc, Xdf, fit.cellWs, fit.ActI; showXmat = showXmat, murine = murine)

function nnls_fit(Xfc, extra, Y::AbstractVector, ActI::AbstractVector)
    cY = inv_exponential.(Y)
    ansType = promote_type(eltype(Xfc), eltype(Y), eltype(ActI))
    Xmat = Matrix{ansType}(undef, size(Xfc, 3), size(Xfc, 1))
    for i = 1:size(Xfc, 1)
        for j = 1:size(Xfc, 3)
            Xmat[j, i] = Xfc[i, :, j]' * ActI
        end
    end
    Xmat[Xmat .< 0.0] .= 0.0
    if !any(size(extra) .<= 0)
        extraM = Matrix(extra)
        @assert size(extraM, 1) == size(Xfc, 3)
        Xmat = hcat(Xmat, extraM)
    end

    cY = convert(Vector{ansType}, cY)
    w = vec(nonneg_lsq(Xmat, cY; alg=:nnls))
    Yr = Xmat * w
    residual = norm(cY - Yr, 2) / length(Y)

    return w, residual
end

function fitRegression_woActI(Xfc, Xdf, Y, ActI)
    extra = Xdf[!, in(["C1q", "Neutralization"]).(names(Xdf))]
    cellWs, residual = nnls_fit(Xfc, extra, Y, ActI)
    return optResult(cellWs, ActI, residual)
end

function fitRegression(Xfc, Xdf, Y; murine::Bool = true)
    extra = Xdf[!, in(["C1q", "Neutralization"]).(names(Xdf))]
    cellWlen = size(Xfc, 1) + size(extra, 2)
    ActIL = length(murine ? murineActI : humanActI)
    init = ones(ActIL)

    func = x -> nnls_fit(Xfc, extra, Y, x)[2]
    opt = optimize(func, init; autodiff = :forward, method = LBFGS(manifold = Optim.Sphere()))
    ActI = opt.minimizer

    return fitRegression_woActI(Xfc, Xdf, Y, ActI)
end


function LOOCrossVal(Xfc, Xdf, Y; murine, ActI::Union{Nothing, Vector} = nothing)
    n = size(Xdf, 1)
    fitResults = Vector{optResult}(undef, n)
    LOOindex = LOOCV(n)
    if ActI === nothing
        for (i, idx) in enumerate(LOOindex)
            fitResults[i] = fitRegression(Xfc[:, :, idx], Xdf[idx, :], Y[idx]; murine = murine)
        end
    else
        for (i, idx) in enumerate(LOOindex)
            fitResults[i] = fitRegression_woActI(Xfc[:, :, idx], Xdf[idx, :], Y[idx], ActI)
        end
    end
    return fitResults
end


function regressionResult(dataType; L0, f, murine::Bool)
    df = murine ? importDepletion(dataType) : importHumanized(dataType)

    Xfc, Xdf, Y = modelPred(df; L0 = L0, f = f, murine = murine)
    res = fitRegression(Xfc, Xdf, Y; murine = murine)
    loo_res = LOOCrossVal(Xfc, Xdf, Y; murine = murine, ActI = res.ActI)

    odf = df[!, in(["Condition", "Background"]).(names(df))]
    odf[!, "Concentration"] .= ("Concentration" in names(df)) ? (df[!, "Concentration"] .* L0) : L0
    odf[!, "Y"] = Y
    odf[!, "Fitted"] = exponential(regressionPred(Xfc, Xdf, res; murine = murine))
    odf[!, "LOOPredict"] = exponential([regressionPred(Xfc[:, :, ii], Xdf[[ii], :], loo_res[ii]; murine = murine)[1] for ii = 1:length(loo_res)])
    if "Label" in names(df)
        odf[!, "Label"] .= df[!, "Label"]
    end
    if "Genotype" in names(df)
        odf[!, "Genotype"] .= df[!, "Genotype"]
    end

    # Prepare for cell type weights in wildtype
    wildtype = copy(importKav(; murine = murine, c1q = ("C1q" in names(df)), IgG2bFucose = (:IgG2bFucose in df.Condition), retdf = true))
    wildtype[!, "Background"] .= "wt"
    wildtype[!, "Target"] .= 0.0
    if !murine
        wildtype[!, "Genotype"] .= "ZZZ"
    end
    if "Neutralization" in names(Xdf)
        wildtype[!, "Neutralization"] .= 0.0
    end
    rename!(wildtype, "IgG" => "Condition")

    wtXfc, wtXdf, wtY = modelPred(wildtype; L0 = L0, f = f, murine = murine)
    Xmat, _ = regressionPred(wtXfc, wtXdf, res.cellWs, res.ActI; showXmat = true, murine = murine)
    Xmat = Xmat .* res.cellWs'
    Xmat[!, "Condition"] = wtXdf[!, "Condition"]
    Cell_df = stack(Xmat, Not("Condition"))
    rename!(Cell_df, "value" => "Weight")
    rename!(Cell_df, "variable" => "Component")

    # ActI interval
    ActI_df = DataFrame(Receptor = (murine ? murineFcgR : humanFcgR), Activity = res.ActI)

    return res, odf, Cell_df, ActI_df
end
