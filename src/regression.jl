import MLBase.LOOCV
import Statistics: mean, quantile, std
import Distributions: cdf, Exponential
import StatsBase: sample, mode
using NonNegLeastSquares
import Base: tanh

exponential(x::Real) = cdf(Exponential(), x)
exponential(X::Array) = cdf.(Exponential(), X)
exponential(X::Matrix, p::Vector) = cdf.(Exponential(), X * p)
inv_exponential(y::Real) = -log(1 - y)

tanh(X::Array) = tanh.(X)
tanh(X::Matrix, p::Vector) = tanh.(X * p)

mutable struct regResult{T}
    cellWs::Array{T}
    R2::T
end


""" plug in reg df, and output binding model results """
function modelPred(df; L0, f, murine::Bool = true, cellTypes = nothing)
    df = copy(df)
    FcRecep = murine ? murineFcgR : humanFcgR
    ActI = murine ? murineActI : humanActI
    if cellTypes == nothing
        cellTypes = murine ? murineCellTypes : humanCellTypes
    end

    if "Concentration" in names(df)
        df[!, "Concentration"] .*= L0
    else
        insertcols!(df, 3, "Concentration" => L0)
    end
    # Xfc = items * cellTypes
    Xfc = Array{Float64}(undef, size(df, 1), length(cellTypes))
    for k = 1:size(Xfc, 1)    # dataframe row
        Kav = Vector{Float64}(df[k, FcRecep])
        Kav = reshape(Kav, 1, :)
        Rtot = copy(importRtot(; murine = murine, genotype = murine ? "NA" : df[k, :Genotype], cellTypes = cellTypes))
        if "Background" in names(df)
            if df[k, "Background"] == "NeuKO"
                Rtot[:, "Neu" .== cellTypes] .= 0.0
            elseif df[k, "Background"] == "ncMOKO"
                Rtot[:, "ncMO" .== cellTypes] .= 0.0
            elseif df[k, "Background"] == "cMOKO"
                Rtot[:, "cMO" .== cellTypes] .= 0.0
            elseif df[k, "Background"] == "EOKO"
                Rtot[:, "EO" .== cellTypes] .= 0.0
            end
        end
        for i = 1:size(Xfc, 2)    # cell type
            Rtotc = Rtot[:, i]
            if dot(Rtotc, Kav) == 0.0
                Xfc[k, i] = 0.0
            else
                pred = polyfc(df[k, :Concentration], KxConst, f, Rtotc, [1.0], Kav).Rmulti_n
                Xfc[k, i] = maximum([dot(ActI, pred), 0.0])
            end
        end
    end
    Xdf = df[!, in(["Condition", "Background", "Concentration", "Genotype", "C1q", "Neutralization", "Target"]).(names(df))]
    Xdf = hcat(Xdf, DataFrame(Xfc, cellTypes))
    if "C1q" in names(df)
        Xdf[!, "C1q"] = Xdf[!, "C1q"] .* Xdf[!, "Concentration"]
    end
    return Xdf
end


function regPred(Xfc::Matrix, cellWeights::Vector; exp_method = true)
    Yr = Xfc * cellWeights
    return exp_method ? exponential(Yr) : tanh(Yr)
end

regPred(Xfc::Matrix, opt::regResult; exp_method) = regPred(Xfc, opt.cellWs; exp_method = exp_method)

function regPred(Xdf::DataFrame, opt::regResult; murine, exp_method, cellTypes = nothing)
    if cellTypes === nothing
        cellTypes = murine ? murineCellTypes : humanCellTypes
    end
    Xmat = Matrix(Xdf[!, in(cellTypes).(names(Xdf))])
    return regPred(Xmat, opt; exp_method = exp_method)
end


function fitRegNNLS(Xdf::DataFrame; murine = true, cellTypes = nothing, exp_method = true)
    if cellTypes === nothing
        cellTypes = murine ? murineCellTypes : humanCellTypes
    end
    Xmat = Matrix(Xdf[!, in(cellTypes).(names(Xdf))])
    Y = Xdf[!, "Target"]
    cY = exp_method ? inv_exponential.(Y) : atanh.(Y)

    w = vec(nonneg_lsq(Xmat, cY; alg = :nnls))  # cell type weight found by NNLS
    Yr = Xmat * w
    R2val = R2(Y, (exp_method ? exponential(Yr) : tanh(Yr)); logscale = false)
    return regResult(w, R2val)
end

function regLOO(Xdf; kwargs...)
    n = size(Xdf, 1)
    fitResults = Vector{regResult}(undef, n)
    LOOindex = LOOCV(n)
    for (i, idx) in enumerate(LOOindex)
        fitResults[i] = fitRegNNLS(Xdf[idx, :]; kwargs...)
    end
    return fitResults
end

function regBootstrap(bootsize, Xdf; kwargs...)
    ress = Vector{regResult}(undef, bootsize)
    n = size(Xdf, 1)
    for b = 1:bootsize
        idx = rand(1:n, n)
        ress[b] = fitRegNNLS(Xdf[idx, :]; kwargs...)
    end
    return ress
end


function wildtypeWeights(res::regResult, df; L0 = 1e-9, f = 4, murine = true)
    # Prepare for cell type weights in wildtype
    wildtype = copy(importKav(; murine = murine, c1q = ("C1q" in names(df)), IgG2bFucose = (:IgG2bFucose in df.Condition), retdf = true))
    wildtype[!, "Background"] .= "wt"
    wildtype[!, "Target"] .= 0.0
    if !murine
        wildtype[!, "Genotype"] .= "ZZZ"
    end
    if "Neutralization" in names(df)
        wildtype[!, "Neutralization"] .= 0.0
    end
    rename!(wildtype, "IgG" => "Condition")

    wtXdf = modelPred(wildtype; L0 = L0, f = f, murine = murine, cellTypes = nothing)
    cellTypes = murine ? murineCellTypes : humanCellTypes

    wtXdf[!, cellTypes] .*= res.cellWs'
    wtXdf = wtXdf[!, vcat(["Condition"], cellTypes)]
    Cell_df = stack(wtXdf, Not("Condition"))
    rename!(Cell_df, "value" => "Weight")
    rename!(Cell_df, "variable" => "Component")
    return Cell_df
end


function regResult(dataType; L0, f, murine::Bool, exp_method = true)
    df = murine ? importDepletion(dataType) : importHumanized(dataType)

    Xdf = modelPred(df; L0 = L0, f = f, murine = murine, cellTypes = nothing)
    res = fitRegNNLS(Xdf; murine = murine, cellTypes = nothing, exp_method = exp_method)
    loo_res = regLOO(Xdf; murine = murine, cellTypes = nothing, exp_method = exp_method)
    boot_res = regBootstrap(10, Xdf; murine = murine, cellTypes = nothing, exp_method = exp_method)

    odf = df[!, in(["Condition", "Background"]).(names(df))]
    odf[!, "Concentration"] .= ("Concentration" in names(df)) ? (df[!, "Concentration"] .* L0) : L0
    odf[!, "Y"] = df[!, "Target"]
    odf[!, "Fitted"] = regPred(Xdf, res; murine = murine, exp_method = exp_method)
    odf[!, "LOOPredict"] = [regPred(Xdf[[ii], :], loo_res[ii]; murine = murine, exp_method = exp_method)[1] for ii = 1:length(loo_res)]

    if "Label" in names(df)
        odf[!, "Label"] .= df[!, "Label"]
    end
    if "Genotype" in names(df)
        odf[!, "Genotype"] .= df[!, "Genotype"]
    end

    return res, odf, loo_res, boot_res
end
