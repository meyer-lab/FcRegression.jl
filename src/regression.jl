import MLBase.LOOCV
import Statistics: mean, quantile, std
import Distributions: cdf, Exponential
import StatsBase: sample, mode
using NonNegLeastSquares
import Base: tanh
using InverseFunctions

exponential(x::Real) = cdf(Exponential(), x)
exponential(X::Array) = cdf.(Exponential(), X)
exponential(X::Matrix, p::Vector) = cdf.(Exponential(), X * p)
InverseFunctions.inverse(::typeof(exponential)) = y -> -log(1 - y)

tanh(X::Array) = tanh.(X)
tanh(X::Matrix, p::Vector) = tanh.(X * p)
InverseFunctions.inverse(::typeof(tanh)) = atanh

mutable struct regResult{T}
    cellWs::Array{T}
    R2::T
end


""" plug in reg df, and output binding model results """
function modelPred(df; L0, f, murine::Bool = true, cellTypes = nothing, ActI = nothing)
    df = deepcopy(df)
    FcRecep = murine ? murineFcgR : humanFcgR
    if ActI == nothing
        ActI = murine ? murineActI : humanActI
    end
    if cellTypes === nothing
        cellTypes = murine ? murineCellTypes : humanCellTypes
    end

    if "Concentration" in names(df)
        df[!, "Concentration"] ./= maximum(df[!, "Concentration"])
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
            for cname in ["Neu", "ncMO", "cMO", "EO"]
                if occursin(cname, df[k, "Background"])
                    Rtot[:, cname .== cellTypes] .= 0.0
                end
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

function regPred(
        Xdf::DataFrame, 
        opt::regResult; 
        murine = true, 
        link::Function = exponential, 
        cellTypes = nothing, 
        ActI = nothing
    )
    if cellTypes === nothing
        cellTypes = murine ? murineCellTypes : humanCellTypes
    end
    Xmat = Matrix(Xdf[!, in(cellTypes).(names(Xdf))])

    return link(Xmat * opt.cellWs)
end


function fitRegNNLS(Xdf::DataFrame; murine = true, cellTypes = nothing, link::Function = exponential)
    if cellTypes === nothing
        cellTypes = murine ? murineCellTypes : humanCellTypes
    end
    Xmat = Matrix(Xdf[!, in(cellTypes).(names(Xdf))])
    Y = Xdf[!, "Target"]
    @assert !(inverse(link) isa NoInverse)
    cY = inverse(link).(Y)

    w = vec(nonneg_lsq(Xmat, cY; alg = :nnls))  # cell type weight found by NNLS
    Yr = Xmat * w
    R2val = R2(Y, link(Yr); logscale = false)
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


function wildtypeWeights(res::regResult, df; L0 = 1e-9, f = 4, murine = true, Kav = nothing, cellTypes = nothing, ActI = nothing)
    # Prepare for cell type weights in wildtype
    Kavd = importKav(; murine = murine, c1q = ("C1q" in names(df)), IgG2bFucose = (:IgG2bFucose in df.Condition), retdf = true)
    if Kav !== nothing
        if murine
            Kav[Kav."IgG" .== "IgG2c", "IgG"] .= "IgG2a"
        end
        # replace the value in
        for igg in Kav."IgG"
            Kavd[Kavd."IgG" .== igg, names(Kav)[2:end]] = Kav[Kav."IgG" .== igg, names(Kav)[2:end]]
        end
    end

    wildtype = copy(Kavd)
    wildtype[!, "Background"] .= "wt"
    wildtype[!, "Target"] .= 0.0
    if !murine
        wildtype[!, "Genotype"] .= "ZZZ"
    end
    if "Neutralization" in names(df)
        wildtype[!, "Neutralization"] .= 0.0
    end
    rename!(wildtype, "IgG" => "Condition")
    
    if cellTypes == nothing
        cellTypes = murine ? murineCellTypes : humanCellTypes
    end
    wtXdf = modelPred(wildtype; L0 = L0, f = f, murine = murine, cellTypes = cellTypes, ActI = ActI)

    wtXdf[!, cellTypes] .*= res.cellWs'
    wtXdf = wtXdf[!, vcat(["Condition"], cellTypes)]
    Cell_df = stack(wtXdf, Not("Condition"), variable_name = "Component", value_name = "Weight")
    return Cell_df
end


function regResult(dataType; L0 = 1e-9, f = 4, murine::Bool = true, link::Function = exponential, Kav = nothing, cellTypes = nothing, ActI = nothing)
    df = murine ? importDepletion(dataType; Kav = Kav) : importHumanized(dataType)

    Xdf = modelPred(df; L0 = L0, f = f, murine = murine, cellTypes = cellTypes, ActI = ActI)
    res = fitRegNNLS(Xdf; murine = murine, cellTypes = cellTypes, link = link)
    loo_res = regLOO(Xdf; murine = murine, cellTypes = cellTypes, link = link)
    boot_res = regBootstrap(10, Xdf; murine = murine, cellTypes = cellTypes, link = link)

    odf = df[!, in(["Condition", "Background"]).(names(df))]
    odf[!, "Concentration"] .= ("Concentration" in names(df)) ? (df[!, "Concentration"] .* L0) : L0
    odf[!, "Y"] = df[!, "Target"]
    odf[!, "Fitted"] = regPred(Xdf, res; murine = murine)
    odf[!, "LOOPredict"] = [regPred(Xdf[[ii], :], loo_res[ii]; murine = murine, link = link)[1] for ii = 1:length(loo_res)]

    if "Label" in names(df)
        odf[!, "Label"] .= df[!, "Label"]
    end
    if "Genotype" in names(df)
        odf[!, "Genotype"] .= df[!, "Genotype"]
    end

    # Cell type weight
    Cell_df = wildtypeWeights(res, df; L0 = L0, f = f, murine = murine, Kav = Kav, cellTypes = cellTypes, ActI = ActI)
    Cell_loo = vcat([wildtypeWeights(loo, df; Kav = Kav, cellTypes = cellTypes, ActI = ActI) for loo in loo_res]...)
    Cell_conf = combine(groupby(Cell_loo, ["Condition", "Component"]), "Weight" => lower => "ymin", "Weight" => upper => "ymax")
    Cell_df = innerjoin(Cell_df, Cell_conf, on = ["Condition", "Component"])

    return res, odf, loo_res, boot_res, Cell_df
end
