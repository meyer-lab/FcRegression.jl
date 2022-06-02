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
inv_exponential(y::Real) = -log(1 - y)
InverseFunctions.inverse(::typeof(exponential)) = y -> -log(1 - y)

tanh(X::Array) = tanh.(X)
tanh(X::Matrix, p::Vector) = tanh.(X * p)
InverseFunctions.inverse(::typeof(tanh)) = atanh

mutable struct regResult{T}
    cellWs::Array{T}
    R2::T
end

mutable struct regParams{T}
    cellWs::Vector{T}
    ActIs::Vector{T}
end


function modelPred(dfr::DataFrameRow; f = 4, ActI = murineActI,
    Kav = importKav(; murine = true, retdf = true, IgG2bFucose = true), 
    Rtot = importRtot(; murine = true, retdf = true))

    ## upper layer deal with L0, murine, cellTypes
    ## cellTypes selects Rtot in upper layer!

    IgGs = String[]
    cps = [1.0]
    if "Condition" in names(dfr)
        IgGs = [dfr."Condition"]
    elseif any(startswith.(names(dfr), "subclass"))
        iggv = Vector(dfr[startswith.(names(dfr), "subclass")])
        cps = Vector(dfr[startswith.(names(dfr), "%_")])
        cps = cps[iggv .!= "PBS"]
        cps ./= sum(cps)
        append!(IgGs, iggv[iggv .!= "PBS"])
    end
    IgGs[IgGs .== "IgG2c"] .= "IgG2a"
    Kav = Kav[in(IgGs).(Kav."IgG"), :]

    if ("Background" in names(dfr)) && (dfr."Background" != "wt")
        rr_names = ["FcgRI", "FcgRIIB", "FcgRIII", "FcgRIV", ["FcgRI", "FcgRIIB", "FcgRIII", "FcgRIV"], 
            ["FcgRI", "FcgRIII"], ["FcgRI", "FcgRIV"]]
        for (ii, rr) in enumerate(["R1", "R2", "R3", "R4", "gc", "R1/3KO", "R1/4KO"])
            if occursin(rr, dfr."Background")
                Kav[:, rr_names[ii]] .= 0.0
            end
        end
        for cname in ["Neu", "ncMO", "cMO", "EO"]
            if occursin(cname, dfr."Background")
                Rtot[:, cname .== names(Rtot)] .= 0.0
            end
        end
    end
    if dfr."Condition" == "IgG1D265A"
        Kav[:, ["FcgRI", "FcgRIIB", "FcgRIII", "FcgRIV"]] .= 0.0
    end
    Kav = Matrix(Kav[:, Not("IgG")])
    Rtot = Matrix(Rtot[:, Not("Receptor")])

    cellActs = zeros(size(Rtot)[2])
    for jj = 1:size(Rtot)[2]
        pred = polyfc(dfr."Concentration", KxConst, f, Rtot[:, jj], cps, Kav).Rmulti_n
        cellActs[jj] = maximum([dot(ActI, pred), 0.0])
    end
    return cellActs
end

function modelPred(df::DataFrame; L0 = 1e-9, murine::Bool = true, cellTypes = nothing, ActI = nothing, kwargs...)
    df = deepcopy(df)
    if ActI === nothing
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

    ansType = promote_type(eltype(df."Target"), eltype(ActI))
    Xfc = Array{ansType}(undef, size(df, 1), length(cellTypes))  
    Threads.@threads for k = 1:size(df, 1)
        Xfc[k, :] = modelPred(df[k, :]; ActI = ActI, 
             Kav = importKav(; murine = murine, retdf = true, IgG2bFucose = true),
            Rtot = importRtot(; murine = murine, retdf = true, cellTypes = cellTypes),
            kwargs...)
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

function regPred(df, opt::regParams; cellTypes = nothing, murine::Bool = true, link::Function = exponential, kwargs...)
    if cellTypes === nothing
        cellTypes = murine ? murineCellTypes : humanCellTypes
    end
    @assert length(cellTypes) == length(opt.cellWs)
    Xdf = modelPred(df; ActI = opt.ActIs, cellTypes = cellTypes, murine = murine, kwargs...)
    Xmat = Matrix(Xdf[!, in(cellTypes).(names(Xdf))])
    return link(Xmat * opt.cellWs)
end


function fitRegNNLS(Xdf::DataFrame; murine = true, cellTypes = nothing, link::Function = exponential, inv_link::Function = inv_exponential)
    if cellTypes === nothing
        cellTypes = murine ? murineCellTypes : humanCellTypes
    end
    Xmat = Matrix(Xdf[!, in(cellTypes).(names(Xdf))])
    Y = Xdf[!, "Target"]
    cY = if inverse(link) isa NoInverse
        inv_link.(Y)
    else
        inverse(link).(Y)
    end

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


function regResult(dataType; L0 = 1e-9, f = 4, murine::Bool = true, link::Function = exponential, Kav = nothing, cellTypes = nothing, ActI = nothing, inv_link::Function = inv_exponential)
    df = murine ? importDepletion(dataType; Kav = Kav) : importHumanized(dataType)

    if (cellTypes == nothing) && (dataType in ["melanoma", "ITP"])
        if dataType == "melanoma"
            cellTypes = ["ncMO", "cMO", "NKs", "Neu", "EO"]
        elseif dataType == "ITP"
            cellTypes = ["ncMO", "cMO", "NKs", "Neu", "EO", "Kupffer", "KupfferHi"]
        end
    end

    Xdf = modelPred(df; L0 = L0, f = f, murine = murine, cellTypes = cellTypes, ActI = ActI)
    res = fitRegNNLS(Xdf; murine = murine, cellTypes = cellTypes, link = link, inv_link = inv_link)
    loo_res = regLOO(Xdf; murine = murine, cellTypes = cellTypes, link = link, inv_link = inv_link)
    boot_res = regBootstrap(10, Xdf; murine = murine, cellTypes = cellTypes, link = link, inv_link = inv_link)

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


@model function regmodel(df, targets; murine = true, L0 = 1e-9, f = 4, cellTypes = nothing)
    ActI_means = murine ? murineActI : humanActI
    ActIs = Vector(undef, length(ActI_means))
    for ii in eachindex(ActI_means)
        ActIs[ii] ~ Normal(ActI_means[ii], 1.0)
    end

    if cellTypes === nothing
        cellTypes = murine ? murineCellTypes : humanCellTypes
    end
    cellWs = Vector(undef, length(cellTypes))
    for ii in eachindex(cellTypes)
        cellWs[ii] ~ Exponential(Float64(length(cellTypes)))
    end

    Xdf = modelPred(df; L0 = L0, f = f, murine = murine, cellTypes = cellTypes, ActI = ActIs)
    Yfit = regPred(Xdf, regResult(cellWs, 0.0); murine = murine, cellTypes = cellTypes, ActI = ActIs)

    stdv = std(Yfit - targets) / 10
    targets ~ MvNormal(Yfit, stdv * I)
    nothing
end

function runRegMCMC(dataType; mcmc_iter = 1_000, ret_opt = false)
    df = importDepletion(dataType)
    m = regmodel(df, df."Target")
    opts = Optim.Options(iterations = 500, show_every = 10, show_trace = true)
    opt = optimize(m, MAP(), LBFGS(; m = 20), opts)
    if ret_opt
        return opt
    end
    c = sample(m, NUTS(), mcmc_iter, init_params = opt.values.array)

    # Put model parameters into a df
    cdf = DataFrame(c)
    cdf = cdf[!, startswith.(names(cdf), "ActIs") .| startswith.(names(cdf), "cellWs")]
    cdf = stack(cdf, variable_name = "Parameter")
    cdf = combine(
        groupby(cdf, "Parameter"),
        "value" => median => "Value",
        "value" => (xs -> quantile(xs, 0.25)) => "xmin",
        "value" => (xs -> quantile(xs, 0.75)) => "xmax",
    )
    cdf."MAP" = opt.values.array
    return c, cdf
end

""" Run a MAP parameter estimation, with LOO/jackknife as errorbar """
function runRegMAP(dataType; mcmc_iter = 1_000, ret_opt = false)
    df = importDepletion(dataType)
    m = regmodel(df, df."Target")
    opts = Optim.Options(iterations = 500, show_every = 10, show_trace = true)
    opt = optimize(m, MAP(), LBFGS(; m = 20), opts)
    cdf = DataFrame(Parameter = String.(names(opt.values)[1]), Value = opt.values.array)

    LOOindex = LOOCV(size(df)[1])
    for (i, idx) in enumerate(LOOindex)
        mv = regmodel(df[idx, :], df[idx, :]."Target")
        optv = optimize(mv, MAP(), LBFGS(; m = 20), opts)
        cdf[!, "CV$i"] = optv.values.array
    end

    cdf = combine(
        groupby(DataFrames.stack(cdf, Not(["Parameter", "Value"])), ["Parameter", "Value"]),
        "value" => median => "Median",
        "value" => (xs -> quantile(xs, 0.25)) => "xmin",
        "value" => (xs -> quantile(xs, 0.75)) => "xmax",
    )
    return cdf
end

function extractRegMCMC(c::Union{Chains, StatisticalModel})
    pnames = [String(s) for s in (c isa Chains ? c.name_map[1] : names(c.values)[1])]
    ext(s::String) = c isa Chains ? median(c[s].data) : c.values[Symbol(s)]

    cellWs = [ext("cellWs[$i]") for i = 1:sum(startswith.(pnames, "cellWs"))]
    ActIs = [ext("ActIs[$i]") for i = 1:sum(startswith.(pnames, "ActIs"))]
    return regParams(cellWs, ActIs)
end

function plotRegMCMC(c::Union{Chains, StatisticalModel, regParams}, df::Union{DataFrame, String}; L0 = 1e-9, f = 4, ptitle = "",  kwargs...)
    if df isa String
        if ptitle === nothing
            ptitle = df
        end
        df = importDepletion(df)
    end
    if c isa Chains
        fits = hcat([regPred(df, extractRegMCMC(c[ii]); L0 = 1e-9, f = 4) for ii = 1:length(c)]...)
        df."Fitted" .= mapslices(median, fits, dims = 2)
        df."ymax" .= mapslices(xs -> quantile(xs, 0.75), fits, dims = 2)
        df."ymin" .= mapslices(xs -> quantile(xs, 0.25), fits, dims = 2)
    else
        if c isa StatisticalModel
            c = extractRegMCMC(c)
        end
        df."Fitted" = regPred(df, c; L0 = L0, f = f, kwargs...)
    end

    setGadflyTheme()
    R2anno = "<i>R</i><sup>2</sup>" * @sprintf("=%.3f", R2(df.Target, df.Fitted; logscale = false))
    pl = plot(
        df,
        x = "Target",
        y = "Fitted",
        ymin = (c isa Chains ? "ymin" : "Fitted"),
        ymax = (c isa Chains ? "ymax" : "Fitted"),
        Geom.point,
        "ymin" in names(df) ? Geom.errorbar : Guide.xlabel(xxlabel),
        color = "Background",
        shape = "Condition",
        Guide.colorkey(),
        Guide.shapekey(),
        Scale.y_continuous(minvalue = 0.0, maxvalue = 1.0),
        Geom.abline(color = "red"),
        Guide.xlabel("Actual effect"),
        Guide.ylabel("Fitted effect"),
        Guide.title("Actual vs fitted effect ($ptitle)"),
        Guide.annotation(compose(context(), text(0.1, 0.8, R2anno), fill("black"), fontsize(10pt), font("Helvetica"))),
        style(point_size = 5px, key_position = :right),
    )
    return pl
end

function plotRegParams(cdf::DataFrame)
    
end