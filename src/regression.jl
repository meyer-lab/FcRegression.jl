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

tanh(X::Array) = tanh.(X)
tanh(X::Matrix, p::Vector) = tanh.(X * p)


mutable struct regParams{T}
    cellWs::Vector{T}
    ActIs::Vector{T}
    isMurine::Bool
end


function modelPred(dfr::DataFrameRow; f = 4, ActI = murineActI, Kav::DataFrame, Rtot = importRtot(; murine = true, retdf = true))
    Kav[Kav."IgG" .== "IgG2c", "IgG"] .= "IgG2a"

    IgGs = String[]
    cps = [1.0]
    if any(startswith.(names(dfr), "subclass"))
        iggv = Vector(dfr[startswith.(names(dfr), "subclass")])
        cps = Vector(dfr[startswith.(names(dfr), "%_")]) .* 1.0
        cps = cps[iggv .!= "PBS"]
        cps ./= sum(cps)
        append!(IgGs, iggv[iggv .!= "PBS"])
    elseif "Condition" in names(dfr)
        IgGs = [dfr."Condition"]
    elseif "Subclass" in names(dfr)
        IgGs = [dfr."Subclass"]
    end
    @assert length(IgGs) > 0
    IgGs[IgGs .== "IgG2c"] .= "IgG2a"
    Kav = Kav[in(IgGs).(Kav."IgG"), :]
    @assert size(Kav)[1] > 0

    if ("Background" in names(dfr)) && (dfr."Background" != "wt")
        rr_names = ["FcgRI", "FcgRIIB", "FcgRIII", "FcgRIV", ["FcgRI", "FcgRIIB", "FcgRIII", "FcgRIV"], ["FcgRI", "FcgRIII"], ["FcgRI", "FcgRIV"]]
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
    if ("Condition" in names(dfr)) && (dfr."Condition" == "IgG1D265A")
        Kav[:, ["FcgRI", "FcgRIIB", "FcgRIII", "FcgRIV"]] .= 0.0
    end
    Kav = Matrix(Kav[:, Not("IgG")])
    Rtot = Matrix(Rtot[:, Not("Receptor")])

    cellActs = Vector(undef, size(Rtot)[2])
    if length(cps) <= 0
        cellActs .= 0.0
        return cellActs
    end
    for jj = 1:size(Rtot)[2]
        pred = polyfc(dfr."Concentration", KxConst, f, Rtot[:, jj], cps, Kav).Rmulti_n
        cellActs[jj] = maximum([dot(ActI, pred), 0.0])
    end
    return cellActs
end

function modelPred(df::DataFrame; L0 = 1e-9, murine::Bool, cellTypes = nothing, ActI = nothing, Kav::DataFrame, kwargs...)
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
    #=if Kav === nothing
        @warn "Kav unprovided to modelPred(); use default "
        Kav = importKav(; murine = murine, retdf = true, IgG2bFucose = true)
    end=#

    ansType = ("Target" in names(df)) ? promote_type(eltype(df."Target"), eltype(ActI)) : eltype(ActI)
    Xfc = Array{ansType}(undef, size(df, 1), length(cellTypes))
    Threads.@threads for k = 1:size(df, 1)
        Xfc[k, :] = modelPred(df[k, :]; ActI = ActI, Kav = Kav, Rtot = importRtot(; murine = murine, retdf = true, cellTypes = cellTypes), kwargs...)
    end

    colls = murine ? murineFcgR : humanFcgR
    colls = vcat(["Target", "Concentration", "Baseline", "Measurement"], colls, cellTypes)
    ldf = df[!, (!).(in(colls).(names(df)))]
    Xdf = hcat(ldf, DataFrame(Xfc, cellTypes))
    if "C1q" in names(df)
        Xdf[!, "C1q"] = df[!, "C1q"] .* df[!, "Concentration"]
    end
    return Xdf
end


function regPred(df, opt::regParams; cellTypes = nothing, link::Function = exponential, kwargs...)
    if cellTypes === nothing
        cellTypes = opt.isMurine ? murineCellTypes : humanCellTypes
    end
    @assert length(cellTypes) == length(opt.cellWs)
    Xdf = modelPred(df; ActI = opt.ActIs, cellTypes = cellTypes, murine = opt.isMurine, kwargs...)
    Xmat = Matrix(Xdf[!, in(cellTypes).(names(Xdf))])
    return link(Xmat * opt.cellWs)
end

function wildtypeWeights(opt::regParams; cellTypes = nothing, kwargs...)
    # Prepare for cell type weights in wildtype
    IgGs = opt.isMurine ? murineIgG[murineIgG .!= "IgG3"] : humanIgG
    df = DataFrame(:Condition => IgGs, :Background .=> "wt")
    if !opt.isMurine
        df."Genotype" .= "ZZZ"
    end
    df = modelPred(df; ActI = opt.ActIs, murine = opt.isMurine, kwargs...)
    if cellTypes === nothing
        cellTypes = opt.isMurine ? murineCellTypes : humanCellTypes
    end
    df[!, cellTypes] .*= opt.cellWs'
    df = stack(df, cellTypes, variable_name = "Component", value_name = "Weight")
    return df[!, ["Condition", "Component", "Weight"]]
end


@model function regmodel(df, targets; murine::Bool, L0 = 1e-9, f = 4, cellTypes = nothing, Kav::DataFrame)
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

    Xdf = modelPred(df; L0 = L0, f = f, murine = murine, cellTypes = cellTypes, ActI = ActIs, Kav = Kav)
    Yfit = regPred(Xdf, regParams(cellWs, ActIs, murine); cellTypes = cellTypes, ActI = ActIs, Kav = Kav)

    stdv = std(Yfit - targets) / 10
    targets ~ MvNormal(Yfit, stdv * I)
    nothing
end

function runRegMCMC(dataType::Union{DataFrame, String}; mcmc_iter = 1_000, kwargs...)
    df = (dataType isa String) ? importDepletion(dataType) : dataType
    m = regmodel(df, df."Target"; kwargs...)
    opts = Optim.Options(iterations = 500, show_every = 10, show_trace = true)
    opt = optimize(m, MAP(), LBFGS(; m = 20), opts)
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
function runRegMAP(dataType::Union{DataFrame, String}; kwargs...)
    df = if (dataType isa String)
        murine ? importDepletion(dataType) : importHumanized(dataType)
    else
        dataType
    end
    m = regmodel(df, df."Target"; kwargs...)
    opts = Optim.Options(iterations = 500, show_every = 10, show_trace = true)
    opt = optimize(m, MAP(), LBFGS(; m = 20), opts)
    cdf = DataFrame(Parameter = String.(names(opt.values)[1]), Value = opt.values.array)

    LOOindex = LOOCV(size(df)[1])
    opts = Optim.Options(iterations = 500, show_trace = false)
    optcv = Vector{StatisticalModel}(undef, size(df)[1])
    for (i, idx) in enumerate(LOOindex)
        mv = regmodel(df[idx, :], df[idx, :]."Target"; kwargs...)
        optcv[i] = optimize(mv, MAP(), LBFGS(; m = 20), opts)
        cdf[!, "CV$i"] = optcv[i].values.array
    end

    cdf = combine(
        groupby(DataFrames.stack(cdf, Not(["Parameter", "Value"])), ["Parameter", "Value"]),
        "value" => median => "Median",
        "value" => (xs -> quantile(xs, 0.25)) => "xmin",
        "value" => (xs -> quantile(xs, 0.75)) => "xmax",
    )
    return opt, optcv, cdf
end

function extractRegMCMC(c::Union{Chains, StatisticalModel})
    pnames = [String(s) for s in (c isa Chains ? c.name_map[1] : names(c.values)[1])]
    ext(s::String) = c isa Chains ? median(c[s].data) : c.values[Symbol(s)]

    cellWs = [ext("cellWs[$i]") for i = 1:sum(startswith.(pnames, "cellWs"))]
    ActIs = [ext("ActIs[$i]") for i = 1:sum(startswith.(pnames, "ActIs"))]
    murine = length(ActIs) <= 4  # assume mice have only 4 receptors
    return regParams(cellWs, ActIs, murine)
end

function plotRegMCMC(
    c::Union{Chains, StatisticalModel, regParams},
    df::Union{DataFrame, String};
    ptitle = "",
    colorL = nothing,
    shapeL = nothing,
    kwargs...,
)
    if df isa String
        if ptitle === nothing
            ptitle = df
        end
        df = importDepletion(df)
    end
    if c isa Chains
        fits = hcat([regPred(df, extractRegMCMC(c[ii]); kwargs...) for ii = 1:length(c)]...)
        df."Fitted" .= mapslices(median, fits, dims = 2)
        df."ymax" .= mapslices(xs -> quantile(xs, 0.75), fits, dims = 2)
        df."ymin" .= mapslices(xs -> quantile(xs, 0.25), fits, dims = 2)
    else
        if c isa StatisticalModel
            c = extractRegMCMC(c)
        end
        df."Fitted" = regPred(df, c; kwargs...)
    end

    if shapeL === nothing
        shapeL = names(df)[1]
    end
    if colorL === nothing
        colorL = names(df)[2]
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
        "ymin" in names(df) ? Geom.errorbar : Guide.xlabel("Measured"),
        color = colorL,
        shape = shapeL,
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

function plotRegParams(c::Union{Chains, Vector{StatisticalModel}}; ptitle::String = "", legend = true, retdf = false, Kav::DataFrame)
    murine = extractRegMCMC(c[1]).isMurine
    df = if c isa Vector
        vcat([wildtypeWeights(extractRegMCMC(cc); murine = murine, Kav = Kav) for cc in c]...)
    else
        vcat([wildtypeWeights(extractRegMCMC(c[ii]); murine = murine, Kav = Kav) for ii = 1:length(c)]...)
    end
    df = combine(
        groupby(df, Not("Weight")),
        "Weight" => median => "Weight",
        "Weight" => (xs -> quantile(xs, 0.25)) => "ymin",
        "Weight" => (xs -> quantile(xs, 0.75)) => "ymax",
    )
    if retdf
        return plotCellTypeEffects(df, ptitle; legend = legend), df
    else
        return plotCellTypeEffects(df, ptitle; legend = legend)
    end
end
