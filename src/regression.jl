import MLBase.LOOCV
import Statistics: mean, quantile, std
import Distributions: cdf, Exponential
import StatsBase: sample, mode
import Base: tanh
import ProgressMeter: @showprogress

exponential(x::Real) = cdf(Exponential(), x)
exponential(X::Array) = cdf.(Exponential(), X)
exponential(X::Matrix, p::Vector) = cdf.(Exponential(), X * p)

tanh(X::Array) = tanh.(X)
tanh(X::Matrix, p::Vector) = tanh.(X * p)


mutable struct regParams{T}
    cellWs::NamedArray{T}
    ActIs::NamedArray{T}
    isMurine::Bool
end

@memoize receptorType(fcgr) = String(split(fcgr, "-")[1])

""" Calculate the Activity value, considering different genotype should have the same ActI """
function assembleActs(ActI::NamedVector, pred::NamedVector)
    return maximum([sum([ActI[receptorType(nn)] * pred[nn] for nn in names(pred)[1]]), 0.0])
end


"""
Make a 3D array that include all Rmulti predictions for each sample, each receptor, 
each cell type prediction (Row x Receptor x Cell).
"""
function regPrepareData(df::DataFrame = importHumanized("ITP");
        L0 = 1e-9,
        f = 4, 
        murine::Bool = false,
        cellTypes = nothing,
        Kav::DataFrame = importKav(; murine = false, retdf = true), 
    )

    # some initial variable check
    df = deepcopy(df)
    if cellTypes === nothing
        cellTypes = murine ? murineCellTypes : humanCellTypes
    end
    if "Concentration" in names(df)
        df[!, "Concentration"] ./= maximum(df[!, "Concentration"])
        df[!, "Concentration"] .*= L0
    else
        insertcols!(df, 3, "Concentration" => L0)
    end

    # create an empty 3D array for storing model predictions
    pred = NamedArray(
                zeros(eltype(df."Target"), nrow(df), size(Kav)[2]-1, length(cellTypes)), 
                (1:nrow(df), names(Kav)[2:end], cellTypes), 
                ("Row", "Receptor", "Cell")
            )

    for ii = 1:nrow(df)
        # import Rtot with the right genotype
        Kar = deepcopy(Kav)
        genotype = ("Genotype" in names(df)) ? df[ii, "Genotype"] : "XXX"
        genotype = genotype[1] * "I" * genotype[3:end]    # FcgRIIB has default genotype is 232I
        Rtott = importRtot(; murine = murine, retdf = true, cellTypes = cellTypes, genotype = genotype)
        Rtott = Rtott[in(names(Kar[!, Not("IgG")])).(Rtott."Receptor"), :]

        # checking all those complicated cases of murine/human
        Kar[Kar."IgG" .== "IgG2c", "IgG"] .= "IgG2a"
        IgGs = String[]
        cps = [1.0]
        if any(startswith.(names(df), "subclass"))
            iggv = Vector(df[ii, startswith.(names(df), "subclass")])
            cps = Vector(df[ii, startswith.(names(df), "%_")]) .* 1.0
            cps = cps[iggv .!= "PBS"]
            cps ./= sum(cps)
            append!(IgGs, iggv[iggv .!= "PBS"])
        elseif "Condition" in names(df)
            IgGs = [df[ii, "Condition"]]
        elseif "Subclass" in names(df)
            IgGs = [df[ii, "Subclass"]]
        end
        @assert length(IgGs) > 0
        IgGs[IgGs .== "IgG2c"] .= "IgG2a"
        Kar = Kar[in(IgGs).(Kar."IgG"), :]
        @assert size(Kar)[1] > 0
    
        if ("Background" in names(df)) && (df[ii, "Background"] != "wt")
            rr_names = ["FcgRI", "FcgRIIB", "FcgRIII", "FcgRIV", ["FcgRI", "FcgRIIB", "FcgRIII", "FcgRIV"], ["FcgRI", "FcgRIII"], ["FcgRI", "FcgRIV"]]
            for (ii, rr) in enumerate(["R1", "R2", "R3", "R4", "gc", "R1/3KO", "R1/4KO"])
                if occursin(rr, df[ii, "Background"])
                    Kar[:, rr_names[ii]] .= 0.0
                end
            end
            for cname in ["Neu", "ncMO", "cMO", "EO"]
                if occursin(cname, df[ii, "Background"])
                    Rtott[:, cname .== names(Rtott)] .= 0.0
                end
            end
        end
        if ("Condition" in names(df)) && (df[ii, "Condition"] == "IgG1D265A")
            Kar[:, ["FcgRI", "FcgRIIB", "FcgRIII", "FcgRIV"]] .= 0.0
        end
        
        # run model and make predictions
        Kar = Matrix(Kar[:, Not("IgG")])
        Rtott = Matrix(Rtott[:, Not("Receptor")])
        for jj = 1:size(Rtott)[2]
            pred[ii, :, jj] = polyfc(df[ii, "Concentration"], KxConst, f, Rtott[:, jj], cps, Kar).Rmulti_n
        end
    end
    return pred
end


function regPred(Rmulti::NamedArray, opt::regParams; link::Function = exponential)
    X = deepcopy(Rmulti[:, in(names(opt.ActIs)[1]).(receptorType.(names(Rmulti)[2])), names(opt.cellWs)[1]])
    Ys = zeros(promote_type(eltype(X), eltype(opt.cellWs), eltype(opt.ActIs)), size(X)[1])
    for ii = 1:size(X)[1]
        Xws = [assembleActs(opt.ActIs, X[ii, :, kk]) for kk = 1:size(X)[3]]
        Ys[ii] = link(sum(Xws .* opt.cellWs))
    end
    return Ys
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


@model function regmodel(Rmulti::NamedArray, targets)
    murine = !any(occursin.("-", names(Rmulti)[2]))      # if any receptor name contains "-", assume this is human
    ActI_means = murine ? murineActI : humanActI
    ActI_means = ActI_means[unique(receptorType.(names(Rmulti)[2]))]
    ActIs = NamedArray(repeat([0.0], length(ActI_means)), names(ActI_means)[1], ("Receptor"))
    for ii in eachindex(ActI_means)
        ActIs[ii] ~ Normal(ActI_means[ii], 1.0)
    end

    cellTypes = names(Rmulti)[3]
    cellWs = NamedArray(repeat([0.0], length(cellTypes)), cellTypes, ("CellType"))
    for ii in eachindex(cellTypes)
        cellWs[ii] ~ Exponential(Float64(length(cellTypes)))
    end

    Yfit = regPred(Rmulti, regParams(cellWs, ActIs, murine))
    stdv = std(Yfit - targets) / 10
    targets ~ MvNormal(Yfit, stdv * I)
    nothing
end


function runRegMCMC(df::DataFrame, fname = nothing; mcmc_iter = 1_000, kwargs...)
    if fname !== nothing
        fname = "cached/" * fname
    end
    if (fname !== nothing) && isfile(fname)
        deserial = deserialize(fname)
        if all(df."Target" .== deserial[3]."Target") && (Dict(kwargs) == deserial[4])
            println("Loading cached regression MCMC results from $fname...")
            return deserial[1:2]
        end
    end

    m = regmodel(regPrepareData(df; kwargs...), df."Target")
    opt = optimizeMAP(df; kwargs...)
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

    if fname !== nothing
        f = serialize(fname, [c, cdf, df, Dict(kwargs)])
    end
    return c, cdf
end


function optimizeMAP(df::DataFrame; repeat = 10, kwargs...)
    m = regmodel(regPrepareData(df; kwargs...), df."Target")
    success = false
    min_val = Inf
    min_opt = nothing
    opts = Optim.Options(iterations = 200, show_every = 10, show_trace = false)
    for r = 1:repeat
        opt = optimize(m, MAP(), LBFGS(; m = 50), opts)
        println(opt.optim_result.minimum, opt.optim_result.ls_success)
        if opt.optim_result.minimum < min_val
            min_val = opt.optim_result.minimum
            min_opt = opt
        end
        if opt.optim_result.ls_success
            success = true
        end
    end
    @assert success "MAP optimization was never successful even after $repeat tries."
    return min_opt
end


""" Run a MAP parameter estimation, with LOO/jackknife as errorbar """
function runRegMAP(df::DataFrame, fname = nothing; bootstrap = 10, kwargs...)
    if fname !== nothing
        fname = "cached/" * fname
    end
    if (fname !== nothing) && isfile(fname)
        deserial = deserialize(fname)
        if (df == deserial[4]) && (Dict(kwargs) == deserial[5])
            println("Loading cached regression MAP results from $fname...")
            return deserial[1:3]
        end
    end

    opt = optimizeMAP(df; kwargs...)
    mapdf = DataFrame(Parameter = String.(names(opt.values)[1]), Value = opt.values.array)

    optcv = Vector{StatisticalModel}(undef, bootstrap)
    @showprogress for b = 1:bootstrap
        idx = rand(1:size(df)[1], size(df)[1])
        optcv[b] = optimizeMAP(df[idx, :]; kwargs...)
        mapdf[!, "Boot$b"] = optcv[b].values.array
    end

    mapdf = combine(
        groupby(DataFrames.stack(mapdf, Not(["Parameter", "Value"])), ["Parameter", "Value"]),
        "value" => median => "Median",
        "value" => (xs -> quantile(xs, 0.25)) => "xmin",
        "value" => (xs -> quantile(xs, 0.75)) => "xmax",
    )
    if fname !== nothing
        f = serialize(fname, [opt, optcv, mapdf, df, Dict(kwargs)])
    end
    return opt, optcv, mapdf
end

function extractRegMCMC(c::Union{Chains, StatisticalModel}; cellTypes = nothing, FcgRs::Union{Nothing, Vector{String}, DataFrame} = nothing)
    pnames = [String(s) for s in (c isa Chains ? c.name_map[1] : names(c.values)[1])]
    ext(s::String) = c isa Chains ? median(c[s].data) : c.values[Symbol(s)]

    cellWs = [ext("cellWs[$i]") for i = 1:sum(startswith.(pnames, "cellWs"))]
    ActIs = [ext("ActIs[$i]") for i = 1:sum(startswith.(pnames, "ActIs"))]
    murine = length(ActIs) <= 4  # assume mice have only 4 receptors
    if any(in(FcgRs).(["FcgRIIA", "FcgRIIIA"]))
        murine = false
    end

    if cellTypes === nothing
        cellTypes = murine ? murineCellTypes : humanCellTypes
    end
    if FcgRs === nothing
        FcgRs = names(murine ? murineActI : humanActI)[1]
    elseif FcgRs isa DataFrame
        FcgRs = unique([receptorType(fcgr) for fcgr in names(FcgRs[!, Not("IgG")])])
    end
    return regParams(NamedArray(cellWs, cellTypes, "CellType"), NamedArray(ActIs, FcgRs, "Receptor"), murine)
end

function plotRegMCMC(
    c::Union{Chains, StatisticalModel, regParams},
    df::Union{DataFrame, String};
    ptitle = "",
    colorL = nothing,
    shapeL = nothing,
    Kav::DataFrame,
    cellTypes = nothing,
    legend = true,
    kwargs...,
)
    extractReg = x -> extractRegMCMC(x; cellTypes = cellTypes, FcgRs = unique([receptorType(fcgr) for fcgr in names(Kav[!, Not("IgG")])]))

    if df isa String
        if ptitle === nothing
            ptitle = df
        end
        df = importDepletion(df)
    end
    if c isa Chains
        fits = hcat([regPred(regPrepareData(df), extractReg(c[ii]); Kav = Kav, kwargs...) for ii = 1:length(c)]...)
        df."Fitted" .= mapslices(median, fits, dims = 2)
        df."ymax" .= mapslices(xs -> quantile(xs, 0.75), fits, dims = 2)
        df."ymin" .= mapslices(xs -> quantile(xs, 0.25), fits, dims = 2)
    else
        if c isa StatisticalModel
            c = extractReg(c)
        end
        df."Fitted" = regPred(regPrepareData(df), c; Kav = Kav, kwargs...)
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
        "ymin" in names(df) ? Geom.errorbar : Geom.point,
        color = colorL,
        shape = shapeL,
        Guide.colorkey(),
        Guide.shapekey(),
        Scale.y_continuous(minvalue = 0.0, maxvalue = 1.0),
        Geom.abline(color = "red"),
        Guide.xlabel("Actual effect"),
        Guide.ylabel("Fitted effect"),
        Guide.title("Actual vs fitted effect\n($ptitle)"),
        Guide.annotation(compose(context(), text(0.1, 0.8, R2anno), fill("black"), fontsize(10pt), font("Helvetica"))),
        style(point_size = 5px, key_position = legend ? :right : :none),
    )
    return pl
end

function plotRegParams(
    c::Union{Chains, Vector{StatisticalModel}};
    ptitle::String = "",
    legend = true,
    retdf = false,
    Kav::DataFrame,
    cellTypes = nothing,
)
    extractReg = x -> extractRegMCMC(x; cellTypes = cellTypes, FcgRs = unique([receptorType(fcgr) for fcgr in names(Kav[!, Not("IgG")])]))
    murine = extractReg(c[1]).isMurine
    df = vcat([wildtypeWeights(extractReg(c[ii]); cellTypes = names(extractReg(c[1]).cellWs)[1], murine = murine, Kav = Kav) for ii = 1:length(c)]...)

    ActI_df = vcat([DataFrame(:Receptor => names(extractReg(c[1]).ActIs)[1], :Weight => extractReg(c[ii]).ActIs) for ii = 1:length(c)]...)

    df = combine(
        groupby(df, Not("Weight")),
        "Weight" => median => "Weight",
        "Weight" => (xs -> quantile(xs, 0.25)) => "ymin",
        "Weight" => (xs -> quantile(xs, 0.75)) => "ymax",
    )
    ActI_df = combine(
        groupby(ActI_df, Not("Weight")),
        "Weight" => median => "Weight",
        "Weight" => (xs -> quantile(xs, 0.25)) => "ymin",
        "Weight" => (xs -> quantile(xs, 0.75)) => "ymax",
    )

    if retdf
        return plotCellTypeEffects(df, ptitle; legend = legend), plotActI(ActI_df, ptitle), df, ActI_df
    else
        return plotCellTypeEffects(df, ptitle; legend = legend), plotActI(ActI_df, ptitle)
    end
end
