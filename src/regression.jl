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


mutable struct regParams
    cellWs::NamedArray
    ActIs::NamedArray
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
function regPrepareData(
    df::DataFrame = importHumanized();
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
        zeros(eltype(Kav[!, 2]), nrow(df), size(Kav)[2] - 1, length(cellTypes)),
        (1:nrow(df), names(Kav)[2:end], cellTypes),
        ("Row", "Receptor", "Cell"),
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
    X = Rmulti[:, in(names(opt.ActIs)[1]).(receptorType.(names(Rmulti)[2])), names(opt.cellWs)[1]]
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
        df."Genotype" .= "ZIZ"
    end
    Rmulti = regPrepareData(df; kwargs...)
    X = Rmulti[:, in(names(opt.ActIs)[1]).(receptorType.(names(Rmulti)[2])), names(opt.cellWs)[1]]

    cellTypes = names(opt.cellWs)[1]
    for cell in cellTypes
        df[!, cell] = [assembleActs(opt.ActIs, X[ii, :, cell]) for ii = 1:size(X)[1]]
    end
    df[!, cellTypes] .*= opt.cellWs'
    df = stack(df, cellTypes, variable_name = "Component", value_name = "Weight")
    return df[!, ["Condition", "Component", "Weight"]]
end


@model function regmodel(Rmulti::NamedArray, targets; link::Function = exponential, fitActI = true)
    murine = !any(occursin.("-", names(Rmulti)[2]))      # if any receptor name contains "-", assume this is human
    rcpNames = unique(receptorType.(names(Rmulti)[2]))
    ActIs = NamedArray(repeat([0.0], length(rcpNames)), rcpNames, ("Receptor"))
    for rcp in names(ActIs)[1]
        if contains(rcp, "RIIB")
            ActIs[rcp] = -1.0                   # RIIB is set to be always -1.0
        elseif fitActI
            ActIs[rcp] ~ Normal(1.0, 0.1)       # only fit activating receptors
        else
            ActIs[rcp] = 1.0                    # not fitting ActI then all 1.0
        end
    end

    cellTypes = names(Rmulti)[3]
    cellWs = NamedArray(repeat([0.0], length(cellTypes)), cellTypes, ("CellType"))
    for ct in names(cellWs)[1]
        cellWs[ct] ~ Exponential(Float64(length(cellTypes)) * 2.0)
    end

    Yfit = regPred(Rmulti, regParams(cellWs, ActIs, murine); link = link)
    stdv = std(Yfit - targets) / 10
    targets ~ MvNormal(Yfit, stdv * I)
    nothing
end


function runRegMCMC(df::DataFrame, fname = nothing; mcmc_iter = 1_000, link::Function = exponential, fitActI = true, kwargs...)
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

    m = regmodel(regPrepareData(df; kwargs...), df."Target"; link = link, fitActI = fitActI)
    opt = optimizeMAP(df; link = link, fitActI = fitActI, kwargs...)
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


function optimizeMAP(df::DataFrame; repeat = 10, link::Function = exponential, fitActI = true, kwargs...)
    m = regmodel(regPrepareData(df; kwargs...), df."Target"; link = link, fitActI = fitActI)
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
function runRegMAP(df::DataFrame, fname = nothing; bootstrap = 10, link::Function = exponential, fitActI = true, kwargs...)
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

    opt = optimizeMAP(df; link = link, fitActI = fitActI, kwargs...)
    mapdf = DataFrame(Parameter = String.(names(opt.values)[1]), Value = opt.values.array)

    optcv = Vector{StatisticalModel}(undef, bootstrap)
    @showprogress for b = 1:bootstrap
        idx = rand(1:size(df)[1], size(df)[1])
        optcv[b] = optimizeMAP(df[idx, :]; link = link, fitActI = fitActI, kwargs...)
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

function extractRegMCMC(c::Union{Chains, StatisticalModel})
    pnames = [String(s) for s in (c isa Chains ? c.name_map[1] : names(c.values)[1])]
    ext(s::String) = c isa Chains ? median(c[s].data) : c.values[Symbol(s)]
    extVar(x::String) = split(x, "[")[2][1:(end - 1)]
    splitVar(xs::Vector{String}, pref::String) = String.(extVar.(xs[startswith.(xs, pref)]))

    murine = false    # only dealing with humanized mice in this work
    #murine = any(contains.(pnames, "RIIA")) || any(contains.(pnames, "RIIIA")) ? false : true
    ActIs = murine ? murineActI : humanActI
    for rcp in splitVar(pnames, "ActIs")
        ActIs[rcp] = ext("ActIs[$rcp]")
    end
    cell_entries = splitVar(pnames, "cellWs")
    cellWs = NamedArray([ext("cellWs[$s]") for s in cell_entries], cell_entries, "Cell")
    return regParams(cellWs, ActIs, murine)
end
