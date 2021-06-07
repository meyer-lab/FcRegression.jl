using Optim

function mixSqLoss(df; logscale = false)
    # square root differences of model prediction and adjusted measurements
    adj = logscale ? log.(df."Adjusted") : df."Adjusted"
    pred = logscale ? log.(df."Predict") : df."Predict"
    return sum((adj .- pred) .^ 2)
end

function fit_experiments(conversion_facs, df)
    # called by fit coordinator, fit experiment conversion factors
    ansType = promote_type(eltype(df."Value"), eltype(conversion_facs))
    if !("Adjusted" in names(df))
        df."Adjusted" = convert.(ansType, df."Value")
    end
    exps = sort(unique(df."Experiment"))
    @assert length(conversion_facs) == length(exps)
    for (ii, exp) in enumerate(exps)
        df[df."Experiment" .== exp, "Adjusted"] .*= conversion_facs[ii]
    end
    return df
end

function fit_valencies(new_vals, df)
    # called by fit coordinator, fit new valencies
    if !("NewValency" in names(df))
        df."NewValency" = convert.(eltype(new_vals), df."Valency")
    end
    vals = sort(unique(df."Valency"))
    @assert length(vals) == length(new_vals)
    for (ii, val) in enumerate(vals)
        df[df."Valency" .== val, "NewValency"] .= new_vals[ii]
    end
    return df
end

function fitMixFunc(x, optDict::Dict, df = loadMixData())
    # order: experiments (log), valencies (log), receptor amounts (log), Kx* (log)
    df = copy(df)
    iidx = 0
    if optDict["Experiment"] > 0
        df = fit_experiments(exp.(x[(iidx + 1):(iidx + optDict["Experiment"])]), df)
        iidx += optDict["Experiment"]
    end
    if optDict["Valency"] > 0
        df = fit_valencies(exp.(x[(iidx + 1):(iidx + optDict["Valency"])]), df)
        iidx += optDict["Valency"]
    end
    # receptor amounts and Kx*
    recepExp = Dict{String, Any}(measuredRecepExp)
    if optDict["Receptor"] > 0
        rcps = sort(collect(measuredRecepExp), by = x->x[1])
        for (rcp, _) in rcps
            recepExp[rcp] = exp(x[iidx + 1])
            iidx += 1
        end
    end
    KxStar = (optDict["KxStar"] > 0) ? exp(x[iidx + 1]) : KxConst
    iidx += optDict["KxStar"]
    @assert length(x) == iidx
    return mixSqLoss(predictMix(df; recepExp = recepExp, KxStar = KxStar))
end

function fitMixX0(optDict::Dict)
    """ Provide the initial point for optimization """
    x0 = Vector{Float64}([])
    if optDict["Experiment"] > 0
        append!(x0, zeros(optDict["Experiment"]))
    end
    if optDict["Valency"] > 0
        append!(x0, log.([4.0, 33.0]))
    end
    if optDict["Receptor"] > 0
        recepExp = copy(measuredRecepExp)
        rcps = sort(collect(recepExp), by = x->x[1])
        append!(x0, [log(recepExp[rcp]) for (rcp, _) in rcps])
    end
    if optDict["KxStar"] > 0
        append!(x0, [log(KxConst)])
    end
    return x0
end

function fitMixX0Bound(optDict::Dict; max = true)
    """ Provide the max/min bound for optimization """
    x0 = Vector{Float64}([])
    if optDict["Experiment"] > 0
        append!(x0, ones(optDict["Experiment"]) .* (max ? 6 : -6))
    end
    if optDict["Valency"] > 0
        append!(x0, log.([4.0, 33.0]) .+ (max ? 1 : -1))
    end
    if optDict["Receptor"] > 0
        recepExp = copy(measuredRecepExp)
        rcps = sort(collect(recepExp), by = x->x[1])
        append!(x0, [log(recepExp[rcp]) .+ (max ? 3 : -3) for (rcp, _) in rcps])
    end
    if optDict["KxStar"] > 0
        append!(x0, [log(KxConst) .+ (max ? 7 : -7)])
    end
    return x0
end


function fitMixMaster()
    df = loadMixData(; avg = false)
    optDict = Dict(
        "Experiment" => length(unique(df."Experiment")),
        "Valency" => 2,
        "Receptor" => length(measuredRecepExp),
        "KxStar" => 0,
    )

    lb = fitMixX0Bound(optDict; max = false)
    ub = fitMixX0Bound(optDict)
    x0 = fitMixX0(optDict)
    res = optimize(x -> fitMixFunc(x, optDict, df), x0)


    
end