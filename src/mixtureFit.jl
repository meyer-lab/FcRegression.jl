""" Fitting mixture measurements """

using Optim

function mixSqLoss(df; logscale = true)
    # square root differences of model prediction and adjusted measurements
    if "Adjusted" in names(df)
        adj = logscale ? log.(df."Adjusted") : df."Adjusted"
    else
        adj = logscale ? log.(df."Value") : df."Value"
    end
    pred = logscale ? log.(df."Predict") : df."Predict"
    return sum((adj .- pred) .^ 2)
end

function fit_experiments(conversion_facs, df)
    # called by fit coordinator, fit experiment conversion factors
    ansType = promote_type(eltype(df."Value"), eltype(conversion_facs))
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

function fitMixFunc(x, optDict::Dict, df = averageMixData(loadMixData()))
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
    df = averageMixData(loadMixData())
    optDict = Dict(
        "Experiment" => 0, #length(unique(df."Experiment")),
        "Valency" => 2,
        "Receptor" => length(measuredRecepExp),
        "KxStar" => 1
    )

    if !("Adjusted" in names(df))
        df."Adjusted" = df."Value"
    end

    lb = fitMixX0Bound(optDict; max = false)
    ub = fitMixX0Bound(optDict)
    x0 = fitMixX0(optDict)
    res = optimize(x -> fitMixFunc(x, optDict, df), lb, ub, x0)
    return lb, ub, x0, res
    
end




function PCAData(; cutoff = 0.9)
    df = loadMixData()
    exps_sets = [["6/23/20", "6/30/20", "7/14/20", "7/23/20", "9/11/20"], ["5/15/20", "5/20/20", "5/28/20", "6/2/20", "9/2/20"]]
    retdf = DataFrame()

    for (i, exps) in enumerate(exps_sets)
        ndf = df[in(exps).(df."Experiment"), :]
        widedf = unstack(ndf, ["Valency", "Cell", "subclass_1", "%_1", "subclass_2", "%_2"], "Experiment", "Value")
        widedf = coalesce.(widedf, 0)

        # Perform PCA
        mat = Matrix(widedf[!, exps])
        mat[mat .< 1.0] .= 1.0
        mat = log.(mat)
        M = fit(PCA, mat; maxoutdim = 2)
        recon = reconstruct(M, MultivariateStats.transform(M, mat))
        error = ((recon .- mat) .^ 2)

        # Impute by SVD
        matmiss = convert(Array{Union{Float64, Missing}}, mat)
        matmiss[error .> quantile(reshape(error, :), [cutoff])] .= missing
        Impute.impute!(matmiss, Impute.SVD())
        widedf[!, exps] .= exp.(matmiss)

        ndf = stack(widedf, exps)
        rename!(ndf, "variable" => "Experiment")
        rename!(ndf, "value" => "Value")

        append!(retdf, ndf)
    end
    return sort!(retdf, ["Valency", "Cell", "subclass_1", "subclass_2", "Experiment", "%_2"])
end


function PCA_dimred()
    df = loadMixData()
    mdf = unstack(df, ["Valency", "Cell", "subclass_1", "%_1", "subclass_2", "%_2"], "Experiment", "Value")
    mat = Matrix(mdf[!, Not(["Valency", "Cell", "subclass_1", "%_1", "subclass_2", "%_2"])])
    Impute.impute!(mat, Impute.SVD())

    # to change Matrix type without Missing
    rmat = zeros(size(mat))
    rmat .= mat

    M = fit(PCA, rmat; maxoutdim = 1)
    mdf = mdf[!, ["Valency", "Cell", "subclass_1", "%_1", "subclass_2", "%_2"]]
    mdf."PCA" = projection(M)[:, 1]
    mdf = predictMix(mdf)
    mdf."PCA" *= mean(mdf."Predict") / mean(mdf."PCA")
    mdf[mdf."PCA" .< 1.0, "PCA"] .= 1.0
    return mdf
end
