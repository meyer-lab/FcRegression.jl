using NonNegLeastSquares


mutable struct optResult{T}
    cellWs::Array{T}
    ActI::Array{T}
    residual::T
end

function modelPred(df; L0, f, murine::Bool = true)
    df = copy(df)
    FcRecep = murine ? murineFcgR : humanFcgR

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
    Xdf = df[!, in(["Condition", "Concentration", "C1q", "Neutralization"]).(names(df))]
    if "C1q" in names(df)
        Xdf[!, "C1q"] = Xdf[!, "C1q"] .* Xdf[!, "Concentration"]
    end

    Y = df[!, "Target"]
    return Xfc, Xdf, Y
end


function regressionPred(Xfc, Xdf, cellWeights, recepActI; showXmat=false)
    ansType = promote_type(eltype(Xfc), eltype(cellWeights), eltype(recepActI))
    extra = Xdf[!, in(["C1q", "Neutralization"]).(names(Xdf))]
    @assert length(cellWeights) == size(Xfc, 1) + size(extra, 2)
    @assert length(recepActI) == size(Xfc, 2)
    @assert size(Xfc, 3) == size(Xdf, 1)

    Xmat = DataFrame(repeat([ansType], length(cellWeights)), vcat(FcRegression.cellTypes, names(extra)))
    Ypred = Array{ansType}(undef, size(Xfc, 3))
    for k = 1:size(Xfc, 3)
        recepEff = Xfc[:, :, k] * recepActI
        recepEff[recepEff .< 0.0] .= 0.0
        if size(extra, 2) > 0
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

regressionPred(Xfc, Xdf, fit::optResult; showXmat=false) = regressionPred(Xfc, Xdf, fit.cellWs, fit.ActI; showXmat=showXmat)

function old_opt(Xfc, extra, Y, ActI)
    cY = inv_exponential.(Y)
    Xmat = Matrix{eltype(Xfc)}(undef, size(Xfc, 3), size(Xfc, 1))
    for i in 1:size(Xfc, 1)
        for j in 1:size(Xfc, 3)
            Xmat[j, i] = Xfc[i, :, j]' * ActI
        end
    end
    Xmat[Xmat .< 0.0] .= 0.0
    if !any(size(extra) .<= 0)
        extraM = Matrix(extra)
        @assert size(extraM, 1) == size(Xfc, 3)
        Xmat = hcat(Xmat, extraM)
    end
    w = nonneg_lsq(Xmat, cY)
    Yr = Xmat * w
    residual = norm(Y-Yr, 2)/length(Y)

    return w, residual
end


function fitRegression2(Xfc, Xdf, Y; murine::Bool=true)
    cY = inv_exponential.(Y)
    extra = Xdf[!, in(["C1q", "Neutralization"]).(names(Xdf))]
    cellWlen = size(Xfc, 1) + size(extra, 2)

    """func = x -> old_opt(Xfc, extra, Y, x)[2]
    ActI_init = Float64.(murine ? murineActI : humanActI)

    inits = ones(length(ActI_init))
    lower = repeat([-2.0], length(ActI_init))
    upper = repeat([10.0], length(ActI_init))

    opt = optimize(func, lower, upper, ActI_init)
    ActI = opt.minimizer"""

    ActI = Float64.(murine ? murineActI : humanActI)
    cellWs, residual = old_opt(Xfc, extra, Y, ActI)
    res = optResult(cellWs, ActI, residual)
    return res
end


function LOOCrossVal2(Xfc, Xdf, Y; murine)
    n = size(Xdf, 1)
    fitResults = Vector{optResult}(undef, n)
    LOOindex = LOOCV(n)
    for (i, idx) in enumerate(LOOindex)
        fitResults[i] = fitRegression2(Xfc[:,:,idx], Xdf[idx,:], Y[idx]; murine = murine)
    end
    return fitResults
end


function bootstrap2(Xfc, Xdf, Y; nsample = 100, murine)
    n = size(Xdf, 1)
    fitResults = Vector{optResult}(undef, nsample)
    for i = 1:nsample
        for j = 1:5
            idx = sample(1:n, n, replace = true)
            fit = try
                fitRegression2(Xfc[:,:,idx], Xdf[idx,:], Y[idx]; murine = murine)
            catch e
                @warn "This bootstrapping set failed at fitRegression"
                nothing
            end
            if fit != nothing
                fitResults[i] = fit
                break
            end
        end
    end
    return fitResults
end


function regressionResult(df; L0, f, murine::Bool)
    Xfc, Xdf, Y = modelPred(df; L0 = L0, f = f, murine = murine)
    res = fitRegression2(Xfc, Xdf, Y; murine = murine)
    loo_res = LOOCrossVal2(Xfc, Xdf, Y; murine = murine)
    btp_res = bootstrap2(Xfc, Xdf, Y; murine = murine)

    odf = df[!, in(["Condition", "Background"]).(names(df))]
    odf[!, "Concentration"] .= ("Concentration" in names(df)) ? (df[!, "Concentration"] .* L0) : L0
    odf[!, "Y"] = Y
    odf[!, "Fitted"] = exponential(regressionPred(Xfc, Xdf, res))
    odf[!, "LOOPredict"] = exponential([regressionPred(Xfc[:, :, ii], Xdf[[ii], :], loo_res[ii])[1] for ii in 1:length(loo_res)])

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
    Xmat, _ = regressionPred(wtXfc, wtXdf, res.cellWs, res.ActI; showXmat=true)
    Xmat = Xmat .* res.cellWs'
    Xmat[!, "Condition"] = wtXdf[!, "Condition"]
    effects = stack(Xmat, Not("Condition"))
    rename!(effects, "value" => "Weight")

    # Assemble bootstrap results
    for bres in btp_res
        Xmatb, _ = regressionPred(wtXfc, wtXdf, bres; showXmat=true)
        Xmatb = Xmatb .* bres.cellWs'
        Xmatb[!, "Condition"] = wtXdf[!, "Condition"]
        effects = innerjoin(effects, stack(Xmatb, Not("Condition")), on = [:Condition, :variable], makeunique=true)
    end
    rename!(effects, "variable" => "Component")

    btp_qtlmat = Matrix(effects[!, Not(["Condition", "Component"])])
    effects = effects[!, ["Condition", "Component", "Weight"]]
    btp_qtl = mapslices(x -> quantile(x, [0.1, 0.5, 0.9]), btp_qtlmat, dims = [2])
    effects[!, :Q10] .= vec(btp_qtl[:, 1])
    effects[!, :Median] .= vec(btp_qtl[:, 2])
    effects[!, :Q90] .= vec(btp_qtl[:, 3])
    effects = effects[effects[!, "Component"] .!= "Neutralization", :]

    # ActI interval
    ActI_btp = hcat([bres.ActI for bres in btp_res]...)
    ActI_qtl = mapslices(x -> quantile(x, [0.1, 0.5, 0.9]), ActI_btp, dims = [2])
    ActI_df = DataFrame(Q10 = ActI_qtl[:, 1], Median = ActI_qtl[:, 2], Q90 = ActI_qtl[:, 3])
    ActI_df[!, "Receptor"] = murine ? murineFcgR : humanFcgR

    return res, odf, effects, ActI_df
end
