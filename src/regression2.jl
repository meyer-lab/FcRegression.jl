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

    Xmat = Array{ansType}(undef, size(Xfc, 3), length(cellWeights))
    Ypred = Array{ansType}(undef, size(Xfc, 3))
    for k = 1:size(Xfc, 3)
        recepEff = Xfc[:, :, k] * recepActI
        recepEff[recepEff .< 0.0] .= 0.0
        if size(extra, 2) > 0
            append!(recepEff, extra[k, :])
        end
        Xmat[k, :] .= recepEff
        Ypred[k] = cellWeights' * recepEff      # exponential?
    end
    if showXmat
        return Xmat, Ypred
    else
        return Ypred
    end
end


"""
function fitRegression2(df; L0, f, murine::Bool=true)
    Xfc, Xdf, Y = modelPred(df; L0 = L0, f = f, murine::Bool = murine)
    cellWlen = size(Xfc, 1) + size(Xdf[!, in(["C1q", "Neutralization"]).(names(Xdf))], 2)
    f = x -> regressionPred(Xfc, Xdf, x[1:cellWlen], x[(cellWlen+1):end])


    optimize(f, lower, upper, method)

end
"""
