function genModelPred(df; L0 = 1e-9, f = 6, murine::Bool)
    df = copy(df)
    FcRecep = murine ? murineFcgR : humanFcgR

    if :Concentration in propertynames(df)
        df[!, :Concentration] .*= L0
    else
        insertcols!(df, 3, :Concentration => L0)
    end

    X = Array{Float64}(undef, length(cellTypes), length(FcRecep), size(df, 1))
    for k = 1:size(X, 3)    # dataframe row
        Kav = convert(Vector{Float64}, df[k, FcRecep])
        Kav = reshape(Kav, 1, :)
        Rtot = importRtot(; murine = murine, genotype = murine ? "NA" : df[k, :Genotype])
        for i = 1:size(X, 1)    # cell type
            Rtotc = Rtot[:, i]
            X[i, :, k] = polyfc(df[k, :Concentration], KxConst, f, Rtotc, [1.0], Kav).Rmulti_n
        end
    end

    Y = df[!, :Target]
    return (X, Y)
end


columnnames = ["x", "y", "z"]
columns = [Symbol(col) => Float64[] for col in columnnames]


function genModelPred2(df; L0 = 1e-9, f = 6, murine::Bool)
    df = copy(df)
    FcRecep = murine ? murineFcgR : humanFcgR

    if :Concentration in propertynames(df)
        df[!, :Concentration] .*= L0
    else
        insertcols!(df, 3, :Concentration => L0)
    end

    colnames = copy(cellTypes)
    columns = [Symbol(col) => Float64[] for col in columnnames]

    X = Array{Float64}(undef, length(cellTypes), length(FcRecep), size(df, 1))

    X = Array{Union{Float64, Vector{Float64}}}(undef, size(df, 1))

    for k = 1:size(X, 3)    # dataframe row
        Kav = convert(Vector{Float64}, df[k, FcRecep])
        Kav = reshape(Kav, 1, :)
        Rtot = importRtot(; murine = murine, genotype = murine ? "NA" : df[k, :Genotype])
        for i = 1:size(X, 1)    # cell type
            Rtotc = Rtot[:, i]
            X[i, :, k] = polyfc(df[k, :Concentration], KxConst, f, Rtotc, [1.0], Kav).Rmulti_n
        end
    end

    Y = df[!, :Target]
    return (X, Y)
end


function genRegPred(X, cellWeights, recepActI)
    ansType = promote_type(eltype(X), eltype(cellWeights), eltype(recepActI))
    @assert length(cellWeights) == size(X, 1)
    @assert length(recepActI) == size(X, 2)
    Ypred = Array{ansType}(undef, size(X, 3))
    for k = 1:size(X, 3)
        recepEff = X[:, :, k] * recepActI
        recepEff[recepEff .< 0.0] .= 0.0
        Ypred[k] = cellWeights' * recepEff      # exponential?
    end
    return Ypred
end
