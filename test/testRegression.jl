function old_regGenData(df; L0, f, murine::Bool, ActI::Union{Vector, Nothing} = nothing)
    df = copy(df)
    FcRecep = murine ? FcRegression.murineFcgR : FcRegression.humanFcgR
    if ActI == nothing
        ActI = murine ? FcRegression.FcRegression.murineActI : FcRegression.humanActI
    else
        @assert length(ActI) == length(murine ? FcRegression.murineActI : FcRegression.humanActI)
    end

    if :Concentration in propertynames(df)
        df[!, :Concentration] .*= L0
    else
        insertcols!(df, 3, :Concentration => L0)
    end

    X = DataFrame(repeat([Float64], length(FcRegression.cellTypes)), FcRegression.cellTypes)
    for i = 1:size(df, 1)
        Kav = convert(Vector{Float64}, df[i, FcRecep])
        Kav = reshape(Kav, 1, :)
        Rtot = FcRegression.importRtot(; murine = murine, genotype = murine ? "NA" : df[i, :Genotype])
        Rmulti = polyfc_ActV(df[i, :Concentration], FcRegression.KxConst, f, Rtot, [1.0], Kav, false)
        item = dropdims(Rmulti, dims = 3) * ActI
        item[item .<= 0.0] .= 0.0
        push!(X, item)
    end

    if :C1q in propertynames(df)
        X[!, :C1q] = df[!, :C1q] .* df[!, :Concentration]
    end
    if :Neutralization in propertynames(df)
        X[!, :Neutralization] = df[!, :Neutralization]
    end

    if :Background in propertynames(df)
        X[df[:, :Background] .== "NeuKO", :Neu] .= 0.0
        X[df[:, :Background] .== "ncMOKO", :ncMO] .= 0.0
    end
    Y = df[!, :Target]

    @assert all(isfinite.(Matrix(X)))
    @assert all(isfinite.(Y))
    return (Matrix{Float64}(X), Y)
end


@testset "Test genRegPred() can take ForwardDiff." begin
    for data in ("ITP", "melanoma", "blood", "bone")
        df = FcRegression.importDepletion(data)
        Xo, Yo = old_regGenData(df; L0 = 1e-10, f = 6, murine = true)

        Xfc, Xdf, Yex = FcRegression.modelPred(df; L0 = 1e-10, f = 6, murine = true)
        extra = Xdf[!, in(["C1q", "Neutralization"]).(names(Xdf))]
        weights = rand(size(Xfc, 1) + size(extra, 2))
        Xmat, _ = FcRegression.regressionPred(Xfc, Xdf, weights, FcRegression.murineActI; showXmat = true)
        @test all(Matrix(Xmat) .≈ Xo)
        @test all(Yex .≈ Yo)

        f = x -> sum(FcRegression.regressionPred(Xfc, Xdf, weights, x))
        @test isa(ForwardDiff.gradient(f, [1.0, 1, -1, 1]), Vector{Float64})
        f = x -> sum(FcRegression.regressionPred(Xfc, Xdf, x, [1.0, 1, -1, 1]))
        @test isa(ForwardDiff.gradient(f, weights), Vector{Float64})
    end
end

@testset "Test regressionResult()" begin
    for data in ("ITP", "blood", "bone", "melanoma", "HIV", "Bcell")
        res, odf, effects, ActI_df = FcRegression.regressionResult(data; L0 = 1e-9, f = 6, murine = true)
        @test "Fitted" in names(odf)
        @test "LOOPredict" in names(odf)
        @test all(effects[!, "Q10"] .<= effects[!, "Median"])
        @test all(effects[!, "Q90"] .>= effects[!, "Median"])
    end
    for data in ("blood", "spleen", "bone")
        df = FcRegression.importHumanized(data)
        res, odf, effects, ActI_df = FcRegression.regressionResult(df; L0 = 1e-9, f = 6, murine = false)
        @test "Fitted" in names(odf)
        @test "LOOPredict" in names(odf)
        @test all(effects[!, "Q10"] .<= effects[!, "Median"])
        @test all(effects[!, "Q90"] .>= effects[!, "Median"])
    end
end
