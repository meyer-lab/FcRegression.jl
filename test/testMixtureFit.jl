using ForwardDiff

@testset "" begin
    df = FcRegression.loadMixData()
    optDict = Dict(
        "Experiment" => length(unique(df."Experiment")),
        "Valency" => 2,
        "Receptor" => length(FcRegression.measuredRecepExp),
        "KxStar" => 0,
    )
    x0 = FcRegression.fitMixX0(optDict)
    diff = ForwardDiff.gradient(x -> FcRegression.fitMixFunc(x, optDict, df), x0)
    @test eltype(diff) == Float64
end