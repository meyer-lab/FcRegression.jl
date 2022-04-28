import Turing: sample, Prior

@testset "Test the MCMC model can make sane Priors" begin
    df = FcRegression.loadMixData()
    m = FcRegression.sfit(df, df."Value")
    c = sample(m, Prior(), 2)
end

@testset "Test predictMurine() perform correctly" begin
    df = FcRegression.importMurineInVitro()
    mKav = FcRegression.murineKavDist(; regularKav = true)
    mRecepExp = deepcopy(FcRegression.InVitroMurineRcpExp)
    ndf0 = FcRegression.predictMurine(df; Kav = mKav, recepExp = mRecepExp)
    @test all(ndf0."Predict" .> 1e-8)

    # remove a single IgG-FcgR affinity
    mKav[mKav."IgG" .== "IgG2a", "FcgRI"] .= 0.0
    ndf = FcRegression.predictMurine(df; Kav = mKav, recepExp = mRecepExp)
    @test all(ndf[(ndf."Receptor" .== "FcgRI") .& (ndf."Subclass" .== "IgG2c"), "Predict"] .<= 1e-8)
    @test all(ndf[Not((ndf."Receptor" .== "FcgRI") .& (ndf."Subclass" .== "IgG2c")), "Predict"] 
        .== ndf0[Not((ndf0."Receptor" .== "FcgRI") .& (ndf0."Subclass" .== "IgG2c")), "Predict"])

    # remove one IgG's affinities
    mKav = FcRegression.murineKavDist(; regularKav = true)
    mKav[mKav."IgG" .== "IgG2b", Not("IgG")] .= 0.0
    ndf = FcRegression.predictMurine(df; Kav = mKav, recepExp = mRecepExp)
    @test all(ndf[(ndf."Subclass" .== "IgG2b"), "Predict"] .<= 1e-8)
    @test all(ndf[Not(ndf."Subclass" .== "IgG2b"), "Predict"] 
        .== ndf[Not(ndf."Subclass" .== "IgG2b"), "Predict"])

    # remove one Recep's abundance
    mKav = FcRegression.murineKavDist(; regularKav = true)
    mRecepExp["FcgRIIB"] = 0.0
    ndf = FcRegression.predictMurine(df; Kav = mKav, recepExp = mRecepExp)
    @test all(ndf[(ndf."Receptor" .== "FcgRIIB"), "Predict"] .<= 1e-8)
    @test all(ndf[Not(ndf."Receptor" .== "FcgRIIB"), "Predict"] 
        .== ndf[Not(ndf."Receptor" .== "FcgRIIB"), "Predict"])
end