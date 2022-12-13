import Turing: sample, Prior, MAP
import StatsBase: mean, std
using Optim
using ForwardDiff
using Random

dfs = [FcRegression.loadMixData(), FcRegression.importRobinett(), FcRegression.importMurineLeukocyte(), FcRegression.importMurineInVitro()]
dats = [:hCHO, :hRob, :mLeuk, :mCHO]
murines = [false, false, true, true]

@testset "Test predMix() perform correctly" begin
    for (ii, df) in enumerate(dfs)
        Kav = FcRegression.importKavDist(; murine = murines[ii], regularKav = true, retdf = true)   # no 0 affinity
        Rtot = FcRegression.importRtotDist(dats[ii]; regular = true, retdf = true)
        ndf0 = FcRegression.predMix(df; Kav = Kav, Rtot = Rtot)
        @test all(ndf0."Predict" .>= 0.0)

        if "Subclass" in names(df)
            # remove one IgG's affinities
            rIgG = rand(unique(df."Subclass"))
            Kav1 = deepcopy(Kav)
            Kav1[Kav1."IgG" .== rIgG, Not("IgG")] .= 0.0
            ndf1 = FcRegression.predMix(df; Kav = Kav1, Rtot = Rtot)
            @test all(ndf1[ndf1."Subclass" .== rIgG, "Predict"] .<= 0.0)
            @test all(ndf1[ndf1."Subclass" .!= rIgG, "Predict"] .>= 0.0)
        end

        if "Receptor" in names(df)
            # remove one FcgR's abundance
            rFcgR = rand(unique(df."Receptor"))
            Rtot2 = deepcopy(Rtot)
            Rtot2[rFcgR] = 0.0
            ndf2 = FcRegression.predMix(df; Kav = Kav, Rtot = Rtot2)
            @assert all(ndf2[ndf2."Receptor" .== rFcgR, "Predict"] .<= 0.0)
            @assert all(ndf2[ndf2."Receptor" .!= rFcgR, "Predict"] .>= 0.0)
        end
    end
end

@testset "Test the MCMC model can make sane Priors" begin
    for (ii, df) in enumerate(dfs)
        m = FcRegression.gmodel(df, df."Value"; dat = dats[ii])
        c = sample(m, Prior(), 2)
    end
end

@testset "Building the MCMC model can make sane MAP" begin
    opts = Optim.Options(iterations = 1000, show_every = 10, show_trace = true)
    for (ii, df) in enumerate(dfs)
        m = FcRegression.gmodel(df, df."Value"; dat = dats[ii], Kavd = nothing)
        rng = MersenneTwister(1234)
        m(rng)
        opt = optimize(m, MAP(), LBFGS(; m = 20), opts)
        p = FcRegression.extractMCMC(opt; dat = dats[ii])
        ndf = FcRegression.predMix(df; Kav = p["Kav"], KxStar = p["KxStar"], Rtot = p["Rtot"], fs = [p["f4"], p["f33"]])
        @test FcRegression.R2(ndf."Value", ndf."Predict") > 0.4
    end
end
