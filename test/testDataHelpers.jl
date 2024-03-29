@testset "dataHelpers.jl tests" begin
    @testset "Testing deepcopy indeed gives us different objects" begin
        @test FcRegression.importRtot(; murine = true, retdf = true) == FcRegression.importRtot(; murine = true, retdf = true)    # same results
        @test FcRegression.importRtot(; murine = true, retdf = true) !== FcRegression.importRtot(; murine = true, retdf = true)   # different references
        @test FcRegression.importRtot(; murine = false) == FcRegression.importRtot(; murine = false)
        @test FcRegression.importRtot(; murine = false) !== FcRegression.importRtot(; murine = false)
        @test FcRegression.importKav(; murine = true, retdf = true) == FcRegression.importKav(; murine = true, retdf = true)
        @test FcRegression.importKav(; murine = true, retdf = true) !== FcRegression.importKav(; murine = true, retdf = true)
        @test FcRegression.importKav(; murine = false) == FcRegression.importKav(; murine = false)
        @test FcRegression.importKav(; murine = false) !== FcRegression.importKav(; murine = false)
    end

    @testset "Testing importRtot()" begin
        @test FcRegression.importRtot(; murine = true, retdf = false) isa Matrix
        @test FcRegression.importRtot(; murine = false, retdf = false) isa Matrix
        @test FcRegression.importRtot(; murine = true, retdf = true) isa DataFrame
        @test FcRegression.importRtot(; murine = false, retdf = true) isa DataFrame
        @test FcRegression.importKav(; murine = true, retdf = false) isa Matrix
        @test size(FcRegression.importKav(murine = false, retdf = false), 2) == length(FcRegression.humanFcgR)
        @test FcRegression.importKav(; murine = true, retdf = true) isa DataFrame
        @test FcRegression.importKav(; murine = false, retdf = true) isa DataFrame
    end

    @testset "Testing human FcgR abundance import for genotypes" begin
        arr_HIV = FcRegression.importRtot(murine = false, retdf = false)
        arr_heterozygote = FcRegression.importRtot(murine = false, genotype = "ZZZ", retdf = false)
        arr_RTF = FcRegression.importRtot(murine = false, genotype = "RTF", retdf = false)
        arr_HZF = FcRegression.importRtot(murine = false, genotype = "HZF", retdf = false)

        ## Set up mapping
        receps = FcRegression.humanFcgR
        d = Dict(receps .=> (1:length(receps)))

        ## Assert correct assigning of the genotypes (can add reverse cases too)
        @test(arr_HIV[d["FcgRIIA-131H"], :] == arr_RTF[d["FcgRIIA-131R"], :] && arr_HIV[d["FcgRIIA-131R"], :] == arr_RTF[d["FcgRIIA-131H"], :])
        @test(arr_HIV[d["FcgRIIIA-158V"], :] == arr_RTF[d["FcgRIIIA-158F"], :] && arr_HIV[d["FcgRIIIA-158F"], :] == arr_RTF[d["FcgRIIIA-158V"], :])

        ## Check Heterozygote Case
        @test(arr_RTF[d["FcgRIIA-131R"], :] == arr_heterozygote[d["FcgRIIA-131H"], :] .+ arr_heterozygote[d["FcgRIIA-131R"], :])
        @test(arr_RTF[d["FcgRIIIA-158F"], :] == arr_heterozygote[d["FcgRIIIA-158V"], :] .+ arr_heterozygote[d["FcgRIIIA-158F"], :])

        ## Check mixed genotype case
        @test(arr_HZF[d["FcgRIIA-131H"], :] == arr_HIV[d["FcgRIIA-131H"], :])
    end

end
