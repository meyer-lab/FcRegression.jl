@testset "dataHelpers.jl tests" begin
    @testset "Testing human FcgR abundance import" begin
        arr_HIV = FcgR.importRtot(murine=false)
        arr_heterozygote = FcgR.importRtot(murine=false, genotype="ZZZ")
        arr_RTF = FcgR.importRtot(murine=false, genotype="RTF")
        arr_HZF = FcgR.importRtot(murine=false, genotype="HZF")
        
        ## Set up mapping
        receps = ["FcgRI", "FcgRIIA-131H", "FcgRIIA-131R", "FcgRIIB-232I", "FcgRIIB-232T", "FcgRIIC-13N", "FcgRIIIA-158V", "FcgRIIIA-158F", "FcgRIIIB"]
        idxs = 1:9
        d = Dict(receps.=> idxs)

        ## Assert correct assigning of the genotypes (can add reverse cases too)
        @test(arr_HIV[d["FcgRIIA-131H"], :] == arr_RTF[d["FcgRIIA-131R"], :] && arr_HIV[d["FcgRIIA-131R"], :] == arr_RTF[d["FcgRIIA-131H"], :])
        @test(arr_HIV[d["FcgRIIB-232I"], :] == arr_RTF[d["FcgRIIB-232T"], :] && arr_HIV[d["FcgRIIB-232T"], :] == arr_RTF[d["FcgRIIB-232I"], :])
        @test(arr_HIV[d["FcgRIIIA-158V"], :] == arr_RTF[d["FcgRIIIA-158F"], :] && arr_HIV[d["FcgRIIIA-158F"], :] == arr_RTF[d["FcgRIIIA-158V"], :])

        ## Check Heterozygote Case
        @test(arr_RTF[d["FcgRIIA-131R"], :] == arr_heterozygote[d["FcgRIIA-131H"], :] .+ arr_heterozygote[d["FcgRIIA-131R"], :])
        @test(arr_HIV[d["FcgRIIB-232I"], :] == arr_heterozygote[d["FcgRIIB-232I"], :] .+ arr_heterozygote[d["FcgRIIB-232T"], :])
        @test(arr_RTF[d["FcgRIIIA-158F"], :] == arr_heterozygote[d["FcgRIIIA-158V"], :] .+ arr_heterozygote[d["FcgRIIIA-158F"], :])

        ## Check mixed genotype case
        @test(arr_HZF[d["FcgRIIA-131H"], :] == arr_HIV[d["FcgRIIA-131H"], :])
        @test(arr_HZF[d["FcgRIIB-232I"], :] == arr_HZF[d["FcgRIIB-232T"], :])
        @test(arr_RTF[d["FcgRIIB-232T"], :] == arr_HZF[d["FcgRIIB-232I"], :] .+ arr_HZF[d["FcgRIIB-232T"], :])
    end
end
