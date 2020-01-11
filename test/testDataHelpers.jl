@testset "dataHelpers.jl tests" begin
    @testset "Testing human FcgR abundance import" begin
        arr_HIV = FcgR.importRtot(murine=false)
        arr_heterozygote = FcgR.importRtot(murine=false, genotype="ZZZ")
        arr_RTF = FcgR.importRtot(murine=false, genotype="RTF")
        arr_HZF = FcgR.importRtot(murine=false, genotype="HZF")

        ## Assert correct assigning of the genotypes (can add reverse cases too)
        @test(arr_HIV[2, :] == arr_RTF[3, :] && arr_HIV[3, :] == arr_RTF[2, :])
        @test(arr_HIV[4, :] == arr_RTF[5, :] && arr_HIV[5, :] == arr_RTF[4, :])
        @test(arr_HIV[7, :] == arr_RTF[8, :] && arr_HIV[8, :] == arr_RTF[7, :])

        ## Check Heterozygote Case
        @test(arr_HIV[2, :] == 2 .* arr_heterozygote[2, :] && arr_HIV[2, :] == 2 .* arr_heterozygote[3, :])
        @test(arr_RTF[3, :] == 2 .* arr_heterozygote[2, :] && arr_RTF[3, :] == 2 .* arr_heterozygote[3, :])
        @test(arr_HIV[4, :] == arr_heterozygote[4, :] .+ arr_heterozygote[5, :])
        @test(arr_RTF[5, :] == arr_heterozygote[4, :] .+ arr_heterozygote[5, :])

        ## Check mixed genotype case
        @test(arr_HZF[2, :] == arr_HIV[2, :])
        @test(arr_HZF[4, :] == arr_HZF[5, :])
        @test(arr_RTF[5, :] == arr_HZF[4, :] .+ arr_HZF[5, :])
    end
end
