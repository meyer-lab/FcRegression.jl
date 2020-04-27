@testset "translation.jl tests" begin
    @testset "Test that more refined brute force improves match." begin
        resultTwo = FcgR.brute_force_optimum(1)
        resultFour = FcgR.brute_force_optimum(3)

        @test resultTwo[1] <= resultFour[1]
    end
end
