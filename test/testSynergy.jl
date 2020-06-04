@testset "synergy.jl tests" begin
    @testset "test calculateIsobologram" begin
        Kav = ones(6, 5) * 1e9
        FcExpr = ones(5) * 1e3

        # All the receptors are the same, so this should be flat
        output = FcgR.calculateIsobologram(1, 5, 16, 1.2e-9, FcExpr, Kav, quantity = :Rbound, nPoints = 33)

        @test length(output) == 33
        @test all(output .≈ output[1])

        # Activity should be zero if valency is 1
        output = FcgR.calculateIsobologram(1, 5, 1, 1.2e-9, FcExpr, Kav, actV = ones(length(FcExpr)), nPoints = 33)
        @test all(output .≈ 0.0)
    end

    @testset "test synergyGrid" begin
        Kav = ones(6, 5) * 1e9
        FcExpr = ones(5) * 1e3

        #the grid should be 6x6
        grid = FcgR.synergyGrid(4, 1.2e-9, FcExpr, Kav)

        @test all(size(grid) .== size(Kav)[1])
        @test all(grid .== transpose(grid))
    end
end
