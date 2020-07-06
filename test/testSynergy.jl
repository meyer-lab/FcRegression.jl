@testset "synergy.jl tests" begin
    @testset "test calculateIsobologram" begin
        Kav = ones(6, 5) * 1e9
        FcExpr = ones(5) * 1e3
        murine = true
        nPoints = 100

        # All the receptors are the same, so this should be flat
        #output = FcgR.calculateIsobologram(1, 5, 16, 1.2e-9, FcExpr, Kav, nPoints = 33)
        ActI = murine ? murineActI : humanActI
        IgGC[1, :] = range(0.0, 1.0; length = nPoints)
        IgGC[5, :] = range(1.0, 0.0; length = nPoints)
        output = FcgR.polyfc_ActV(1.2e-9, FcgR.KxConst, 16, FcExpr, IgGC, Kav, ActI)

        @test length(output) == 100
        @test all(output .≈ output[1])

        # Activity should be zero if valency is 1
        #output = FcgR.calculateIsobologram(1, 5, 1, 1.2e-9, FcExpr, Kav, actV = ones(length(FcExpr)), nPoints = 33)
        output = FcgR.polyfc_ActV(1.2e-9, FcgR.KxConst, 1, FcExpr, IgGC, Kav, actI = ones(length(FcExpr)))
        @test all(output .≈ 0.0)
    end

    @testset "test synergyGrid" begin
        Kav = ones(6, 5) * 1e9
        FcExpr = ones(5) * 1e3

        #the grid should be 6x6
        grid = FcgR.synergyGrid(1.2e-9, 4, FcExpr, Kav)

        @test all(size(grid) .== size(Kav)[1])
        #@test all(grid .== transpose(grid))
    end
    end
