
@testset "synergy.jl tests" begin
	@testset "test calculateIsobologram" begin
		Kav = ones(6, 5) * 1e9
		FcExpr = ones(5) * 1e3

		# All the receptors are the same, so this should be flat
		output = calculateIsobologram(1, 5, 16, 1.2e-9, FcExpr, Kav, quantity="Rbound", nPoints=33)

		@test length(output) == 33
		@test all(output .≈ output[1])

		# Activity should be zero if valency is 1
		output = calculateIsobologram(1, 5, 1, 1.2e-9, FcExpr, Kav, actV=ones(length(FcExpr)), nPoints=33)
		@test all(output .≈ 0.0)
	end

	@testset "test calcSynergy" begin
		# A flat curve should have zero synergy
		@test abs(FcgR.calcSynergy(ones(7))) < 1e-14

		# Line should have zero synergy
		@test abs(FcgR.calcSynergy(range(0.0, stop=1.0, length=7))) < 1e-14
	end
    
    @testset "test synergyGrid" begin
        Kav = ones(6, 5) * 1e9
        FcExpr = ones(5) * 1e3
        #the grid should be 5x5
        grid = synergyGrid(4, random.random(), FcExpr, Kav) #not sure if want to keep random functionality
        
    end
end
