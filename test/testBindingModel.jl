@testset "fcBindingModel.jl tests" begin
	L0 = 1.0e-8
	KxStar = 1.0e-10
	f = 16
	Rtot = [100.0, 1000.0, 10.0, 10000.0]
	IgGC = [1.0]
	Kav = ones((1, 4))*1.0e6
	ActI = [1.0, 1.0, -1.0, 0.0]

	@testset "Can successfully assemble the parameters." begin
	    polyfc(L0, KxStar, f, Rtot, IgGC, Kav, ActI)
	    @time polyfc(L0, KxStar, f, Rtot, IgGC, Kav, ActI)
	
	    for i in 1:20
	    	@profile polyfc(L0, KxStar, f, Rtot, IgGC, Kav, ActI)
	    end
	
	    Profile.print(noisefloor=2.0)
	end
	
	@testset "Can use forwardDiff on parameters." begin
		func = x -> polyfc(L0, KxStar, f, x, IgGC, Kav, ActI)["ActV"]
		out = ForwardDiff.gradient(func, Rtot)
		@test eltype(out) == Float64
	
		func = x -> polyfc(L0, KxStar, x, Rtot, IgGC, Kav, ActI)["ActV"]
		out = ForwardDiff.derivative(func, Float64(f))
		@test typeof(out) == Float64
	
		func = x -> polyfc(x, KxStar, f, Rtot, IgGC, Kav, ActI)["ActV"]
		out = ForwardDiff.derivative(func, L0)
		@test typeof(out) == Float64
	
		func = x -> polyfc(L0, x, f, Rtot, IgGC, Kav, ActI)["ActV"]
		out = ForwardDiff.derivative(func, KxStar)
		@test typeof(out) == Float64
	end
end