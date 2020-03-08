@testset "fcBindingModel.jl tests" begin
    L0 = 1.0e-8
    KxStar = 1.0e-10
    f = 4
    Rtot = [100.0, 1000.0, 10.0, 10000.0]
    IgGC = [1.0]
    Kav = ones((1, 4)) * 1.0e6
    ActI = [1.0, 1.0, -1.0, 0.0]

    @testset "Can successfully assemble the parameters and get a sane result." begin
        out = polyfc(L0, KxStar, f, Rtot, IgGC, Kav, ActI)
        # Mass balance
        @test all(out.Rbound_n .<= Rtot)
        @test all(out.Req .<= Rtot)
        @test all(isapprox.(out.Rbound_n .+ out.Req, Rtot))
    end

    @testset "Can use forwardDiff on parameters." begin
        func = x -> polyfc(L0, KxStar, f, x, IgGC, Kav, ActI).ActV
        out = ForwardDiff.gradient(func, Rtot)
        @test eltype(out) == Float64
        @test length(out) == length(Rtot)

        func = x -> polyfc(L0, KxStar, x, Rtot, IgGC, Kav, ActI).ActV
        out = ForwardDiff.derivative(func, Float64(f))
        @test typeof(out) == Float64
        @test out > 0.0

        func = x -> polyfc(x, KxStar, f, Rtot, IgGC, Kav, ActI).ActV
        out = ForwardDiff.derivative(func, L0)
        @test typeof(out) == Float64

        func = x -> polyfc(L0, x, f, Rtot, IgGC, Kav, ActI).ActV
        out = ForwardDiff.derivative(func, KxStar)
        @test typeof(out) == Float64
    end

    @testset "Test monovalent case." begin
        out = polyfc(L0, KxStar, 1, Rtot, [1], Kav, ActI)

        # Note f is not used
        comp = vec(Kav .* L0 .* Rtot' ./ (1 .+ (Kav .* L0)))

        @test all(out.Lbound .≈ sum(comp))
        @test all(out.Rbound_n .≈ comp)
        @test all(out.Rbound_n .+ out.Req .≈ Rtot)
        @test out.Rmulti == 0.0
        @test out.ActV == 0.0
    end
end
