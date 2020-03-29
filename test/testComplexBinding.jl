import Combinatorics.multinomial

function polyfc_via_polyc(L0::Real, KxStar::Real, f::Number, Rtot::Vector, LigC::Vector, Kav::Matrix)
    LigC /= sum(LigC)
    Cplx = FcgR.vec_comb(length(LigC), f, repeat([f], length(LigC)))'
    Ctheta = exp.(Cplx * log.(LigC)) .* [multinomial(x...) for x in eachrow(Cplx)]
    @assert sum(Ctheta) ≈ 1.0   "Ctheta is $(Ctheta) with sum $(sum(Ctheta)) != 1.0"

    return FcgR.polyc(L0, KxStar, Rtot, Cplx, Ctheta, Kav)
end


@testset "complexBinding.jl tests" begin
    @testset "Give the same results as fcBindingModel" begin
        for i in 1:10
            L0  = rand() * 10.0^rand(-15:-5)
            KxStar = rand() * 10.0^rand(-15:-5)
            f = rand(2:20)
            nl = rand(1:10)
            nr = rand(1:10)
            Rtot = floor.(100 .+ rand(nr) .* (10 .^ rand(4:6, nr)))
            LigC = rand(nl) .* (10 .^ rand(1:2, nl))
            Kav = rand(nl, nr) .* (10 .^ rand(3:7, nl, nr))

            old_res = FcgR.polyfc(L0, KxStar, f, Rtot, LigC, Kav).Lbound
            new_res = polyfc_via_polyc(L0, KxStar, f, Rtot, LigC, Kav)
            @test old_res ≈ new_res
        end
    end

    @testset "complexBinding can take ForwardDiff" begin
        L0 = 1.0e-8
        KxStar = 1.2e-10
        Rtot = [100.0, 1000.0, 10.0, 10000.0]
        Cplx = [1 0 3; 2 2 0; 1 1 2; 4 0 0]
        Ctheta = rand(4)
        Ctheta = Ctheta / sum(Ctheta)
        Kav = rand(3, 4) * 1.0e7
        Ltheta = L0 .* Ctheta

        func = x -> FcgR.polyc(x, KxStar, Rtot, Cplx, Ctheta, Kav)
        out = ForwardDiff.derivative(func, L0)
        @test typeof(out) == Float64

        func = x -> FcgR.polyc(L0, x, Rtot, Cplx, Ctheta, Kav)
        out = ForwardDiff.derivative(func, KxStar)
        @test typeof(out) == Float64

        func = x -> FcgR.polyc(L0, KxStar, x, Cplx, Ctheta, Kav)
        out = ForwardDiff.gradient(func, Rtot)
        @test eltype(out) == Float64
        @test length(out) == length(Rtot)

        func = x -> FcgR.polyc(L0, KxStar, Rtot, Cplx, x, Kav)
        out = ForwardDiff.gradient(func, Ctheta)
        @test eltype(out) == Float64
        @test length(out) == length(Ctheta)

        func = x -> FcgR.polycm(KxStar, Rtot, Cplx, x, Kav)
        out = ForwardDiff.gradient(func, Ltheta)
        @test eltype(out) == Float64
        @test length(out) == length(Ltheta)
    end
end
