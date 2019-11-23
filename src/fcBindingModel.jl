using NLsolve
import LinearAlgebra.diagind

function Req_func!(F, J, x, L0::Float64, f, Rtot, Av, KxStar::Float64)
    Phisum = sum(x .* Av)
    if !(F == nothing)
        F .= x + L0 * f / KxStar .* (x .* Av) .* (1 + Phisum)^(f-1) - vec(Rtot)
    end
    if !(J == nothing)
        J .= L0 * f / KxStar * (f-1) * (1 + Phisum)^(f-2)
        J .*= x .* Av
        J .*= transpose(Av)
        J[diagind(J)] .= 1 .+ L0 * f / KxStar .* Av * (1+Phisum)^(f-2) .* (1+Phisum .+ (f-1)*Av)
    end
end

function Req_Regression(L0::Float64, KxStar::Float64, f, Rtot, IgGC, Kav)
    Av = transpose(Kav) * IgGC * KxStar
    fj! = (F, J, x) -> Req_func!(F, J, x, L0, f, Rtot, Av, KxStar)
    solve_res = nlsolve(only_fj!( fj! ), Rtot)

    if solve_res.f_converged == false
        @warn "Req_Regression fails to converge"
    end

    return solve_res.zero
end

function polyfc(L0::Float64, KxStar::Float64, f, Rtot::Vector, IgGC::Vector, Kav::AbstractMatrix, ActI=nothing)
    # Data consistency check
    (ni, nr) = size(Kav)
    @assert ni == length(IgGC)
    @assert nr == length(Rtot)
    IgGC /= sum(IgGC)
    Rtot = convert(Vector{Float64}, Rtot)

    Req = Req_Regression(L0, KxStar, f, Rtot, IgGC, Kav)

    Phi = ones(Float64, ni, nr + 1) .* IgGC
    Phi[:, 1:nr] .*= Kav .* transpose(Rtot) .* KxStar
    Phisum = sum(Phi[:, 1:nr])
    Phisum_n = sum(Phi[:, 1:nr], dims=1)

    w = Dict()
    w["Lbound"] = L0 / KxStar * ((1 + Phisum)^f - 1)
    w["Rbound"] = L0 / KxStar * f * Phisum * (1 + Phisum)^(f - 1)
    w["Rbound_n"] = L0 / KxStar * f .* Phisum_n * (1 + Phisum)^(f - 1)
    w["Rmulti"] = L0 / KxStar * f * Phisum * ((1 + Phisum)^(f - 1) - 1)
    w["Rmulti_n"] = L0 / KxStar * f .* Phisum_n * ((1 + Phisum)^(f - 1) - 1)
    w["nXlink"] = L0 / KxStar * (1 + (1 + Phisum)^(f - 1) * ((f - 1) * Phisum - 1))
    w["Req"] = Req
    w["vtot"] = L0 / KxStar * (1 + Phisum)^f

    if typeof(f) == Int
        w["vieq"] = L0 / KxStar .* [binomial(f,i) for i in 0:f] .* Phisum .^ (0:f)
    end

    if !(ActI == nothing)
        ActI = vec(ActI)
        @assert nr == length(ActI)
        w["ActV"] = max((w["Rmulti_n"] * ActI)[1, 1], 0)
    end
    return w
end

polyfcm = (KxStar, f, Rtot, IgG, Kav, ActI=nothing) -> polyfc(sum(IgG), KxStar, f, Rtot, IgG/sum(IgG), Kav, ActI)
