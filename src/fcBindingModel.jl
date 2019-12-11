using NLsolve
import LinearAlgebra.diagind
import LinearAlgebra.dot

function Req_func!(F, J, x, L0, f, Rtot, Av, KxStar)
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

function Req_Regression(L0, KxStar, f, Rtot, IgGC, Kav)
    ansType = promote_type(typeof(L0), typeof(KxStar), typeof(f), eltype(Rtot), eltype(IgGC))

    Av = transpose(Kav) * IgGC * KxStar
    fj! = (F, J, x) -> Req_func!(F, J, x, L0, f, Rtot, Av, KxStar)
    solve_res = nlsolve(only_fj!( fj! ), convert(Vector{ansType}, Rtot))

    if solve_res.f_converged == false
        @warn "Req_Regression fails to converge"
    end

    @assert all(solve_res.zero .<= Rtot)
    return solve_res.zero
end

function polyfc(L0, KxStar, f, Rtot::Vector, IgGC::Vector, Kav::AbstractMatrix, ActI=nothing)
    # Data consistency check
    (ni, nr) = size(Kav)
    @assert ni == length(IgGC)
    @assert nr == length(Rtot)
    IgGC /= sum(IgGC)

    Req = Req_Regression(L0, KxStar, f, Rtot, IgGC, Kav)

    ansType = promote_type(typeof(L0), typeof(KxStar), typeof(f), eltype(Rtot), eltype(IgGC))
    Phi = ones(ansType, ni, nr+1) .* IgGC
    Phi[:, 1:nr] .*= Kav .* transpose(Req) .* KxStar
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
        w["ActV"] = max(dot(w["Rmulti_n"], ActI), 0.0)
    end
    return w
end

polyfcm = (KxStar, f, Rtot, IgG, Kav, ActI=nothing) -> polyfc(sum(IgG), KxStar, f, Rtot, IgG/sum(IgG), Kav, ActI)

function polyfc_ActV(L0, KxStar, f, Rtot::Array, IgGC::Array, Kav::AbstractMatrix, ActI::Vector)
    """
    Input:
    Rtot: nr * nct matrix, nct = # cell types
    IgGC: ni * nset matrix, nset = # IgGC combinations
    Output:
    Matrix of size nct * nset filled with ActV
    """
    nct = size(Rtot, 2)
    nset = size(IgGC, 2)
    ansType = promote_type(typeof(L0), typeof(KxStar), typeof(f), eltype(Rtot), eltype(IgGC))
    res = Matrix{ansType}(undef, nct, nset)
    for ict in 1:nct
        for iset in 1:nset
            res[ict, iset] = polyfc(L0, KxStar, f, Rtot[:,ict], IgGC[:,iset], Kav, ActI)["ActV"]
        end
    end
    return res
end
