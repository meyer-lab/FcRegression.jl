using NLsolve
using Optim
using LinearAlgebra


function Req_Regression(L0::Real, KxStar::Real, f::Number, Rtot::Vector, IgGC, Kav)
    ansType = promote_type(typeof(L0), typeof(KxStar), typeof(f), eltype(Rtot), eltype(IgGC))

    Av = transpose(Kav) * IgGC * KxStar
    f! = (F, x) -> F .= x + L0 * f / KxStar .* (x .* Av) .* (1 + sum(x .* Av))^(f - 1) - Rtot
    j! = (J, x) -> J[diagind(J)] .= 1 .+ L0 * f / KxStar .* Av * (1 + sum(x .* Av))^(f - 2) .* (1 + sum(x .* Av) .+ (f - 1) * Av)

    x0 = convert(Vector{ansType}, Rtot)
    df = OnceDifferentiable(f!, j!, x0, copy(x0), Diagonal(x0))

    local solve_res
    try
        solve_res = nlsolve(df, x0, method = :newton, iterations = 5000)
        @assert solve_res.f_converged == true
        @assert all(solve_res.zero .<= Rtot)
    catch e
        println("Req solving failed")
        show(Base.@locals)
        println()
        rethrow(e)
    end

    return solve_res.zero
end

function polyfc(L0::Real, KxStar::Real, f::Number, Rtot::Vector, IgGC::Vector, Kav::AbstractMatrix, ActI = nothing)
    # Data consistency check
    (ni, nr) = size(Kav)
    @assert ni == length(IgGC)
    @assert nr == length(Rtot)
    IgGC /= sum(IgGC)

    Req = Req_Regression(L0, KxStar, f, Rtot, IgGC, Kav)

    ansType = promote_type(typeof(L0), typeof(KxStar), typeof(f), eltype(Rtot), eltype(IgGC))
    Phi = ones(ansType, ni, nr + 1) .* IgGC
    Phi[:, 1:nr] .*= Kav .* transpose(Req) .* KxStar
    Phisum = sum(Phi[:, 1:nr])
    Phisum_n = sum(Phi[:, 1:nr], dims = 1)

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
        w["vieq"] = L0 / KxStar .* [binomial(f, i) for i = 0:f] .* Phisum .^ (0:f)
    end

    if ActI != nothing
        ActI = vec(ActI)
        @assert nr == length(ActI)
        w["ActV"] = max(dot(w["Rmulti_n"], ActI), 0.0)
    end
    return w
end

polyfcm = (KxStar, f, Rtot, IgG, Kav, ActI = nothing) -> polyfc(sum(IgG), KxStar, f, Rtot, IgG / sum(IgG), Kav, ActI)

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
    for ict = 1:nct
        for iset = 1:nset
            res[ict, iset] = polyfc(L0, KxStar, f, Rtot[:, ict], IgGC[:, iset], Kav, ActI)["ActV"]
        end
    end
    return res
end
