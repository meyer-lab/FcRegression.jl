using NLsolve


mutable struct fcOutput{T}
    Lbound::T
    Rbound::T
    Rmulti::T
    ActV::T
    Req::Vector{T}
    Rbound_n::Vector{T}
end


function Req_Regression(L0::Real, KxStar::Real, f::Number, Rtot::Vector, IgGC, Kav)
    ansType = promote_type(typeof(L0), typeof(KxStar), typeof(f), eltype(Rtot), eltype(IgGC))

    Av = transpose(Kav) * IgGC * KxStar
    f! = (F, x) -> F .= x + L0 * f / KxStar .* (x .* Av) .* (1 + sum(x .* Av))^(f - 1) - Rtot

    x0 = convert(Vector{ansType}, Rtot / 1.1)

    local solve_res
    try
        solve_res = nlsolve(f!, x0, method = :newton, autodiff = :forward, iterations = 5000)
        @assert solve_res.f_converged == true
        @assert all(solve_res.zero .<= Rtot .+ 1.0e-12)
        @assert all(-1.0e-12 .<= solve_res.zero)
    catch e
        println("Req solving failed")
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

    w = fcOutput{ansType}(
        L0 / KxStar * ((1 + Phisum)^f - 1),
        L0 / KxStar * f * Phisum * (1 + Phisum)^(f - 1),
        L0 / KxStar * f * Phisum * ((1 + Phisum)^(f - 1) - 1),
        NaN,
        Req,
        vec(L0 / KxStar * f .* Phisum_n * (1 + Phisum)^(f - 1)),
    )

    if ActI != nothing
        ActI = vec(ActI)
        @assert nr == length(ActI)
        Rmulti_n = L0 / KxStar * f .* Phisum_n * ((1 + Phisum)^(f - 1) - 1)
        w.ActV = max(dot(Rmulti_n, ActI), 0.0)
    end
    return w
end

polyfcm = (KxStar, f, Rtot, IgG, Kav, ActI = nothing) -> polyfc(sum(IgG) / f, KxStar, f, Rtot, IgG / sum(IgG), Kav, ActI)

function polyfc_ActV(L0, KxStar, f, Rtot::Array, IgGC::Array, Kav::AbstractMatrix, ActI::Vector)
    """
    Input:
    Rtot: nr * nct matrix, nct = # cell types
    IgGC: ni * nset matrix, nset = # IgGC combinations
    Output:
    Matrix of size nct * nset filled with ActV
    """
    ansType = promote_type(typeof(L0), typeof(KxStar), typeof(f), eltype(Rtot), eltype(IgGC))
    res = Matrix{ansType}(undef, size(Rtot, 2), size(IgGC, 2))
    for ict = 1:size(res, 1)
        for iset = 1:size(res, 2)
            res[ict, iset] = polyfc(L0, KxStar, f, Rtot[:, ict], IgGC[:, iset], Kav, ActI).ActV
        end
    end
    return res
end
