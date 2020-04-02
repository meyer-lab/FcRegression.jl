using NLsolve

function Req_func(Req::Vector, L0::Real, KxStar::Real, Rtot::Vector, Cplx::AbstractMatrix, Ctheta::Vector, Kav::Matrix)
    Psi = Kav .* transpose(Req) .* KxStar
    PsiRS = sum(Psi, dims = 2) .+ 1.0
    PsiNorm = (Psi ./ PsiRS)
    Rbound = L0 / KxStar * dropdims(sum(Ctheta .* Cplx * PsiNorm .* exp.(Cplx * log1p.(PsiRS .- 1)), dims = 1), dims = 1)
    return Req .+ Rbound .- Rtot
end

function Req_Regression(L0::Real, KxStar::Real, Rtot::Vector, Cplx::AbstractMatrix, Ctheta::Vector, Kav::Matrix)
    f! = (F, x) -> F .= Req_func(x, L0, KxStar, Rtot, Cplx, Ctheta, Kav)
    solve_res = nlsolve(f!, Rtot .* 0.9, method = :newton, autodiff = :forward, iterations = 500)
    @assert solve_res.f_converged == true
    @assert all(solve_res.zero .<= Rtot)
    @assert all(solve_res.zero .> 0.0)
    return solve_res.zero
end

function polyc(L0::Real, KxStar::Real, Rtot::Vector, Cplx::AbstractMatrix, Ctheta::Vector, Kav::Matrix)
    (nl, nr) = size(Kav)
    ncplx = size(Cplx, 1)
    @assert size(Cplx, 2) == nl
    @assert ncplx == length(Ctheta)
    @assert nr == length(Rtot)
    Ctheta /= sum(Ctheta)

    ansType = promote_type(typeof(L0), typeof(KxStar), eltype(Rtot), eltype(Ctheta))
    Rtot = ansType.(Rtot)
    Req = Req_Regression(L0, KxStar, Rtot, Cplx, Ctheta, Kav)

    Psi = Kav .* transpose(Req) .* KxStar
    PsiRS = sum(Psi, dims = 2) .+ 1.0
    Lbound = L0 / KxStar * sum(Ctheta .* expm1.(Cplx * log1p.(PsiRS .- 1)))
    Rbound = L0 / KxStar * sum(Ctheta .*  (Cplx * (1 .- 1 ./PsiRS)) .* exp.(Cplx * log.(PsiRS)))
    return Lbound, Rbound
end

polycm = (KxStar, Rtot, Cplx, Ltheta, Kav) -> polyc(sum(Ltheta), KxStar, Rtot, Cplx, Ltheta ./ sum(Ltheta), Kav)
