using LsqFit


shiftsigmoid(X::Array, l::Real, k::Real) = l ./ (1 .+ exp.(-X .+ k))


function ssigmoidEff(X::AbstractMatrix, L::AbstractVector, K::AbstractVector)
    @assert size(X, 2) == length(L)
    @assert length(L) == length(K)
    Xs = copy(X)
    for i in 1:size(X, 2)
        Xs[:, i] .= shiftsigmoid(X[:, i], L[i], K[i])
    end
    return vec(sum(Xs, dims=2))
end

function ssigmoidEff(X::AbstractMatrix, LK::AbstractVector)
    @assert size(X, 2)*2 == length(LK)
    return ssigmoidEff(X, LK[1:size(X, 2)], p[(size(X, 2)+1):end])
end

function fitSsigmoid(df)
    func = (X, p) -> ssigmoidEff(X, p[1:5], p[6:10])
    (X, Y) = regGenData(df; L0=1e-9, f=4, murine=true, retdf=false)
    Y .= inv_exponential.(Y)
    fit = curve_fit(func, X, Y, ones(10))
    return fit
end
