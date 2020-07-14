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
    len = size(X, 2)
    @assert len*2 == length(LK)
    return ssigmoidEff(X, LK[1:len], LK[(len+1):end])
end

function fitSsigmoid(df, murine=true; L0=1e-9, f=4)
    (X, Y) = regGenData(df; L0=L0, f=f, murine=murine, retdf=false)
    len = size(X, 2)
    func = (X, p) -> ssigmoidEff(X, p[1:len], p[(len+1):(2*len)])
    Yi = inv_exponential.(Y)
    fit = curve_fit(func, X, Yi, ones(2*len))
    predicted = exponential.(ssigmoidEff(X, fit.param))
    return fit, Y, predicted
end

function figureSsigmoid()
    for dataType in ["ITP", "blood", "bone", "melanoma", "HIV", "Bcell"]
        fit, Y, predicted = fitSsigmoid(importDepletion(dataType), true)
        pl = plot(x=Y, y=predicted, Geom.point,
            Guide.xlabel("Measured"), Guide.ylabel("Predicted"), Guide.title("Sigmoid Fit for " * dataType))
        draw(SVG("figure_fitSig_M" * dataType * ".svg", 400px, 400px), plotGrid((1, 1), [pl]))
    end
end
