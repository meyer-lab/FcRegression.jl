"""Calculate the single, mixed drug, and additive responses for one IgG pair"""
function calcSynergy(IgGXidx::Int64, IgGYidx::Int64, L0, f, FcExpr = nothing;
        murine, fit::Union{optResult, Nothing}, c1q = false, neutralization = false, nPoints = 100)
    Kav_df = importKav(; murine = murine, IgG2bFucose = murine, c1q = c1q, retdf = true)
    Kav = Matrix{Float64}(Kav_df[!, murine ? murineFcgR : humanFcgR])
    if FcExpr == nothing
        FcExpr = FcRegression.importRtot(;murine = murine)
    end

    IgGC = zeros(Float64, size(Kav, 1), nPoints)
    IgGC[IgGYidx, :] .= eps()
    IgGC[IgGXidx, :] .= 1
    D1 = polyfc_ActV(L0, KxConst, f, FcExpr, IgGC, Kav, Mix = false)  # size: celltype * nPoints
    D1df = c1q ? DataFrame(C1q = IgGC' * Kav_df[!, :C1q] .* range(0, stop = 1, length = size(IgGC, 2)) .* L0) : nothing

    IgGC[IgGXidx, :] .= eps()
    IgGC[IgGYidx, :] .= 1
    D2 = polyfc_ActV(L0, KxConst, f, FcExpr, IgGC, Kav, Mix = false)  # size: celltype * nPoints
    D2df = c1q ? DataFrame(C1q = reverse(IgGC' * Kav_df[!, :C1q] .* range(0, stop = 1, length = size(IgGC, 2))) .* L0) : nothing

    IgGC[IgGXidx, :] = range(0.0, 1.0; length = nPoints)
    IgGC[IgGYidx, :] = range(1.0, 0.0; length = nPoints)
    combine = polyfc_ActV(L0, KxConst, f, FcExpr, IgGC, Kav)  # size: celltype * nPoints
    combinedf = c1q ? DataFrame(C1q = IgGC' * Kav_df[!, :C1q] .* L0) : nothing

    if fit != nothing  # using disease model
        if neutralization
            fit = optResult(fit.cellWs, fit.ActI, fit.residual)
            fit.cellWs = fit.cellWs[1:(end-1)]
        end

        @assert size(combine, 1) + (c1q ? 1 : 0) == length(fit.cellWs)
        @assert size(D1, 1) + (c1q ? 1 : 0) == length(fit.cellWs)
        @assert size(D2, 1) + (c1q ? 1 : 0) == length(fit.cellWs)

        additive = exponential(regressionPred(D1 + reverse(D2; dims = 3), (c1q ? D1df.+D2df : nothing), fit))
        combine = exponential(regressionPred(combine, combinedf, fit))
        D1 = exponential(regressionPred(D1, D1df, fit))
        D2 = reverse(exponential(regressionPred(D2, D2df, fit)))
    else  # single cell activity
        ActI = murine ? murineActI : humanActI
        @assert ndims(FcExpr) == 1
        D1 = dropdims(D1, dims=1)
        D1 = D1' * ActI
        D2 = dropdims(D2, dims=1)
        D2 = reverse(D2' * ActI)
        combine = dropdims(combine, dims=1)
        combine = combine' * ActI
        additive = D1 + D2
        D1[D1 .<= 0.0] .= 0.0
        D2[D2 .<= 0.0] .= 0.0
        combine[combine .<= 0.0] .= 0.0
        additive[additive .<= 0.0] .= 0.0
    end
    return D1, D2, additive, combine
end


"""Calculate the IgG mixture at the point of maximum synergy or antagonism for a pair of IgGs"""
function maxSynergy(IgGXidx::Int64, IgGYidx::Int64, L0, f; fit = nothing, c1q = false, neutralization = false, nPoints = 100)

    D1, D2, additive, output = calcSynergy(IgGXidx, IgGYidx, L0, f; fit = fit, c1q = c1q, neutralization = neutralization, nPoints = nPoints)
    sampleAxis = range(0, stop = 1, length = length(output))

    # Subtract a line
    diff = output - additive
    maxValue, maxIndex = findmax(abs.(diff))

    return 1 - sampleAxis[maxIndex], sampleAxis[maxIndex], diff[maxIndex]
end


""" Calculate the synergy metric for all pairs of IgG. """
function synergyGrid(L0, f, FcExpr, Kav; murine, fit = nothing, c1q = false, neutralization = false)
    M = zeros(size(Kav)[1], size(Kav)[1])
    nPoints = 100
    for i = 1:size(Kav)[1]
        for j = 1:(i - 1)
            D1, D2, additive, output = calcSynergy(i, j, L0, f; murine = murine, fit = fit, c1q = c1q, neutralization = neutralization, nPoints = nPoints)
            synergy = sum((output - additive) / nPoints)
            M[i, j] = synergy
        end
    end

    return M
end
