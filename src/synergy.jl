"""Calculate the single, mixed drug, and additive responses for one IgG pair"""
function calcSynergy(IgGXidx::Int64, IgGYidx::Int64, L0, f, FcExpr, Kav; fit = nothing, ActI = nothing, c1q = false, nPoints = 100)
    IgGC = zeros(Float64, size(Kav, 1), nPoints)
    IgGC[IgGYidx, :] .= eps()
    IgGC[IgGXidx, :] .= 1
    D1 = polyfc_ActV(L0, KxConst, f, FcExpr, IgGC, Kav, ActI, Mix = false)  # size: celltype * nPoints

    IgGC[IgGXidx, :] .= eps()
    IgGC[IgGYidx, :] .= 1
    D2 = polyfc_ActV(L0, KxConst, f, FcExpr, IgGC, Kav, ActI, Mix = false)  # size: celltype * nPoints

    IgGC[IgGXidx, :] = range(0.0, 1.0; length = nPoints)
    IgGC[IgGYidx, :] = range(1.0, 0.0; length = nPoints)
    output = polyfc_ActV(L0, KxConst, f, FcExpr, IgGC, Kav, ActI)  # size: celltype * nPoints
            
    if c1q
        output = vcat(output, Kav_df[!, :C1q]' * IgGC)
        D1 = vcat(D1, Kav_df[!, :C1q]' * IgGC)
        D2 = vcat(D2, Kav_df[!, :C1q]' * IgGC)
    end
            
    if fit != nothing #using disease model
        @assert size(output, 1) == length(fit.x)
        @assert size(D1, 1) == length(fit.x)
        @assert size(D2, 1) == length(fit.x)
        additive = exponential(Matrix((D1 + reverse(D2; dims = 2))'), fit)
        output = exponential(Matrix(output'), fit)
        D1 = exponential(Matrix(D1'), fit)
        D2 = reverse(exponential(Matrix(D2'), fit))
    else #single cell activity
        D2 = reverse(D2; dims = 2)
        additive = D1 + D2
    end
    return D1, D2, additive, output
end

"""Calculate the IgG mixture at the point of maximum synergy or antagonism for a pair of IgGs"""
function maxSynergy(IgGXidx::Int64, IgGYidx::Int64, L0, f, FcExpr, Kav; fit = nothing, ActI = nothing, c1q = false, nPoints = 100)
    
    D1, D2, additive, output = calcSynergy(IgGXidx, IgGYidx, L0, f, FcExpr, Kav; fit = fit, ActI = ActI, c1q = c1q, nPoints = nPoints)
    sampleAxis = range(0, stop = 1, length = length(output))

    # Subtract a line
    diff = output - additive
    maxValue, maxIndex = findmax(abs.(diff))

    return 1 - sampleAxis[maxIndex], sampleAxis[maxIndex], diff[maxIndex]
end


""" Calculate the synergy metric for all pairs of IgG. """
function synergyGrid(L0, f, FcExpr, Kav; fit = nothing, ActI = nothing, c1q = false)
    M = zeros(size(Kav)[1], size(Kav)[1])
    nPoints = 100

    for i = 1:size(Kav)[1]
        for j = 1:(i - 1)
            D1, D2, additive, output = calcSynergy(i, j, L0, f, FcExpr, Kav; fit = fit, ActI = ActI, c1q = c1q, nPoints = nPoints)
            synergy = sum((output - additive) / nPoints)
            M[i, j] = synergy
        end
    end

    return M
end
