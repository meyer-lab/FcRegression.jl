using QuadGK

""" Calculates the isobologram between two IgGs under the defined conditions. """
function calculateIsobologramPoint(
    pointt::Float64,
    IgGXidx::Int64,
    IgGYidx::Int64,
    valency,
    ICconc::Float64,
    FcExpr,
    Kav;
    quantity = nothing,
    actV = nothing,
    Mix = true,
)
    @assert length(FcExpr) == size(Kav, 2)

    if actV != nothing
        quantity = :ActV
        @assert length(actV) == length(FcExpr)
    elseif quantity == nothing
        quantity = :Lbound
    end

    @assert 0.0 <= pointt <= 1.0

    IgGC = zeros(size(Kav, 1))
    IgGC[IgGXidx] += pointt
    if Mix
        IgGC[IgGYidx] += 1.0 - pointt
    else
        IgGC[IgGYidx] += eps()
    end

    w = polyfc(ICconc, KxConst, valency, FcExpr, IgGC, Kav, actV)

    return getproperty(w, quantity)
end


function calculateIsobologram(IgGXidx::Int64, IgGYidx::Int64, valency, ICconc::Float64, FcExpr, Kav; quantity = nothing, actV = nothing, Mix = true, nPoints = 100)
    IgGXconc = range(0.0, stop = 1.0, length = nPoints)
    output = zeros(length(IgGXconc))

    for ii = 1:length(IgGXconc)
        output[ii] = calculateIsobologramPoint(IgGXconc[ii], IgGXidx, IgGYidx, valency, ICconc, FcExpr, Kav; quantity = quantity, actV = actV, Mix = Mix)
    end

    return output
end


""" Calculate the synergy index from an isobologram curve. """
function calcSynergy(IgGXidx::Int64, IgGYidx::Int64, valency, ICconc::Float64, FcExpr, Kav; quantity = nothing, actV = nothing)

    function calcFunc(xx)
        return calculateIsobologramPoint(xx, IgGXidx, IgGYidx, valency, ICconc, FcExpr, Kav; quantity = quantity, actV = actV)
    end

    additive = (calcFunc(0.0) + calcFunc(1.0)) / 2.0

    synergy, err = quadgk(calcFunc, 0.0, 1.0, rtol = 1e-8)

    return (synergy - additive) / additive
end


"""Calculate the IgG mixture at the point of maximum synergy or antagonism for a pair of IgGs"""
function maxSynergy(IgGXidx, IgGYidx, valency, ICconc, FcExpr, Kav; quantity = nothing, actV = nothing)
    curve = calculateIsobologram(IgGXidx, IgGYidx, valency, ICconc, FcExpr, Kav, quantity, actV)
    sampleAxis = range(0, stop = 1, length = length(curve))

    # Subtract a line
    diff = curve - range(curve[1], stop = curve[end], length = length(curve))
    maxValue, maxIndex = findmax(abs.(diff))

    return 1 - sampleAxis[maxIndex], sampleAxis[maxIndex], diff[maxIndex]
end


""" Calculate the synergy metric for all pairs of IgG. """
function synergyGrid(valency, ICconc, FcExpr, Kav; quantity = nothing, actV = nothing)
    M = zeros(size(Kav)[1], size(Kav)[1])

    for i = 1:size(Kav)[1]
        for j = 1:(i - 1)
            M[i, j] = calcSynergy(i, j, valency, ICconc, FcExpr, Kav; quantity = quantity, actV = actV)
        end

        M[:, i] = M[i, :]
    end

    return M
end
