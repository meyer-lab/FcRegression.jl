using FastGaussQuadrature
using LinearAlgebra

""" Calculates the isobologram between two IgGs under the defined conditions. """
function calculateIsobologram(IgGXidx, IgGYidx, valency, ICconc, FcExpr, Kav; quantity=nothing, actV=nothing, nPoints=33)
    @assert length(FcExpr) == size(Kav, 2)

    if actV != nothing
        quantity = "ActV"
        @assert length(actV) == length(FcExpr)
    elseif quantity == nothing
        quantity = "Lbound"
    end

    IgGYconc = range(0.0, stop=1.0, length=nPoints)
    output = zeros(length(IgGYconc))

    for idx in 1:length(IgGYconc)
        IgGC = zeros(size(Kav, 1))
        IgGC[IgGYidx] = IgGYconc[idx]
        IgGC[IgGXidx] = 1 - IgGYconc[idx]

        output[idx] = fcBindingModel.polyfc(ICconc, KxConst, valency, FcExpr, IgGC, Kav, actV)[quantity]
    end

    return output
end


""" Calculate the synergy index from an isobologram curve. """
function calcSynergy(curve)
	# TODO: This isn't _quite_ right as the edge points don't go to the end
	nodes, weights = gausslegendre(length(curve))

	synergy = dot(weights, curve) / 2.0
	synergy -= (curve[1] + curve[end]) / 2.0

	return synergy
end

""" Calculate the synergy metric for all pairs of IgG. """
function synergyGrid(valency, ICconc, FcExpr, mouse=false, quantity=nothing, actV=nothing)
    if mouse
        Kav = murine_affinities()
    else
        Kav = human_affinities()
    end

    M = zeros(size(Kav)[1],size(Kav)[1])
    for i in 1:size(Kav)[1]
        for j in 1:i
            I = calculateIsobologram(i, j, valency, ICconc, FcExpr, mouse, quantity, actV, nPoints=17)
            M[i,j] = calcSynergy(I)
            M[j,i] = M[i,j]
        end
    end
    return M 
end
