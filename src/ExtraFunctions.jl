"""Calculate EC50"""
function EC50(IgGXidx::Int64, IgGYidx::Int64, L0 = 1e-9, f = 4, FcExpr = nothing; fit = nothing, Rbound = true, murine = true)

    D1, D2, additive, output = calcSynergy(IgGXidx, IgGYidx, L0, f, FcExpr; murine = murine, fit = fit, Rbound = Rbound, nPoints = 100)
    sampleAxis = range(0, stop = 1, length = length(output))

    EC50value = 0.5 * maximum(output)
    diff = output .- EC50value
    EC50index = findmin(abs.(diff))[2]
    Xpercent = sampleAxis[EC50index]

    return Xpercent
end

""" Calculate the EC50 for all pairs of IgG """
function EC50Grid(L0, f, FcExpr, Kav, RecepKav; murine = true, fit = nothing, Rbound = false)
    M = zeros(size(Kav)[1], size(Kav)[1])
    Affinity = zeros(size(Kav)[1], size(Kav)[1])
    Idx = Array{Int64}(undef, size(Kav)[1], size(Kav)[1])
    for i = 1:size(Kav)[1]
        for j = 1:(i - 1)
            xPercent = EC50(i, j, L0, f, FcExpr; murine = murine, fit = fit, Rbound = Rbound)
            if xPercent > 0.5
                EC = xPercent
                Aff = RecepKav[i]
                Idx[i, j] = i
            else
                EC = 1 - xPercent
                Aff = RecepKav[j]
                Idx[i, j] = j
            end
            M[i, j] = EC
            Affinity[i, j] = Aff
        end
    end

    #display(RecepKav)
    return M, Affinity, Idx
end

function plotEC50(L0, f, Cellidx, Recepidx; murine = true, dataType = nothing, fit = nothing, Rbound = false)
    Kav_df = importKav(; murine = true, IgG2bFucose = true, c1q = false, retdf = true)
    Kav = Matrix{Float64}(Kav_df[!, murineFcgR])

    FcExpr = zeros(length(murineFcgR))
    FcExpr[Recepidx] = importRtot(; murine = true)[Recepidx, Cellidx]
    ylabel = "Activity"
    RecepKav = Kav[:, Recepidx]
    if fit == nothing
        title = "ActI not Fit"
    else
        title = "$dataType"
    end
    if Rbound
        title = "$title Rbound"
    end

    M, AffM, Index = EC50Grid(L0, f, FcExpr, Kav, RecepKav; murine = murine, fit = fit, Rbound = Rbound)

    flat = collect(Iterators.flatten(AffM))
    Affinity = zeros(length(receptorNamesB1))
    Affinity[1:4] = flat[2:5]
    Affinity[5:7] = flat[8:10]
    Affinity[8:9] = flat[14:15]
    Affinity[10] = flat[20]

    flat = collect(Iterators.flatten(Index))
    idx = Array{Int64}(undef, length(receptorNamesB1))
    idx[1:4] = flat[2:5]
    idx[5:7] = flat[8:10]
    idx[8:9] = flat[14:15]
    idx[10] = flat[20]
    Ig = Vector(undef, 10)
    for i = 1:10
        Ig[i] = (murineIgGFucose[idx[i]])
    end

    h = collect(Iterators.flatten(M))
    S = zeros(length(receptorNamesB1))
    S[1:4] = h[2:5]
    S[5:7] = h[8:10]
    S[8:9] = h[14:15]
    S[10] = h[20]
    S = DataFrame(Tables.table(S', header = receptorNamesB1))
    S = stack(S)
    #display(S)

    pl = plot(
        S,
        y = :value,
        x = Affinity,
        Geom.point,
        color = :variable,
        shape = Ig,
        Guide.colorkey(),
        Scale.x_log10,
        Guide.xlabel("Kav"),
        Guide.ylabel("Predicted EC50 (% Mixture Concentration"),
        Guide.title("EC50 vs Kav: $(murineCellTypes[Cellidx]) $(murineFcgR[Recepidx]) $dataType"),
        style(point_size = 5px, key_position = :right),
    )
    return pl
end
