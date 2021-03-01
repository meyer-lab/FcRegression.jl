"""Calculate EC50"""
function EC50(IgGXidx::Int64, IgGYidx::Int64, L0 = 1e-9, f = 4, FcExpr = nothing; fit = nothing, Rbound = true, murine = true)

    D1, D2, additive, output = calcSynergy(IgGXidx, IgGYidx, L0, f, FcExpr; murine = murine, fit = fit, Rbound = Rbound, nPoints = 100)
    sampleAxis = range(0, stop = 1, length = length(output))

    EC50value = 0.5 * (maximum(output) - minimum(output))
    diff = output .- EC50value
    Value, EC50index = findmin(abs.(diff))
    Xpercent = sampleAxis[EC50index]

    return Xpercent, Value
end

""" Calculate the EC50 for all pairs of IgG """
function EC50Grid(L0, f, FcExpr, Kav, RecepKav; murine = true, fit = nothing, Rbound = true)
    M = zeros(size(Kav)[1], size(Kav)[1])
    Affinity = zeros(size(Kav)[1], size(Kav)[1])
    Idx = Array{Int64}(undef, size(Kav)[1], size(Kav)[1])
    for i = 1:size(Kav)[1]
        for j = 1:(i - 1)
            xPercent, binding = EC50(i, j, L0, f, FcExpr; murine = murine, fit = fit, Rbound = Rbound)
            if xPercent > 0.5
                Aff = RecepKav[i]
                Idx[i, j] = i
            else
                Aff = RecepKav[j]
                Idx[i, j] = j
            end
            M[i, j] = binding
            Affinity[i, j] = Aff
        end
    end
    return M, Affinity, Idx
end

""" Calculate the synergy metric for all pairs of IgG. """
function maxSynergyGrid(L0, f, FcExpr, Kav; murine, fit = nothing, Rbound = false, c1q = false, neutralization = false)
    M = zeros(size(Kav)[1], size(Kav)[1])
    nPoints = 100
    for i = 1:size(Kav)[1]
        for j = 1:(i - 1)
            
            D1, D2, additive, output = calcSynergy(i, j, L0, f, FcExpr; murine = murine, fit = fit, Rbound = Rbound, c1q = c1q, neutralization = neutralization, nPoints = nPoints,)
            sampleAxis = range(0, stop = 1, length = length(output))

            # Subtract a line
            diff = output - additive
            maxValue, maxIndex = findmax(abs.(diff))

            M[i, j] = maxIndex
        end
    end

    return M
end

function plotEC50(L0, f, Cellidx, Recepidx; murine = true, dataType = nothing, fit = nothing, Rbound = true)
    Kav_df = importKav(; murine = murine, IgG2bFucose = murine, c1q = false, retdf = true)
    Kav = Matrix{Float64}(Kav_df[!, murineFcgR])

    FcExpr = zeros(length(murineFcgR))
    FcExpr[Recepidx] = importRtot(; murine = murine)[Recepidx, Cellidx]
    ylabel = "Activity"
    RecepKav = Kav[:, Recepidx]
    if fit === nothing
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
        Guide.ylabel("Predicted EC50 Rbound"),
        Guide.title("EC50 vs Kav: $(murineCellTypes[Cellidx]) $(murineFcgR[Recepidx]) $dataType"),
        style(point_size = 5px, key_position = :right),
    )
    return pl
end

function maxSynergyBar(L0, f; Cellidx = nothing, Recepidx = nothing, murine = true, dataType = nothing, fit = nothing, Rbound = false)

    Receps = murine ? murineFcgR : humanFcgR
    Kav_df = importKav(; murine = murine, IgG2bFucose = murine, retdf = true)
    Kav = Matrix{Float64}(Kav_df[!, Receps])

    if Rbound
        title = "$title Rbound"
    end

    if Recepidx !== nothing # look at only one receptor
        FcExpr = zeros(length(Receps))
        FcExpr[Recepidx] = importRtot(murine = murine)[Recepidx, Cellidx]
        ylabel = "Activity"
    elseif Cellidx !== nothing # look at only one cell FcExpr
        FcExpr = importRtot(murine = murine)[:, Cellidx]
        ylabel = "Activity"
    else
        FcExpr = importRtot(; murine = murine)
    end

    if fit !== nothing  # use disease model
        if Recepidx !== nothing # look at only one receptor
            title = "$(cellTypes[Cellidx]) $(murine ? murineFcgR[Recepidx] : humanFcgR[Recepidx]) $dataType"
        elseif Cellidx !== nothing # look at only one cell FcExpr
            title = "$(cellTypes[Cellidx]) $dataType"
        else
            ylabel = "Depletion"
            title = "$dataType"
            ymax = 1.0
        end
    elseif Recepidx !== nothing  # bind to one receptor
        title = "$(murine ? murineFcgR[Recepidx] : humanFcgR[Recepidx]), $(cellTypes[Cellidx])"
    elseif Cellidx !== nothing  # bind to one cell type
        title = "$(cellTypes[Cellidx])"
    else
        @error "Not allowed combination of fit/Cellidx/Recepidx."
    end

    M = maxSynergyGrid(L0, f, FcExpr, Kav; murine = murine, fit = fit, Rbound = Rbound)

    h = collect(Iterators.flatten(M))
    if murine
        S = zeros(length(receptorNamesB1))
        S[1:4] = h[2:5]
        S[5:7] = h[8:10]
        S[8:9] = h[14:15]
        S[10] = h[20]
        S = convert(DataFrame, S')
        rename!(S, receptorNamesB1)
    else
        S = zeros(length(humanreceptorNamesB1))
        S[1:3] = h[2:4]
        S[4:5] = h[7:8]
        S[6] = h[12]
        S = convert(DataFrame, S')
        rename!(S, humanreceptorNamesB1)
    end

    S = stack(S)

    pl = plot(
        S,
        y = :value,
        x = :variable,
        color = :variable,
        Geom.bar(position = :dodge),
        style(key_position = :none),
        Guide.xlabel("Mixture", orientation = :vertical),
        Guide.xlabel("Synergy", orientation = :horizontal),
        Guide.title("Synergy vs Mixture ($title)"),
    )
    return pl
end
