"""Calculate the single, mixed drug, and additive responses for one IgG pair"""
function calcSynergy(IgGXidx::Int64, IgGYidx::Int64, L0, f, FcExpr = nothing; fit::Union{optResult, Nothing} = nothing, Rbound = false, nPoints = 100)
    Kav_df = importKav(; murine = true, IgG2bFucose = true, retdf = true)
    Kav = Matrix{Float64}(Kav_df[!, murineFcgR])
    if FcExpr == nothing
        FcExpr = importRtot(; murine = true)
    end

    IgGC = zeros(Float64, size(Kav, 1), nPoints)
    IgGC[IgGYidx, :] .= eps()
    IgGC[IgGXidx, :] .= 1
    D1 = polyfc_ActV(L0, KxConst, f, FcExpr, IgGC, Kav, Rbound; Mix = false)  # size: celltype * FcRecep * nPoints

    IgGC[IgGXidx, :] .= eps()
    IgGC[IgGYidx, :] .= 1
    D2 = polyfc_ActV(L0, KxConst, f, FcExpr, IgGC, Kav, Rbound; Mix = false)  # size: celltype * FcRecep * nPoints

    IgGC[IgGXidx, :] = range(0.0, 1.0; length = nPoints)
    IgGC[IgGYidx, :] = range(1.0, 0.0; length = nPoints)
    combine = polyfc_ActV(L0, KxConst, f, FcExpr, IgGC, Kav, Rbound;)  # size: celltype * nPoints

    if fit != nothing  # using disease model
        if ndims(FcExpr) == 1
            if !Rbound #Activity of single cell or receptor with disease type
                ActI = fit.ActI
                D1 = dropdims(D1, dims = 1)
                D2 = dropdims(D2, dims = 1)
                combine = dropdims(combine, dims = 1)
                D1 = D1' * ActI
                D2 = (D2' * ActI)
                combine = combine' * ActI
            end
            D2 = reverse(D2) #Binding only
            D1[D1 .<= 0.0] .= 0.0
            D2[D2 .<= 0.0] .= 0.0
            additive = D1 + D2
            combine[combine .<= 0.0] .= 0.0
            additive[additive .<= 0.0] .= 0.0
        else #Depletion
            @assert size(combine, 1) == length(fit.cellWs)
            @assert size(D1, 1) == length(fit.cellWs)
            @assert size(D2, 1) == length(fit.cellWs)

            additive = exponential(regressionPred(D1 + reverse(D2; dims = 3), fit))
            combine = exponential(regressionPred(combine, fit))
            D1 = exponential(regressionPred(D1, fit))
            D2 = reverse(exponential(regressionPred(D2, fit)))
        end
    else
        @assert ndims(FcExpr) == 1
        if !Rbound # Activity of single cell without fit activity
            ActI = murineActI
            D1 = dropdims(D1, dims = 1)
            D2 = dropdims(D2, dims = 1)
            combine = dropdims(combine, dims = 1)
            D1 = D1' * ActI
            D2 = (D2' * ActI)
            combine = combine' * ActI
        end
        D2 = reverse(D2) # Binding only
        D1[D1 .<= 0.0] .= 0.0
        D2[D2 .<= 0.0] .= 0.0
        additive = D1 + D2
        combine[combine .<= 0.0] .= 0.0
        additive[additive .<= 0.0] .= 0.0
    end
    return D1, D2, additive, combine
end


"""Calculate the IgG mixture at the point of maximum synergy or antagonism for a pair of IgGs"""
function maxSynergy(IgGXidx::Int64, IgGYidx::Int64, L0, f, FcExpr; fit = nothing, Rbound = false, nPoints = 100)

    D1, D2, additive, output = calcSynergy(IgGXidx, IgGYidx, L0, f, FcExpr; fit = fit, Rbound = Rbound, nPoints = nPoints)
    sampleAxis = range(0, stop = 1, length = length(output))

    # Subtract a line
    diff = output - additive
    maxValue, maxIndex = findmax(abs.(diff))

    return 1 - sampleAxis[maxIndex], sampleAxis[maxIndex], diff[maxIndex]
end


""" Calculate the synergy metric for all pairs of IgG. """
function synergyGrid(L0, f, FcExpr, Kav; fit = nothing, Rbound = false)
    M = zeros(size(Kav)[1], size(Kav)[1])
    nPoints = 100
    for i = 1:size(Kav)[1]
        for j = 1:(i - 1)
            D1, D2, additive, output = calcSynergy(i, j, L0, f, FcExpr; fit = fit, Rbound = Rbound, nPoints = nPoints)
            synergy = sum((output - additive) / nPoints)
            M[i, j] = synergy
        end
    end

    return M
end

"""Calculate EC50"""
function EC50(
    IgGXidx::Int64,
    IgGYidx::Int64,
    L0 = 1e-9,
    f = 4,
    FcExpr = nothing;
    fit = nothing,
    Rbound = true,
    )

    D1, D2, additive, output = calcSynergy(IgGXidx, IgGYidx, L0, f, FcExpr; fit = fit, Rbound = Rbound, nPoints = 100)
    sampleAxis = range(0, stop = 1, length = length(output))

    EC50value = 0.5*maximum(output)
    diff = output .- EC50value
    EC50index = findmin(abs.(diff))[2]
    Xpercent = sampleAxis[EC50index]
    return Xpercent
end

""" Calculate the EC50 for all pairs of IgG """
function EC50Grid(L0, f, FcExpr, Kav, RecepKav; fit = nothing, Rbound = false)
    M = zeros(size(Kav)[1], size(Kav)[1])
    Affinity = zeros(size(Kav)[1], size(Kav)[1])
    Idx = Array{Int64}(undef, size(Kav)[1], size(Kav)[1])
    for i = 1:size(Kav)[1]
        for j = 1:(i - 1)
            xPercent = EC50(i, j, L0, f, FcExpr; fit = fit, Rbound = Rbound)
            if xPercent > 0.5
                EC = xPercent
                Aff = RecepKav[i]
                Idx[i,j] = i
            else
                EC = 1 - xPercent
                Aff = RecepKav[j]
                Idx[i,j] = j
            end
            M[i, j] = EC
            Affinity[i, j] = Aff
        end
    end

    return M, Affinity, Idx
end


function plotDepletionSynergy(
    IgGXidx::Int64,
    IgGYidx::Int64;
    L0 = 1e-9,
    f = 4,
    dataType = nothing,
    fit = nothing,
    Cellidx = nothing,
    Recepidx = nothing,
    Rbound = false,
)

    if IgGXidx > length(murineIgGFucose) || IgGYidx > length(murineIgGFucose)
        IgGXidx, IgGYidx = 1, 2
    end

    if Cellidx == nothing && Recepidx != nothing
        @error "Must specify Cellidx AND Recepidx"
    end
    if dataType != nothing && fit == nothing
        @error "Fit must be provided with dataType"
    end
    @assert IgGXidx != IgGYidx

    Xname = murineIgGFucose[IgGXidx]
    Yname = murineIgGFucose[IgGYidx]
    Receps = murineFcgR
    nPoints = 100
    ymax = nothing

    if Recepidx != nothing # look at only one receptor
        FcExpr = zeros(length(Receps))
        FcExpr[Recepidx] = importRtot(murine = true)[Recepidx, Cellidx]
        ylabel = "Activity"
    elseif Cellidx != nothing # look at only one cell FcExpr
        FcExpr = importRtot(murine = true)[:, Cellidx]
        ylabel = "Activity"
    else
        FcExpr = importRtot(; murine = true)
    end

    if fit != nothing  # use disease model
        if Recepidx != nothing # look at only one receptor
            title = "$(cellTypes[Cellidx]) $murineFcgR[Recepidx] $dataType"
        elseif Cellidx != nothing # look at only one cell FcExpr
            title = "$(cellTypes[Cellidx]) $dataType"
        else
            ylabel = "Depletion"
            title = "$dataType"
            ymax = 1.0
        end
    elseif Recepidx != nothing  # bind to one receptor
        title = "$murineFcgR[Recepidx], $cellTypes[Cellidx]"
    elseif Cellidx != nothing  # bind to one cell type
        title = "$(cellTypes[Cellidx])"
    else
        @error "Not allowed combination of fit/Cellidx/Recepidx."
    end

    D1, D2, additive, output = calcSynergy(IgGXidx, IgGYidx, L0, f, FcExpr; fit = fit, Rbound = Rbound)
    if Rbound
        ylabel = "Binding"
        title = "$title Rbound"
    end
    if ymax == nothing
        ymax = 1.1 * maximum(additive)
    end

    x = range(0.0, 1.0; length = nPoints)
    @assert length(x) == length(D1)
    @assert length(x) == length(D2)
    @assert length(x) == length(additive)
    @assert length(x) == length(output)

    pl = plot(
        layer(x = x, y = D1, Geom.line, Theme(default_color = colorant"blue", line_width = 1px)),
        layer(x = x, y = D2, Geom.line, Theme(default_color = colorant"orange", line_width = 1px)),
        layer(x = x, y = output, Geom.line, Theme(default_color = colorant"green", line_width = 2px)),
        layer(x = x, y = additive, Geom.line, Theme(default_color = colorant"red", line_width = 3px)),
        Scale.x_continuous(labels = n -> "$Xname $(n*100)%\n$Yname $(100-n*100)%"),
        Guide.xticks(orientation = :horizontal),
        Guide.xlabel("L0 = $L0", orientation = :horizontal),
        Guide.ylabel("Predicted $ylabel", orientation = :vertical),
        Coord.cartesian(ymin = 0, ymax = ymax),
        Guide.manual_color_key("", ["Predicted", "Additive", "$Xname only", "$Yname only"], ["green", "red", "blue", "orange"]),
        Guide.title("$Xname-$Yname ($title)"),
        style(key_position = :inside),
    )
    return pl
end
