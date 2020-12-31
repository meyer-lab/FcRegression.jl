"""Calculate the single, mixed drug, and additive responses for one IgG pair"""
function calcSynergy(
    IgGXidx::Int64,
    IgGYidx::Int64,
    L0,
    f,
    FcExpr = nothing;
    murine,
    fit::Union{optResult, Nothing} = nothing,
    Rbound = false,
    c1q = false,
    neutralization = false,
    nPoints = 100,
)
    Kav_df = importKav(; murine = murine, IgG2bFucose = murine, c1q = c1q, retdf = true)
    Kav = Matrix{Float64}(Kav_df[!, murine ? murineFcgR : humanFcgR])
    if FcExpr == nothing
        FcExpr = importRtot(; murine = murine)
    end

    IgGC = zeros(Float64, size(Kav, 1), nPoints)
    IgGC[IgGYidx, :] .= eps()
    IgGC[IgGXidx, :] .= 1
    D1 = polyfc_ActV(L0, KxConst, f, FcExpr, IgGC, Kav, Rbound; Mix = false)  # size: celltype * FcRecep * nPoints
    D1df = c1q ? DataFrame(C1q = IgGC' * Kav_df[!, :C1q] .* range(0, stop = 1, length = nPoints) .* L0) : nothing

    IgGC[IgGXidx, :] .= eps()
    IgGC[IgGYidx, :] .= 1
    D2 = polyfc_ActV(L0, KxConst, f, FcExpr, IgGC, Kav, Rbound; Mix = false)  # size: celltype * FcRecep * nPoints
    D2df = c1q ? DataFrame(C1q = IgGC' * Kav_df[!, :C1q] .* range(0, stop = 1, length = nPoints) .* L0) : nothing

    IgGC[IgGXidx, :] = range(0.0, 1.0; length = nPoints)
    IgGC[IgGYidx, :] = range(1.0, 0.0; length = nPoints)
    combine = polyfc_ActV(L0, KxConst, f, FcExpr, IgGC, Kav, Rbound;)  # size: celltype * nPoints
    combinedf = c1q ? DataFrame(C1q = IgGC' * Kav_df[!, :C1q] .* L0) : nothing

    if fit != nothing  # using disease model
        if neutralization
            fit = optResult(fit.cellWs, fit.ActI, fit.residual)
            fit.cellWs = fit.cellWs[1:(end - 1)]
        end
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
            @assert size(combine, 1) + (c1q ? 1 : 0) == length(fit.cellWs)
            @assert size(D1, 1) + (c1q ? 1 : 0) == length(fit.cellWs)
            @assert size(D2, 1) + (c1q ? 1 : 0) == length(fit.cellWs)

            additive = exponential(regressionPred(D1 + reverse(D2; dims = 3), (c1q ? D1df .+ D2df[nPoints:-1:1, :] : nothing), fit))
            combine = exponential(regressionPred(combine, combinedf, fit))
            D1 = exponential(regressionPred(D1, D1df, fit))
            D2 = reverse(exponential(regressionPred(D2, D2df, fit)))
        end
    else
        @assert ndims(FcExpr) == 1
        if !Rbound # Activity of single cell without fit activity
            ActI = murine ? murineActI : humanActI
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
function maxSynergy(IgGXidx::Int64, IgGYidx::Int64, L0, f, FcExpr; fit = nothing, Rbound = false, c1q = false, neutralization = false, nPoints = 100)

    D1, D2, additive, output =
        calcSynergy(IgGXidx, IgGYidx, L0, f, FcExpr; fit = fit, Rbound = Rbound, c1q = c1q, neutralization = neutralization, nPoints = nPoints)
    sampleAxis = range(0, stop = 1, length = length(output))

    # Subtract a line
    diff = output - additive
    maxValue, maxIndex = findmax(abs.(diff))

    return 1 - sampleAxis[maxIndex], sampleAxis[maxIndex], diff[maxIndex]
end


""" Calculate the synergy metric for all pairs of IgG. """
function synergyGrid(L0, f, FcExpr, Kav; murine, fit = nothing, Rbound = false, c1q = false, neutralization = false)
    M = zeros(size(Kav)[1], size(Kav)[1])
    nPoints = 100
    for i = 1:size(Kav)[1]
        for j = 1:(i - 1)
            D1, D2, additive, output = calcSynergy(
                i,
                j,
                L0,
                f,
                FcExpr;
                murine = murine,
                fit = fit,
                Rbound = Rbound,
                c1q = c1q,
                neutralization = neutralization,
                nPoints = nPoints,
            )
            synergy = sum((output - additive) / nPoints)
            M[i, j] = synergy
        end
    end

    return M
end

function plotDepletionSynergy(
    IgGXidx::Int64,
    IgGYidx::Int64;
    L0 = 1e-9,
    f = 4,
    murine = true,
    dataType = nothing,
    fit = nothing,
    Cellidx = nothing,
    c1q = false,
    neutralization = false,
    Recepidx = nothing,
    Rbound = false,
)

    if murine
        if IgGXidx > length(murineIgGFucose) || IgGYidx > length(murineIgGFucose)
            IgGXidx, IgGYidx = 1, 2
        end
    else
        if IgGXidx > length(humanIgG) || IgGYidx > length(humanIgG)
            IgGXidx, IgGYidx = 1, 2
        end
    end
    cellTypes = murine ? murineCellTypes : humanCellTypes

    if Cellidx == nothing && Recepidx != nothing
        @error "Must specify Cellidx AND Recepidx"
    end
    if dataType != nothing && fit == nothing
        @error "Fit must be provided with dataType"
    end
    @assert IgGXidx != IgGYidx

    Xname = murine ? murineIgGFucose[IgGXidx] : humanIgG[IgGXidx]
    Yname = murine ? murineIgGFucose[IgGYidx] : humanIgG[IgGYidx]
    Kav_df = importKav(; murine = murine, IgG2bFucose = murine, c1q = c1q, retdf = true)
    Kav = Matrix{Float64}(Kav_df[!, murine ? murineFcgR : humanFcgR])
    Receps = murine ? murineFcgR : humanFcgR
    nPoints = 100
    ymax = nothing

    if Recepidx != nothing # look at only one receptor
        FcExpr = zeros(length(Receps))
        FcExpr[Recepidx] = importRtot(murine = murine)[Recepidx, Cellidx]
        ylabel = "Activity"
    elseif Cellidx != nothing # look at only one cell FcExpr
        FcExpr = importRtot(murine = murine)[:, Cellidx]
        ylabel = "Activity"
    else
        FcExpr = importRtot(; murine = murine)
    end

    if fit != nothing  # use disease model
        if Recepidx != nothing # look at only one receptor
            title = "$(cellTypes[Cellidx]) $(murine ? murineFcgR[Recepidx] : humanFcgR[Recepidx]) $dataType"
        elseif Cellidx != nothing # look at only one cell FcExpr
            title = "$(cellTypes[Cellidx]) $dataType"
        else
            ylabel = "Depletion"
            title = "$dataType"
            ymax = 1.0
        end
    elseif Recepidx != nothing  # bind to one receptor
        title = "$(murine ? murineFcgR[Recepidx] : humanFcgR[Recepidx]), $(cellTypes[Cellidx])"
    elseif Cellidx != nothing  # bind to one cell type
        title = "$(cellTypes[Cellidx])"
    else
        @error "Not allowed combination of fit/Cellidx/Recepidx."
    end

    D1, D2, additive, output =
        calcSynergy(IgGXidx, IgGYidx, L0, f, FcExpr; murine = murine, fit = fit, Rbound = Rbound, c1q = c1q, neutralization = neutralization)
    if Rbound
        ylabel = "Binding"
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
