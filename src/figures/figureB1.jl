""" This file builds the depletion manuscript, Figure 1. """

""" Plot an example isobologram. """
function plotCellIsobologram(IgGXidx::Int64, IgGYidx::Int64, Cellidx::Int64; L0 = 1e-9, f = 4, murine = true, c1q = false, ex = false)
    Cell = cellTypes[Cellidx]
    Xname = murine ? murineIgGFucose[IgGXidx] : humanIgG[IgGXidx]
    Yname = murine ? murineIgGFucose[IgGYidx] : humanIgG[IgGYidx]
    Kav_df = importKav(; murine = murine, IgG2bFucose = true, c1q = c1q, retdf = true)
    Kav = Matrix{Float64}(Kav_df[!, murine ? murineFcgR : humanFcgR])
    ActI = murine ? murineActI : humanActI
    FcExpr = importRtot(murine = murine)[:, Cellidx]
    if ex
        FcExpr = zeros(length(humanFcgR))
        FcExpr[7] = importRtot(murine = murine)[7, Cellidx]
        ActI = nothing
    end

    output = calculateIsobologram(IgGXidx, IgGYidx, f, L0, FcExpr, Kav, actV = ActI)
    D1 = calculateIsobologram(IgGXidx, IgGYidx, f, L0, FcExpr, Kav, actV = ActI, Mix = false)
    D2 = reverse(calculateIsobologram(IgGYidx, IgGXidx, f, L0, FcExpr, Kav, actV = ActI, Mix = false))

    if ex
        title = "Receptor Binding"
    else
        title = "Activity"
    end

    X = range(0, stop = 1, length = length(output))

    pl = plot(
        layer(x = X, y = D1, Geom.line, Theme(default_color = colorant"blue", line_width = 1px)),
        layer(x = X, y = D2, Geom.line, Theme(default_color = colorant"orange", line_width = 1px)),
        layer(x = X, y = output, Geom.line, Theme(default_color = colorant"green", line_width = 2px)),
        layer(x = X, y = D1 + D2, Geom.line, Theme(default_color = colorant"red", line_width = 3px)),
        Scale.x_continuous(labels = n -> "$Xname $(n*100)%\n$Yname $(100-n*100)%"),
        Guide.xticks(orientation = :horizontal),
        Guide.ylabel("$Cell Predicted $title", orientation = :vertical),
        Guide.manual_color_key("", ["Predicted", "Additive", "$Xname only", "$Yname only"], ["green", "red", "blue", "orange"]),
        Guide.title("$title vs IgG Composition"),
        style(key_position = :inside),
    )
    return pl
end

const receptorNamesB1 =
    Symbol.([
        "IgG1/2a",
        "IgG1/2b",
        "IgG1/3",
        "IgG1/2b-Fucose",
        "IgG2a/2b",
        "IgG2a/3",
        "IgG2a/2b-Fucose",
        "IgG2b/3",
        "IgG2b/2b-Fucose",
        "IgG3/2b-Fucose",
    ])

const humanreceptorNamesB1 = Symbol.(["IgG1/2", "IgG1/3", "IgG1/4", "IgG2/3", "IgG2/4", "IgG3/4"])


"""Figure shows the affect of increasing immune complex concentration on synergies for each IgG combination"""
function PlotSynGraph()
    Kav = importKav(murine = true, IgG2bFucose = true)
    FcgR = importRtot(murine = true)[:, 2] #2 = mean cMO
    IC = exp10.(range(-12, stop = -6, length = 20))
    S = zeros((length(IC), 10))

    for (ii, value) in enumerate(IC)
        A = synergyGrid(4, value, FcgR, Kav)
        h = collect(Iterators.flatten(A))
        S[ii, 1:5] = h[2:6]
        S[ii, 5:8] = h[8:11]
        S[ii, 8:9] = h[14:15]
        S[ii, 10] = h[20]
    end

    S = convert(DataFrame, S)
    rename!(S, receptorNamesB1)
    S = stack(S)

    pl = plot(
        S,
        x = repeat(IC; outer = [10]),
        y = :value,
        Geom.line,
        color = :variable,
        Guide.colorkey(title = "IgG Combination"),
        Scale.x_log10,
        Guide.xlabel("IC Concentration"),
        Guide.ylabel("Synergy"),
        Guide.title("Effect of IC Concentration on Synergy"),
        Guide.xticks(),
    )
    return pl
end

""" Figure shows how immune complex valency affects synergy """
function PlotSynValency()
    Kav = importKav(murine = true, IgG2bFucose = true)
    FcgR = importRtot(murine = true)[:, 2] #2 = mean cMO
    IC = 10e-9
    Valency = range(1, stop = 24)
    S = zeros((length(Valency), 10))

    for (ii, value) in enumerate(Valency)
        A = synergyGrid(value, IC, FcgR, Kav)
        h = collect(Iterators.flatten(A))
        S[ii, 1:5] = h[2:6]
        S[ii, 5:8] = h[8:11]
        S[ii, 8:9] = h[14:15]
        S[ii, 10] = h[20]
    end

    S = convert(DataFrame, S)
    rename!(S, receptorNamesB1)
    S = stack(S)

    pl = plot(
        S,
        x = repeat(Valency; outer = [10]),
        y = :value,
        Geom.point,
        color = :variable,
        Guide.colorkey(title = "IgG Combination"),
        Guide.xlabel("IC Valency"),
        Guide.ylabel("Synergy"),
        Guide.title("Effect of IC Valency on Synergy"),
    )
    return pl
end

function PlotSynvFcrExpr()
    """ Figure shows how Fc receptor expression affects synergy """
    Kav = importKav(murine = true, IgG2bFucose = true)
    FcgR = importRtot(murine = true)[:, 2] #2 = mean cMO
    IC = 10e-9
    multiplier = exp10.(range(-2, stop = 0, length = 50))
    S = zeros((50, 10))

    for (ii, value) in enumerate(multiplier)
        A = synergyGrid(4, IC, (FcgR * value), Kav)
        h = collect(Iterators.flatten(A))
        S[ii, 1:5] = h[2:6]
        S[ii, 5:8] = h[8:11]
        S[ii, 8:9] = h[14:15]
        S[ii, 10] = h[20]
    end

    S = convert(DataFrame, S)
    rename!(S, receptorNamesB1)
    S = stack(S)

    pl = plot(
        S,
        x = repeat(multiplier; outer = [10]),
        y = :value,
        Geom.line,
        Scale.x_log10,
        color = :variable,
        Guide.colorkey(title = "IgG Combination"),
        Guide.xlabel("Fcr Expression"),
        Guide.ylabel("Synergy"),
        Guide.title("Effect of Fc Expression on Synergy"),
    )
    return pl
end


function figureS(Cellidx; L0 = 1e-9, f = 4, murine = true)
    setGadflyTheme()
    p1 = plotCellIsobologram(1, 2, Cellidx; L0 = L0, f = f, murine = murine)
    p2 = plotCellIsobologram(1, 3, Cellidx; L0 = L0, f = f, murine = murine)
    p3 = plotCellIsobologram(1, 4, Cellidx; L0 = L0, f = f, murine = murine)
    p4 = plotCellIsobologram(1, 5, Cellidx; L0 = L0, f = f, murine = murine)
    p5 = plotCellIsobologram(2, 3, Cellidx; L0 = L0, f = f, murine = murine)
    p6 = plotCellIsobologram(2, 4, Cellidx; L0 = L0, f = f, murine = murine)
    p7 = plotCellIsobologram(2, 5, Cellidx; L0 = L0, f = f, murine = murine)
    p8 = plotCellIsobologram(3, 4, Cellidx; L0 = L0, f = f, murine = murine)
    p9 = plotCellIsobologram(3, 5, Cellidx; L0 = L0, f = f, murine = murine)
    p10 = plotCellIsobologram(4, 5, Cellidx; L0 = L0, f = f, murine = murine)

    return p1, p2, p3, p4, p5, p6, p7, p8, p9, p10
end
