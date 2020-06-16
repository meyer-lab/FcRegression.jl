""" This file builds the depletion manuscript, Figure 1. """

""" Plot an example isobologram. """
function plotCellIsobologram(IgGXidx::Int64, IgGYidx::Int64, Cellidx::Int64; L0 = 1e-9, f = 4, murine = true, c1q = false, ex = false)
    CellName = ["ncMO", "cMO", "NKs", "Neu", "EO"]
    Cell = CellName[Cellidx]
    Xname = murine ? murineIgG[IgGXidx] : humanIgG[IgGXidx]
    Yname = murine ? murineIgG[IgGYidx] : humanIgG[IgGYidx]
    Kav_df = importKav(; murine = murine, c1q = c1q, retdf = true)
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
        layer(x = X, y = output, Geom.line, Theme(default_color = colorant"green")),
        layer(x = X, y = D1, Geom.line, Theme(default_color = colorant"blue")),
        layer(x = X, y = D2, Geom.line, Theme(default_color = colorant"yellow")),
        layer(x = X, y = D1 + D2, Geom.line, Theme(default_color = colorant"red")),
        Scale.x_continuous(labels = n -> "$Xname $(n*100)%\n$Yname $(100-n*100)%"),
        Guide.ylabel("$Cell Predicted $title"),
        Guide.manual_color_key("", ["Predicted", "Additive", "$Xname only", "$Yname only"], ["green", "red", "blue", "yellow"]),
        Guide.title("$title vs IgG Composition"),
        Theme(key_position = :inside),
    )
    return pl
end

function receptorNamesB1()
    return [
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
    ]
end


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
    rename!(S, Symbol.(receptorNamesB1()))
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
    rename!(S, Symbol.(receptorNamesB1()))
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
    rename!(S, Symbol.(receptorNamesB1()))
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

function figureB1()
    p1 = plotCellIsobologram(2, 3, 2; L0 = 1e-8, f = 24, murine = false, ex = true)
    p2 = plotCellIsobologram(2, 4, 2)
    p3 = PlotSynGraph()
    p4 = PlotSynValency()
    p5 = PlotSynvFcrExpr()

    draw(SVG("figureB1.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p5 p5]))
end

function figureS(Cellidx; L0 = 1e-9, f = 4, murine = true)
    p1 = plotCellIsobologram(1, 2, Cellidx; L0 = L0, f = f, murine = murine)
    p2 = plotCellIsobologram(1, 3, Cellidx; L0 = L0, f = f, murine = murine)
    p3 = plotCellIsobologram(1, 4, Cellidx; L0 = L0, f = f, murine = murine)
    p4 = plotCellIsobologram(2, 3, Cellidx; L0 = L0, f = f, murine = murine)
    p5 = plotCellIsobologram(2, 4, Cellidx; L0 = L0, f = f, murine = murine)
    p6 = plotCellIsobologram(3, 4, Cellidx; L0 = L0, f = f, murine = murine)

    return p1, p2, p3, p4, p5, p6
end
