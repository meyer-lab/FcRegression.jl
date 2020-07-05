""" This file builds the depletion manuscript, Figure 1. """

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
    p1 = plotDepletionSynergy(1, 2; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)
    p2 = plotCellIsobologram(1, 3; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)
    p3 = plotCellIsobologram(1, 4; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)
    p4 = plotCellIsobologram(1, 5; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)
    p5 = plotCellIsobologram(2, 3; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)
    p6 = plotCellIsobologram(2, 4; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)
    p7 = plotCellIsobologram(2, 5; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)
    p8 = plotCellIsobologram(3, 4; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)
    p9 = plotCellIsobologram(3, 5; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)
    p10 = plotCellIsobologram(4, 5; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)

    return p1, p2, p3, p4, p5, p6, p7, p8, p9, p10
end
