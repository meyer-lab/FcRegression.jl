""" This file builds the depletion manuscript, Figure 1. """

""" Plot an example isobologram. """
function plotIsobologram()
    Kav = importKav(murine = false)
    # TODO: Should import actual receptor abundance
    FcExpr = zeros(length(humanFcgR))
    FcExpr[7] = importRtot(murine = false)[7, 2]

    output = calculateIsobologram(2, 3, 24, 1.0e-8, FcExpr, Kav)

    X = range(0, stop = 1, length = length(output))

    pl = plot(
        x = X, 
        y = output, 
        Geom.line,
        layer(x = [0, 1], y = [output[1], output[33]], Geom.line),
        Scale.y_continuous(minvalue = -1, maxvalue = maximum(output) * 1.1),
        #Guide.annotation()
        Guide.xlabel("Percent hIgG3"),
        Guide.ylabel("hFcgRIIIA-158V Binding"),
        Guide.title("Receptor Binding vs IgG Composition"),
        Guide.xticks(),
        Guide.annotation(compose(context(), text(0, -1000, "100% hIgG2"), text(.8, -1000, "100% hIgG3"))),
    )
    return pl
end

""" Plot an example isobologram. """
function plotIsobologramTwo()
    Kav = importKav(murine = true, IgG2bFucose = true)
    FcExpr = importRtot(murine = true)[:, 2]

    output = calculateIsobologram(2, 4, 4, 1.0e-9, FcExpr, Kav, actV = murineActI)
    output /= maximum(output)

    X = range(0, stop = 1, length = length(output))

    pl = plot(
        x = X,
        y = output,
        Geom.line,
        layer(x = [0, 1], y = [output[1], output[33]], Geom.line),
        Scale.y_continuous(minvalue = -0.02, maxvalue = maximum(output) * 1.1),
        Guide.xlabel("Percent mIgG2bFucose"),
        Guide.ylabel("cMO Predicted Activity"),
        Guide.title("Activity vs IgG Composition"),
        Guide.xticks(),
        Guide.annotation(compose(context(), text(0, -.1, "100% mIgG2a"), text(.65, -.1, "100% mIgG2bFucose"))),
    )
    return pl
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
    
    names = ["IgG1/2a", "IgG1/2b", "IgG1/3", "IgG1/2b-Fucose", "IgG2a/2b", "IgG2a/3", "IgG2a/2b-Fucose", "IgG2b/3", "IgG2b/2b-Fucose", "IgG3/2b-Fucose"]
    S = convert(DataFrame, S)
    rename!(S, Symbol.(names))
    S = stack(S)

    pl = plot(
        S,
        x = repeat(IC; outer=[10]),
        y = :value,
        Geom.point,
        color = :variable,
        Guide.colorkey(title = "IgG Combination"),
        Scale.x_log10,
        Guide.xlabel("Percent mIgG2bFucose"),
        Guide.ylabel("Synergy"),
        Guide.title("IC Concentration"),
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
    
    names = ["IgG1/2a", "IgG1/2b", "IgG1/3", "IgG1/2b-Fucose", "IgG2a/2b", "IgG2a/3", "IgG2a/2b-Fucose", "IgG2b/3", "IgG2b/2b-Fucose", "IgG3/2b-Fucose"]
    S = convert(DataFrame, S)
    rename!(S, Symbol.(names))
    S = stack(S)

    pl = plot(
        S,
        x = repeat(Valency; outer=[10]),
        y = :value,
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
    IC = 1e-9
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
    
    names = ["IgG1/2a", "IgG1/2b", "IgG1/3", "IgG1/2b-Fucose", "IgG2a/2b", "IgG2a/3", "IgG2a/2b-Fucose", "IgG2b/3", "IgG2b/2b-Fucose", "IgG3/2b-Fucose"]
    S = convert(DataFrame, S)
    rename!(S, Symbol.(names))
    S = stack(S)

    pl = plot(
        S,
        x = repeat(multiplier; outer=[10]),
        y = :value,
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
    p1 = plotIsobologram()
    p2 = plotIsobologramTwo()
    p3 = PlotSynGraph()
    p4 = PlotSynValency()
    p5 = PlotSynvFcrExpr()

    draw(SVG("figureB1.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p5 p5]))
end
