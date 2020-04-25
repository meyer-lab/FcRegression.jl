""" This file builds the depletion manuscript, Figure 1. """

""" Plot an example isobologram. """
function plotIsobologram(IgGXidx::Int64, IgGYidx::Int64; murine = false, c1q = false)
    Xname = murine ? murineIgG[IgGXidx] : humanIgG[IgGXidx]
    Yname = murine ? murineIgG[IgGYidx] : humanIgG[IgGYidx]
    Kav = importKav(murine = murine)
    # TODO: Should import actual receptor abundance
    FcExpr = zeros(length(humanFcgR))
    FcExpr[7] = importRtot(murine = murine)[7, 2]

    output = calculateIsobologram(IgGXidx, IgGYidx, 24, 1.0e-8, FcExpr, Kav)

    X = range(0, stop = 1, length = length(output))

    pl = plot(
        layer(x = X, y = output, Geom.line, Theme(default_color = colorant"green")),
        layer(x = [0, 1], y = [output[1], output[end]], Geom.line, Theme(default_color = colorant"red")),
        Scale.x_continuous(labels = n -> "$Xname $(n*100)%\n$Yname $(100-n*100)%"),
        Scale.y_continuous(minvalue = 0.0, maxvalue = 1.0),
        Guide.xlabel("Percent hIgG3"),
        Guide.ylabel("hFcgRIIIA-158V Binding"),
        Guide.manual_color_key("", ["Predicted", "Linear Addition"], ["green", "red"]),
        Guide.title("Receptor Binding vs IgG Composition"),
        Theme(key_position = :inside),
    )
    return pl
end

""" Plot an example isobologram. """
function plotIsobologramTwo(IgGXidx::Int64, IgGYidx::Int64; murine = true, c1q = false)
    Xname = murine ? murineIgG[IgGXidx] : humanIgG[IgGXidx]
    Yname = murine ? murineIgG[IgGYidx] : humanIgG[IgGYidx]
    Kav = importKav(murine = murine, IgG2bFucose = true)
    FcExpr = importRtot(murine = murine)[:, 2]

    output = calculateIsobologram(IgGXidx, IgGYidx, 4, 1.0e-9, FcExpr, Kav, actV = murineActI)
    output /= maximum(output)

    X = range(0, stop = 1, length = length(output))

    pl = plot(
        layer(x = X, y = output, Geom.line, Theme(default_color = colorant"green")),
        layer(x = [0, 1], y = [output[1], output[end]], Geom.line, Theme(default_color = colorant"red")),
        Scale.x_continuous(labels = n -> "$Xname $(n*100)%\n$Yname $(100-n*100)%"),
        Scale.y_continuous(minvalue = 0.0, maxvalue = 1.0),
        Guide.xlabel("Percent mIgG2bFucose"),
        Guide.ylabel("cMO Predicted Activity"),
        Guide.manual_color_key("", ["Predicted", "Linear Addition"], ["green", "red"]),
        Guide.title("Activity vs IgG Composition"),
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
    p1 = plotIsobologram(2, 3)
    p2 = plotIsobologramTwo(2, 4)
    p3 = PlotSynGraph()
    p4 = PlotSynValency()
    p5 = PlotSynvFcrExpr()

    draw(SVG("figureB1.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p5 p5]))
end
