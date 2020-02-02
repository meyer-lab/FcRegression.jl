""" This file builds the depletion manuscript, Figure 1. """

""" Plot an example isobologram. """
function plotIsobologram()
    Kav = importKav(murine = false)
    # TODO: Should import actual receptor abundance
    FcExpr = zeros(length(humanFcgR))
    FcExpr[7] = importRtot(murine = false)[7, 2]

    output = calculateIsobologram(2, 3, 24, 1.0e-8, FcExpr, Kav)

    X = range(0, stop = 1, length = length(output))

    pl = plot(X, output, title = "Receptor Binding vs IgG Composition", xticks = false, legend = false)
    plot!(pl, [0, 1], [output[1], output[33]])
    annotate!(pl, [(0, 0, text("100% hIgG2", 8, :right, rotation = 45)), (1.0, 0, text("100% hIgG3", 8, :right, rotation = 45))])
    ylabel!(pl, "hFcgRIIIA-158V Binding")
    xlabel!(pl, "Percent hIgG3")
    ylims!(pl, (-1, maximum(output) * 1.1))

    return pl
end

""" Plot an example isobologram. """
function plotIsobologramTwo()
    Kav = importKav(murine = true, IgG2bFucose = true)
    # TODO: Should import actual receptor abundance
    FcExpr = importRtot(murine = true)[:, 2]

    output = calculateIsobologram(2, 4, 4, 1.0e-9, FcExpr, Kav, actV = murineActI)
    output /= maximum(output)

    X = range(0, stop = 1, length = length(output))

    pl = plot(X, output, title = "Activity vs IgG Composition", xticks = false, legend = false)
    plot!(pl, [0, 1], [output[1], output[33]])
    annotate!(pl, [(0, 0, text("100% mIgG2a", 8, :right, rotation = 45)), (1.0, 0, text("100% mIgG2bFucose", 8, :right, rotation = 45))])
    ylabel!(pl, "cMO Predicted Activity")
    xlabel!(pl, "Percent mIgG2bFucose")
    ylims!(pl, (-0.02, maximum(output) * 1.1))

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

    pl = plot(
        IC,
        S,
        xaxis = :log,
        title = "Effect of Concentration on Synergy",
        label = ["IgG1/2a" "IgG1/2b" "IgG1/3" "IgG1/2b-Fucose" "IgG2a/2b" "IgG2a/3" "IgG2a/2b-Fucose" "IgG2b/3" "IgG2b/2b-Fucose" "IgG3/2b-Fucose"],
        legend = :topleft,
    )
    xlabel!(pl, "IC Concentration")
    ylabel!(pl, "Synergy")

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

    pl = plot(
        Valency,
        S,
        title = "Effect of IC Valency on Synergy",
        label = ["IgG1/2a" "IgG1/2b" "IgG1/3" "IgG1/2b-Fucose" "IgG2a/2b" "IgG2a/3" "IgG2a/2b-Fucose" "IgG2b/3" "IgG2b/2b-Fucose" "IgG3/2b-Fucose"],
        legend = :topleft,
    )
    xlabel!(pl, "IC Valency")
    ylabel!(pl, "Synergy")

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

    pl = plot(
        multiplier,
        S,
        xaxis = :log,
        title = "Effect of Fc Expression on Synergy",
        label = ["IgG1/2a" "IgG1/2b" "IgG1/3" "IgG1/2b-Fucose" "IgG2a/2b" "IgG2a/3" "IgG2a/2b-Fucose" "IgG2b/3" "IgG2b/2b-Fucose" "IgG3/2b-Fucose"],
        legend = :topleft,
    )
    xlabel!(pl, "Fcr Expression")
    ylabel!(pl, "Synergy")

    return pl
end

function figureB1()
    p1 = plotIsobologram()
    p2 = plotIsobologramTwo()
    p3 = PlotSynGraph()
    p4 = PlotSynValency()
    p5 = PlotSynvFcrExpr()
    p = plot(p1, p1, p2, p3, p4, p5, layout = (3, 2), size = (2100, 1200), dpi = 300)

    savefig(p, "figureB1.pdf")
end
