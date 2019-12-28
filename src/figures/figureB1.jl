
""" Plot an example isobologram. """
function plotIsobologram()
    Kav = importKav(murine = false)
    # TODO: Should import actual receptor abundance
    FcExpr = zeros(6)
    FcExpr[5] = 1000.0

    output = calculateIsobologram(2, 3, 24, 1.0e-8, FcExpr, Kav)

    X = range(0, stop = 1, length = length(output))

    pl = plot(X, output, title = "Receptor Binding vs IgG Composition", xticks = false, legend = false, dpi = 72)
    plot!(pl, [0, 1], [output[1], output[33]])
    annotate!(pl, [(0, 0, text("100% hIgG2", 8, :right, rotation = 45)), (1.0, 0, text("100% hIgG3", 8, :right, rotation = 45))])
    ylabel!(pl, "hFcgRIIIA-158V Binding")
    xlabel!(pl, "Percent hIgG3")
    ylims!(pl, (-1, maximum(output) * 1.1))

    return pl
end

""" Plot an example isobologram. """
function plotIsobologramTwo()
    Kav = importKav(murine = true)
    # TODO: Should import actual receptor abundance
    FcExpr = [2571.0, 12886.0, 12563.0, 2371.0]

    output = calculateIsobologram(2, 3, 4, 1.0e-9, FcExpr, Kav, actV = murineActI)
    output /= maximum(output)

    X = range(0, stop = 1, length = length(output))

    pl = plot(X, output, title = "Activity vs IgG Composition", xticks = false, legend = false, dpi = 72)
    plot!(pl, [0, 1], [output[1], output[33]])
    annotate!(pl, [(0, 0, text("100% mIgG2a", 8, :right, rotation = 45)), (1.0, 0, text("100% mIgG2b", 8, :right, rotation = 45))])
    ylabel!(pl, "cMO Predicted Activity")
    xlabel!(pl, "Percent mIgG2b")
    ylims!(pl, (-0.02, maximum(output) * 1.1))

    return pl
end

"""Figure shows the affect of increasing immune complex concentration on synergies for each IgG combination"""
function PlotSynGraph()
    Kav = importKav(murine = true)
    df = importRtot()
    FcgR = df[:, 2] #2 = mean cMO
    IC = exp10.(range(-12, stop = -6, length = 20))
    S = zeros((length(IC), 10))

    for (ii, value) in enumerate(IC)
        A = synergyGrid(4, value, FcgR, Kav)
        h = collect(Iterators.flatten(A))
        S[ii, 1:4] = h[1:4]
        S[ii, 5:7] = h[6:8]
        S[ii, 8:9] = h[11:12]
        S[ii, 10] = h[16]
    end

    pl = plot(
        IC,
        S,
        xaxis = :log,
        title = "Effect of Concentration on Synergy",
        label = ["IgG1/2a" "IgG1/2b" "IgG1/3" "IgG2a/2b" "IgG2a/3" "IgG2b/3"],
        legend = :topleft,
        dpi = 72,
    )
    xlabel!(pl, "IC Concentration")
    ylabel!(pl, "Synergy")

    return pl
end

function figureB1()
    l = @layout [a; b c]

    p1 = plotIsobologram()
    p2 = plotIsobologramTwo()
    p3 = PlotSynGraph()
    p = plot(p1, p2, p3, layout = l)

    savefig(p, "figureB1.pdf")
end
