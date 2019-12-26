using FcgR
using Plots

function plotIsobologram()
    """ Plot an example isobologram. """

    Kav = importKav(murine = false)
    FcExpr = zeros(6)
    FcExpr[5] = 1000.0

    output = calculateIsobologram(2, 3, 24, 1.0e-8, FcExpr, Kav)

    X = range(0, stop = 1, length = length(output))

    plot(X, output, title = "Receptor Binding vs IgG Composition", xticks = false, legend = false, dpi = 72)
    plot!([0, 1], [output[1], output[33]])
    annotate!([(0, 0, text("100% hIgG2", 8, :right, rotation = 45)), (1.0, 0, text("100% hIgG3", 8, :right, rotation = 45))])
    ylabel!("hFcgRIIIA-158V Binding")
    xlabel!("Percent hIgG3")
    ylims!((-1, maximum(output) * 1.1))
end
