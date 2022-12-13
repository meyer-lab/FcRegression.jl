""" Figure 1: show the mixture IC binding data """

""" 
General function to make subplots for every cell types and IgG pairs
splot() is a function that take dataframe with only a single cell type and IgG pair
    and output a plot
"""
function plotMixSubplots(splot::Function, df = loadMixData(); widths = [], kwargs...)
    setGadflyTheme()

    cells = unique(df."Receptor")
    pairs = unique(df[df."subclass_2" .!= "None", ["subclass_1", "subclass_2"]])
    lcells = length(cells)
    lpairs = size(pairs, 1)
    pls = Vector(undef, lcells * lpairs)

    for (i, pairrow) in enumerate(eachrow(pairs))
        for (j, cell) in enumerate(cells)
            IgGXname, IgGYname = pairrow."subclass_1", pairrow."subclass_2"
            ndf = df[(df."Receptor" .== cell) .& (df."subclass_1" .== IgGXname) .& (df."subclass_2" .== IgGYname), :]
            with_legend = i == size(pairs)[1]
            pls[(j - 1) * lpairs + (i - 1) + 1] = splot(ndf; legend = with_legend, y_label = (i == 1), kwargs...)
        end
    end
    return plotGrid((lcells, lpairs), pls; sublabels = false, widths = widths)
end


function figureS1(; figsize = (15inch, 13inch), widths = [3, 3, 3, 3, 3, 3.5], kwargs...)
    setGadflyTheme()
    draw(
        PDF("output/figureS1.pdf", figsize[1], figsize[2]),
        plotMixSubplots(splot_origData, averageMixData(); widths = widths, match_y = false, kwargs...),
    )
end
