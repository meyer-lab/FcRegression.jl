function plotActualvFit()
    (odf, effects, btp_std) = CVResults("blood")
    pl = plot(odf[!, :Y], odf[!, :Fitted], seriestype=:scatter, smooth = true, legend = false)
    xlabel!(pl, "Actual effect")
    ylabel!(pl, "Fitted effect")
    title!(pl, "Actual effect vs fitted effect for Blood CD20")
    return pl
end

function plotActualvPredict()
    (odf, effects, btp_std) = CVResults("blood")
    pl = plot(odf[!, :Y], odf[!, :LOOPredict], seriestype=:scatter, smooth = true, legend = false)
    xlabel!(pl, "Actual effect")
    ylabel!(pl, "LOO predicted effect")
    title!(pl, "Actual effect vs LOO predicted for Blood CD20")
    return pl
end

function plotCellTypeEffects()
    dataType = "blood"
    ## blood data has different concentration and can't use this -- is not called for this figure
    (odf, effects, btp_std) = CVResults("blood")
    wtLineNo = odf[!, :Background] .== "wt"
    IgGcategory = odf[wtLineNo, :Condition]
    itemName = [String(i) * "_" * String(c) for c in cellTypes for i in IgGcategory]
    values = effects[wtLineNo, :]
    stdevs = btp_std[wtLineNo, :]

    pl = bar(itemName, vec(values), xrotation=40, yerr = vec(stdevs))
    title!(pl, "Weights of IgGx + celltype in wt for Blood CD20")
    return pl
end

function figureB5()
    p1 = plotActualvFit()
    p2 = plotActualvPredict()
    p = plot(p1, p2, layout = (1, 2), size = (1200, 1200), dpi = 300)
    savefig(p, "figureB5.pdf")
end
