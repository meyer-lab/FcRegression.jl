function plotActualvFit()
    (odf, effects, btp_std) = CVResults("melanoma")
    pl = plot(odf[!, :Y], odf[!, :Fitted], seriestype=:scatter, smooth = true, legend = false)
    xlabel!(pl, "Actual effect")
    ylabel!(pl, "Fitted effect")
    title!(pl, "Actual effect vs fitted effect for melanoma")
    return pl
end

function plotActualvPredict()
    (odf, effects, btp_std) = CVResults("melanoma")
    pl = plot(odf[!, :Y], odf[!, :LOOPredict], seriestype=:scatter, smooth = true, legend = false)
    xlabel!(pl, "Actual effect")
    ylabel!(pl, "LOO predicted effect")
    title!(pl, "Actual effect vs LOO predicted for melanoma")
    return pl
end

function plotCellTypeEffects()
    dataType = "melanoma"
    ## blood data has different concentration and can't use this
    (odf, effects, btp_std) = CVResults("melanoma")
    wtLineNo = odf[!, :Background] .== "wt"
    IgGcategory = odf[wtLineNo, :Condition]
    itemName = [String(i) * "_" * String(c) for c in cellTypes for i in IgGcategory]
    values = effects[wtLineNo, :]
    stdevs = btp_std[wtLineNo, :]

    pl = bar(itemName, vec(values), xrotation=40, yerr = vec(stdevs))
    title!(pl, "Weights of IgGx + celltype in wt for melanoma")
    return pl
end

function figureB3()
    p1 = plotActualvFit()
    p2 = plotActualvPredict()
    p3 = plotCellTypeEffects()
    p = plot(p1, p2, p3, p3, layout = (2, 2), size = (1200, 1200), dpi = 300)
    savefig(p, "figureB3.pdf")
end
