using Plots

function plotActualvFit(odf, dataType)
    pl = plot(odf[!, :Y], odf[!, :Fitted], seriestype = :scatter, smooth = true, legend = false)
    xlabel!(pl, "Actual effect")
    ylabel!(pl, "Fitted effect")
    title!(pl, "Actual effect vs fitted effect for $dataType")
    return pl
end


function plotActualvPredict(odf, dataType)
    pl = plot(odf[!, :Y], odf[!, :LOOPredict], seriestype = :scatter, smooth = true, legend = false)
    xlabel!(pl, "Actual effect")
    ylabel!(pl, "LOO predicted effect")
    title!(pl, "Actual effect vs LOO predicted for $dataType")
    return pl
end


function plotCellTypeEffects(dataType)
    ## blood data has different concentration and can't use this
    (fit_w, odf, effects, btp_std) = CVResults(dataType)
    wtLineNo = odf[!, :Background] .== "wt"
    IgGcategory = odf[wtLineNo, :Condition]
    itemName = [String(i) * "_" * String(c) for c in cellTypes for i in IgGcategory]
    values = effects[wtLineNo, :]
    stdevs = btp_std[wtLineNo, :]

    pl = bar(itemName, vec(values), xrotation = 40, yerr = vec(stdevs), legend = false)
    ylabel!(pl, "Regressed effect")
    title!(pl, "Weights of IgGx + celltype in wt for $dataType")
    return pl
end


function plotDepletionSynergy(IgGXidx::Int64, IgGYidx::Int64, L0, f, weights::Vector; murine = true, nPoints = 40)

    Xname = murine ? murineIgG[IgGXidx] : humanIgG[IgGXidx]
    Yname = murine ? murineIgG[IgGYidx] : humanIgG[IgGYidx]
    Kav = importKav(; murine = murine)
    FcExpr = importRtot(; murine = murine)
    ActI = murine ? murineActI : humanActI

    IgGC = zeros(Float64, size(Kav, 1), nPoints)
    IgGC[IgGXidx, :] = range(0.0, 1.0; length = nPoints)
    IgGC[IgGYidx, :] = range(1.0, 0.0; length = nPoints)
    @assert size(FcExpr, 2) == length(weights)
    X = polyfc_ActV(L0, KxConst, f, FcExpr, IgGC, Kav, ActI) # size: celltype * nPoints
    output = X' * weights

    pl = plot(IgGC[IgGXidx, :], output, xticks = false, legend = false)
    plot!(pl, [0, 1], [output[1], output[end]])
    annotate!(pl, [(0, minimum(output), text("100% $Xname", 12, :right)), (1.0, minimum(output), text("100% $Yname", 12, :right))])
    title!(pl, "Total predicted effects vs IgG Composition")
    return pl
end


function figureW(dataType)
    (fit_w, odf, effects, btp_std) = CVResults(dataType)
    p1 = plotActualvFit(odf, dataType)
    p2 = plotActualvPredict(odf, dataType)
    p3 = plotCellTypeEffects(dataType)
    p4 = plotDepletionSynergy(2, 3, 1e-9, 4, fit_w)
    p = plot(p1, p2, p3, p4, layout = (2, 2), size = (1200, 1200), dpi = 300)
    savefig(p, "figureW_$(dataType).pdf")
end
