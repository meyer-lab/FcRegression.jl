"""Plot figure B4, CD 20 Bone model"""
function figureB4()
    p1 = plotActualvFit("bone")
    p2 = plotActualvPredict("bone")
    p3 = plotCellTypeEffects("bone")
    p = plot(p1, p2, p3, p3, layout = (2, 2), size = (1200, 1200), dpi = 300)
    savefig(p, "figureB4.pdf")
end
