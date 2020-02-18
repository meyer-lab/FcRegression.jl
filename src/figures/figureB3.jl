"""Plot Figure B3, melanoma model"""
function figureB3()
    p1 = plotActualvFit("melanoma")
    p2 = plotActualvPredict("melanoma")
    p3 = plotCellTypeEffects("melanoma")
    p = plot(p1, p2, p3, p3, layout = (2, 2), size = (1200, 1200), dpi = 300)
    savefig(p, "figureB3.pdf")
end
