""" Figure 3: fitted murine affinities """

function figure3()
    df = importMurineInVitro()
    ndf = predictMurine(df)
    pl1 = plotPredvsMeasured(ndf; xx = "Value", yy = "Predict", 
        color = "Receptor", shape = "Subclass", 
        clip2one = false, R2pos = (-1.5, 1), 
        title = "Raw murine prediction without fitting")
    pl2 = MAPmurineLikelihood()

    pp = plotGrid((1, 2), [pl1 ,pl2])

    #draw(SVG("figure2.svg", 6inch, 3inch), pp)
    draw(PDF("figure3.pdf", 6inch, 3inch), pp)

end