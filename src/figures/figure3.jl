""" Figure 3: fitted murine affinities """

function figure3()
    df = importMurineInVitro()
    ndf = predictMurine(df)
    pl1 = plotPredvsMeasured(ndf; xx = "Value", yy = "Predict", 
        color = "Receptor", shape = "Subclass", 
        clip2one = false, R2pos = (-1.5, 1), 
        title = "Raw murine prediction without fitting")
    pl2 = MAPmurineLikelihood()
    _, pl3 = predictLeukocyte(; average = true, title = "Leukocyte binding raw predictions")

    pp = plotGrid((1, 3), [pl1, pl2, pl3])
    draw(PDF("figure3.pdf", 11inch, 3inch), pp)
end