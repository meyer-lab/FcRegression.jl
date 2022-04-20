""" Figure 3: fitted murine affinities """

function figure3()
    df = importMurineInVitro()
    ndf = predictMurine(df)
    pl = plotPredvsMeasured(ndf; xx = "Value", yy = "Predict", 
        color = "Receptor", shape = "Subclass", 
        clip2one = false, R2pos = (-1.5, 1), 
        title = "Raw murine prediction without fitting")
end