""" Figure 2: we can accurately account for mixed ICs """

function figure2()
    raw_predict = predMix(
        averageMixData(loadMixData());
        Kav = importKav(; murine = false, retdf = true),
        Rtot = importRtotDist(:hCHO; regular = true, retdf = true),
    )

    raw_pred_pl = plotPredvsMeasured(
        raw_predict;
        xx = "Value",
        xxlabel = "Measured",
        color = "Receptor",
        shape = "Valency",
        title = "Raw model predictions without fitting",
        R2pos = (0, -2.7),
    )

    c = rungMCMC("humanNUTSfit_0505.dat"; dat = :hCHO, mcmc_iter = 1_000)
    df = loadMixData()

    pl1 = plotMCMCPredict(c, df; dat = :hCHO, title = "All predictions with \nsingle hIgG fitted parameters", R2pos = (0, -2.5))
    pl2 = plotMCMCPredict(
        c,
        df[(df."%_1" .!= 1.0) .& (df."%_2" .!= 1.0), :];
        dat = :hCHO,
        title = "Mixture predictions with \nsingle IgG fitted parameters",
        R2pos = (0, -2.5),
    )
    pl_igg = plotAffinityViolin(c; murine = false, y_range = (5, 8))
    rob1, rob2 = validateFittedKav(c; murine = false)

    pp = plotGrid((4, 3), [nothing, nothing, raw_pred_pl, pl1, pl2, pl_igg[1], pl_igg[2], pl_igg[3], pl_igg[4], rob1, rob2])
    draw(PDF("figure2.pdf", 12inch, 12inch), pp)
end
