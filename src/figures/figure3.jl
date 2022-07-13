""" Figure 3: Explain model, fit data, and validate with Robinett """

function figure3()
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

    df = loadMixData()
    Kav_old = importKavDist(; murine = false, regularKav = true, retdf = true)
    c_noKav = rungMCMC("humanfit_0701_noKav.dat"; dat = :hCHO, Kavd = Kav_old)
    pl_noKav = plotMCMCPredict(c_noKav, df; dat = :hCHO, Kav = Kav_old, R2pos = (0, -2), title = "Predictions with all but affinity fitting")

    c = rungMCMC("humanKavfit_0701.dat"; dat = :hCHO, mcmc_iter = 1_000)

    pl1 = plotMCMCPredict(c, df; dat = :hCHO, title = "All predictions with \nsingle hIgG fitted parameters", R2pos = (0, -2.5))
    pl2 = plotMCMCPredict(
        c,
        df[(df."%_1" .!= 1.0) .& (df."%_2" .!= 1.0), :];
        dat = :hCHO,
        title = "Mixture predictions with \nsingle IgG fitted parameters",
        R2pos = (0, -2.5),
    )
    rob1, rob2 = validateFittedKav(c; murine = false)

    pp = plotGrid((3, 3), [nothing, raw_pred_pl, nothing, pl_noKav, pl1, pl2, rob1, rob2])
    draw(PDF("figure3.pdf", 12inch, 10inch), pp)
end
