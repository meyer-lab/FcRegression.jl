""" Figure 3: fitted murine affinities """
function figure3()
    df = importMurineLeukocyte(; average = true)
    Kav_old = importKavDist(; murine = true, regularKav = true, retdf = true)

    ndf1 = predMix(df; Kav = Kav_old, Rtot = nothing, fs = [4, 33])
    pl1 = plotPredvsMeasured(
        ndf1;
        xx = "Value",
        yy = "Predict",
        color = "ImCell",
        shape = "Subclass",
        R2pos = (0, -1.5),
        title = "Raw murine leukocyte prediction\nwith documented affinities",
    )

    c_noKav = rungMCMC("leukNUTSfit_0524_KavOnly_v1.7.dat"; dat = :mLeuk, Kavd = Kav_old)
    pl_noKav = plotMCMCPredict(
        c_noKav,
        df;
        dat = :mLeuk,
        Kav = Kav_old,
        R2pos = (0, -2),
        title = "Murine leukocyte prediction\nwith all but affinity fitting",
    )

    c = rungMCMC("leukNUTSfit_0524_v1.7.dat"; dat = :mLeuk)
    pl2 = plotMCMCPredict(c, df; dat = :mLeuk, Kav = nothing, R2pos = (0, -2), title = "Murine leukocyte prediction\nwith updated affinities")

    apls = plotAffinityViolin(c; murine = true, y_range = (4, 9))
    vpl1, vpl2 = validateFittedKav(c; murine = true)

    pp = plotGrid((3, 3), [pl1, pl_noKav, pl2, apls[1], apls[2], apls[3], vpl1, vpl2])
    draw(PDF("figure3.pdf", 12inch, 9inch), pp)
end
