""" Figure 3: fitted murine affinities """
function figure3()
    df = importMurineLeukocyte(; average = true)
    #Kav_old = importKav(; murine = true, retdf = true)
    #Kav_old[Kav_old."IgG" .== "IgG2a", "IgG"] .= "IgG2c";
    #Kav_old = Kav_old[Kav_old."IgG" .!= "IgG3", :];
    #sort!(Kav_old, "IgG");

    Kav_old = importKavDist(; murine = true, regularKav = true, retdf = true)

    ndf1 = predictLeukocyte(df; Kav = Kav_old)
    pl1 = plotPredvsMeasured(
        ndf1;
        xx = "Value",
        yy = "Predict",
        color = "Cell",
        shape = "Subclass",
        R2pos = (0, -1.5),
        title = "Raw murine leukocyte prediction\nwith documented affinities",
    )

    c_noKav = fitLeukocyteMCMC("leukNUTSfit_0511_KavOnly_v1.8.dat"; Kavd = Kav_old)
    pl_noKav = plotMCMCPredict(c_noKav, df; murine = true, CHO = false, Kav = Kav_old,
        R2pos = (0, -2), title = "Murine leukocyte prediction\nwith all but affinity fitting")

    c = fitLeukocyteMCMC("leukNUTSfit_0511_v1.8.dat")
    pl2 = plotMCMCPredict(c, df; murine = true, CHO = false, Kav = nothing,
        R2pos = (0, -2), title = "Murine leukocyte prediction\nwith updated affinities")

    apls = plotAffinityViolin(c; murine = true)
    vpl1, vpl2 = validateMurineInVitro(c)

    pp = plotGrid((3, 3), [pl1, pl_noKav, pl2, apls[1], apls[2], apls[3], vpl1, vpl2])
    draw(PDF("figure3.pdf", 12inch, 9inch), pp)
end
