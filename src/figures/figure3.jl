""" Figure 3: fitted murine affinities """
function figure3()
    df = importMurineLeukocyte(; average = true)
    #Kav_old = importKav(; murine = true, retdf = true)
    #Kav_old[Kav_old."IgG" .== "IgG2a", "IgG"] .= "IgG2c";
    #Kav_old = Kav_old[Kav_old."IgG" .!= "IgG3", :];
    #sort!(Kav_old, "IgG");

    Kav_old = murineKavDist(; regularKav = true, retdf = true)

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

    c_noKav = fitLeukocyteMCMC(nothing; Kavd = Kav_old)
    ndf_noKav = predictLeukocyte(c_noKav, df; Kavd = Kav_old)
    pl_noKav = plotPredvsMeasured(
        ndf_noKav;
        xx = "Value",
        yy = "Predict",
        color = "Cell",
        shape = "Subclass",
        R2pos = (0, -2),
        title = "Murine leukocyte prediction\nwith all but affinity fitting",
    )

    c = fitLeukocyteMCMC("leukNUTSfit_0509.dat")
    ndf2 = predictLeukocyte(c, df)
    pl2 = plotPredvsMeasured(
        ndf2;
        xx = "Value",
        yy = "Predict",
        color = "Cell",
        shape = "Subclass",
        R2pos = (0, -2),
        title = "Murine leukocyte prediction\nwith updated affinities",
    )

    apls = plotAffinityViolin(c; murine = true)
    vpl1, vpl2 = validateMurineInVitro(c)

    pp = plotGrid((3, 3), [pl1, pl_noKav, pl2, apls[1], apls[2], apls[3], vpl1, vpl2])
    draw(PDF("figure3.pdf", 12inch, 9inch), pp)
end
