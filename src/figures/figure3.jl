""" Figure 3: fitted murine affinities """

function plot_murine_MCMC_affinity(c::Chains = runMurineMCMC())
    Kav_priors = murineKavDist()
    pls = Vector{Union{Gadfly.Plot, Context}}(undef, 3)
    Kav_posts = extractLeukocyteKav(c; retdf = true)
    for (ii, igg) in enumerate(Kav_priors[!, "IgG"])
        priors = reshape(Matrix(Kav_priors[Kav_priors."IgG" .== igg, Not("IgG")]), :)
        posts = DataFrame(hcat([reshape(Kav_posts[Kav_posts."IgG" .== igg, Not("IgG")][1, i], :) for i = 1:4]...), names(Kav_posts)[2:end])
        pls[ii] = dist_violin_plot(posts, priors; y_range = (4, 9), title = "m$igg Affinities Distributions")
    end
    return pls
end

function figure3_v1()
    # fetch human KxStar here
    KxStar = median(runMCMC("humanNUTSfit_0505.dat")["KxStar"])

    df = importMurineInVitro()
    ndf = predictMurine(df; Kav = importKav(; murine = true, retdf = true))
    pl1 = plotPredvsMeasured(
        ndf;
        xx = "Value",
        yy = "Predict",
        color = "Receptor",
        shape = "Subclass",
        R2pos = (0, -1),
        title = "Raw murine prediction\nwith documented affinities",
    )

    c = runMurineMCMC("murineNUTSdepfit_0509.dat"; KxStar = KxStar)
    pl2 = plot_murineMCMC_predict(c, df; title = "Murine prediction with fitted parameters", KxStar = KxStar, R2pos = (0, -0.3))

    apls = plot_murine_MCMC_affinity(c)
    #leuk_old, leuk_new = validateLeukocyte(c; KxStar = KxStar)

    pp = plotGrid((3, 3), [pl1, pl2, apls[1], apls[2], apls[3]])
    draw(PDF("figure3_v1.pdf", 12inch, 9inch), pp)
end

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

    c = fitLeukocyteMCMC("leukNUTSfit_0509_2.dat")
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

    apls = plot_murine_MCMC_affinity(c)
    vpl1, vpl2 = validateMurineInVitro(c)

    pp = plotGrid((3, 3), [pl1, pl_noKav, pl2, apls[1], apls[2], apls[3], vpl1, vpl2])
    draw(PDF("figure3_v2.pdf", 12inch, 9inch), pp)
end
