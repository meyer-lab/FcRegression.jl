""" Figure 3: fitted murine affinities """

function plot_murine_ADVI_affinity(c::Chains = runMurineMCMC())
    Kav_priors = murineKavDist()
    pls = Vector{Union{Gadfly.Plot, Context}}(undef, 3)
    for (ii, igg) in enumerate(Kav_priors[!, "IgG"])
        priors = reshape(Matrix(Kav_priors[Kav_priors."IgG" .== igg, Not("IgG")]), :)
        posts = DataFrame(hcat([reshape(Matrix(c["Kav[$i]"]), :) for i = (ii * 4 - 3):(ii * 4)]...), names(Kav_priors)[2:end])
        pls[ii] = dist_violin_plot(posts, priors; y_range = (4, 9), title = "m$igg Affinities Distributions")
    end
    return pls
end

function figure3()
    # fetch human KxStar here

    df = importMurineInVitro()
    ndf = predictMurine(df)
    pl1 = plotPredvsMeasured(
        ndf;
        xx = "Value",
        yy = "Predict",
        color = "Receptor",
        shape = "Subclass",
        clip2one = false,
        R2pos = (0, -1),
        title = "Raw murine prediction without fitting",
    )
    #pl2, _ = MAPmurineLikelihood()

    c = runMurineMCMC("murine_NUTS_0502_5.62332e-16.dat")
    apls = plot_murine_ADVI_affinity(c)

    ndf2 = predictMurine(c, df; KxStar = 5.62332e-16)
    pl3 = plotPredvsMeasured(
        ndf2;
        xx = "Value",
        yy = "Predict",
        color = "Receptor",
        shape = "Subclass",
        clip2one = false,
        R2pos = (0, -1),
        title = "Murine prediction after fitting",
    )

    leuk = importMurineLeukocyte()
    leuk_old, leuk_new = validateMurine(c, leuk; KxStar = 5.62332e-16)

    pp = plotGrid((3, 3), [pl1, pl3, apls[1], apls[2], apls[3], leuk_old, leuk_new])
    draw(PDF("figure3.pdf", 12inch, 9inch), pp)
end
