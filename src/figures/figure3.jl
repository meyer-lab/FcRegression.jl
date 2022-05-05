""" Figure 3: fitted murine affinities """

function plot_murine_MCMC_affinity(c::Chains = runMurineMCMC())
    Kav_priors = murineKavDist()
    pls = Vector{Union{Gadfly.Plot, Context}}(undef, 3)
    for (ii, igg) in enumerate(Kav_priors[!, "IgG"])
        priors = reshape(Matrix(Kav_priors[Kav_priors."IgG" .== igg, Not("IgG")]), :)
        posts = DataFrame(hcat([reshape(Matrix(c["Kav[$i]"]), :) for i = (ii*4-3):(ii*4)]...), names(Kav_priors)[2:end])
        pls[ii] = dist_violin_plot(posts, priors; y_range = (4, 9), 
            title = "m$igg Affinities Distributions")
    end
    return pls
end

function figure3()
    # fetch human KxStar here
    KxStar = median(runMCMC("humanNUTSfit_0504.dat")["KxStar"])
    # KxStar = 2.5072552083386547e-14

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
        title = "Rawndf =  murine prediction without fitting",
    )

    c = runMurineMCMC("murineNUTSdepfit_0504.dat"; KxStar = KxStar)
    pl2 = plot_murineMCMC_predict(c, df; title = "Murine prediction with fitted parameters", KxStar = KxStar)

    apls = plot_murine_MCMC_affinity(c)
    leuk_old, leuk_new = validateMurine(c; KxStar = KxStar)

    pp = plotGrid((3, 3), [pl1, pl2, apls[1], apls[2], apls[3], leuk_old, leuk_new])
    draw(PDF("figure3.pdf", 12inch, 9inch), pp)
end
