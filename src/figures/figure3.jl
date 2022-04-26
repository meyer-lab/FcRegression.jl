""" Figure 3: fitted murine affinities """

function plot_murine_ADVI_affinity(q = runMurineMCMC())
    s = rand(q, 100000)[5:16, :]
    Kav_priors = murineKavDist()
    pls = Vector{Union{Gadfly.Plot, Context}}(undef, 3)
    for (ii, igg) in enumerate(Kav_priors[!, "IgG"])
        priors = reshape(Matrix(Kav_priors[Kav_priors."IgG" .== igg, Not("IgG")]), :)
        posts = DataFrame(Matrix(s[(ii*4-3):(ii*4), :]'), names(Kav_priors)[2:end])
        pls[ii] = dist_violin_plot(posts, priors; title = "m$igg Affinities Distributions")
    end
    return pls
end

function figure3()
    df = importMurineInVitro()
    ndf = predictMurine(df)
    pl1 = plotPredvsMeasured(
        ndf;
        xx = "Value",
        yy = "Predict",
        color = "Receptor",
        shape = "Subclass",
        clip2one = false,
        R2pos = (-1.5, 1),
        title = "Raw murine prediction without fitting",
    )
    pl2 = MAPmurineLikelihood()
    _, pl3 = predictLeukocyte(; average = true, title = "Leukocyte binding raw predictions")

    pp = plotGrid((1, 3), [pl1, pl2, pl3])
    draw(PDF("figure3.pdf", 11inch, 3inch), pp)
end
