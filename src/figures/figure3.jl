""" Figure 3: fitted murine affinities """

function plot_murine_ADVI_affinity(q = runMurineMCMC())
    s = rand(q, 100000)[5:16, :]
    Kav_priors = murineKavDist()
    pls = Vector{Union{Gadfly.Plot, Context}}(undef, 3)
    for (ii, igg) in enumerate(Kav_priors[!, "IgG"])
        priors = reshape(Matrix(Kav_priors[Kav_priors."IgG" .== igg, Not("IgG")]), :)
        posts = DataFrame(Matrix(s[(ii * 4 - 3):(ii * 4), :]'), names(Kav_priors)[2:end])
        pls[ii] = dist_violin_plot(posts, priors; y_range = (5, 10), title = "m$igg Affinities Distributions")
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
    pl2, _ = MAPmurineLikelihood()

    q = runMurineMCMC()
    apls = plot_murine_ADVI_affinity(q)
    _, pl3 = predictLeukocyte(; average = true, title = "Leukocyte binding raw predictions")

    qvals = median(rand(q, 100000), dims = 2)
    Kav = importKav(; murine = true, retdf = true)
    Kav[Kav."IgG" .!= "IgG3", Not("IgG")] = reshape(qvals[5:16], 3, 4)
    _, pl4 = predictLeukocyte(; average = true, Kav = Kav, KxStar = qvals[17], title = "Leukocyte binding prediction\nwith updated affinities")

    pp = plotGrid((3, 3), [nothing, pl1, pl2, apls[1], apls[2], apls[3], pl3, pl4])
    draw(PDF("figure3.pdf", 11inch, 9inch), pp)
end
