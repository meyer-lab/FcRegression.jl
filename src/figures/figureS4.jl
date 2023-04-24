""" Supplemental figure to Fig. 5 """


function figureS4(; ssize = (6inch, 3inch), widths = [3 4])
    setGadflyTheme()

    # FcgRIIIB: IgG2 ~= 60000, IgG4 ~= 80000

    c = inferCD16b(; mcmc_iter = 5000)
    af1, af2 = plotCD16bCHOaff(c[4000:end])

    draw(
        PDF("output/figureS4.pdf", ssize...),
        plotGrid(
            (1, 2),
            [af1, af2];
            widths = widths,
        )
    )
end