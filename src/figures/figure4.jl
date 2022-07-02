""" Figure 4: Affinity updates """

"""
Ideas:
Fitted affinity vs. measured binding

"""

function figure4()
    c = rungMCMC("humanKavfit_0701.dat"; dat = :hCHO, mcmc_iter = 1_000)
    pl_igg = plotAffinityViolin(c; murine = false, y_range = (5, 8))

    pp = plotGrid((1, 4), [pl_igg[1], pl_igg[2], pl_igg[3], pl_igg[4]])
    draw(PDF("figure4.pdf", 12inch, 3inch), pp)
end
