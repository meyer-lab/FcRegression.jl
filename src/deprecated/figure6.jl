function regLOCellOut(df, Kav)
    # Leave one cell type out
    cellTs = FcRegression.humanCellTypes
    R2s = zeros(length(cellTs))
    LOOindex = LOOCV(length(cellTs))

    for (i, idx) in enumerate(LOOindex)
        m = FcRegression.regmodel(df, df."Target"; cellTypes = cellTs[idx], murine = false, Kav = Kav)
        opts = Optim.Options(iterations = 1000, show_every = 50, show_trace = false)
        opt = optimize(m, MAP(), LBFGS(; m = 20), opts)
        param = FcRegression.extractRegMCMC(
            opt;
            cellTypes = cellTs[idx],
            FcgRs = unique([String(split(fcgr, "-")[1]) for fcgr in names(Kav[!, Not("IgG")])]),
        )
        pred = FcRegression.regPred(df, param; cellTypes = cellTs[idx], murine = false, Kav = Kav)
        R2s[i] = FcRegression.R2(df."Target", pred; logscale = false)
    end
    return DataFrame("Model" => ["No " * c for c in cellTs], "R2" => R2s)
end

function regOnly1Cell(df, Kav)
    # Keep only one cell types
    cellTs = FcRegression.humanCellTypes
    R2s = zeros(length(cellTs))

    for i = 1:length(cellTs)
        m = FcRegression.regmodel(df, df."Target"; cellTypes = [cellTs[i]], murine = false, Kav = Kav)
        opts = Optim.Options(iterations = 1000, show_every = 50, show_trace = false)
        opt = optimize(m, MAP(), LBFGS(; m = 20), opts)
        param = FcRegression.extractRegMCMC(
            opt;
            cellTypes = [cellTs[i]],
            FcgRs = unique([String(split(fcgr, "-")[1]) for fcgr in names(Kav[!, Not("IgG")])]),
        )
        pred = FcRegression.regPred(df, param; cellTypes = [cellTs[i]], murine = false, Kav = Kav)
        R2s[i] = FcRegression.R2(df."Target", pred; logscale = false)
    end
    return DataFrame("Model" => ["Only " * c for c in cellTs], "R2" => R2s)
end

function regLOReceptorOut(Kav, rcps = ["FcgRI", "FcgRIIA", "FcgRIIB", "FcgRIIIA", "FcgRIIIB"])
    Rtot = FcRegression.importRtot(; murine = false, retdf = true)
    Rtot = Rtot[in(rcps).([split(r, "-")[1] for r in Rtot."Receptor"]), :]

    Kav = Kav[!, [true; in(rcps).([split(r, "-")[1] for r in names(Kav[!, 2:end])])]]

    Kav0 = FcRegression.importKav(; murine = false)
    Kav1 = FcRegression.extractNewHumanKav()
    Kav0 = Kav0[!, Not(["FcgRIIB-232T", "FcgRIIC-13N"])]
    Kav1 = Kav1[!, Not(["FcgRIIB-232T", "FcgRIIC-13N"])]

end
