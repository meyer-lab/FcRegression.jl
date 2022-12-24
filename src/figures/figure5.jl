""" Plot in vivo regression results. Return predictions, cell weights, and receptor weights for MAP and MCMC. """
function plot_regressions(df; Kav::DataFrame, murine = false, cellTypes = nothing, ptitle = "")
    opt, optcv, cdf = runRegMAP(df; murine = murine, Kav = Kav, cellTypes = cellTypes)
    c0, ccdf0 = FcRegression.runRegMCMC(df; murine = murine, Kav = Kav0, mcmc_iter = 200, cellTypes = cellTypes)

    pl_map = plotRegMCMC(opt, deepcopy(df); ptitle = ptitle * "[MAP]", Kav = Kav, cellTypes = cellTypes, colorL = "Genotype", shapeL = "Condition")
    pl_mcmc = plotRegMCMC(c, deepcopy(df); ptitle = ptitle * "[MCMC]", Kav = Kav, cellTypes = cellTypes, colorL = "Genotype", shapeL = "Condition")
    cell_map, act_map = plotRegParams(optcv; ptitle = ptitle * "[MAP]", legend = true, Kav = Kav, cellTypes = cellTypes)
    cell_mcmc, act_mcmc = plotRegParams(c; ptitle = ptitle * "[MCMC]", legend = true, Kav = Kav, cellTypes = cellTypes)
    return [pl_map, cell_map, act_map, pl_mcmc, cell_mcmc, act_mcmc]
end


function predictLbound(
    Kav = FcRegression.extractNewHumanKav(),
    Rtot = FcRegression.importRtot(; murine = false, retdf = true);
    L0 = 1e-9,
    fs = [4, 33],
    KxStar = KxConst,
)
    Rtot = Rtot[in(names(Kav[!, Not("IgG")])).(Rtot."Receptor"), :]

    """ Predict Lbound of each cell type based on Kav """
    df = DataFrame((IgG=x, Cell=y, Valency=z) for x in Kav."IgG" for y in names(Rtot)[2:end] for z in fs)
    df."Lbound" .= 0.0

    for igg in unique(df."IgG")
        kav = Matrix(Kav[Kav."IgG" .== igg, 2:end])
        for cn in unique(df."Cell")
            for f in unique(df."Valency")
                df[(df."IgG".==igg) .& (df."Cell" .== cn) .& (df."Valency" .== f), "Lbound"] .= polyfc(L0, KxStar, f, Rtot[!, cn], [1.0], kav).Lbound
            end
        end
    end
    return df
end


function plotLbound(Rtot = importRtot(; murine = false, retdf = true); 
        title = "", 
        cellTypes = ["ncMO", "cMO", "Neu"], 
        kwargs...)

    Kav0 = FcRegression.importKav(; murine = false)
    Kav1 = FcRegression.extractNewHumanKav()
    df0 = FcRegression.predictLbound(Kav0, Rtot; kwargs...)
    df0."Affinity" .= "Documented"
    df1 = FcRegression.predictLbound(Kav1, Rtot; kwargs...)
    df1."Affinity" .= "Updated"
    df = vcat(df0, df1)

    df = df[in(cellTypes).(df."Cell"), :]
    df."Valency" .= Symbol.(df."Valency")

    return [
            plot(
                df[df."IgG" .== igg, :],
                x = "Valency",
                xgroup = "Cell",
                y = "Lbound",
                color = "Affinity",
                Geom.subplot_grid(Geom.bar(position = :dodge)),
                Guide.title("Predicted bound $igg"),
                Guide.xlabel(nothing),
                Guide.ylabel(igg == "IgG1" ? "Predict binding" : nothing),
                Scale.color_discrete_manual(colorAffinity...),
                style(
                    bar_spacing = 0.0pt,
                    plot_padding = [0.0pt, 0.0pt, 0.0pt, 0.0pt],
                    key_position = igg == "IgG4" ? :right : :none,
                    major_label_font_size = 8pt,
                    minor_label_font_size = 8pt,
                ),
            ) for igg in unique(df."IgG")
        ]
end




function figure5(ssize = (8.5inch, 9.5inch); cellTypes = ["ncMO", "cMO", "Neu"], mcmc_iter = 10000, suffix = "1213C_", kwargs...)
    setGadflyTheme()
    df = FcRegression.importHumanized("ITP")

    Kav0 = FcRegression.extractNewHumanKav(; old = true)
    Kav1 = FcRegression.extractNewHumanKav(; old = false)

    c1, ccdf1 = FcRegression.runRegMCMC(df, "regMCMC_$(suffix)1.dat"; murine = false, Kav = Kav1, mcmc_iter = mcmc_iter, cellTypes = cellTypes)
    c0, ccdf0 = FcRegression.runRegMCMC(df, "regMCMC_$(suffix)0.dat"; murine = false, Kav = Kav0, mcmc_iter = mcmc_iter, cellTypes = cellTypes)

    c0 = c0[(mcmc_iter ÷ 2):end]
    c1 = c1[(mcmc_iter ÷ 2):end]

    lbounds = plotLbound(; cellTypes = cellTypes)

    pl_map0 = FcRegression.plotRegMCMC(
        c0,
        deepcopy(df);
        ptitle = "documented affinities",
        colorL = "Genotype",
        shapeL = "Condition",
        legend = false,
        Kav = Kav0,
        cellTypes = cellTypes,
    )
    cell_map0, act_map0 = FcRegression.plotRegParams(c0; ptitle = "documented affinities", legend = false, Kav = Kav0, cellTypes = cellTypes)

    pl_map1 = FcRegression.plotRegMCMC(
        c1,
        deepcopy(df);
        ptitle = "updated affinities",
        colorL = "Genotype",
        shapeL = "Condition",
        legend = false,
        Kav = Kav1,
        cellTypes = cellTypes,
    )
    cell_map1, act_map1 = FcRegression.plotRegParams(c1; ptitle = "updated affinities", legend = false, Kav = Kav1, cellTypes = cellTypes)

    pl_mapL = FcRegression.plotRegMCMC(
        c1,
        deepcopy(df);
        ptitle = "updated affinities",
        colorL = "Genotype",
        shapeL = "Condition",
        legend = true,
        Kav = Kav1,
        cellTypes = cellTypes,
    )
    cell_mapL, act_mapL = FcRegression.plotRegParams(c1; ptitle = "updated affinities", legend = true, Kav = Kav1, cellTypes = cellTypes)

    pl = FcRegression.plotGrid(
        (4, 4),
        [lbounds[1], lbounds[2], lbounds[3], lbounds[4], 
        nothing, nothing, pl_map0, pl_map1, 
        pl_mapL, nothing, cell_map0, cell_map1,
        cell_mapL, nothing, act_map0, act_map1];
        sublabels = "abcde fg  hi  jk",
        widths = [1.1 1 1 1.4; 1 0.3 1 1; 1 0.3 1 1; 1 0.3 1 1],
        heights = [1.3, 1.5, 1.5, 1.5],
        kwargs...,
    )
    draw(PDF("output/figure5_$suffix.pdf", ssize[1], ssize[2]), pl)
end



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
