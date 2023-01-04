

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


function plotEffectorPred(; Kav = FcRegression.extractNewHumanKav(), title = "", legend = true, kwargs...)
    df = importEffectorBind(; avg = true)
    pred = predictLbound(Kav; kwargs...)
    rename!(pred, "IgG" => "Subclass")

    jdf = innerjoin(df, pred, on = ["Valency", "Subclass", "Cell"])
    jdf."Valency" .= Symbol.(jdf."Valency")
    jdf."Lbound" ./= geocmean(jdf."Lbound") / geocmean(jdf."Value")
    jdf[jdf."xmin" .<= 1.0, "xmin"] .= 10.0
    jdf[jdf."Value" .<= 1.0, "Value"] .= 10.0

    r2 = R2(jdf."Value", jdf."Lbound"; logscale = true)

    return plot(
        jdf,
        x = "Value",
        y = "Lbound",
        xmin = "xmin",
        xmax = "xmax",
        color = "Subclass",
        shape = "Cell",
        Geom.point,
        Geom.errorbar,
        Guide.xlabel("Measurements", orientation = :horizontal),
        Guide.ylabel("Predicted binding", orientation = :vertical),
        Guide.title(title),
        Scale.x_log10,
        Scale.y_log10,
        Coord.cartesian(xmin = 1),
        Scale.color_discrete_manual(FcRegression.colorSubclass...),
        Geom.abline(color = "black"),
        Guide.annotation(
            compose(
                context(),
                text(3.5, 1, "<i>R</i><sup>2</sup> = " * @sprintf("%.4f", r2)),
                stroke("black"),
                fill("black"),
                font("Helvetica-Bold"),
            ),
        ),
        style(errorbar_cap_length = 1px, key_position = legend ? :right : :none),
    )
end

function figure5(ssize = (8.5inch, 5inch); cellTypes = ["ncMO", "cMO", "Neu"], kwargs...)
    setGadflyTheme()

    lbounds = plotLbound(; cellTypes = cellTypes)
    oldPred = plotEffectorPred(; Kav = extractNewHumanKav(; old = true), 
        title = "Documented Affinity", legend = false)
    newPred = plotEffectorPred(; Kav = extractNewHumanKav(; old = false), 
        title = "Updated Affinity", legend = true)
    
    c = FcRegression.rungMCMC("humanKavfit_0701.dat"; dat = :hCHO, mcmc_iter = 1_000);
    pms = FcRegression.extractMCMC(c; dat = :hCHO)
    newPred2 = plotEffectorPred(; Kav = extractNewHumanKav(; old = false), 
        title = "Updated Affinity", legend = true, KxConst = pms["KxStar"])

    pl = FcRegression.plotGrid(
        (2, 4),
        [lbounds[1], lbounds[2], lbounds[3], lbounds[4],
        oldPred, newPred, nothing, nothing];
        sublabels = "abcdef  ",
        widths = [1.1 1 1 1.4; 1 1 0.1 0.1],
        heights = [1.3, 1.5],
        kwargs...,
    )
    draw(PDF("output/figure5.pdf", ssize[1], ssize[2]), pl)
end
