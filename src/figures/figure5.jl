function plotLbound(Rtot = importRtot(; murine = false, retdf = true); 
        title = "", 
        cellTypes = ["ncMO", "cMO", "Neu"], 
        kwargs...)

    Kav0 = importKav(; murine = false)
    Kav1 = extractNewHumanKav()
    df0 = predictLbound(Kav0, Rtot; kwargs...)
    df0."Affinity" .= "Documented"
    df1 = predictLbound(Kav1, Rtot; kwargs...)
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
                    key_position = :none,
                    major_label_font_size = 8pt,
                    minor_label_font_size = 8pt,
                ),
            ) for igg in unique(df."IgG")
        ]
end

function plotEffectorMeasured()
    df = importEffectorBind(; avg = true)
    df."Valency" .= Symbol.(df."Valency")
    df[df."Value" .< 0.0, "Value"] .= 0.0
    df[df."xmin" .< 0.0, "xmin"] .= 0.0

    return [
            plot(
                df[df."Subclass" .== igg, :],
                x = "Valency",
                xgroup = "Cell",
                y = "Value",
                color = [colorant"hsl(350, 40%, 45%)"], #"Valency",
                ymin = "xmin",
                ymax = "xmax",
                Scale.x_discrete,
                Geom.subplot_grid(Geom.errorbar, Geom.bar(position = :dodge),),
                Guide.title("Measured bound $igg"),
                Guide.xlabel(nothing),
                Guide.ylabel(igg == "IgG1" ? "ΔMFI" : nothing),
                style(
                    bar_spacing = 3px,
                    plot_padding = [0.0pt, 0.0pt, 0.0pt, 0.0pt],
                    key_position = igg == "IgG4" ? :right : :none,
                    major_label_font_size = 8pt,
                    minor_label_font_size = 8pt,
                    stroke_color = c -> "black", 
                    errorbar_cap_length = 6px,
                ),
            ) for igg in unique(df."Subclass")
        ]
end


function plotEffectorPred(; Kav = extractNewHumanKav(), title = "", legend = true, kwargs...)
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
        Guide.xlabel("Measurements (MFI)", orientation = :horizontal),
        Guide.ylabel("Predicted binding", orientation = :vertical),
        Guide.title(title),
        Scale.x_log10,
        Scale.y_log10,
        Coord.cartesian(xmin = 1),
        Guide.xticks(orientation = :horizontal),
        Scale.color_discrete_manual(colorSubclass...),
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

function figure5(ssize = (8.5inch, 7inch); cellTypes = ["ncMO", "cMO", "Neu"], kwargs...)
    setGadflyTheme()

    measured = plotEffectorMeasured()
    lbounds = plotLbound(; cellTypes = cellTypes)

    c = FcRegression.rungMCMC("humanKavfit_0701.dat"; dat = :hCHO, mcmc_iter = 1_000);
    pms = FcRegression.extractMCMC(c; dat = :hCHO)

    oldPred = plotEffectorPred(; Kav = extractNewHumanKav(; old = true), 
        title = "Documented Affinity", legend = false, KxStar = pms["KxStar"])  # R2 = 0.6671
    newPred = plotEffectorPred(; Kav = extractNewHumanKav(; old = false), 
        title = "Updated Affinity", legend = true, KxStar = pms["KxStar"])  # R2 = 0.6585

    pl = FcRegression.plotGrid(
        (3, 4),
        [measured[1], measured[2], measured[3], measured[4], 
        lbounds[1], lbounds[2], lbounds[3], lbounds[4],
        oldPred, newPred, nothing, nothing];
        sublabels = "abcdefghij  ",
        widths = [1.1 1 1 1; 1.15 1 1 1; 0.8 1 0.8 0.1],
        heights = [1.3, 1.3, 1.5],
        kwargs...,
    )
    draw(PDF("output/figure5.pdf", ssize[1], ssize[2]), pl)
end
