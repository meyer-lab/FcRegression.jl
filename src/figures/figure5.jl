function plotLbound(Rtot = importRtot(; murine = false, retdf = true, genotype = "ZIZ"); cellTypes = ["ncMO", "cMO", "Neu"], kwargs...)

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

function plotEffectorMeasured(df = importEffectorBind(; avg = true))
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
            Geom.subplot_grid(Geom.errorbar, Geom.bar(position = :dodge)),
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

    # keep only IgG2
    jdf = jdf[in(["IgG2"]).(jdf."Subclass"), :]

    r2 = R2(jdf."Value", jdf."Lbound"; logscale = false)
    println("R2: $r2")

    return plot(
        jdf,
        x = "Value",
        y = "Lbound",
        xmin = "xmin",
        xmax = "xmax",
        color = "Valency",
        shape = "Cell",
        Geom.point,
        Geom.errorbar,
        Guide.xlabel("Measurements (MFI)", orientation = :horizontal),
        Guide.ylabel("Predicted binding", orientation = :vertical),
        Guide.title(title),
        Guide.xticks(orientation = :horizontal),
        Scale.color_discrete_manual(colorValency...),
        Geom.abline(color = "black"),
        Guide.annotation(
            compose(context(), text(3, 1, "<i>R</i><sup>2</sup> = " * @sprintf("%.4f", r2)), stroke("black"), fill("black"), font("Helvetica-Bold")),
        ),
        style(errorbar_cap_length = 1px, key_position = legend ? :right : :none),
    )
end

function figure5(ssize = (8.5inch, 7.5inch); cellTypes = ["ncMO", "cMO", "Neu"], kwargs...)
    setGadflyTheme()

    measured = plotEffectorMeasured()
    lbounds = plotLbound(; cellTypes = cellTypes)

    p0 = plotEffectorPredict(
        importEffectorBind(; avg = true),
        predictLbound(importKav(; murine = false)); 
        title = "Leukocyte binding with documented affinities"
    )
    p1 = plotEffectorPredict(
        importEffectorBind(; avg = true),
        predictLbound(extractNewHumanKav()); 
        title = "Leukocyte binding with updated affinities"
    )
    
    pl = plotGrid(
        (3, 4),
        [measured[1], measured[2], measured[3], measured[4], 
        lbounds[1], lbounds[2], lbounds[3], lbounds[4], 
        p0, p1, nothing, nothing];
        sublabels = "abcdefghij  ",
        widths = [1.1 1 1 1; 1.15 1 1 1; 1 1 .3 .3],
        heights = [2.25, 2.25, 3],
        kwargs...,
    )
    draw(PDF("output/figure5.pdf", ssize[1], ssize[2]), pl)
end
