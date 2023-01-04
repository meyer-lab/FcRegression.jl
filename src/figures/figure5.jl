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

function figure5(ssize = (8.5inch, 2.5inch); cellTypes = ["ncMO", "cMO", "Neu"], kwargs...)
    setGadflyTheme()

    lbounds = plotLbound(; cellTypes = cellTypes)
    pl = FcRegression.plotGrid(
        (1, 4),
        [lbounds[1], lbounds[2], lbounds[3], lbounds[4]];
        sublabels = "abcd",
        widths = [1.1 1 1 1.4],
        heights = [1.3],
        kwargs...,
    )
    draw(PDF("output/figure5.pdf", ssize[1], ssize[2]), pl)
end
