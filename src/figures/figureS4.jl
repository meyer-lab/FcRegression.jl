function figureS4(ssize = (10inch, 5inch); cellTypes = ["ncMO", "cMO", "Neu"], kwargs...)
    setGadflyTheme()

    df0 = predictLbound(extractNewHumanKav(; old = true); specificRcp = true)
    df0."Affinity" .= :Documented
    df1 = predictLbound(extractNewHumanKav(; old = false); specificRcp = true)
    df1."Affinity" .= :Updated
    df = vcat(df0, df1)
    df."Valency" .= Symbol.(df."Valency")
    @assert all(in(unique(df."Cell")).(cellTypes))
    df."Receptor" = replace.(df."Receptor", "gR" => "γR")

    Rtot = importRtot(; murine = false, retdf = true)
    existRcp = all(Matrix(Rtot[!, 2:end]) .> 0.0, dims = 2)

    pls = Matrix(undef, length(unique(df."IgG")), length(cellTypes))

    for (c, cell) in enumerate(cellTypes)
        for (i, igg) in enumerate(unique(df."IgG"))
            pls[i, c] = plot(
                df[(df."Cell" .== cell) .& (df."IgG" .== igg), :],
                x = "Valency",
                xgroup = "Affinity",
                y = "Lbound",
                color = "Receptor",
                Geom.subplot_grid(Geom.bar(position = :dodge)),
                Guide.title("Predicted bound of $igg to..."),
                Guide.xlabel(nothing),
                Guide.ylabel(igg == "IgG1" ? "Predict binding to $cell" : nothing),
                Scale.color_discrete_manual(colorReceptor[BitVector([existRcp...])]...),
                style(
                    bar_spacing = 0.0pt,
                    plot_padding = [0.0pt, 0.0pt, 0.0pt, 0.0pt],
                    key_position = igg == "IgG4" ? :right : :none,
                    major_label_font_size = 8pt,
                    minor_label_font_size = 8pt,
                ),
            )
        end
    end

    pl = plotGrid(
        (size(pls)[2], size(pls)[1]),
        pls;
        widths = [1.1 1 1 1.4],
        sublabels = false,
        kwargs...,
    )
    draw(PDF("output/figureS4.pdf", ssize[1], ssize[2]), pl)
end
