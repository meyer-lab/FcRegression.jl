""" Figure 2: Explore mixture binding data with PCA """

using ColorSchemes
igg_color_designation = Dict([humanIgG[i] => Scale.color_discrete().f(4)[i] for i = 1:length(humanIgG)])
igg_pair_color(iggA, iggB; tot = 5) = reverse([i for i in ColorScheme(range(igg_color_designation[iggA], igg_color_designation[iggB], length = tot))])

function plot_PCA_score(df; title = "Score", xx = "PC 1", yy = "PC 2")
    df[!, "Valency"] .= Symbol.(df[!, "Valency"])
    layers = []
    for val in unique(df."Valency")
        for pair in unique(df."Subclass Pair")
            ddf = df[(df."Subclass Pair" .== pair) .& (df."Valency" .== val), :]
            sort!(ddf, ["%_2"])
            arrdf = DataFrame(xstart = Float64[], ystart = Float64[], xend = Float64[], yend = Float64[], Subclass = String[])
            for ii = 1:(nrow(ddf) - 1)
                push!(arrdf, [ddf[ii, xx], ddf[ii, yy], ddf[ii + 1, xx], ddf[ii + 1, yy], "Mixed"])
            end
            append!(layers, layer(arrdf, x = :xstart, y = :ystart, xend = :xend, yend = :yend, color = [colorant"black"], Geom.segment))
            #color=igg_pair_color(ddf."subclass_1"[1], ddf."subclass_2"[1]; tot=nrow(arrdf))
        end
    end

    df."Subclass" = copy(df."subclass_1")
    df[df."%_2" .== 1.0, "Subclass"] .= df[df."%_2" .== 1.0, "subclass_2"]
    df[(df."%_1" .< 1.0) .& (df."%_2" .< 1.0), "Subclass"] .= "Mixed"
    sdf = df[df."Subclass" .!= "Mixed", :]
    append!(layers, layer(df, x = xx, y = yy, color = [colorant"black"], size = [1mm], Geom.point))
    return plot(
        sdf,
        layers...,
        x = xx,
        y = yy,
        color = "Subclass",
        size = [3mm],
        Geom.point,
        Guide.title(title),
        Guide.xticks(ticks = [-15, 0, 15], orientation = :horizontal),
        Guide.yticks(ticks = [-5, 0, 5]),
    )
end

function figure2()
    setGadflyTheme()

    score, loading, vars_expl = mixtureDataPCA()
    vars = plot(
        DataFrame(Components = 1:length(vars_expl), R2X = vars_expl),
        x = "Components",
        y = "R2X",
        label = [@sprintf("%.2f %%", i * 100) for i in vars_expl],
        Geom.point,
        Geom.line,
        Geom.label,
        Scale.x_discrete,
        Scale.y_continuous(minvalue = 0.5),
        Guide.title("Variance Explained by PCA"),
        Guide.xlabel("Number of components"),
        Guide.ylabel("R2X"),
    )

    ## TODO: add percent variance explained on each PC

    SP4 = plot_PCA_score(score[score."Valency" .== 4, :]; title = "PCA Score, 4-valent ICs", xx = "PC 1", yy = "PC 2")
    SP33 = plot_PCA_score(score[score."Valency" .== 33, :]; title = "PCA Score, 33-valent ICs", xx = "PC 1", yy = "PC 2")
    SP4_13 = plot_PCA_score(score[score."Valency" .== 4, :]; title = "PCA Score, 4-valent ICs", xx = "PC 1", yy = "PC 3")
    SP33_13 = plot_PCA_score(score[score."Valency" .== 33, :]; title = "PCA Score, 33-valent ICs", xx = "PC 1", yy = "PC 3")
    LP = plot(
        loading,
        x = "PC 1",
        y = "PC 2",
        color = "Receptor",
        label = "Receptor",
        Geom.point,
        Geom.label,
        Guide.title("PCA Loadings"),
        Scale.x_continuous(minvalue = -1.0, maxvalue = 1.0),
        Scale.y_continuous(minvalue = -1.0, maxvalue = 1.0),
    )
    LP_13 = plot(
        loading,
        x = "PC 1",
        y = "PC 3",
        color = "Receptor",
        label = "Receptor",
        Geom.point,
        Geom.label,
        Guide.title("PCA Loadings"),
        Scale.x_continuous(minvalue = -1.0, maxvalue = 1.0),
        Scale.y_continuous(minvalue = -1.0, maxvalue = 1.0),
    )

    pl = plotGrid((2, 4), [vars, SP4, SP33, LP, nothing, SP4_13, SP33_13, LP_13]; sublabels="abcd efg")
    return draw(PDF("figure2.pdf", 13inch, 6inch), pl)
end