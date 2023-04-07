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
            append!(
                layers, 
                layer(arrdf, x = :xstart, y = :ystart, xend = :xend, yend = :yend, color = [colorant"black"], Geom.segment, 
                style(line_style=[(val == Symbol("4")) ? :solid : :dot]))
            )
        end
    end

    df."Subclass" = copy(df."subclass_1")
    df[df."%_2" .== 1.0, "Subclass"] .= df[df."%_2" .== 1.0, "subclass_2"]
    df[(df."%_1" .< 1.0) .& (df."%_2" .< 1.0), "Subclass"] .= "Mixed"
    sdf = df[df."Subclass" .!= "Mixed", :]
    append!(layers, layer(df, x = xx, y = yy, color = [colorant"black"], size = [1mm], Geom.point))
    all_subclass = unique(vcat(df."subclass_1", df."subclass_2"))
    return plot(
        sdf,
        layers...,
        x = xx,
        y = yy,
        color = "Subclass",
        size = [2.5mm],
        shape = "Valency",
        Scale.color_discrete_manual(colorSubclass[[in("IgG" * string(i), all_subclass) for i in range(1,4)]]...),
        Geom.point,
        Guide.title(title),
        Guide.xticks(ticks = [-15, 0, 15], orientation = :horizontal),
        Guide.yticks(ticks = [-5, 0, 5]),
    )
end

function plot_PCA_line(score)
    pc_line_style = Dict("PC 1" => :solid, "PC 2" => :dot)
    pc_point_style = Dict("PC 1" => Shape.circle, "PC 2" => Shape.square)
    score."Valency" = Symbol.(score."Valency")

    pls = Vector{Union{Gadfly.Plot, Context}}(undef, length(unique(score."Subclass Pair")))
    for (ip, pair) in enumerate(unique(score."Subclass Pair"))
        pcs = stack(score[score."Subclass Pair" .== pair, :], 6:8, variable_name="PC", value_name="Value")
        IgGX, IgGY = split(pair, "-")
        pls[ip] = plot(
            [layer(
                pcs[(pcs."Valency" .== val) .& (pcs."PC" .== pc), :], 
                x="%_2", 
                y="Value", 
                color=[colorValency[iv]], 
                shape=[pc_point_style[pc]],
                Geom.point, 
                Geom.line,
                style(line_style = [pc_line_style[pc]]))
            for (iv, val) in enumerate([Symbol("4"), Symbol("33")]) for pc in ["PC 1", "PC 2"]]...,
            Scale.x_continuous(labels = n -> "$IgGX $(trunc(Int, n*100))%\n$IgGY $(trunc(Int, 100-n*100))%"),
            Scale.y_continuous(minvalue = -15, maxvalue = 15),
            Guide.xlabel(nothing),
            Guide.ylabel("Scores"),
            Guide.title("$pair mixtures"),
            Guide.manual_color_key("Valency", ["4", "33"], colorValency),
            Guide.manual_color_key("Component", ["PC 1", "PC 2"], [colorant"black", colorant"black"], shape=[Shape.circle, Shape.square]),
        )
    end
    return pls
end

function figure2(ssize = (13inch, 6inch); widths = [3, 3, 3, 3.2])
    setGadflyTheme()

    score, loading, vars_expl = mixtureDataPCA()
    loading.Receptor = replace.(loading.Receptor, "gR" => "γR")
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
        Guide.ylabel("Variance Explained"),
    )

    ## TODO: add percent variance explained on each PC

    SP4 = plot_PCA_score(score[score."Valency" .== 4, :]; title = "PCA Score, 4-valent ICs", xx = "PC 1", yy = "PC 2")
    SP33 = plot_PCA_score(score[score."Valency" .== 33, :]; title = "PCA Score, 33-valent ICs", xx = "PC 1", yy = "PC 2")
    SP4_13 = plot_PCA_score(score[score."Valency" .== 4, :]; title = "PCA Score, 4-valent ICs", xx = "PC 1", yy = "PC 3")
    SP33_13 = plot_PCA_score(score[score."Valency" .== 33, :]; title = "PCA Score, 33-valent ICs", xx = "PC 1", yy = "PC 3")
    
    SPs = [plot_PCA_score(score[score."Subclass Pair" .== pair, :]; title = "PCA Score, $pair") for pair in unique(score."Subclass Pair")]
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
        Scale.color_discrete_manual(colorReceptor...),
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
        Scale.color_discrete_manual(colorReceptor...),
    )

    pl = plotGrid((2, 4), [vars, SPs..., LP]; sublabels = true, widths = widths)
    #pl = plotGrid((1, 4), [vars, SP4, SP33, LP]; sublabels = "abcd", widths = widths)
    draw(PDF("output/figure2.pdf", ssize[1], ssize[2]), pl)
    return SPs
end
