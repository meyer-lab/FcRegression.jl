""" Figure S3/S4: predicted binding with two affinities """

function splot_predData(df; legend = true, ll = 100, match_y = false,
        Kav::DataFrame,
        Rtot = importRtotDist(:hCHO; retdf = true, regular = true))
    cell = unique(df."Receptor")[1]
    IgGX = unique(df."subclass_1")[1]
    IgGY = unique(df."subclass_2")[1]
    palette = [Scale.color_discrete().f(3)[1], Scale.color_discrete().f(3)[3]]
    cell_name = replace(cell, "FcgR" => "FcγR")

    ymax = Dict(
        "FcgRI" => 4.0e4,
        "FcgRIIA-131H" => 2.5e5,
        "FcgRIIA-131R" => 5.0e4,
        "FcgRIIB-232I" => 1000,
        "FcgRIIIA-158F" => 1.5e5,
        "FcgRIIIA-158V" => 2.5e5,
    )

    gdf = predMix(
        DataFrame(
            Dict(
                "Valency" => repeat([4, 33], ll), 
                "Receptor" => cell,
                "subclass_1" => IgGX,
                "%_1" => range(1.0, 0.0, ll*2),
                "subclass_2" => IgGY,
                "%_2" => range(0.0, 1.0, ll*2),
            )
        ); 
        Kav = Kav, 
        Rtot = Rtot
    )
    gdf."Valency" = Symbol.(gdf."Valency")

    setGadflyTheme()
    return plot(
        gdf,
        x = "%_1",
        y = "Predict",
        color = "Valency",
        Geom.line,
        Scale.x_continuous(labels = n -> "$IgGX $(n*100)%\n$IgGY $(100-n*100)%"),
        Scale.y_continuous(; minvalue = 0.0, maxvalue = match_y ? ymax[cell] : maximum(gdf."Predict")),
        Scale.color_discrete_manual(palette[1], palette[2]),
        Guide.xlabel("", orientation = :horizontal),
        Guide.ylabel("AU", orientation = :vertical),
        Guide.xticks(orientation = :horizontal),
        Guide.title("$IgGX-$IgGY bind to $cell_name"),
        style(key_position = legend ? :right : :none),
    )
end

function figureS3(; figsize = (15inch, 13inch), widths = [3, 3, 3, 3, 3, 3.5], kwargs...)
    setGadflyTheme()
    draw(
        SVG("output/figureS3.svg", figsize[1], figsize[2]), 
        FcRegression.plotMixSubplots(FcRegression.splot_predData, FcRegression.averageMixData(); 
            widths = widths, 
            match_y = true, 
            Kav = FcRegression.importKav(; murine = false, retdf = true),
            kwargs...
        )
    )
end

function figureS4(; figsize = (15inch, 13inch), widths = [3, 3, 3, 3, 3, 3.5], kwargs...)
    setGadflyTheme()
    draw(
        SVG("output/figureS4.svg", figsize[1], figsize[2]), 
        FcRegression.plotMixSubplots(FcRegression.splot_predData, FcRegression.averageMixData(); 
            widths = widths, 
            match_y = true, 
            Kav = FcRegression.extractNewHumanKav(),
            kwargs...)
    )
end