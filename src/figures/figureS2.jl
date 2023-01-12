""" Figure S2/S3: predicted binding with two affinities """

function splot_predData(
    df;
    legend = true,
    ll = 100,
    match_y = false,
    y_label = true,
    Kav::DataFrame,
    Rtot = importRtotDist(:hCHO; retdf = true, regular = true),
    yticks = :auto,
)
    cell = unique(df."Receptor")[1]
    IgGX = unique(df."subclass_1")[1]
    IgGY = unique(df."subclass_2")[1]
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
                "%_1" => range(1.0, 0.0, ll * 2),
                "subclass_2" => IgGY,
                "%_2" => range(0.0, 1.0, ll * 2),
            ),
        );
        Kav = Kav,
        Rtot = Rtot,
    )
    gdf."Valency" = Symbol.(gdf."Valency")

    setGadflyTheme()
    return plot(
        gdf,
        x = "%_1",
        y = "Predict",
        color = "Valency",
        Geom.line,
        Scale.x_continuous(labels = n -> "$IgGX $(trunc(Int, n*100))%\n$IgGY $(trunc(Int, 100-n*100))%"),
        Scale.y_continuous(; minvalue = 0.0, maxvalue = match_y ? ymax[cell] : maximum(gdf."Predict")),
        Scale.color_discrete_manual(colorValency...),
        Guide.xlabel("", orientation = :horizontal),
        Guide.ylabel(y_label ? "Predicted CHO binding" : nothing, orientation = :vertical),
        Guide.xticks(orientation = :horizontal),
        Guide.yticks(ticks = yticks, orientation = :horizontal),
        Guide.title("$IgGX-$IgGY to $cell_name"),
        style(key_position = legend ? :right : :none),
    )
end

function figureS2(; figsize = (14inch, 13inch), widths = [3.4, 3, 3, 3, 3, 3.7], kwargs...)
    setGadflyTheme()
    draw(
        PDF("output/figureS2.pdf", figsize[1], figsize[2]),
        plotMixSubplots(
            splot_predData,
            averageMixData();
            widths = widths,
            match_y = true,
            Kav = importKav(; murine = false, retdf = true),
            kwargs...,
        ),
    )
end

function figureS3(; figsize = (14inch, 13inch), widths = [3.4, 3, 3, 3, 3, 3.7], kwargs...)
    setGadflyTheme()
    draw(
        PDF("output/figureS3.pdf", figsize[1], figsize[2]),
        plotMixSubplots(splot_predData, averageMixData(); widths = widths, match_y = true, Kav = extractNewHumanKav(), kwargs...),
    )
end
