""" Figure S2: Anja's receptor quantification """

function figureS2()
    setGadflyTheme()

    df1 = CSV.File(joinpath(dataDir, "murine-FcgR-abundance.csv"), comment = "#") |> DataFrame
    df1 = combine(groupby(df1, ["Cells", "Receptor"]), names(df1, :Count) .=> geocmean, names(df1, :Count) .=> std)
    rename!(df1, :"Count_geocmean" => "Abundance")
    df1[!, "ymin"] .= df1[!, "Abundance"] .- df1[!, "Count_std"]
    df1[df1[!, "ymin"] .< 1.0, "ymin"] .= 1.0
    df1[!, "ymax"] .= df1[!, "Abundance"] .+ df1[!, "Count_std"]
    pl1 = plot(
        df1,
        x = "Cells",
        y = "Abundance",
        ymin = "ymin",
        ymax = "ymax",
        color = "Receptor",
        Geom.bar(position = :dodge),
        Geom.errorbar,
        Stat.dodge,
        Scale.x_discrete,
        Scale.y_log10,
        style(bar_spacing = 1.5mm),
        Guide.title("Murine immune cell FcR expression"),
    )

    df2 = CSV.File(joinpath(dataDir, "human-FcgR-abundance.csv"), comment = "#") |> DataFrame
    df2 = combine(groupby(df2, ["Cells", "Receptor"]), names(df2, :Count) .=> geocmean, names(df2, :Count) .=> std)
    rename!(df2, :"Count_geocmean" => "Abundance")
    df2[!, "ymin"] .= df2[!, "Abundance"] .- df2[!, "Count_std"]
    df2[df2[!, "ymin"] .< 1.0, "ymin"] .= 1.0
    df2[!, "ymax"] .= df2[!, "Abundance"] .+ df2[!, "Count_std"]
    pl2 = plot(
        df2,
        x = "Cells",
        y = "Abundance",
        ymin = "ymin",
        ymax = "ymax",
        color = "Receptor",
        Geom.bar(position = :dodge),
        Geom.errorbar,
        Stat.dodge,
        Scale.x_discrete,
        Scale.y_log10,
        style(bar_spacing = 1.5mm),
        Guide.title("Human immune cell FcR expression"),
    )
    draw(PDF("output/figureS2.pdf", 1000px, 400px), plotGrid((1, 2), [pl1, pl2]))
end
