function figureS2()
    df1 = CSV.File(joinpath(dataDir, "murine-FcgR-abundance.csv"), comment = "#") |> DataFrame!
    df1 = combine(groupby(df1, [:Cells, :Receptor]), names(df1, :Count) .=> geocmean)
    rename!(df1, :Count_geocmean => :Abundance)
    pl1 = plot(df1, x = :Cells, y = :Abundance, color = :Receptor, Geom.bar(position = :dodge))
    df2 = CSV.File(joinpath(dataDir, "human-FcgR-abundance.csv"), comment = "#") |> DataFrame!
    df2 = combine(groupby(df2, [:Cells, :Receptor]), names(df2, :Count) .=> geocmean)
    rename!(df2, :Count_geocmean => :Abundance)
    pl2 = plot(df2, x = :Cells, y = :Abundance, color = :Receptor, Geom.bar(position = :dodge))
    draw(SVG("figureS2.svg", 800px, 400px), plotGrid((1, 2), [pl1, pl2]))
end
