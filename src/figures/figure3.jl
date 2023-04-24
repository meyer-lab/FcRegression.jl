""" Figure 3: Explain model, fit data, and validate with Robinett """

function splot_predVorig(
    df,
    p2o_scale = 1.0;
    legend = true,
    ll = 100,
    y_label = true,
    Kav::DataFrame,
    Rtot = importRtotDist(:hCHO; retdf = true, regular = true),
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
    df."Valency" = Symbol.(df."Valency")
    gdf."Predict" ./= p2o_scale

    setGadflyTheme()

    return plot(
        layer(gdf, x = "%_1", y = "Predict", color = "Valency", Geom.line),
        layer(df, x = "%_1", y = "Value", ymin = "xmin", ymax = "xmax", color = "Valency", Geom.line, Geom.errorbar, style(line_style = [:dot])),
        Scale.x_continuous(labels = n -> "$IgGX $(trunc(Int, n*100))%\n$IgGY $(trunc(Int, 100-n*100))%"),
        Scale.y_continuous(; minvalue = 0.0, labels = n -> "$(round(n, digits=1))\n $(trunc(Int, n*p2o_scale))"),
        Scale.color_discrete_manual(colorValency..., colorValency...),
        Guide.xlabel(nothing),
        Guide.ylabel(y_label ? "Predicted binding - Normalized RFU" : nothing, orientation = :vertical),
        Guide.title("$IgGX-$IgGY to $cell_name, updated"),
        style(key_position = legend ? :right : :none),
    )
end


function figure3(ssize = (10inch, 8inch); kwargs...)
    raw_predict = predMix(
        averageMixData(loadMixData());
        Kav = importKav(; murine = false, retdf = true),
        Rtot = importRtotDist(:hCHO; regular = true, retdf = true),
    )

    raw_pred_pl = plotPredvsMeasured(
        raw_predict;
        xx = "Value",
        xxlabel = "Measured",
        color = "Receptor",
        shape = "Valency",
        title = "Raw model predictions without fitting",
        R2pos = (0, -2.7),
        legend = false,
    )

    df = loadMixData()
    Kav_old = importKavDist(; murine = false, regularKav = true, retdf = true)
    c_noKav = rungMCMC("fit_mixture_no_Kav.dat"; dat = :hCHO, Kavd = Kav_old)
    pl_noKav = plotMCMCPredict(c_noKav, df; dat = :hCHO, Kav = Kav_old, R2pos = (0, -2), 
        title = "Predictions with all but affinity fitting", legend = false)

    c = rungMCMC("fit_mixture.dat"; dat = :hCHO, mcmc_iter = 1_000)
    pl1 = plotMCMCPredict(
        c,
        df;
        dat = :hCHO,
        title = "All predictions with affinities\ninferred from single IgG measurements",
        R2pos = (0, -2.5),
        legend = false,
    )
    pl2 = plotMCMCPredict(
        c,
        df[(df."%_1" .!= 1.0) .& (df."%_2" .!= 1.0), :];
        dat = :hCHO,
        title = "Mixture predictions with affinities\ninferred from single IgG measurements",
        R2pos = (0, -2.5),
        legend = false,
    )
    # Robinett fitting
    rob1, rob2 = validateFittedKav(c, "fit_robinett.dat"; murine = false, legend = false)
    pl_legend = plotMCMCPredict(c, df; dat = :hCHO, R2pos = (0, -2.5), legend = true)   # just to insert the legend


    # Visually compare a few examples of old and new predictions
    df = averageMixData()
    df_igg24_1 = df[(df."Receptor" .== "FcgRI") .& (df."subclass_1" .== "IgG2") .& (df."subclass_2" .== "IgG4"), :]
    igg24_old = splot_predData(
        df_igg24_1;
        legend = false,
        ll = 100,
        y_label = true,
        Kav = importKav(; murine = false, retdf = true),
        yticks = [0, 10000, 20000, 30000, 40000, 50000],
    )
    igg24_new = splot_predVorig(df_igg24_1, 20000; legend = false, ll = 100, y_label = true, Kav = extractNewHumanKav())

    df_igg34_2b = df[(df."Receptor" .== "FcgRIIB-232I") .& (df."subclass_1" .== "IgG3") .& (df."subclass_2" .== "IgG4"), :]
    igg34_old = splot_predData(
        df_igg34_2b;
        legend = false,
        ll = 100,
        y_label = true,
        Kav = importKav(; murine = false, retdf = true),
        yticks = [0, 400, 800, 1200],
    )
    igg34_new = splot_predVorig(df_igg34_2b, 400.0; legend = false, ll = 100, y_label = true, Kav = extractNewHumanKav())

    # Put everything together
    pp = plotGrid(
        (3, 4),
        [
            nothing pl1 igg24_old
            nothing pl2 igg24_new
            raw_pred_pl rob1 igg34_old
            pl_noKav rob2 igg34_new
        ];
        sublabels = "abcdefghijkl",
        widths = [1 1 1 1; 1 1 1 1; 1 1 1 1],
        heights = [1, 1, 0.8],
        kwargs...,
    )
    draw(PDF("output/figure3.pdf", ssize[1], ssize[2]), pp)
end
