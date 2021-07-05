""" Figure 2: we can accurately account for mixed ICs """

fitRecepExp_avg = Dict(
    "FcgRIIA-131H"  => 121276.42463581044,
    "FcgRIIIA-158F" => 140477.31559988233,
    "FcgRI"         => 5053.07323314973,
    "FcgRIIA-131R"  => 388164.6494841347,
    "FcgRIIB-232I"  => 110809.66312386713,
    "FcgRIIIA-158V" => 130955.9353147472
)

fitRecepExp = Dict(
    "FcgRIIA-131H"  => 50100.85,
    "FcgRIIIA-158F" => 66406.39,
    "FcgRI"         => 5053.07,
    "FcgRIIA-131R"  => 134438.67,
    "FcgRIIB-232I"  => 43598.72,
    "FcgRIIIA-158V" => 48764.04
)

fitValencies_avg = [10.05236308028924, 29.164651398831566]
fitValencies = [10.873127313836179, 58.23798729531203]

fitKx_avg = 5.753975202649048e-16
fitKx = 5.753975202780001e-16

conversion_facs = [0.3402021565455237, 1.2248156708972613, 1.361209344637733, 0.9090341348129763, 0.6854599604408271, 1.0552201581533605, 1.5124684176685206, 0.6962631533301593, 0.19556733467114384, 0.5662981145220682, 0.45755828877894, 1.0953797563645433, 2.230368949209136, 0.6922007199858611, 1.2057509006401947, 0.43086556742756954, 0.2052997458136844]

function plotPredvsMeasured(df; xx = "Adjusted", yy = "Predict", xxlabel = "Actual", yylabel = "Predicted", color = "Valency", shape = "Cell")
    setGadflyTheme()

    df[!, "Valency"] .= Symbol.(df[!, "Valency"])
    df[(df[!, xx]) .< 1.0, xx] .= 1.0
    df[(df[!, yy]) .< 1.0, yy] .= 1.0

    xmins = "StdDev" in names(df) ? (df[!, xx] .- df[!, "StdDev"]) : df[!, xx]
    xmaxs = "StdDev" in names(df) ? (df[!, xx] .+ df[!, "StdDev"]) : xmins
    xmins[xmins .< 1.0] .= 1.0
    xmaxs[xmaxs .< 1.0] .= 1.0

    r2 = R2((df[!, xx]), (df[!, yy]))
    return plot(
        df,
        x = xx,
        y = yy,
        xmin = xmins,
        xmax = xmaxs,
        color = color,
        shape = shape,
        Geom.point,
        "StdDev" in names(df) ? Geom.errorbar : Guide.xlabel(xxlabel),
        Guide.ylabel(yylabel, orientation = :vertical),
        Guide.title("R^2: $r2"),
        Scale.x_log10,
        Scale.y_log10,
        Scale.color_discrete_manual(
            Scale.color_discrete().f(10)[1],
            Scale.color_discrete().f(10)[3],
            Scale.color_discrete().f(10)[2],
            Scale.color_discrete().f(10)[4:end]...,
        ),
        Geom.abline(color = "green"),
    )
end


function figure2(IgGx_Only = false, avg = true)
    data = avg ? averageData(loadMixData()) : loadMixData()

    data."NewValency" = convert.(Float64, data."Valency")
    data[data[!, "NewValency"] .== 4.0, "NewValency"] .= avg ? fitValencies_avg[1] : fitValencies[1]
    data[(data[!, "NewValency"]) .== 33.0, "NewValency"] .= avg ? fitValencies_avg[2] : fitValencies[2]

    if IgGx_Only
        data = data[data[!, "%_1"] .!= 10 / 100, :]
        data = data[data[!, "%_1"] .!= 33 / 100, :]
        data = data[data[!, "%_1"] .!= 66 / 100, :]
        data = data[data[!, "%_1"] .!= 90 / 100, :]
    end

    """if adjusted
        df = MixtureFit(data; logscale = true)["df"]
        xvar = "Adjusted"
    else"""
    df = predictMix(data; recepExp = avg ? fitRecepExp_avg : fitRecepExp, KxStar = avg ? fitKx_avg : fitKx)
    df[!, :Adjusted] = df[!, "Value"]
    if !avg
        df = fit_experiments(conversion_facs, df)
    end

    draw(SVG("figure2.svg", 1300px, 600px), plotGrid((1, 2), [nothing, plotPredvsMeasured(df)]))
end

function figure2c()
    pl = plotPredvsMeasured(PCA_dimred(); xx = "PCA", yy = "Predict", xxlabel = "Actual on imputed PC1")
    draw(SVG("figure2c.svg", 700px, 600px), pl)
end
