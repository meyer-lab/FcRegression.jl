""" Figure 2: we can accurately account for mixed ICs """

function plotPredvsMeasured(df; xx = "Adjusted", yy = "Predict", 
    xxlabel = "Actual", yylabel = "Predicted", color = "Valency", shape = "Cell")
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
        "StdDev" in names(df) ? Geom.errorbar : 
        Guide.xlabel(xxlabel),
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


function figure2(adjusted = false, IgGx_Only = false, avg = true)
    data = avg ? averageData(loadMixData()) : loadMixData()

    data[(data[!, "Valency"]) .== 4, "Valency"] .= 11
    data[(data[!, "Valency"]) .== 33, "Valency"] .= 56

    if IgGx_Only
        data = data[data[!, "%_1"] .!= 10 / 100, :]
        data = data[data[!, "%_1"] .!= 33 / 100, :]
        data = data[data[!, "%_1"] .!= 66 / 100, :]
        data = data[data[!, "%_1"] .!= 90 / 100, :]
    end

    fitrecepExp = Dict(
        "FcgRIIA-131H"  => 59413.0,
        "FcgRIIIA-158F" => 73689.0,
        "FcgRI"         => 5053.07,
        "FcgRIIA-131R"  => 1.51025e5,
        "FcgRIIB-232I"  => 69149.7,
        "FcgRIIIA-158V" => 71607.7
    )

    if adjusted
        df = MixtureFit(data; logscale = true)["df"]
        xvar = "Adjusted"
    else
        df = predictMix(data; recepExp = fitrecepExp)
        xvar = "Value"
    end

    draw(SVG("figure2.svg", 1300px, 600px), plotGrid((1, 2), [nothing, plotPredvsMeasured(df; xx = xvar)]))
end

function figure2c()
    pl = plotPredvsMeasured(PCA_dimred(); xx="PCA", yy="Predict", xxlabel="Actual on imputed PC1")
    draw(SVG("figure2c.svg", 700px, 600px), pl)
end
