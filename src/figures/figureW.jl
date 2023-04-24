function figureW(dataType::String; legend = true, murine::Bool = true, title = nothing, kwargs...)
    res, odf, loo_res, boot_res, Cell_df = regResult(dataType; murine = murine, kwargs...)
    if murine
        df = importDepletion(dataType)
        if dataType == "HIV"
            color = "Label"
        elseif dataType == "Bcell"
            color = "Condition"
        else
            color = "Background"
        end
        shape = "Condition"
    else
        df = importHumanized(dataType)
        color = "Genotype"
        shape = (dataType == "ITP") ? "Condition" : "Concentration"
    end
    @assert all(in(names(odf)).([color, shape]))

    setGadflyTheme()
    ptitle = "$dataType"
    if title !== nothing
        ptitle *= " $title"
    end
    p1 = plotActualvFit(odf, color, shape, ptitle; legend = legend)
    p2 = plotActualvPredict(odf, color, shape, ptitle; legend = legend)
    p3 = plotCellTypeEffects(Cell_df, ptitle; legend = legend)
    return p1, p2, p3
end


function plotActualvFit(odf, colorL::Union{Symbol, String}, shapeL::Union{Symbol, String}, ptitle = ""; legend = true)
    R2anno = "<i>R</i><sup>2</sup>" * @sprintf("=%.3f", R2(odf.Y, odf.Fitted; logscale = false))
    pl = plot(
        odf,
        x = "Y",
        y = "Fitted",
        Geom.point,
        color = colorL,
        shape = shapeL,
        Guide.colorkey(),
        Guide.shapekey(),
        Scale.y_continuous(minvalue = 0.0, maxvalue = 1.0),
        Geom.abline(color = "red"),
        Guide.xlabel("Actual effect"),
        Guide.ylabel("Fitted effect"),
        Guide.title("Actual vs fitted effect ($ptitle)"),
        Guide.annotation(compose(context(), text(0.1, 0.8, R2anno), fill("black"), fontsize(10pt), font("Helvetica"))),
        style(point_size = 5px, key_position = legend ? :right : :none),
    )
    return pl
end


function plotActualvPredict(odf, colorL::Union{Symbol, String}, shapeL::Union{Symbol, String}, ptitle = ""; legend = true)
    R2anno = "<i>R</i><sup>2</sup>" * @sprintf("=%.3f", R2(odf.Y, odf.LOOPredict; logscale = false))
    pl = plot(
        odf,
        x = "Y",
        y = "LOOPredict",
        Geom.point,
        color = colorL,
        shape = shapeL,
        Guide.colorkey(),
        Guide.shapekey(),
        Geom.abline(color = "red"),
        Guide.xlabel("Actual effect"),
        Guide.ylabel("LOO predicted effect"),
        Guide.title("Actual vs LOO prediction ($ptitle)"),
        Guide.annotation(compose(context(), text(0.1, 0.8, R2anno), fill("black"), fontsize(10pt), font("Helvetica"))),
        style(point_size = 5px, key_position = legend ? :right : :none),
    )
    return pl
end


function plotCellTypeEffects(Cell_df, ptitle = ""; legend = true, maxy = nothing)
    setGadflyTheme()
    pl = plot(
        Cell_df,
        x = "Condition",
        y = "Weight",
        ymin = "ymin",
        ymax = "ymax",
        color = "Component",
        Guide.colorkey(),
        Geom.errorbar,
        Stat.dodge(axis = :x),
        Geom.bar(position = :dodge),
        Scale.x_discrete(levels = unique(Cell_df.Condition)),
        Scale.y_continuous(minvalue = 0.0, maxvalue = maxy),
        Scale.color_discrete_manual(colorant"#008f48", colorant"#ffc984", colorant"#de76b8", colorant"#9fae4f", colorant"#ff968f"),
        Guide.title("Predicted cell type weights\n($ptitle)"),
        style(key_position = legend ? :right : :none, stroke_color = c -> "black", errorbar_cap_length = 4px),
    )
    return pl
end

function plotActI(ActI_df, ptitle = ""; legend = true)
    setGadflyTheme()
    ActI_df."Receptor" = replace.(ActI_df."Receptor", "FcgR" => "FcγR")
    pl = plot(
        ActI_df,
        x = "Receptor",
        y = "Weight",
        ymin = "ymin",
        ymax = "ymax",
        Geom.errorbar,
        Stat.dodge(axis = :x),
        Geom.bar(position = :dodge),
        Scale.x_discrete(levels = unique(ActI_df.Receptor)),
        Scale.y_continuous(),
        Guide.title("Predicted receptor weights\n($ptitle)"),
        style(stroke_color = c -> "black", errorbar_cap_length = 4px, bar_spacing = 10px),
    )
    return pl
end
