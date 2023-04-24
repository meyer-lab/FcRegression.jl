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

function plotRegMCMC(
    c::Union{Chains, StatisticalModel, regParams},
    df::Union{DataFrame, String};
    ptitle = "",
    colorL = nothing,
    shapeL = nothing,
    Kav::DataFrame,
    legend = true,
    link::Function = exponential,
    kwargs...,
)
    if df isa String
        if ptitle === nothing
            ptitle = df
        end
        df = importDepletion(df)
    end
    if c isa Chains
        prep = regPrepareData(df; Kav = Kav, murine = extractRegMCMC(c[1]).isMurine)
        fits = hcat([regPred(prep, extractRegMCMC(c[ii]); link = link) for ii = 1:length(c)]...)
        df."Fitted" .= mapslices(median, fits, dims = 2)
        df."ymax" .= mapslices(xs -> quantile(xs, 0.75), fits, dims = 2)
        df."ymin" .= mapslices(xs -> quantile(xs, 0.25), fits, dims = 2)
    else
        if c isa StatisticalModel
            c = extractRegMCMC(c)
        end
        # c isa regParams at this point
        df."Fitted" = regPred(regPrepareData(df; Kav = Kav, c.isMurine), c; link = link)
    end

    if shapeL === nothing
        shapeL = names(df)[1]
    end
    if colorL === nothing
        colorL = names(df)[2]
    end

    setGadflyTheme()
    R2anno = "<i>R</i><sup>2</sup>" * @sprintf("=%.3f", R2(df.Target, df.Fitted; logscale = false))
    pl = plot(
        df,
        x = "Target",
        y = "Fitted",
        ymin = (c isa Chains ? "ymin" : "Fitted"),
        ymax = (c isa Chains ? "ymax" : "Fitted"),
        Geom.point,
        "ymin" in names(df) ? Geom.errorbar : Geom.point,
        color = colorL,
        shape = shapeL,
        Guide.colorkey(),
        Guide.shapekey(),
        Scale.y_continuous(minvalue = 0.0, maxvalue = 1.0),
        Scale.color_discrete_manual(
            colorant"#4b40de",
            colorant"#00a8ff",
            colorant"#00f0ff",
            colorant"#00c398",
            colorant"#008f18",
            colorant"#9fb338",
            colorant"#fed57a",
            colorant"#f18c3e",
            colorant"#de2c2c",
        ),
        Geom.abline(color = "red"),
        Guide.xlabel("Actual effect"),
        Guide.ylabel("Fitted effect"),
        Guide.title("Actual vs fitted effect\n($ptitle)"),
        Guide.annotation(compose(context(), text(0.1, 0.8, R2anno), fill("black"), fontsize(10pt), font("Helvetica"))),
        style(point_size = 5px, key_position = legend ? :right : :none),
    )
    return pl
end

function plotRegParams(
    c::Union{Chains, Vector{StatisticalModel}};
    ptitle::String = "",
    legend = true,
    retdf = false,
    Kav::DataFrame,
    cellTypes = nothing,
    cell_max_y = nothing,
)
    murine = extractRegMCMC(c[1]).isMurine
    df = vcat(
        [
            wildtypeWeights(extractRegMCMC(c[ii]); cellTypes = names(extractRegMCMC(c[1]).cellWs)[1], murine = murine, Kav = Kav) for ii = 1:length(c)
        ]...,
    )

    ActI_df = vcat([DataFrame(:Receptor => names(extractRegMCMC(c[1]).ActIs)[1], :Weight => extractRegMCMC(c[ii]).ActIs) for ii = 1:length(c)]...)

    df = combine(
        groupby(df, Not("Weight")),
        "Weight" => median => "Weight",
        "Weight" => (xs -> quantile(xs, 0.25)) => "ymin",
        "Weight" => (xs -> quantile(xs, 0.75)) => "ymax",
    )
    ActI_df = combine(
        groupby(ActI_df, Not("Weight")),
        "Weight" => median => "Weight",
        "Weight" => (xs -> quantile(xs, 0.25)) => "ymin",
        "Weight" => (xs -> quantile(xs, 0.75)) => "ymax",
    )

    if retdf
        return plotCellTypeEffects(df, ptitle; legend = legend, maxy = cell_max_y), plotActI(ActI_df, ptitle), df, ActI_df
    else
        return plotCellTypeEffects(df, ptitle; legend = legend, maxy = cell_max_y), plotActI(ActI_df, ptitle)
    end
end
