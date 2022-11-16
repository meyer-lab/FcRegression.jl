""" Figure 1: Explain mixture binding experiment and explore data """
function plotDFwithGreekGamma(df::DataFrame)
    df = deepcopy(df)
    df."Receptor" = replace.(df."Receptor", "FcgR" => "FcγR")
    return df
end

""" Original measurements with middle 50% as error bar """
function splot_origData(df; match_y = true, legend = true)
    cell = unique(df."Receptor")[1]
    IgGX = unique(df."subclass_1")[1]
    IgGY = unique(df."subclass_2")[1]
    palette = [Scale.color_discrete().f(3)[1], Scale.color_discrete().f(3)[3]]
    cell_name = replace(cell, "FcgR" => "FcγR")

    ymax = Dict(
        "FcgRI" => 6,
        "FcgRIIA-131H" => 20,
        "FcgRIIA-131R" => 15,
        "FcgRIIB-232I" => 8,
        "FcgRIIIA-158F" => 20,
        "FcgRIIIA-158V" => 15,
    )
    return plot(
        df,
        x = "%_1",
        y = "Value",
        ymin = "xmin",
        ymax = "xmax",
        color = "Valency",
        Geom.point,
        Geom.line,
        Geom.errorbar,
        Scale.x_continuous(labels = n -> "$IgGX $(n*100)%\n$IgGY $(100-n*100)%"),
        Scale.y_continuous(; minvalue = 0.0, maxvalue = match_y ? ymax[cell] : maximum(df."xmax")),
        Scale.color_discrete_manual(palette[1], palette[2]),
        Guide.xlabel("", orientation = :horizontal),
        Guide.ylabel("RFU", orientation = :vertical),
        Guide.xticks(orientation = :horizontal),
        Guide.title("$IgGX-$IgGY bind to $cell_name"),
        style(key_position = legend ? :right : :none),
    )
end

function bindVSaff(hKav = importKav(; murine = false, retdf = true); affinity_name = "Documented")
    setGadflyTheme()

    # Binding data, keep single IgG subclass only
    df = loadMixData()
    df = df[(df."%_1" .== 1.0) .| (df."%_2" .== 1.0), :]
    df."Subclass" = [r."%_1" >= 1 ? r."subclass_1" : r."subclass_2" for r in eachrow(df)]
    df = df[!, ["Valency", "Receptor", "Subclass", "Value"]]
    df = combine(
        groupby(df, ["Valency", "Receptor", "Subclass"]),
        "Value" => StatsBase.median => "Value",
        "Value" => geocmean => "Geomean",
        "Value" => lower => "xmin",
        "Value" => upper => "xmax",
    )
    df."Affinity" = [hKav[hKav."IgG" .== r."Subclass", r."Receptor"][1] for r in eachrow(df)]
    df[df."Affinity" .< 1e3, "Affinity"] .= 1e3
    df[!, "Valency"] .= Symbol.(df[!, "Valency"])
    pearson_cor = cor(log.(df."Affinity"), log.(df."Value"))

    pl1 = plot(
        plotDFwithGreekGamma(df),
        x = "Affinity",
        y = "Value",
        ymin = "xmin",
        ymax = "xmax",
        color = "Receptor",
        shape = "Subclass",
        Geom.point,
        Geom.errorbar,
        Scale.x_log10,
        Scale.y_log10,
        Guide.title("$affinity_name affinity vs. single IgG binding"),
        Guide.xlabel("$affinity_name Affinity (M<sup>-1</sup>)"),
        Guide.ylabel("Binding quantification"),
        Guide.annotation(
            compose(context(), text(6, -3, "<i>ρ</i> = " * @sprintf("%.4f", pearson_cor)), stroke("black"), fill("black"), font("Helvetica-Bold")),
        ),
        style(key_position = :none),
    )

    val_ratio = combine(groupby(df, ["Receptor", "Subclass"])) do df
        (Ratio = df[df."Valency" .== Symbol("33"), "Value"][1] / df[df."Valency" .== Symbol("4"), "Value"][1],)
    end

    val_ratio."Affinity" = [hKav[hKav."IgG" .== r."Subclass", r."Receptor"][1] for r in eachrow(val_ratio)]
    val_ratio."Affinity"[val_ratio."Affinity" .< 1000] .= 1000
    ratio_cor = cor(log.(val_ratio."Affinity"), log.(val_ratio."Ratio"))

    pl2 = plot(
        plotDFwithGreekGamma(val_ratio),
        x = "Affinity",
        y = "Ratio",
        color = "Receptor",
        shape = "Subclass",
        Geom.point,
        Scale.x_log10,
        Scale.y_log10,
        Guide.title("$affinity_name affinity vs. intervalency ratio"),
        Guide.xlabel("$affinity_name Affinity (M<sup>-1</sup>)"),
        Guide.ylabel("33- to 4-valent binding ratio"),
        Guide.annotation(
            compose(context(), text(6.0, 1.5, "<i>ρ</i> = " * @sprintf("%.4f", ratio_cor)), stroke("black"), fill("black"), font("Helvetica-Bold")),
        ),
    )
    return pl1, pl2
end

#=
function mixtureANOVA()
    ## GLM is not included in the package. Results will be saved separately after run
    import GLM: lm, formula, ftest

    df = loadMixData()
    df."Condition" = string.(df."Valency") .* df."Receptor" .* df."subclass_1" .* " " .* 
            string.(df."%_1") .* ", " .* df."subclass_2" .* " " .* string.(df."%_2")    
    df."Dose" = df."subclass_1" .* " " .* string.(df."%_1") .* ", " .* df."subclass_2" .* 
            " " .* string.(df."%_2")

    nullmodel = lm(@formula(Value ~ 1), df)
    model = lm(@formula(Value ~ 1 + Condition), df)
    return ftest(nullmodel.model, model.model)

    vmodel = lm(@formula(Value ~ 1 + Valency), df)
    vrmodel = lm(@formula(Value ~ 1 + Valency + Receptor), df)
    vrdmodel = lm(@formula(Value ~ 1 + Valency + Receptor + Dose), df)
    ftest(nullmodel.model, vmodel.model, vrmodel.model, vrdmodel.model)
end
=#

function figure1(; kwargs...)
    setGadflyTheme()

    p1, p2 = bindVSaff()

    # Specific IgG pair - receptor interaction
    df = averageMixData()
    igg12_1 = splot_origData(
        df[(df."Receptor" .== "FcgRI") .& (df."subclass_1" .== "IgG1") .& (df."subclass_2" .== "IgG2"), :];
        match_y = false,
        legend = false,
    )
    igg14_1 = splot_origData(
        df[(df."Receptor" .== "FcgRIIIA-158F") .& (df."subclass_1" .== "IgG1") .& (df."subclass_2" .== "IgG4"), :];
        match_y = false,
        legend = true,
    )

    pl = plotGrid((2, 3), [nothing, p1, p2, nothing, igg12_1, igg14_1]; sublabels = "acdbef", widths = [4, 3, 3.8], kwargs...)
    return draw(PDF("output/figure1.pdf", 10inch, 6inch), pl)
end
