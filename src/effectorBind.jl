function importEffectorBind(; avg = false)
    df = CSV.File(joinpath(dataDir, "human_donor_leukocyte_binding.csv"), comment = "#") |> DataFrame
    # remove days with many negative values
    exp_cols = Vector{Bool}([1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1])
    df = df[!, exp_cols]
    df = stack(df, Not(["Valency", "Subclass", "Cell"]), variable_name = "Experiment", value_name = "Value")
    df = dropmissing(df)
    df[!, "Value"] = convert.(Float64, df[!, "Value"])

    df[(df[!, "Value"]) .< 1.0, "Value"] .= 1.0
    baseline = combine(groupby(df, "Experiment"), "Value" => geomean => "Baseline")
    df = innerjoin(df, baseline, on = "Experiment")
    df[!, "Value"] ./= df[!, "Baseline"]    # normalize fluorescence by daily geomean
    df = df[!, Not("Baseline")]
    df = df[!, Not("Experiment")]

    if avg
        df = combine(groupby(df, Not("Value")), "Value" => StatsBase.median => "Value", "Value" => lower => "xmin", "Value" => upper => "xmax")
    end
    return df
end


function predictLbound(
    Kav = extractNewHumanKav(),
    Rtot = importRtot(; murine = false, retdf = true, genotype = "ZIZ");
    specificRcp = false,
    L0 = 1e-9,
    fs = [4, 33],
    KxStar = KxConst,
)
    @assert size(Kav[!, Not("IgG")])[2] == size(Rtot)[1]
    vals = [(f > 11) ? 33 : 4 for f in fs]

    """ Predict Lbound of each cell type based on Kav """
    df = if specificRcp
        DataFrame((IgG = x, Cell = y, Valency = z, Receptor = w) for x in Kav."IgG" for y in names(Rtot)[2:end] for z in vals for w in Rtot."Receptor")
    else
        DataFrame((IgG = x, Cell = y, Valency = z) for x in Kav."IgG" for y in names(Rtot)[2:end] for z in vals)
    end

    colname = specificRcp ? "Rmulti" : "Lbound"
    for igg in unique(df."IgG")
        kav = Matrix(Kav[Kav."IgG" .== igg, 2:end])
        for cn in unique(df."Cell")
            for fi = 1:length(fs)
                res = polyfc(L0, KxStar, fs[fi], Rtot[!, cn], [1.0], kav)
                res = specificRcp ? res.Rmulti_n : res.Lbound
                if !(colname in names(df))
                    df[!, colname] .= specificRcp ? res[1] : res
                end
                df[(df."IgG" .== igg) .& (df."Cell" .== cn) .& (df."Valency" .== vals[fi]), colname] .= res
            end
        end
    end
    return df
end


function plotEffectorPredict(measured = importEffectorBind(; avg = false), pred = predictLbound(); title = nothing, minvalues = 1e-2)
    setGadflyTheme()

    df = innerjoin(measured, pred, on = ["Valency" => "Valency", "Subclass" => "IgG", "Cell" => "Cell"])
    df[df."Lbound" .< minvalues, "Lbound"] .= minvalues
    df[df."Value" .< minvalues, "Value"] .= minvalues
    if "xmin" in names(df)
        df[df."xmin" .< minvalues, "xmin"] .= minvalues
        df[df."xmax" .< minvalues, "xmax"] .= minvalues
    end

    r2 = R2((df[!, "Value"]), (df[!, "Lbound"]); logscale = true)
    println("R2 = $r2")
    linek, lineb = bestfitline(df."Value", df."Lbound"; logscale = true)

    return plot(
        df,
        x = "Value",
        y = "Lbound",
        color = "Subclass",
        shape = "Cell",
        Geom.point,
        intercept = [10^lineb],
        slope = [linek],
        Geom.abline(color = "black"),
        Scale.x_log10,
        Scale.y_log10,
        Scale.color_discrete_manual(colorSubclass...),
        xmin = ("xmin" in names(df) ? "xmin" : "Value"),
        xmax = ("xmax" in names(df) ? "xmax" : "Value"),
        "xmin" in names(df) ? Geom.errorbar : Guide.xlabel(nothing),
        Guide.xlabel("Measured binding"),
        Guide.ylabel("Predicted binding"),
        Guide.title(title),
        Guide.annotation(
            compose(context(), text(-2, 5, "<i>R</i><sup>2</sup> = " * @sprintf("%.4f", r2)), stroke("black"), fill("black"), font("Helvetica-Bold")),
        ),
    )
end
