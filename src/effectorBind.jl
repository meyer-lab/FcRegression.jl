function importEffectorBind(; avg = false)
    df = CSV.File(joinpath(dataDir, "dec2022_h7B4-IC_hPBL_binding.csv"), comment = "#") |> DataFrame
    df = stack(df, Not(["Valency", "Subclass", "Cell"]), variable_name = "Experiment", value_name = "Value")
    df = dropmissing(df)
    df[!, "Value"] = convert.(Float64, df[!, "Value"])

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

    for igg in unique(df."IgG")
        kav = Matrix(Kav[Kav."IgG" .== igg, 2:end])
        for cn in unique(df."Cell")
            for fi = 1:length(fs)
                res = polyfc(L0, KxStar, fs[fi], Rtot[!, cn], [1.0], kav)
                res = specificRcp ? res.Rmulti_n : res.Lbound
                if "Lbound" in names(df)
                    df[(df."IgG" .== igg) .& (df."Cell" .== cn) .& (df."Valency" .== vals[fi]), "Lbound"] .= res
                else
                    df[!, "Lbound"] .= res
                end
            end
        end
    end
    return df
end
