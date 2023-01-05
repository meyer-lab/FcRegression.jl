function importEffectorBind(; avg = false)
    df = CSV.File(joinpath(dataDir, "dec2022_h7B4-IC_hPBL_binding.csv"), comment = "#") |> DataFrame
    df = stack(df, Not(["Valency", "Subclass", "Cell"]), variable_name = "Experiment", value_name = "Value")
    df = dropmissing(df)
    df[!, "Value"] = convert.(Float64, df[!, "Value"])

    df = df[!, Not("Experiment")]

    if avg
        df = combine(
            groupby(df, Not("Value")), 
            "Value" => StatsBase.median => "Value", 
            "Value" => lower => "xmin", 
            "Value" => upper => "xmax")
    end
    return df
end


function predictLbound(
    Kav = extractNewHumanKav(),
    Rtot = importRtot(; murine = false, retdf = true);
    specificRcp = false,
    L0 = 1e-9,
    fs = [4, 33],
    KxStar = KxConst,
)
    Rtot = Rtot[in(names(Kav[!, Not("IgG")])).(Rtot."Receptor"), :]
    existRcp = all(Matrix(Rtot[!, 2:end]) .> 0.0, dims = 2)
    Kav = Kav[!, BitVector([1; existRcp...])]
    Rtot = Rtot[BitVector([existRcp...]), :]

    """ Predict Lbound of each cell type based on Kav """
    df = if specificRcp
        DataFrame((IgG=x, Cell=y, Valency=z, Receptor=w) for x in Kav."IgG" for y in names(Rtot)[2:end] for z in fs for w in Rtot."Receptor")
    else
        DataFrame((IgG=x, Cell=y, Valency=z) for x in Kav."IgG" for y in names(Rtot)[2:end] for z in fs)
    end
    df."Lbound" .= 0.0

    for igg in unique(df."IgG")
        kav = Matrix(Kav[Kav."IgG" .== igg, 2:end])
        for cn in unique(df."Cell")
            for f in unique(df."Valency")
                if specificRcp
                    # Specific receptor prediction, don't change rcp order before this
                    df[(df."IgG".==igg) .& (df."Cell" .== cn) .& (df."Valency" .== f), "Lbound"] = polyfc(L0, KxStar, f, Rtot[!, cn], [1.0], kav).Rmulti_n
                else
                    df[(df."IgG".==igg) .& (df."Cell" .== cn) .& (df."Valency" .== f), "Lbound"] .= polyfc(L0, KxStar, f, Rtot[!, cn], [1.0], kav).Lbound
                end
            end
        end
    end
    return df
end