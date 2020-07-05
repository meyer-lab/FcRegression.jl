using DataFrames
import Statistics: mean, std

function plotTemporal()
    setGadflyTheme()
    plot_hIgGtemp()
    plot_Murphy1(1e-11, 4)
    plot_Murphy2(1e-11, 4)
end

function decay(x, x0, N0, thalf)
    tau = thalf / log(2)
    y = N0 .* exp.(-(x .- x0) ./ tau)
    y[x .< x0] .= 0.0
    return y
end

function hIgGtemp(x)
    # data from Almagro, Front. Immunol., 2018
    thalf = [21, 21, 7, 21] # days
    #N0 = [12.2/146, 3.7/146, 0.7/165, 2.4/146]
    N0 = [9 / 146, 3 / 146, 1 / 165, 0.5 / 146]
    # data estimated according to Collins, Front. Immunol., 2013
    x0 = [21.0, 22.0, 17.7, 27.1] .* 7 ./ 24
    m = zeros(length(x), length(humanIgG))
    for i = 1:length(humanIgG)
        m[:, i] .= decay(x, x0[i], N0[i], thalf[i])
    end
    df = convert(DataFrame, m)
    df = rename(df, humanIgG)
    df.x = x
    return df
end

function hCTtemp(x)
    df = hIgGtemp(x)
    m = zeros(length(x), length(cellTypes))
    for i = 1:length(x)
        IgG = Vector(df[i, humanIgG]) .* 1e-8
        if sum(IgG) > 0
            m[i, :] .= reshape(polyfc_ActV(sum(IgG), KxConst, 4, importRtot(murine = false), IgG / sum(IgG), importKav(murine = false), humanActI), :)
        end
    end
    df = convert(DataFrame, m)
    df = rename(df, cellTypes)
    for col in eachcol(df)
        col .= col ./ maximum(col)
    end
    df.x = x
    return df
end

function plot_hIgGtemp()
    setGadflyTheme()
    x = [0:0.1:20;]
    IgGlvl = hIgGtemp(x)
    CTlvl = hCTtemp(x)
    pl1 = plot(
        IgGlvl,
        x = :x,
        y = Col.value(humanIgG...),
        color = Col.index(humanIgG...),
        Geom.line,
        Guide.xlabel("time (days)"),
        Guide.ylabel("Amount of IgG"),
        Guide.colorkey(title = "Isotypes"),
    )
    pl2 = plot(
        CTlvl,
        x = :x,
        y = Col.value(cellTypes...),
        color = Col.index(cellTypes...),
        Geom.line,
        Guide.xlabel("time (days)"),
        Guide.ylabel("Effector response (au)", orientation = :vertical),
        Guide.colorkey(title = "Effectors"),
    )
    draw(SVG("IgGtemp.svg", 5inch, 6inch), vstack(pl1, pl2))
end


function plot_Murphy1(L0, f)
    df = CSV.File(joinpath(dataDir, "murphy-jmedvir-2009-tab1.csv"), comment = "#") |> DataFrame!
    res = polyfcm_ActV(Matrix(df)' .* L0, KxConst, f, importRtot(murine = false, retdf = false), importKav(murine = false, retdf = false), humanActI)
    #return res
    avgs = mean(res; dims = [2])
    std = mapslices(x -> quantile(x, [0.1, 0.9]), res, dims = [2])
    pl = plot(x = cellTypes, y = avgs, ymin = std[:, 1], ymax = std[:, 2], Geom.bar, Geom.yerrorbar)
    draw(SVG("murphy1.svg", 4inch, 4inch), plotGrid((1, 1), [pl]))
end


function plot_Murphy2(L0, f)
    df =
        CSV.File(joinpath(dataDir, "murphy-jmedvir-2009-tab2.csv"), comment = "#"; types = [String, Float64, Int64, Int64, Int64, Int64]) |>
        DataFrame!
    df[!, :Subject] .= Symbol.(df[!, :Subject])
    res = polyfcm_ActV(
        Matrix(df[:, FcgR.humanIgG])' .* L0,
        KxConst,
        f,
        importRtot(murine = false, retdf = false),
        importKav(murine = false, retdf = false),
        humanActI,
    )
    for i = 1:size(res, 1)
        df[!, cellTypes[i]] = res[i, :] ./ maximum(res[i, :])
    end

    df = stack(df, cellTypes, [:Subject, :Time])
    #return df
    pls = Vector{Plot}(undef, 4)
    for (i, s) in enumerate(Symbol.(["0101", "0106", "0205", "0207"]))
        pls[i] = plot(
            df[df[!, :Subject] .== s, :],
            x = :Time,
            y = :value,
            color = :variable,
            Geom.line,
            Guide.title("Subject " * String(s)),
            Guide.colorkey(title = "Cell Types"),
        )
    end
    draw(SVG("murphy2.svg", 7inch, 7inch), plotGrid((2, 2), pls))
end
