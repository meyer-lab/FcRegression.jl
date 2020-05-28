using DataFrames

function decay(x, x0, N0, thalf)
    tau = thalf / log(2)
    y = N0 .* exp.(-(x.-x0)./tau)
    y[x .< x0] .= 0.0
    return y
end

function hIgGtemp(x)
    # data from Almagro, Front. Immunol., 2018
    thalf = [21, 21, 7, 21] # days
    #N0 = [12.2/146, 3.7/146, 0.7/165, 2.4/146]
    N0 = [9/146, 3/146, 1/165, 0.5/146]
    # data estimated according to Collins, Front. Immunol., 2013
    x0 = [21.0, 22.0, 17.7, 27.1] .* 7 ./ 24
    m = zeros(length(x), length(humanIgG))
    for i in 1:length(humanIgG)
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
    for i in 1:length(x)
        IgG = Vector(df[i, humanIgG]) .* 1e-8
        if sum(IgG) > 0
            m[i, :] .= reshape(polyfc_ActV(sum(IgG), KxConst, 4,
                importRtot(murine=false), IgG/sum(IgG), importKav(murine=false), humanActI), :)
        end
    end
    df = convert(DataFrame, m)
    df = rename(df, cellTypes)
    for col = eachcol(df)
        col .= col ./ maximum(col)
    end
    df.x = x
    return df
end

function plot_hIgGtemp()
    x = [0:0.1:20;];
    IgGlvl = hIgGtemp(x)
    CTlvl = hCTtemp(x)
    pl1 = plot(IgGlvl, x=:x, y=Col.value(humanIgG...), color=Col.index(humanIgG...), Geom.line,
        Guide.xlabel("time (days)"), Guide.ylabel("Amount of IgG"), Guide.colorkey(title="isotypes"))
    pl2 = plot(CTlvl, x=:x, y=Col.value(cellTypes...), color=Col.index(cellTypes...), Geom.line,
        Guide.xlabel("time (days)"), Guide.ylabel("Effector response (au)"), Guide.colorkey(title="Effector cells"))
    draw(SVG("IgGtemp.svg"), vstack(pl1, pl2))
end
