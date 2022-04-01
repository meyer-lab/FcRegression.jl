function cornerplot(df::AbstractDataFrame; nbins = 20, hexbins = 50, plotsize = 20cm, logscale = true)
    setGadflyTheme()

    nsamps, ndims = size(df)
    @assert ndims <= nsamps
    if logscale
        ranges = [(log10(minimum(df[!, n])), log10(maximum(df[!, n]))) for n in names(df)]
    else
        ranges = [(minimum(df[!, n]), maximum(df[!, n])) for n in names(df)]
    end

    # set up the plots
    set_default_plot_size(plotsize, plotsize)
    subplots = Array{Context}(undef, ndims, ndims)

    for i = 1:ndims
        # for last diagonal, add axis label if given
        xlabel = (i == ndims) ? String(names(df)[i]) : nothing
        ylabel = (i == 1) ? String(names(df)[i])  : nothing

        # plot the diagonal histogram
        subplots[i, i] = render(plot(layer(x = df[!, i], Geom.histogram(bincount = nbins),),
                                     Guide.xlabel(xlabel),
                                     Guide.ylabel(ylabel, orientation = :vertical),
                                     (logscale ? Scale.x_log10 : Scale.x_continuous),
                                     Coord.Cartesian(xmin = ranges[i][1], xmax = ranges[i][end],),))
        
        for j = (i + 1):ndims
            xlabel = (j == ndims) ? String(names(df)[i]) : nothing
            ylabel = (i == 1) ? String(names(df)[j]) : nothing

            subplots[j, i] = render(plot(layer(x=df[!, i], y=df[!, j], Geom.density2d,),
                                        Guide.xlabel(xlabel),
                                        Guide.ylabel(ylabel, orientation = :vertical),
                                        (logscale ? Scale.x_log10 : Scale.x_continuous),
                                        (logscale ? Scale.y_log10 : Scale.y_continuous),
                                        Guide.yticks(label=(ylabel != nothing)),
                                        Coord.Cartesian(xmin = ranges[i][1], xmax = ranges[i][end],
                                                        ymin = ranges[j][1], ymax = ranges[j][end]),
                                        style(key_position = :none)))

            # make subplots above diagonal empty
            subplots[i, j] = Compose.context()

        end # for j = (i + 1):ndims
    end # for i = 1:ndims
    return gridstack(subplots)
end


function MCMC_cornerplot(c = runMCMC())
    df = DataFrame(c)
    df = df[500:1000, [contains(n, "[") | contains(n, "f") | contains(n, "Star") for n in names(df)]]
    pl = cornerplot(ndf; plotsize = 200cm)
    draw(PDF("figure2a.pdf", 200cm, 200cm), pl)
end