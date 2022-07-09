""" Figure 5: Predicted effector cell binding """

function EffectorBinding(IgGs::BitVector = BitVector([0,1,0,1]))
    Rtot = importRtot(; murine = false, retdf = true, genotype = "HIV")
    
    Kav_old = importKav(; murine = false, retdf = true)
    @assert length(IgGs) == length(Kav_old."IgG")
    @assert sum(IgGs) == 2
    Kav_old = Kav_old[IgGs, in(["IgG", Rtot."Receptor"...]).(names(Kav_old))]
    Kav_new = extractNewHumanKav(; replace = true)
    Kav_new = Kav_new[IgGs, in(["IgG", Rtot."Receptor"...]).(names(Kav_new))]

    IgGC = zeros(2, 101)
    IgGC[1, :] = range(0, 1, 101)
    IgGC[2, :] = range(1, 0, 101)
    IgGnames = Kav_old."IgG"

    Rmulti_old = polyfc_ActV(1e-9, KxConst, 4, Matrix(Rtot[!, Not("Receptor")]), IgGC, Matrix(Kav_old[!, Not("IgG")]), false; Mix = true)
    Rmulti_new = polyfc_ActV(1e-9, KxConst, 4, Matrix(Rtot[!, Not("Receptor")]), IgGC, Matrix(Kav_new[!, Not("IgG")]), false; Mix = true)

    pls = Vector(undef, size(Rmulti_old, 1))
    setGadflyTheme()
    for ii = 1:size(Rmulti_old, 1)
        receps = Rtot[Rtot[!, ii+1] .> 2.0, "Receptor"]
        olddf = DataFrame(transpose(Rmulti_old[ii, Rtot[!, ii+1] .> 2.0, :]), receps)
        olddf[!, IgGnames[1]] = IgGC[1, :]
        olddf[!, IgGnames[2]] = IgGC[2, :]
        olddf = stack(olddf, receps, variable_name = "Receptor", value_name = "Old")

        newdf = DataFrame(transpose(Rmulti_new[ii, Rtot[!, ii+1] .> 2.0, :]), receps)
        newdf[!, IgGnames[1]] = IgGC[1, :]
        newdf = stack(newdf, receps, variable_name = "Receptor", value_name = "New")

        df = innerjoin(olddf, newdf, on = [IgGnames[1], "Receptor"])
        #df = stack(df, ["Old", "New"], variable_name = "Affinities", value_name = "Prediction")

        pls[ii] = plot(
            df,
            layer(x = IgGnames[1], y = "Old", Geom.line, style(line_style = [:dash])),
            layer(x = IgGnames[1], y = "New", Geom.line, style(line_style = [:solid])), 
            color = "Receptor",
            Scale.x_continuous(labels = n -> "$(IgGnames[1]) $(n*100)%\n$(IgGnames[2]) $(100-n*100)%"),
            Scale.y_log10,
            Guide.xlabel("$(IgGnames[1]) to $(IgGnames[2]) mixture"),
            Guide.ylabel("Predicted crosslinks"),
            Guide.title("$(names(Rtot)[ii+1]) Rmulti, $(IgGnames[1])-$(IgGnames[2]) mixture"),
        )
    end
    return pls
end


function figure5()
    pls = Matrix(undef, 5, 6)
    IgGpairs = [BitVector(bv) for bv in [[1,1,0,0], [1,0,1,0], [1,0,0,1], [0,1,1,0], [0,1,0,1], [0,0,1,1]]]
    for ii = 1:6
        pls[:, ii] = EffectorBinding(IgGpairs[ii])
    end
    pl = plotGrid((6, 5), pls; sublabels = false)
    draw(PDF("figure5.pdf", 16inch, 18inch), pl)
end