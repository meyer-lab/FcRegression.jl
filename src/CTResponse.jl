function regGenData(df; L0, f, murine::Bool, retdf = false)
    df = copy(df)
    FcRecep = murine ? murineFcgR : humanFcgR
    ActI = murine ? murineActI : humanActI

    if :Concentration in names(df)
        df[!, :Concentration] .*= L0
    else
        insertcols!(df, 3, :Concentration => L0)
    end

    X = DataFrame(repeat([Float64], length(cellTypes)), cellTypes)
    for i = 1:size(df, 1)
        Kav = convert(Vector{Float64}, df[i, FcRecep])
        Kav = reshape(Kav, 1, :)
        Rtot = importRtot(; murine = murine, genotype = murine ? "NA" : df[i, :Genotype])
        push!(X, polyfc_ActV(df[i, :Concentration], KxConst, f, Rtot, [1.0], Kav, ActI))
    end

    if :C1q in names(df)
        X[!, :C1q] = df[!, :C1q] .* df[!, :Concentration]
    end
    if :Neutralization in names(df)
        X[!, :Neutralization] = df[!, :Neutralization]
    end

    if :Background in names(df)
        X[df[:, :Background] .== "NeuKO", :Neu] .= 0.0
        X[df[:, :Background] .== "ncMOKO", :ncMO] .= 0.0
    end
    Y = df[!, :Target]

    @assert all(isfinite.(Matrix(X)))
    @assert all(isfinite.(Y))
    if retdf
        return (X, Y)
    else
        return (Matrix{Float64}(X), Y)
    end
end


function xx(L0, f, IgGC, Rtot, Kav, murine::Bool)
    
end




function plotDepletionSynergy(IgGXidx::Int64, IgGYidx::Int64, weights::Vector; L0, f, murine::Bool, c1q = false)
    Xname = murine ? murineIgG[IgGXidx] : humanIgG[IgGXidx]
    Yname = murine ? murineIgG[IgGYidx] : humanIgG[IgGYidx]
    Kav_df = importKav(; murine = murine, c1q = c1q, retdf = true)
    Kav = Matrix{Float64}(Kav_df[!, murine ? murineFcgR : humanFcgR])
    FcExpr = importRtot(; murine = murine)
    ActI = murine ? murineActI : humanActI

    nPoints = 100
    IgGC = zeros(Float64, size(Kav, 1), nPoints)
    IgGC[IgGXidx, :] = range(0.0, 1.0; length = nPoints)
    IgGC[IgGYidx, :] = range(1.0, 0.0; length = nPoints)
    X = polyfc_ActV(L0, KxConst, f, FcExpr, IgGC, Kav, ActI)  # size: celltype * nPoints
    if c1q
        X = vcat(X, Kav_df[!, :C1q]' * IgGC)
    end
    @assert size(X, 1) == length(weights)
    output = exponential(Matrix(X'), weights)

    pl = plot(
        layer(x = IgGC[IgGXidx, :], y = output, Geom.line, Theme(default_color = colorant"green")),
        layer(x = [0, 1], y = [output[1], output[end]], Geom.line, Theme(default_color = colorant"red")),
        Scale.x_continuous(labels = n -> "$Xname $(n*100)%\n$Yname $(100-n*100)%"),
        Scale.y_continuous(minvalue = 0.0, maxvalue = 1.0),
        Guide.ylabel("Predicted Depletion"),
        Guide.manual_color_key("", ["Predicted", "Linear Addition"], ["green", "red"]),
        Guide.title("Total predicted effects vs $Xname-$Yname Composition"),
        Theme(key_position = :inside),
    )
    return pl
end
