function modelPred(dfr::DataFrameRow; f = 4, ActI::NamedVector = humanActI, Kav::DataFrame, Rtot = importRtot(; murine = false, retdf = true))
    Kav[Kav."IgG" .== "IgG2c", "IgG"] .= "IgG2a"

    IgGs = String[]
    cps = [1.0]
    if any(startswith.(names(dfr), "subclass"))
        iggv = Vector(dfr[startswith.(names(dfr), "subclass")])
        cps = Vector(dfr[startswith.(names(dfr), "%_")]) .* 1.0
        cps = cps[iggv .!= "PBS"]
        cps ./= sum(cps)
        append!(IgGs, iggv[iggv .!= "PBS"])
    elseif "Condition" in names(dfr)
        IgGs = [dfr."Condition"]
    elseif "Subclass" in names(dfr)
        IgGs = [dfr."Subclass"]
    end
    @assert length(IgGs) > 0
    IgGs[IgGs .== "IgG2c"] .= "IgG2a"
    Kav = Kav[in(IgGs).(Kav."IgG"), :]
    @assert size(Kav)[1] > 0

    if ("Background" in names(dfr)) && (dfr."Background" != "wt")
        rr_names = ["FcgRI", "FcgRIIB", "FcgRIII", "FcgRIV", ["FcgRI", "FcgRIIB", "FcgRIII", "FcgRIV"], ["FcgRI", "FcgRIII"], ["FcgRI", "FcgRIV"]]
        for (ii, rr) in enumerate(["R1", "R2", "R3", "R4", "gc", "R1/3KO", "R1/4KO"])
            if occursin(rr, dfr."Background")
                Kav[:, rr_names[ii]] .= 0.0
            end
        end
        for cname in ["Neu", "ncMO", "cMO", "EO"]
            if occursin(cname, dfr."Background")
                Rtot[:, cname .== names(Rtot)] .= 0.0
            end
        end
    end
    if ("Condition" in names(dfr)) && (dfr."Condition" == "IgG1D265A")
        Kav[:, ["FcgRI", "FcgRIIB", "FcgRIII", "FcgRIV"]] .= 0.0
    end
    Kav = Matrix(Kav[:, Not("IgG")])
    pred = NamedArray(repeat([0.0], length(Rtot."Receptor")), Rtot."Receptor", "Receptor")
    Rtot = Matrix(Rtot[:, Not("Receptor")])

    cellActs = Vector(undef, size(Rtot)[2])
    if length(cps) <= 0
        cellActs .= 0.0
        return cellActs
    end

    #println("Rtot", Rtot, dfr)
    for jj = 1:size(Rtot)[2]
        pred[:] = polyfc(dfr."Concentration", KxConst, f, Rtot[:, jj], cps, Kav).Rmulti_n
        cellActs[jj] = maximum([assembleActs(ActI, pred), 0.0])
    end
    return cellActs
end

function modelPred(df::DataFrame; L0 = 1e-9, murine::Bool = false, cellTypes = nothing, ActI = nothing, Kav::DataFrame, kwargs...)
    df = deepcopy(df)
    if ActI === nothing
        ActI = murine ? murineActI : humanActI
    end
    if cellTypes === nothing
        cellTypes = murine ? murineCellTypes : humanCellTypes
    end
    if "Concentration" in names(df)
        df[!, "Concentration"] ./= maximum(df[!, "Concentration"])
        df[!, "Concentration"] .*= L0
    else
        insertcols!(df, 3, "Concentration" => L0)
    end

    ansType = ("Target" in names(df)) ? promote_type(eltype(df."Target"), eltype(ActI)) : eltype(ActI)
    Xfc = Array{ansType}(undef, size(df, 1), length(cellTypes))
    Threads.@threads for k = 1:size(df, 1)
        genotype = ("Genotype" in names(df)) ? df[k, "Genotype"] : "XXX"
        genotype = genotype[1] * "I" * genotype[3:end]    # FcgRIIB has default genotype is 232I
        Rtott = importRtot(; murine = murine, retdf = true, cellTypes = cellTypes, genotype = genotype)
        #println("Rtott before", Rtott)
        #println("Kav", Kav)

        Rtott = Rtott[in(names(Kav[!, Not("IgG")])).(Rtott."Receptor"), :]

        #println("Rtott after", Rtott)
        Xfc[k, :] = modelPred(df[k, :]; ActI = ActI, Kav = Kav, Rtot = Rtott, kwargs...)
    end

    colls = murine ? murineFcgR : humanFcgR
    colls = vcat(["Target", "Concentration", "Baseline", "Measurement"], colls, cellTypes)
    ldf = df[!, (!).(in(colls).(names(df)))]
    Xdf = hcat(ldf, DataFrame(Xfc, cellTypes))
    if "C1q" in names(df)
        Xdf[!, "C1q"] = df[!, "C1q"] .* df[!, "Concentration"]
    end
    return Xdf
end


function regPred(df::DataFrame, opt::regParams; link::Function = exponential, kwargs...)
    cellTypes = names(opt.cellWs)[1]
    Xdf = modelPred(df; ActI = opt.ActIs, cellTypes = cellTypes, murine = opt.isMurine, kwargs...)
    Xmat = Matrix(Xdf[!, in(cellTypes).(names(Xdf))])
    return link(Xmat * opt.cellWs)
end


@model function regmodel(df::DataFrame, targets; murine::Bool = false, L0 = 1e-9, f = 4, cellTypes = nothing, Kav::DataFrame)
    ActI_means = murine ? murineActI : humanActI
    ActI_means = ActI_means[unique([receptorType(ss) for ss in names(Kav[!, Not("IgG")])])]
    ActIs = NamedArray(repeat([0.0], length(ActI_means)), names(ActI_means)[1], ("Receptor"))
    for ii in eachindex(ActI_means)
        ActIs[ii] ~ Normal(ActI_means[ii], 1.0)
    end

    if cellTypes === nothing
        cellTypes = murine ? murineCellTypes : humanCellTypes
    end
    cellWs = NamedArray(repeat([0.0], length(cellTypes)), cellTypes, ("CellType"))
    for ii in eachindex(cellTypes)
        cellWs[ii] ~ Exponential(Float64(length(cellTypes)))
    end

    Xdf = modelPred(df; L0 = L0, f = f, murine = murine, cellTypes = cellTypes, ActI = ActIs, Kav = Kav)
    Yfit = regPred(Xdf, regParams(cellWs, ActIs, murine); Kav = Kav)

    stdv = std(Yfit - targets) / 10
    targets ~ MvNormal(Yfit, stdv * I)
    nothing
end
