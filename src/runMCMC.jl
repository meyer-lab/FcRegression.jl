import Turing: ~, sample, NUTS, @model, Chains, rand, MAP, MLE

""" A generic function to extract fitted parameters from an MCMC """
function extractMCMC(c::Chains; murine::Bool, CHO = false)
    pnames = [String(s) for s in c.name_map[1]]
    out = Dict{String, Union{Number, Dict, DataFrame, Nothing}}()
    if "Rtot[1]" in pnames
        local Rtotd
        Rtot = [median(c["Rtot[$i]"].data) for i = 1:sum(startswith.(pnames, "Rtot"))]
        if murine
            if CHO
                Rtotd = Dict([murineFcgR[ii] => Rtot[ii] for ii = 1:length(Rtot)])
            else # Leukocyte
                cells = unique(importMurineLeukocyte()."ImCell")
                Rtotd = importCellRtotDist()
                Rtotd = Rtotd[!, ["Receptor"; names(Rtotd)[in(cells).(names(Rtotd))]]]
                Rtotd[!, Not("Receptor")] = typeof(Rtot[1, 1]).(reshape(Rtot, size(Rtotd)[1], :))
            end
        else # human
            Rtotd = Dict([humanFcgRiv[ii] => Rtot[ii] for ii = 1:length(humanFcgRiv)])
        end
        out["Rtot"] = Rtotd
    end
    if "Kav[1]" in pnames
        Kavd = importKavDist(; murine = murine, regularKav = true, retdf = true)
        Kav = [median(c["Kav[$i]"].data) for i = 1:sum(startswith.(pnames, "Kav"))]
        Kavd[!, Not("IgG")] = typeof(Kav[1, 1]).(reshape(Kav, size(Kavd)[1], :))
        out["Kav"] = Kavd
    else
        out["Kav"] = nothing
    end
    for var in ["f4", "f33", "KxStar"]
        if var in pnames
            out[var] = median(c[var].data)
        else
            out[var] = nothing
        end
    end
    return out
end


""" A general prediction function which can handle all data """
function predMix(dfr::DataFrameRow; Kav::AbstractDataFrame, Rtot = nothing, fs::Vector = [4, 33], KxStar = KxConst)
    val = dfr."Valency" < 12 ? fs[1] : fs[2]

    local recepExp
    if Rtot isa Dict
        recepExp = [Rtot[dfr."Receptor"]]
    elseif Rtot isa Vector
        recepExp = Rtot[names(Kav)[2:end] .== dfr."Receptor"]
    elseif "ImCell" in names(dfr)   # must have receptor amount already looked up
        recepExp = Vector(dfr[names(Kav)[2:end]])
    else
        @error "Failed at predMix(): cannot look up * recepExp *"
    end
    if "Expression" in names(dfr)
        recepExp .*= dfr."Expression" / 100
    end

    local aff
    ratio = [1.0]
    # affinity already added to df
    if "Affinity" in names(dfr)
        aff = reshape([dfr."Affinity"], 1, 1)
    elseif all(in(names(dfr)).(["affinity_1", "affinity_2"]))
        aff = reshape([dfr."affinity_1", dfr."affinity_2"], 2, 1)
        ratio = [dfr."%_1", dfr."%_2"]
    # look up affinity
    elseif all(in(names(dfr)).(["Subclass", "ImCell"]))  # ImCell needs all receptors
        aff = Matrix(Kav[Kav."IgG" .== dfr."Subclass", Not("IgG")])
    elseif all(in(names(dfr)).(["Subclass", "Receptor"]))  # probably slow, avoid this
        aff = reshape([Kav[Kav."IgG" .== dfr."Subclass", dfr."Receptor"][1]], 1, 1)
    elseif all(in(names(dfr)).(["subclass_1", "subclass_2", "Receptor"]))  # probably slow, avoid this
        aff = reshape([Kav[Kav."IgG" .== dfr."subclass_1", dfr."Receptor"][1], 
                       Kav[Kav."IgG" .== dfr."subclass_2", dfr."Receptor"][1]], 2, 1)
        ratio = [dfr."%_1", dfr."%_2"]
    else
        @error "Failed at predMix(): cannot look up * aff *"
    end

    res = try
        polyfc(1e-9, KxStar, val, recepExp, ratio, aff).Lbound
    catch e
        println("Failed solving at predMix():\n KxStar = $KxStar\n f = $val\n Rtot = $(recepExp)\n IgGC = $ratio\n Kav = $aff\n")
        rethrow(e)
    end
    return res
end

function predMix(df::AbstractDataFrame; Kav::AbstractDataFrame, Rtot = nothing, kwargs...)
    dft = deepcopy(df)
    if "IgG2a" in Kav."IgG"
        Kav[Kav."IgG" .== "IgG2a", "IgG"] .= "IgG2c"
    end
    if all(in(names(dft)).(["subclass_1", "subclass_2"]))
        @assert all(in(Kav."IgG").(unique([dft."subclass_1"; dft."subclass_2"]))) "Expected IgG subclass names not in df, possibly wrong murine/human setup."
        dft = innerjoin(dft, stack(Kav, Not("IgG"), variable_name = "recepp", value_name = "affinity_1"), 
            on = ["Receptor" => "recepp", "subclass_1" => "IgG"])
        dft = innerjoin(dft, stack(Kav, Not("IgG"), variable_name = "recepp", value_name = "affinity_2"), 
            on = ["Receptor" => "recepp", "subclass_2" => "IgG"])
    elseif "Subclass" in names(dft)
        @assert all(in(Kav."IgG").(unique(dft."Subclass"))) "Expected IgG subclass names not in df, possibly wrong murine/human setup."
    end

    if "ImCell" in names(dft)
        cellTypes = unique(df."ImCell")
        if !(Rtot isa AbstractDataFrame)  # default mode only works for mice
            Rtot = importRtot(; murine = true, retdf = true, cellTypes = cellTypes)
        end
        @assert all(in(cellTypes).(names(Rtot)[2:end]))
        Rtot = Rtot[!, ["Receptor"; names(Rtot)[in(cellTypes).(names(Rtot))]]]
        Rtot = stack(Rtot, Not("Receptor"), variable_name = "ImCell", value_name = "Abundance")
        Rtot = dropmissing(unstack(Rtot, "ImCell", "Receptor", "Abundance"))
        dft = innerjoin(dft, Rtot, on = "ImCell")
        @assert df."Value" == dft."Value"
    end

    preds = Vector(undef, size(dft)[1])
    Threads.@threads for i = 1:size(dft)[1]
        preds[i] = predMix(dft[i, :]; Kav = Kav, Rtot = Rtot, kwargs...)
    end
    dft."Predict" = typeof(preds[1]).(preds)

    # one conversion factor per valency
    dft[dft."Predict" .< 0.0, "Predict"] .= 0.0
    for val in unique(dft."Valency")
        if any(dft[dft."Valency" .== val, "Predict"] .> 0.0)
            rows = (dft."Valency" .== val) .& (dft."Predict" .> 0.0)
            dft[(dft."Valency" .== val), "Predict"] ./= geomean(dft[rows, "Predict"]) / geomean(dft[rows, "Value"])
        end
    end
    return dft[!, [names(df); "Predict"]]
end

