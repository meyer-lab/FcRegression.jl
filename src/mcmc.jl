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


""" Making one single violin plot. Called by plotAffinityViolin() """
function plot_distribution_violins(df::AbstractDataFrame, dist_list::Vector{T}; title = "", x_name = "Receptor", y_range = (5, 8)) where {T <: Distribution}
    setGadflyTheme()
    @assert size(df)[2] == length(dist_list)

    odf = DataFrame([names(df)[i] => rand(dist_list[i], 100_000) for i = 1:size(df)[2] if dist_list[i].μ > 0])
    df = stack(df, variable_name = x_name, value_name = "Value")
    odf = stack(odf, variable_name = x_name, value_name = "Value")

    return plot(
        layer(df, x = x_name, y = "Value", Geom.violin, Theme(default_color = colorant"firebrick4")),
        layer(odf, x = x_name, y = "Value", Geom.violin, Theme(default_color = colorant"navajowhite2")),
        Coord.cartesian(ymin = y_range[1], ymax = y_range[2]),
        Scale.y_log10,
        Guide.manual_color_key("Legend", ["Prior", "Posterior"], ["navajowhite2", "firebrick4"]),
        Guide.ylabel("<i>K</i><sub>a</sub> (M<sup>-1</sup>)"),
        Guide.title(title),
    )
end

""" Make violin plot of affinities from a chain """
function plotAffinityViolin(c::Chains; murine::Bool)
    pref = murine ? "m" : "h"
    y_range = murine ? (4, 9) : (5, 8)
    Kav_priors = importKavDist(; murine = murine, retdf = true)

    len = length(Matrix(Kav_priors[!, Not("IgG")]))
    Kav_posts = deepcopy(Kav_priors)
    Kav = [c["Kav[$i]"].data for i = 1:len]
    Kav_posts[!, Not("IgG")] = typeof(Kav[1, 1]).(reshape(Kav, size(Kav_posts)[1], :))
    
    pls = Vector{Union{Gadfly.Plot, Context}}(undef, size(Kav_priors)[1])
    for (ii, igg) in enumerate(Kav_priors[!, "IgG"])
        priors = reshape(Matrix(Kav_priors[Kav_priors."IgG" .== igg, Not("IgG")]), :)
        posts = DataFrame(hcat([reshape(Kav_posts[Kav_posts."IgG" .== igg, Not("IgG")][1, i], :) for i = 1:(size(Kav_posts)[2]-1)]...), names(Kav_posts)[2:end])
        pls[ii] = plot_distribution_violins(posts, priors; y_range = y_range, title = "$pref$igg Affinities Distributions")
    end
    return pls
end


""" Making a single subplot for priors and posteriors. Called by plotMCMCdists() """
function plot_dist_histogram(dat::Array{Float64}, dist::Distribution, name::String = ""; bincount = 20)
    dat = reshape(dat, :)
    xxs = exp.(LinRange(dist.μ - 4 * dist.σ, dist.μ + 4 * dist.σ, 100))
    yys = [pdf(dist, xx) for xx in xxs]
    yys = yys ./ maximum(yys)
    pl = plot(
        layer(x = xxs, y = yys, Geom.line, color = [colorant"red"], order = 1),
        layer(DataFrame("Value" => dat), x = "Value", Geom.histogram(bincount = bincount, density = true)),
        Scale.x_log10(),
        Guide.xticks(orientation = :horizontal),
        Guide.xlabel("Value"),
        Guide.ylabel(nothing),
        Guide.title(name),
    )
    return pl
end


""" Plot prior and posterior distributions from a chain. Will generate three figure files """
function plotMCMCdists(c::Chains, fname::String = ""; murine::Bool)
    setGadflyTheme()
    pnames = [String(s) for s in c.name_map[1]]

    local Kav_dist, IgGs
    Kav_dist = importKavDist(; murine = murine, retdf = false)
    IgGs = importKavDist(; murine = murine, retdf = true)."IgG"
    FcgRs = murine ? murineFcgR : humanFcgRiv
    pref = murine ? "m" : "h"
    ligg, lfcr = length(IgGs), length(FcgRs)
    
    # Plot Kav's
    if "Kav[1]" in pnames
        Kav_pls = Matrix{Plot}(undef, ligg, lfcr)
        for ii in eachindex(Kav_pls)
            IgGname = IgGs[(ii - 1) % ligg + 1]
            FcRname = FcgRs[(ii - 1) ÷ ligg + 1]
            name = "$pref$IgGname to $pref$FcRname"
            Kav_pls[ii] = plot_dist_histogram(c["Kav[$ii]"].data, Kav_dist[ii], name)
        end
        Kav_plot = plotGrid((ligg, lfcr), permutedims(Kav_pls, (2, 1)); sublabels = false)
        draw(PDF("MCMC_Kav_$fname.pdf", 11inch, 8inch), Kav_plot)
    end

    # Plot Rtot's
    ## TODO: simplify this
    if murine
        cellTypes = unique(importMurineLeukocyte(; average = true)."ImCell")
        lcell = length(cellTypes)
        Rtotd = importCellRtotDist(; retdf = true)
        Rtotd = Matrix(Rtotd[!, names(Rtotd)[in(cellTypes).(names(Rtotd))]])
        Rtot_pls = Matrix{Plot}(undef, lcell, lfcr)
        for ii in eachindex(Rtot_pls)
            cellname = cellTypes[(ii - 1) % lcell + 1]
            FcRname = murineFcgR[(ii - 1) ÷ lcell + 1]
            name = "m$FcRname on $cellname"
            Rtot_pls[ii] = plot_dist_histogram(c["Rtot[$ii]"].data, Rtotd[ii], name)
        end
        Rtot_plot = plotGrid((lcell, lfcr), permutedims(Rtot_pls, (2, 1)); sublabels = false)
        draw(PDF("MCMC_Rtot_$fname.pdf", 11inch, 14inch), Rtot_plot)
    else # human
        Rtot_pls = Vector{Plot}(undef, lfcr)
        Rtot_dist = importInVitroRtotDist()
        for ii in eachindex(Rtot_pls)
            FcRname = humanFcgRiv[ii]
            Rtot_pls[ii] = plot_dist_histogram(c["Rtot[$ii]"].data, Rtot_dist[ii], FcRname)
        end
        Rtot_plot = plotGrid((1, lfcr), Rtot_pls; sublabels = false)
        draw(PDF("MCMC_Rtot_$fname.pdf", 16inch, 4inch), Rtot_plot)
    end

    # Plot f4, f33, KxStar
    ## TODO: add case where not all three parameters appear
    other_pls = Vector{Plot}(undef, 3)
    other_pls[1] = plot_dist_histogram(c["f4"].data, f4Dist, "f = 4 effective valency")
    other_pls[2] = plot_dist_histogram(c["f33"].data, f33Dist, "f = 33 effective valency")
    other_pls[3] = plot_dist_histogram(c["KxStar"].data, KxStarDist, "KxStar")
    other_plot = plotGrid((1, 3), other_pls; sublabels = false)
    draw(PDF("MCMC_others_$fname.pdf", 9inch, 3inch), other_plot)
end


function plotMCMCPredict(c::Chains, df::AbstractDataFrame; murine::Bool, CHO = false,
        Kav::Union{Nothing, AbstractDataFrame} = nothing, kwargs...)
    p = extractMCMC(c; murine = murine, CHO = CHO)
    Kavd = (Kav !== nothing) ? Kav : p["Kav"]
    if (p["Kav"] !== nothing) && (Kav isa AbstractDataFrame)
        @warn "Using provided Kav, even though Kav's are fitted in MCMC chain"
    end

    if !("xmin" in names(df))
        df = averageMixData(df)
    end
    ndf = predMix(df; Kav = Kavd, KxStar = p["KxStar"], Rtot = p["Rtot"], fs = [p["f4"], p["f33"]])
    return plotPredvsMeasured(ndf; xx = "Value", yy = "Predict", color = ((murine && (!CHO)) ? "ImCell" : "Receptor"), 
        shape = (murine ? "Subclass" : "Valency"), kwargs...)
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