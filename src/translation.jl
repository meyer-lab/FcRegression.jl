import Statistics

function vec_comb(length, rsum, resid)
    if length <= 1
        return fill(rsum, (1, 1))
    end
    enum = Array{Int64}(undef, length, 0)
    for i = max(0, rsum - sum(resid[2:end])):min(resid[1], rsum)
        other_ele = vec_comb(length - 1, rsum - i, resid[2:end])
        first_ele = ones(Int64, 1, size(other_ele, 2)) * i
        enum = hcat(enum, vcat(first_ele, other_ele))
    end
    return enum
end

function mat_comb(nrow, ncol, rsum, cresid)
    if nrow <= 1
        return reshape(cresid, (1, ncol, 1))
    end
    enum = Array{Int64}(undef, nrow, ncol, 0)
    first_rows = vec_comb(ncol, rsum, cresid)
    for first_row in eachcol(first_rows)
        other_row = mat_comb(nrow - 1, ncol, rsum, cresid - first_row)
        first_row = repeat(reshape(first_row, 1, ncol), 1, 1, size(other_row, 3))
        enum = cat(enum, vcat(first_row, other_row), dims = 3)
    end
    return enum
end

function gen_IgGC(length, divisor = 5)
    return vec_comb(length, divisor, repeat([divisor], length)) ./ divisor
end

function gen_discrete_Omega(side, divisor = 3)
    return mat_comb(side, side, divisor, repeat([divisor], side)) ./ divisor
end


function murine2human(Omega, mR, mKav, hR, hKav; L0 = 1e-9, f = 4)
    murineIgGC = gen_IgGC(4)
    murine = polyfc_ActV(L0, KxConst, f, mR, murineIgGC, mKav, murineActI)
    human = polyfc_ActV(L0, KxConst, f, hR, Omega * murineIgGC, hKav, humanActI)
    return murine, human
end


function J_pearson(Omega, mR, mKav, hR, hKav; L0 = 1e-9, f = 4)
    murine, human = murine2human(Omega, mR, mKav, hR, hKav; L0 = L0, f = f)
    pcor = [Statistics.cor(murine[i, :], human[i, :]) for i = 1:size(murine, 1)]
    pcor[isnan.(pcor)] .= 0
    return Statistics.mean(pcor)
end


function brute_force_discrete(divisor = 3)
    mR = importRtot(murine = true)
    mKav = importKav(murine = true, IgG2bFucose = true)
    hR = importRtot(murine = false)
    hKav = importKav(murine = false)
    Omegas = gen_discrete_Omega(4, divisor)

    pcors = zeros(size(Omegas, 3))
    Threads.@threads for i = 1:size(Omegas, 3)
        pcors[i] = J_pearson(Omegas[:, :, i], mR, mKav, hR, hKav)
    end
    return Omegas, pcors
end

function brute_force_optimum(divisor = 3)
    Omegas, pcors = brute_force_discrete(divisor)
    mxval, mxind = findmax(pcors)
    return mxval, Omegas[:, :, mxind]
end

function plot_pearson(divisor = 3)
    _, pcors = brute_force_discrete(divisor)
    return plot(
        x = pcors,
        Coord.Cartesian(xmin = -0.6, xmax = 0.6),
        Guide.xlabel("Average Pearson correlations of predicted cell type responses"),
        Geom.histogram(bincount = 8),
    )
end

function plot_translation(Omega; L0 = 1e-9, f = 4)
    mR = importRtot(murine = true)
    mKav = importKav(murine = true, IgG2bFucose = true)
    hR = importRtot(murine = false)
    hKav = importKav(murine = false)

    murine, human = murine2human(Omega, mR, mKav, hR, hKav; L0 = L0, f = f)
    murine = rename!(DataFrame(murine'), cellTypes)
    human = rename!(DataFrame(human'), cellTypes)
    comb = stack(murine)
    rename!(comb, :value => :murine)
    rename!(comb, :variable => :celltypes)
    comb[!, :human] .= stack(human).value

    pls = Vector(undef, length(cellTypes))
    for i = 1:length(cellTypes)
        pls[i] = plot(comb[comb[!, :celltypes] .== cellTypes[i], :], x = :murine, y = :human, Guide.title("$(cellTypes[i])"))
    end
    pl = hstack(pls...)
    title(pl, "Pearson correlation of predicted responses between human and murine")
    return pl
end
