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


function human2murine(Omega, hR, hKav, mR, mKav; L0 = 1e-9, f = 4)
    humanIgGC = gen_IgGC(4)
    human = polyfc_ActV(L0, KxConst, f, hR, humanIgGC, hKav, humanActI)
    murine = polyfc_ActV(L0, KxConst, f, mR, Omega * humanIgGC, mKav, murineActI)
    return human, murine
end


function J_pearson(Omega, hR, hKav, mR, mKav; L0 = 1e-9, f = 4)
    human, murine = human2murine(Omega, hR, hKav, mR, mKav; L0 = L0, f = f)
    pcor = [Statistics.cor(human[i, :], murine[i, :]) for i = 1:size(human, 1)]
    pcor[isnan.(pcor)] .= 0
    return Statistics.mean(pcor)
end


function brute_force_discrete(divisor = 3)
    hR = importRtot(murine = false)
    hKav = importKav(murine = false)
    mR = importRtot(murine = true)
    mKav = importKav(murine = true)

    Omegas = gen_discrete_Omega(4, divisor)

    pcors = zeros(size(Omegas, 3))
    Threads.@threads for i = 1:size(Omegas, 3)
        pcors[i] = J_pearson(Omegas[:, :, i], hR, hKav, mR, mKav)
    end

    mxval, mxind = findmax(pcors)
    return (mxval, Omegas[:, :, mxind])
end

function plot_translation(Omega; L0 = 1e-9, f = 4)
    hR = importRtot(murine = false)
    hKav = importKav(murine = false)
    mR = importRtot(murine = true)
    mKav = importKav(murine = true)
    human, murine = human2murine(Omega, hR, hKav, mR, mKav; L0 = L0, f = f)
    human = rename!(DataFrame(human'), cellTypes)
    murine = rename!(DataFrame(murine'), cellTypes)
    comb = stack(human)
    rename!(comb, :value => :human)
    rename!(comb, :variable => :celltypes)
    comb[!, :murine] .= stack(murine).value
    return plot(comb, xgroup=:celltypes, x=:human, y=:murine, Geom.subplot_grid(Geom.point, free_x_axis=true))
end
