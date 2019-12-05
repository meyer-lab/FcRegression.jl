using fcBindingModel
using Optim
import Statistics.cor

function vec_comb(length, rsum, resid)
    if length <= 1
        return fill(rsum, (1, 1))
    end
    enum = Array{Int64}(undef, length, 0)
    for i = max(0, rsum-sum(resid[2:end])):min(resid[1], rsum)
        other_ele = vec_comb(length-1, rsum-i, resid[2:end])
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
        other_row = mat_comb(nrow-1, ncol, rsum, cresid-first_row)
        first_row = repeat(reshape(first_row, 1, ncol), 1, 1, size(other_row, 3))
        enum = cat(enum, vcat(first_row, other_row), dims=3)
    end
    return enum
end

function gen_IgGC(length, divisor=5)
    return vec_comb(length, divisor, repeat([divisor],length)) ./ divisor
end

function gen_discrete_Omega(side, divisor=3)
    return mat_comb(side, side, divisor, repeat([divisor], side)) ./ divisor
end


function J_pearson(Omega;
    L0 = 1e-9,
    KxStar = KxConst,
    f = 4)
    humanIgGC = gen_IgGC(4)
    human = polyfc_ActV(L0, KxStar, f, importRtot(murine=false), humanIgGC,
        importKav(murine=false), humanActI)
    murine = polyfc_ActV(L0, KxStar, f, importRtot(murine=true), Omega * humanIgGC,
        importKav(murine=true), murineActI)
    pcor = [cor(human[i,:], murine[i,:]) for i in 1:size(human,1)]
    pcor[isnan.(pcor)] .= 0
    return sum(pcor)/length(pcor)
end


function brute_force_discrete(divisor=3)
    Omegas = gen_discrete_Omega(4, divisor)
    pcors = [J_pearson(Omegas[:,:,i]) for i in 1:size(Omegas,3)]
    mxval, mxind = findmax(pcors)
    return (mxval, Omegas[:,:,mxind])
end
