module FcRegression
using LinearAlgebra
import Distances
using polyBindingModel
using Optim
using StatsBase
import StatsBase: geomean, std
using DataFrames

include("figures/figureCommon.jl")

include("fcBindingModel.jl")
include("dataHelpers.jl")
include("mixture.jl")
include("mixtureFit.jl")
include("regression.jl")
include("synergy.jl")

include("figures/figureW.jl")

include("figures/figure1.jl")
include("figures/figure2.jl")
include("robinett.jl")
include("figures/figure3.jl")
include("figures/figure4.jl")
include("figures/figureS2.jl")

include("figures/extra/figureS.jl")
include("figures/extra/figureCelltypeSynergy.jl")
include("figures/extra/figureDiseaseSpecific.jl")
include("figures/extra/figureDiseaseSynergy.jl")


function figureAll()
    setGadflyTheme()

    figure1()
    figure2()
    figure3()
    figure4()
    figure5()
    
    figureS2()
end

export polyfc, polyfc_ActV

end # module
