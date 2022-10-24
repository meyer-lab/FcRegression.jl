module FcRegression
using LinearAlgebra
import Distances
using polyBindingModel
using Optim
using StatsBase
using DataFrames
import Printf: @sprintf
import Cairo, Fontconfig
using NamedArrays


include("figures/figureCommon.jl")

include("fcBindingModel.jl")
include("dataHelpers.jl")
include("mixture.jl")
include("runMCMC.jl")
include("plotMCMC.jl")
include("dataDepletion.jl")
include("regression.jl")

include("figures/figureW.jl")

include("figures/figureA2a.jl")
include("figures/figureA3.jl")
include("figures/figureA4.jl")
include("figures/figureA5.jl")
include("figures/figureA6.jl")
include("figures/figureS1.jl")
include("figures/figureS2.jl")

include("figures/figure1.jl")
include("figures/figure2.jl")
include("figures/figure3.jl")
include("figures/figure4.jl")
include("figures/figure5.jl")


function figureAll()
    setGadflyTheme()

    figure1()
    figure2()
    figure3()
    figure4()

    figureS1()
    figureS2()
end

export polyfc, polyfc_ActV

end # module
