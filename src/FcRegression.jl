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
include("CD16b.jl")
include("plotMCMC.jl")
include("dataDepletion.jl")
include("effectorBind.jl")
include("regression.jl")

include("figures/figureW.jl")
include("figures/figure1.jl")
include("figures/figure2.jl")
include("figures/figure3.jl")
include("figures/figure4.jl")
include("figures/figure5.jl")
include("figures/figure6.jl")

include("figures/figureS1.jl")
include("figures/figureS2.jl")
include("figures/figureS4.jl")

function figureAll()
    setGadflyTheme()

    figure1()
    figure2()
    figure3()
    figure4()
    figure5()
    figure6()

    figureS1()
    figureS2()
    figureS3()
    figureS4()
end

export importKav, importRtot, importKavDist
end # module
