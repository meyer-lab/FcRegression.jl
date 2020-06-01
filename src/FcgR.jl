module FcgR
using LinearAlgebra
using ForwardDiff
import Distances
using Gadfly
using Compose
using polyBindingModel

include("fcBindingModel.jl")
include("dataHelpers.jl")
include("regression.jl")
include("synergy.jl")
include("translation.jl")
include("systemsSerology.jl")

include("figures/figureCommon.jl")
include("figures/figureW.jl")
include("figures/figure1.jl")
include("figures/figure2.jl")
include("figures/figureB1.jl")
include("figures/figureB2.jl")

function figureAll()
    setGadflyTheme()
    figureB1()
    figureB2()
    figureB3()
    figureB4()
    figureB5()
    figureB6()
    figureB7()
    figureB11()
    figureB12()
    figureB13()
    figureB14()
    figureB15()
    figureB16()
    figureB17()
    figureS1()
    figureS2()
    figureS3()
    figureS4()
    figureS5()
end

export polyfc, polyfc_ActV

end # module
