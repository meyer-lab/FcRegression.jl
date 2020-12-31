module FcRegression
using LinearAlgebra
import Distances
using polyBindingModel
using Optim

include("figures/figureCommon.jl")

include("fcBindingModel.jl")
include("dataHelpers.jl")
include("mixture.jl")
include("regression.jl")
include("synergy.jl")

include("figures/figureW.jl")

include("figures/figure1.jl")
include("figures/figure2.jl")
include("figures/figure3.jl")


function figureAll()
    setGadflyTheme()
    
    figure1()
    figure2()
    figure3()
end

export polyfc, polyfc_ActV

end # module
