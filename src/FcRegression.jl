module FcRegression
using LinearAlgebra
import Distances
using polyBindingModel
using Optim

include("fcBindingModel.jl")
include("dataHelpers.jl")
include("mixture.jl")

include("figures/figureCommon.jl")
include("figures/figure1.jl")
include("figures/figure2.jl")


function figureAll()
    setGadflyTheme()
    
    figure1()
    figure2()
end

export polyfc, polyfc_ActV

end # module
