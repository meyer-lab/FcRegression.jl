module FcgR

include("fcBindingModel.jl")
include("dataHelpers.jl")
include("regression.jl")
include("synergy.jl")
include("translation.jl")
include("figures/figureB1.jl")

export calculateIsobologram, polyfc, polyfc_ActV

end # module
