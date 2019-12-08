module FcgR

include("fcBindingModel.jl")
include("dataHelpers.jl")
include("regression.jl")
include("synergy.jl")

export calculateIsobologram, polyfc, polyfc_ActV

end # module
