module BarotropicAquaplanet

include("TracerAdvection/TracerAdvection.jl")
include("RossbyHaurwitzWave/RossbyHaurwitzWave.jl")

using .TracerAdvection
using .RossbyHaurwitzWave

end # module BarotropicAquaplanet
