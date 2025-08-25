module BarotropicAquaplanet

include("TracerAdvection/TracerAdvection.jl")
include("RossbyHaurwitzWave/RossbyHaurwitzWave.jl")
include("BickleyJet/BickleyJet.jl"  )

using .TracerAdvection
using .RossbyHaurwitzWave
using .BickleyJet

end # module BarotropicAquaplanet
