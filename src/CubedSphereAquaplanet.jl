module CubedSphereAquaplanet

include("BarotropicAquaplanet/BarotropicAquaplanet.jl")
include("BaroclinicAquaplanet/BaroclinicAquaplanet.jl")

using .BarotropicAquaplanet
using .BaroclinicAquaplanet

end # module CubedSphereAquaplanet
