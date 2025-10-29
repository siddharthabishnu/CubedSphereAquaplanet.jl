module BaroclinicInstability

export BaroclinicInstabilityParameters
export BaroclinicInstabilityGrid
export BaroclinicInstabilityBoundaryConditions
export BaroclinicInstabilityClosure
export BaroclinicInstabilityInitialConditions!
export BaroclinicInstabilitySimulation
export BaroclinicInstabilityOutputs!

using Printf
using JLD2

using Oceananigans
using Oceananigans.Grids
using Oceananigans.Operators
using Oceananigans.BoundaryConditions

include("baroclinic_instability_parameters.jl")
include("baroclinic_instability_grid.jl")
include("baroclinic_instability_boundary_conditions.jl")
include("baroclinic_instability_closure.jl")
include("baroclinic_instability_initial_conditions.jl")
include("baroclinic_instability_outputs.jl")
include("baroclinic_instability_simulation.jl")

end # module BaroclinicInstability
