module TracerAdvection

export TracerAdvectionParameters
export TracerAdvectionGrid
export TracerAdvectionBoundaryConditions
export TracerAdvectionClosure
export TracerAdvectionInitialConditions!
export TracerAdvectionSimulation
export TracerAdvectionOutputs!

using Printf
using JLD2

using Oceananigans
using Oceananigans.Grids
using Oceananigans.Operators
using Oceananigans.BoundaryConditions

include("tracer_advection_parameters.jl")
include("tracer_advection_grid.jl")
include("tracer_advection_prescribed_velocities.jl")
include("tracer_advection_initial_conditions.jl")
include("tracer_advection_outputs.jl")
include("tracer_advection_simulation.jl")

end # module TracerAdvection
