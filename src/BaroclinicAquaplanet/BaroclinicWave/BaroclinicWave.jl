module BaroclinicWave

export BaroclinicWaveParameters
export BaroclinicWaveGrid
export BaroclinicWaveBoundaryConditions
export BaroclinicWaveClosure
export BaroclinicWaveInitialConditions!
export BaroclinicWaveSimulation
export BaroclinicWaveOutputs!

using Printf
using JLD2

using Oceananigans
using Oceananigans.Grids
using Oceananigans.Operators
using Oceananigans.BoundaryConditions

include("baroclinic_wave_parameters.jl")
include("baroclinic_wave_grid.jl")
include("baroclinic_wave_boundary_conditions.jl")
include("baroclinic_wave_closure.jl")
include("baroclinic_wave_initial_conditions.jl")
include("baroclinic_wave_outputs.jl")
include("baroclinic_wave_simulation.jl")

end # module BaroclinicWave
