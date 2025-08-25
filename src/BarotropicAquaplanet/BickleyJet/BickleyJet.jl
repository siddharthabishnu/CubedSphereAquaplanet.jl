module BickleyJet

export BickleyJetParameters
export BickleyJetGrid
export BickleyJetBoundaryConditions
export BickleyJetClosure
export BickleyJetInitialConditions!
export BickleyJetSimulation
export BickleyJetOutputs!

using Printf
using JLD2

using Oceananigans
using Oceananigans.Grids
using Oceananigans.Operators
using Oceananigans.BoundaryConditions

include("bickley_jet_parameters.jl")
include("bickley_jet_grid.jl")
include("bickley_jet_initial_conditions.jl")
include("bickley_jet_outputs.jl")
include("bickley_jet_simulation.jl")

end # module BickleyJet
