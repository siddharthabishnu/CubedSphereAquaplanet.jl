using CubedSphereAquaplanet.BarotropicAquaplanet.BickleyJet
using Oceananigans
using Oceananigans.Units: minutes, hours, days
using Oceananigans.Utils: getregion
using Printf

arch = CPU()

unit_sphere = true
bickley_jet_parameters = BickleyJetParameters(; unit_sphere)
bickley_jet_grid = BickleyJetGrid(bickley_jet_parameters; arch)

minimum_grid_spacing = filter(!iszero, getregion(bickley_jet_grid, 1).Δxᶠᶠᵃ) |> minimum
c = bickley_jet_parameters.c

CourantNumber = 0.2
Δt = CourantNumber * minimum_grid_spacing / c # CFL for Bickley jet

stop_time = unit_sphere ? 45 : 45days
Ntime = ceil(Int, stop_time / Δt)
# Redefine the stop time.
stop_time = Ntime * Δt

animation_time = 7.5 # seconds
framerate = 10
Nframes = animation_time * framerate # excluding the initial condition frame
output_interval = stop_time / Nframes # Specify output_interval to be the simulation time per frame.
# Redefine the simulation time per frame.
output_interval = ceil(Int, output_interval / Δt) * Δt
# Redefine the stop time.
stop_time = ceil(Int, stop_time / output_interval) * output_interval
# Redefine the number of time steps.
Ntime = round(Int, stop_time / Δt)
# Redefine the number of frames.
Nframes = round(Int, stop_time / output_interval) # excluding the initial condition frame
# Redefine the animation time.
animation_time = Nframes / framerate
checkpointer_interval_by_output_interval = 15
checkpointer_interval = checkpointer_interval_by_output_interval * output_interval

bickley_jet_simulation = BickleyJetSimulation(arch;
                                              parameters = bickley_jet_parameters,
                                              grid = bickley_jet_grid,
                                              Δt,
                                              stop_time,
                                              Ntime,
                                              checkpointer_interval,
                                              output_interval)
run!(bickley_jet_simulation)

include("bickley_jet_visualization.jl")
Nplots = Nframes
plot_iteration_interval = floor(Int, Nframes / Nplots)
prettytimes = [prettytime((i - 1) * output_interval) for i in 1:Nframes+1]
iPlot_Start = Nframes - 18
iPlot_Δ = 3
make_panelwise_visualization_plots_with_halos = false
make_panelwise_visualization_plots = true
geo_heatmap_type = "heatsphere"
make_geo_heatmap_visualization_plots = true
plot_frames = false
make_panelwise_visualization_animation_with_halos = false
make_panelwise_visualization_animation = true
make_geo_heatmap_visualization_animation = true
BickleyJetVisualization!(bickley_jet_grid, Nplots, Δt, plot_iteration_interval, prettytimes, framerate;
                         iPlot_Start,
                         iPlot_Δ,
                         make_panelwise_visualization_plots_with_halos,
                         make_panelwise_visualization_plots,
                         geo_heatmap_type,
                         make_geo_heatmap_visualization_plots,
                         plot_frames,
                         make_panelwise_visualization_animation_with_halos,
                         make_panelwise_visualization_animation,
                         make_geo_heatmap_visualization_animation)
