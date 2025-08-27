using CubedSphereAquaplanet.BarotropicAquaplanet.RossbyHaurwitzWave
using Oceananigans
using Oceananigans.Units: minutes, hours, days
using Oceananigans.Utils: getregion
using Printf

arch = CPU()

rossby_haurwitz_wave_parameters = RossbyHaurwitzWaveParameters()
rossby_haurwitz_wave_grid = RossbyHaurwitzWaveGrid(rossby_haurwitz_wave_parameters; arch)

# Compute the time required for a 45° rotation.
n = rossby_haurwitz_wave_parameters.n
ω = rossby_haurwitz_wave_parameters.ω
Ω = rossby_haurwitz_wave_parameters.Ω
angular_velocity = (n * (3 + n) * ω - 2Ω) / ((1 + n) * (2 + n))
stop_time = deg2rad(360) / abs(angular_velocity)

minimum_grid_spacing = filter(!iszero, getregion(rossby_haurwitz_wave_grid, 1).Δxᶠᶠᵃ) |> minimum
c = sqrt(rossby_haurwitz_wave_parameters.g * rossby_haurwitz_wave_parameters.Lz)

CourantNumber = 0.2
Δt = CourantNumber * minimum_grid_spacing / c # CFL for Rossby-Haurwitz wave

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
checkpointer_interval_by_output_interval = 25
checkpointer_interval = checkpointer_interval_by_output_interval * output_interval

rossby_haurwitz_wave_simulation = RossbyHaurwitzWaveSimulation(arch;
                                                               parameters = rossby_haurwitz_wave_parameters,
                                                               grid = rossby_haurwitz_wave_grid,
                                                               Δt, stop_time, Ntime, checkpointer_interval,
                                                               output_interval)
run!(rossby_haurwitz_wave_simulation)

include("rossby_haurwitz_wave_visualization.jl")
Nplots = 5
plot_iteration_interval = floor(Int, Nframes / Nplots)
prettytimes = [prettytime((i - 1) * output_interval) for i in 1:Nframes+1]
make_panelwise_visualization_plots_with_halos = false
make_panelwise_visualization_plots = true
geo_heatmap_type = "heatsphere"
make_geo_heatmap_visualization_plots = true
plot_frames = false
make_panelwise_visualization_animation_with_halos = false
make_panelwise_visualization_animation = true
make_geo_heatmap_visualization_animation = true
RossbyHaurwitzWaveVisualization!(rossby_haurwitz_wave_grid, Nplots, Δt,
                                 plot_iteration_interval, prettytimes, framerate;
                                 make_panelwise_visualization_plots_with_halos,
                                 make_panelwise_visualization_plots,
                                 geo_heatmap_type,
                                 make_geo_heatmap_visualization_plots,
                                 plot_frames,
                                 make_panelwise_visualization_animation_with_halos,
                                 make_panelwise_visualization_animation,
                                 make_geo_heatmap_visualization_animation)
