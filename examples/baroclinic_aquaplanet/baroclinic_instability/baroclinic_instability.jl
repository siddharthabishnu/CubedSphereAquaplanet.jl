using CubedSphereAquaplanet.BaroclinicAquaplanet.BaroclinicInstability
using Oceananigans
using Oceananigans.Units: minutes, hours, days
using Oceananigans.Utils: getregion
using Printf

arch = CPU()

baroclinic_instability_parameters = BaroclinicInstabilityParameters()
baroclinic_instability_grid = BaroclinicInstabilityGrid(baroclinic_instability_parameters; arch)

Δt = 5minutes
# Estimate barotropic time step via the CFL condition using the Courant number, minimum grid spacing, and gravity wave 
# speed.
CourantNumber = 0.2
c = sqrt(baroclinic_instability_parameters.g * baroclinic_instability_parameters.Lz)
Δx = minimum_xspacing(baroclinic_instability_grid, Face(), Face(), Face())
Δy = minimum_yspacing(baroclinic_instability_grid, Face(), Face(), Face())
minimum_substeps = ceil(Int, 2c * Δt / (CourantNumber * min(Δx, Δy)))
print("The minimum number of substeps required to satisfy the CFL condition for Courant number $CourantNumber is " 
      * "$minimum_substeps.\n")
substeps = max(minimum_substeps, 50)

stop_time = 200days
Ntime = ceil(Int, stop_time / Δt)
# Redefine the stop time.
stop_time = Ntime * Δt
animation_time = 8 # seconds
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
checkpointer_interval_by_output_interval = 16
checkpointer_interval = checkpointer_interval_by_output_interval * output_interval

baroclinic_instability_simulation = BaroclinicInstabilitySimulation(arch;
                                                                    parameters = baroclinic_instability_parameters,
                                                                    grid = baroclinic_instability_grid,
                                                                    substeps,
                                                                    Δt,
                                                                    stop_time,
                                                                    Ntime,
                                                                    checkpointer_interval,
                                                                    output_interval)
run!(baroclinic_instability_simulation)

include("baroclinic_instability_visualization.jl")
Nplots = Nframes
plot_iteration_interval = floor(Int, Nframes / Nplots)
prettytimes = [prettytime((i - 1) * output_interval) for i in 1:Nframes+1]
iPlot_Start = Nframes - 6
iPlot_Δ = 1
plot_states = Dict(:u => false, :v => false, :w => true, :η => true, :T => true, :S => true, :ζ => true)
make_panelwise_visualization_plots_with_halos = false
make_panelwise_visualization_plots = true
geo_heatmap_type = "heatsphere"
make_geo_heatmap_visualization_plots = true
plot_frames = false
make_panelwise_visualization_animation_with_halos = false
make_panelwise_visualization_animation = true
make_geo_heatmap_visualization_animation = true
BaroclinicInstabilityVisualization!(baroclinic_instability_grid, Nplots, Δt, plot_iteration_interval, prettytimes,
                                    framerate;
                                    iPlot_Start,
                                    iPlot_Δ,
                                    plot_states,
                                    make_panelwise_visualization_plots_with_halos,
                                    make_panelwise_visualization_plots,
                                    geo_heatmap_type,
                                    make_geo_heatmap_visualization_plots,
                                    plot_frames,
                                    make_panelwise_visualization_animation_with_halos,
                                    make_panelwise_visualization_animation,
                                    make_geo_heatmap_visualization_animation)
