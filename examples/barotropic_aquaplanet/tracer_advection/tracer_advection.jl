using CubedSphereAquaplanet.BarotropicAquaplanet.TracerAdvection
using Oceananigans
using Oceananigans.Units: minutes, hours, days
using Oceananigans.Utils: getregion
using Printf

arch = CPU()

tracer_advection_parameters = TracerAdvectionParameters()
tracer_advection_grid = TracerAdvectionGrid(tracer_advection_parameters; arch)

# Estimate time step from the minimum grid spacing based on the CFL condition.
Δx = minimum_xspacing(getregion(tracer_advection_grid, 1))
Δy = minimum_yspacing(getregion(tracer_advection_grid, 1))
CourantNumber = 0.2
Δt = CourantNumber * min(Δx, Δy) / tracer_advection_parameters.U # CFL for tracer advection

stop_time = 2π * tracer_advection_parameters.R / tracer_advection_parameters.U
Ntime = ceil(Int, stop_time / Δt)
# Redefine the stop time.
stop_time = Ntime * Δt

animation_time = 7.5 # seconds
framerate = 10
Nframes = animation_time * framerate # excluding the initial condition frame
tracer_interval = stop_time / Nframes # Specify tracer_interval to be the simulation time per frame.
# Redefine the simulation time per frame.
tracer_interval = ceil(Int, tracer_interval / Δt) * Δt
# Redefine the stop time.
stop_time = ceil(Int, stop_time / tracer_interval) * tracer_interval
# Redefine the number of time steps.
Ntime = round(Int, stop_time / Δt)
# Redefine the number of frames.
Nframes = round(Int, stop_time / tracer_interval) # excluding the initial condition frame
# Redefine the animation time.
animation_time = Nframes / framerate
checkpointer_interval_by_tracer_interval = 100
checkpointer_interval = checkpointer_interval_by_tracer_interval * tracer_interval

tracer_advection_simulation = TracerAdvectionSimulation(arch;
                                                        Δt, stop_time, Ntime, checkpointer_interval,
                                                        tracer_interval)
run!(tracer_advection_simulation)

include("tracer_advection_visualization.jl")
Nplots = 5
plot_iteration_interval = floor(Int, Nframes / Nplots)
prettytimes = [prettytime((i - 1) * tracer_interval) for i in 1:Nframes+1]
make_panelwise_visualization_plots_with_halos = false
make_panelwise_visualization_plots = true
make_geo_heatlatlon_visualization_plots = true
make_panelwise_visualization_animation_with_halos = false
make_panelwise_visualization_animation = true
make_geo_heatlatlon_visualization_animation = true
TracerAdvectionVisualization!(tracer_advection_grid, tracer_advection_parameters.θ₀, Nplots, Δt,
                              plot_iteration_interval, prettytimes, framerate;
                              make_panelwise_visualization_plots_with_halos,
                              make_panelwise_visualization_plots,
                              make_geo_heatlatlon_visualization_plots,
                              make_panelwise_visualization_animation_with_halos,
                              make_panelwise_visualization_animation,
                              make_geo_heatlatlon_visualization_animation)
