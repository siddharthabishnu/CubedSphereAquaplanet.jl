using SeawaterPolynomials.TEOS10: TEOS10EquationOfState

using Oceananigans.Units: minutes, hours, days
using Oceananigans.Utils: getregion

function BaroclinicWaveSimulation(arch;
                                  parameters = BaroclinicWaveParameters(),
                                  grid = BaroclinicWaveGrid(parameters; arch),
                                  momentum_advection = WENOVectorInvariant(vorticity_order = 9),
                                  tracer_advection = WENO(order = 9),
                                  substeps = 50,
                                  free_surface = SplitExplicitFreeSurface(grid;
                                                                          substeps,
                                                                          gravitational_acceleration = parameters.g),
                                  coriolis = HydrostaticSphericalCoriolis(rotation_rate = parameters.Ω),
                                  boundary_conditions = BaroclinicWaveBoundaryConditions(parameters),
                                  closure = BaroclinicWaveClosure(),
                                  tracers = (:T, :S),
                                  buoyancy = SeawaterBuoyancy(equation_of_state = TEOS10EquationOfState()),
                                  Δt = 5minutes,
                                  stop_time = 200days,
                                  Ntime = round(Int, stop_time / Δt),
                                  align_time_step = false,
                                  progress_message_iteration_interval = 100,
                                  checkpointer_interval = 40days,
                                  output_interval = 2.5days)
    #####
    ##### Model setup
    #####

    @info "Building model..."
    baroclinic_wave_model = (
        HydrostaticFreeSurfaceModel(; grid,
                                      momentum_advection,
                                      tracer_advection,
                                      free_surface,
                                      coriolis,
                                      boundary_conditions,
                                      closure,
                                      tracers,
                                      buoyancy))

    @info "Initializing model..."
    BaroclinicWaveInitialConditions!(parameters, baroclinic_wave_model)

    #####
    ##### Simulation setup
    #####
    
    @info "Stop time = $(prettytime(stop_time))"
    @info "Number of time steps = $Ntime"

    baroclinic_wave_simulation = Simulation(baroclinic_wave_model; Δt, stop_time, align_time_step)

    #####
    ##### Callbacks
    #####

    start_time = time_ns()
    wall_time = start_time

    # Print a progress message
    progress_format = Printf.Format(
        "Iteration: %04d, time: %s, Δt: %s, max|u|: %.3f, max|v|: %.3f, max|w|: %.3f, max|η|: %.3f, max|T|: %.3f, " *
        "max|S|: %.3f, wall time for last %03d iterations: %s, wall time: %s\n"
    )

    function progress_message(simulation)
        Printf.format(stdout, progress_format,
                      iteration(simulation),
                      prettytime(simulation),
                      prettytime(simulation.Δt),
                      maximum(abs, simulation.model.velocities.u),
                      maximum(abs, simulation.model.velocities.v),
                      maximum(abs, simulation.model.velocities.w),
                      maximum(abs, simulation.model.free_surface.η),
                      maximum(abs, simulation.model.tracers.T),
                      maximum(abs, simulation.model.tracers.S),
                      progress_message_iteration_interval,
                      prettytime(1e-9 * (time_ns() - wall_time)),
                      prettytime(1e-9 * (time_ns() - start_time)))

        wall_time = time_ns()
    end

    baroclinic_wave_simulation.callbacks[:progress] = Callback(
        progress_message, IterationInterval(progress_message_iteration_interval))

    #####
    ##### Build checkpointer and output writer
    #####

    @info "Building checkpointer and output writer..."
    BaroclinicWaveOutputs!(baroclinic_wave_simulation; checkpointer_interval, output_interval)

    return baroclinic_wave_simulation
end
