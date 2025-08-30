using Oceananigans.Units: minutes, hours, days
using Oceananigans.Utils: getregion

function TracerAdvectionSimulation(
    arch; 
    parameters = TracerAdvectionParameters(), 
    grid = TracerAdvectionGrid(parameters; arch),
    momentum_advection = nothing,
    tracer_advection = WENO(order=9),
    free_surface = ExplicitFreeSurface(; gravitational_acceleration = parameters.g),
    tracers = :θ,
    buoyancy = nothing,
    Courant_number = 0.2,
    # Estimate time step from the minimum grid spacing based on the CFL condition
    Δt = Courant_number * min(minimum_xspacing(getregion(grid, 1)), 
                              minimum_yspacing(getregion(grid, 1))) / parameters.U, # CFL for tracer advection
    stop_time = 2π * parameters.R / parameters.U,
    Ntime = round(Int, stop_time / Δt),
    align_time_step = false,
    progress_message_iteration_interval = 100,
    checkpointer_interval = round(Int, Ntime * Δt / 5),
    tracer_interval = round(Int, Ntime * Δt / 75))
    
    #####
    ##### Model setup
    #####

    @info "Building model..."    
    u, v = TracerAdvectionPrescribedVelocities!(parameters, grid)
    
    tracer_advection_model = (
        HydrostaticFreeSurfaceModel(; grid,
                                      velocities = PrescribedVelocityFields(; u, v),
                                      momentum_advection,
                                      tracer_advection,
                                      free_surface,
                                      tracers,
                                      buoyancy))

    @info "Initializing model..."
    TracerAdvectionInitialConditions!(parameters, tracer_advection_model)

    #####
    ##### Simulation setup
    #####
    
    @info "Stop time = $(prettytime(stop_time))"
    @info "Number of time steps = $Ntime"

    tracer_advection_simulation = Simulation(tracer_advection_model; Δt, stop_time, align_time_step)

    #####
    ##### Callbacks
    #####

    start_time = time_ns()
    wall_time = start_time

    # Print a progress message
    progress_format = Printf.Format(
        "Iteration: %04d, time: %s, Δt: %s, max|θ|: %.3f, wall time for last %03d iterations: %s, wall time: %s\n")

    function progress_message(simulation)
        Printf.format(stdout, progress_format,
                      iteration(simulation),
                      prettytime(simulation),
                      prettytime(simulation.Δt),
                      maximum(abs, simulation.model.tracers.θ),
                      progress_message_iteration_interval,
                      prettytime(1e-9 * (time_ns() - wall_time)),
                      prettytime(1e-9 * (time_ns() - start_time)))

        wall_time = time_ns()
    end

    tracer_advection_simulation.callbacks[:progress] = Callback(progress_message,
                                                                IterationInterval(progress_message_iteration_interval))

    #####
    ##### Build checkpointer and output writer
    #####

    @info "Building checkpointer and output writer..."
    TracerAdvectionOutputs!(tracer_advection_simulation; checkpointer_interval, tracer_interval)

    return tracer_advection_simulation
end
