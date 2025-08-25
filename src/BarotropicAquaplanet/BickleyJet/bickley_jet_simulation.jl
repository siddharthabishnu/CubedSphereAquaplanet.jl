using Oceananigans.Advection: EnstrophyConserving
using Oceananigans.Units: minutes, hours, days
using Oceananigans.Utils: getregion

function BickleyJetSimulation(
    arch;
    parameters = BickleyJetParameters(),
    grid = BickleyJetGrid(parameters; arch),
    momentum_advection = WENOVectorInvariant(vorticity_order = 9),
    tracer_advection = WENO(order=9),
    free_surface = ExplicitFreeSurface(; gravitational_acceleration = parameters.g),
    coriolis = HydrostaticSphericalCoriolis(scheme = EnstrophyConserving(), rotation_rate = parameters.Ω),
    tracers = :c,
    buoyancy = nothing,
    Courant_number = 0.2,
    # Estimate time step from the minimum grid spacing based on the CFL condition
    Δt = Courant_number * min(minimum_xspacing(getregion(grid, 1)), 
                              minimum_yspacing(getregion(grid, 1))) / parameters.c,
    stop_time = deg2rad(360) / abs((parameters.n * (3 + parameters.n) * parameters.ω - 2parameters.Ω) 
                                   / ((1 + parameters.n) * (2 + parameters.n))),
    Ntime = round(Int, stop_time / Δt),
    align_time_step = false,
    progress_message_iteration_interval = 100,
    checkpointer_interval = round(Int, Ntime * Δt / 3),
    output_interval = round(Int, Ntime * Δt / 75))
    
    #####
    ##### Model setup
    #####

    @info "Building model..."
    bickley_jet_model = (
        HydrostaticFreeSurfaceModel(; grid,
                                      momentum_advection,
                                      tracer_advection,
                                      free_surface,
                                      coriolis,
                                      tracers,
                                      buoyancy))

    @info "Initializing model..."
    BickleyJetInitialConditions!(parameters, bickley_jet_model)

    #####
    ##### Simulation setup
    #####
    
    @info "Stop time = $(prettytime(stop_time))"
    @info "Number of time steps = $Ntime"

    bickley_jet_simulation = Simulation(bickley_jet_model; Δt, stop_time, align_time_step)

    #####
    ##### Callbacks
    #####

    start_time = time_ns()
    wall_time = start_time

    # Print a progress message
    progress_format = Printf.Format(
        "Iteration: %04d, time: %s, Δt: %s, max|u|: %.3f, max|v|: %.3f, max|η|: %.3f, max|c|: %.3f, " *
        "wall time for last %03d iterations: %s, wall time: %s\n"
    )

    function progress_message(simulation)
        Printf.format(stdout, progress_format,
                      iteration(simulation),
                      prettytime(simulation),
                      prettytime(simulation.Δt),
                      maximum(abs, simulation.model.velocities.u),
                      maximum(abs, simulation.model.velocities.v),
                      maximum(abs, simulation.model.free_surface.η),
                      maximum(abs, simulation.model.tracers.c),
                      progress_message_iteration_interval,
                      prettytime(1e-9 * (time_ns() - wall_time)),
                      prettytime(1e-9 * (time_ns() - start_time)))

        wall_time = time_ns()
    end

    bickley_jet_simulation.callbacks[:progress] = Callback(
        progress_message, IterationInterval(progress_message_iteration_interval))

    #####
    ##### Build checkpointer and output writer
    #####

    @info "Building checkpointer and output writer..."
    BickleyJetOutputs!(bickley_jet_simulation; checkpointer_interval, output_interval)

    return bickley_jet_simulation
end
