function TracerAdvectionOutputs!(tracer_advection_simulation;
                                 output_directory = "tracer_advection_output",
                                 verbose = false, 
                                 overwrite_existing = true,
                                 checkpointer_interval = 86400/5,
                                 tracer_interval = 86400/75)
    tracer_advection_model = tracer_advection_simulation.model
    tracer_advection_grid = tracer_advection_model.grid

    checkpointer_filename = "tracer_advection_checkpointer"
    tracer_advection_simulation.output_writers[:checkpointer] =
        Checkpointer(tracer_advection_model;
                     dir = output_directory,
                     schedule = TimeInterval(checkpointer_interval),
                     prefix = checkpointer_filename,
                     overwrite_existing)

    outputs = (; θ = tracer_advection_model.tracers.θ)
    output_filename = "tracer_advection_output"
    tracer_advection_simulation.output_writers[:output] =
        JLD2Writer(tracer_advection_model, outputs;
                   dir = output_directory,
                   schedule = TimeInterval(tracer_interval),
                   filename = output_filename,
                   verbose,
                   overwrite_existing)
end
