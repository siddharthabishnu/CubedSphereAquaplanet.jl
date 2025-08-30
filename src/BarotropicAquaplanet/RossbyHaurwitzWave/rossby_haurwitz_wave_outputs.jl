function RossbyHaurwitzWaveOutputs!(rossby_haurwitz_wave_simulation;
                                    output_directory = "rossby_haurwitz_wave_output",
                                    verbose = false, 
                                    overwrite_existing = true,
                                    checkpointer_interval = 1.296e6/5,
                                    output_interval = 1.296e6/75)
    rossby_haurwitz_wave_model = rossby_haurwitz_wave_simulation.model
    rossby_haurwitz_wave_grid = rossby_haurwitz_wave_model.grid

    checkpointer_filename = "rossby_haurwitz_wave_checkpointer"
    rossby_haurwitz_wave_simulation.output_writers[:checkpointer] =
        Checkpointer(rossby_haurwitz_wave_model;
                     dir = output_directory,
                     schedule = TimeInterval(checkpointer_interval),
                     prefix = checkpointer_filename,
                     overwrite_existing)

    outputs = (; u = rossby_haurwitz_wave_model.velocities.u,
                 v = rossby_haurwitz_wave_model.velocities.v,
                 η = rossby_haurwitz_wave_model.free_surface.η)
    output_filename = "rossby_haurwitz_wave_output"
    rossby_haurwitz_wave_simulation.output_writers[:output] =
        JLD2Writer(rossby_haurwitz_wave_model, outputs;
                   dir = output_directory,
                   schedule = TimeInterval(output_interval),
                   filename = output_filename,
                   verbose,
                   overwrite_existing)
end
