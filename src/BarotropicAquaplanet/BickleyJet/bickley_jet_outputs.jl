function BickleyJetOutputs!(bickley_jet_simulation;
                            output_directory = "bickley_jet_output",
                            verbose = false, 
                            overwrite_existing = true,
                            checkpointer_interval = 5days, 
                            output_interval = 12hours)
    bickley_jet_model = bickley_jet_simulation.model
    bickley_jet_grid = bickley_jet_model.grid

    checkpointer_filename = "bickley_jet_checkpointer"
    bickley_jet_simulation.output_writers[:checkpointer] = Checkpointer(
        bickley_jet_model;
        dir = output_directory,
        schedule = TimeInterval(checkpointer_interval), 
        prefix = checkpointer_filename,
        overwrite_existing)

    outputs = (; u = bickley_jet_model.velocities.u, v = bickley_jet_model.velocities.v,
                 η = bickley_jet_model.free_surface.η, c = bickley_jet_model.tracers.c)
    output_filename = "bickley_jet_output"
    bickley_jet_simulation.output_writers[:output] = JLD2Writer(
        bickley_jet_model, outputs;
        dir = output_directory,
        schedule = TimeInterval(output_interval),
        filename = output_filename,
        verbose,
        overwrite_existing)
end
