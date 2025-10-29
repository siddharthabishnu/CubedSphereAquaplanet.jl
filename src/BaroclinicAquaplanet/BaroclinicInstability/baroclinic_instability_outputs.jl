function BaroclinicInstabilityOutputs!(baroclinic_instability_simulation;
                                       output_directory = "baroclinic_instability_output",
                                       verbose = false, 
                                       overwrite_existing = true,
                                       checkpointer_interval = 40days, 
                                       output_interval = 2.5days)
    baroclinic_instability_model = baroclinic_instability_simulation.model
    baroclinic_instability_grid = baroclinic_instability_model.grid

    checkpointer_filename = "baroclinic_instability_checkpointer"
    baroclinic_instability_simulation.output_writers[:checkpointer] =
        Checkpointer(baroclinic_instability_model;
                     dir = output_directory,
                     schedule = TimeInterval(checkpointer_interval), 
                     prefix = checkpointer_filename,
                     overwrite_existing)

    outputs = (; u = baroclinic_instability_model.velocities.u,
                 v = baroclinic_instability_model.velocities.v,
                 T = baroclinic_instability_model.tracers.T,
                 S = baroclinic_instability_model.tracers.S)
    output_filename = "baroclinic_instability_surface_prognostic_fields_output"
    baroclinic_instability_simulation.output_writers[:output] =
        JLD2Writer(baroclinic_instability_model, outputs;
                   dir = output_directory,
                   schedule = TimeInterval(output_interval),
                   filename = output_filename,
                   indices = (:, :, baroclinic_instability_grid.Nz),
                   verbose,
                   overwrite_existing)
        
    outputs = (; w = baroclinic_instability_model.velocities.w,
                 η = baroclinic_instability_model.free_surface.η)
    output_filename = "baroclinic_instability_surface_diagnostic_fields_output"
    baroclinic_instability_simulation.output_writers[:output] =
        JLD2Writer(baroclinic_instability_model, outputs;
                   dir = output_directory,
                   schedule = TimeInterval(output_interval),
                   filename = output_filename,
                   indices = (:, :, baroclinic_instability_grid.Nz + 1),
                   verbose,
                   overwrite_existing)
        
end
