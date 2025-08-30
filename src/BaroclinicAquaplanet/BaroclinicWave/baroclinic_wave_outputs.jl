function BaroclinicWaveOutputs!(baroclinic_wave_simulation;
                                output_directory = "baroclinic_wave_output",
                                verbose = false, 
                                overwrite_existing = true,
                                checkpointer_interval = 40days, 
                                output_interval = 2.5days)
    baroclinic_wave_model = baroclinic_wave_simulation.model
    baroclinic_wave_grid = baroclinic_wave_model.grid

    checkpointer_filename = "baroclinic_wave_checkpointer"
    baroclinic_wave_simulation.output_writers[:checkpointer] =
        Checkpointer(baroclinic_wave_model;
                     dir = output_directory,
                     schedule = TimeInterval(checkpointer_interval), 
                     prefix = checkpointer_filename,
                     overwrite_existing)

    outputs = (; u = baroclinic_wave_model.velocities.u,
                 v = baroclinic_wave_model.velocities.v,
                 T = baroclinic_wave_model.tracers.T,
                 S = baroclinic_wave_model.tracers.S)
    output_filename = "baroclinic_wave_surface_prognostic_fields_output"
    baroclinic_wave_simulation.output_writers[:output] =
        JLD2Writer(baroclinic_wave_model, outputs;
                   dir = output_directory,
                   schedule = TimeInterval(output_interval),
                   filename = output_filename,
                   indices = (:, :, baroclinic_wave_grid.Nz),
                   verbose,
                   overwrite_existing)
        
    outputs = (; w = baroclinic_wave_model.velocities.w,
                 η = baroclinic_wave_model.free_surface.η)
    output_filename = "baroclinic_wave_surface_diagnostic_fields_output"
    baroclinic_wave_simulation.output_writers[:output] =
        JLD2Writer(baroclinic_wave_model, outputs;
                   dir = output_directory,
                   schedule = TimeInterval(output_interval),
                   filename = output_filename,
                   indices = (:, :, baroclinic_wave_grid.Nz + 1),
                   verbose,
                   overwrite_existing)
        
end
