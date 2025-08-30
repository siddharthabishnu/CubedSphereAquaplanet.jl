function RossbyHaurwitzWaveGrid(rossby_haurwitz_wave_parameters;
                                arch = CPU(), FT::DataType = Oceananigans.defaults.FloatType)
    #####
    ##### Grid generation
    #####

    @info "Generating grid..."

    rossby_haurwitz_wave_grid = ConformalCubedSphereGrid(
        arch, FT;
        panel_size = (rossby_haurwitz_wave_parameters.Nx, rossby_haurwitz_wave_parameters.Ny,
                      rossby_haurwitz_wave_parameters.Nz),
        z = (-rossby_haurwitz_wave_parameters.Lz, 0),
        radius = rossby_haurwitz_wave_parameters.R,
        horizontal_direction_halo = rossby_haurwitz_wave_parameters.H,
        non_uniform_conformal_mapping = rossby_haurwitz_wave_parameters.non_uniform_conformal_mapping,
        partition = CubedSpherePartition(; R = 1))

    return rossby_haurwitz_wave_grid
end
