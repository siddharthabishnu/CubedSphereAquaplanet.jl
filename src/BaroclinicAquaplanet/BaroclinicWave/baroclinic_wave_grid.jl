function BaroclinicWaveGrid(baroclinic_wave_parameters; arch = CPU(), FT::DataType = Oceananigans.defaults.FloatType)
    #####
    ##### Grid generation
    #####

    @info "Generating grid..."

    baroclinic_wave_grid = ConformalCubedSphereGrid(
        arch, FT;
        panel_size = (baroclinic_wave_parameters.Nx, baroclinic_wave_parameters.Ny, baroclinic_wave_parameters.Nz),
        z = (-baroclinic_wave_parameters.Lz, 0),
        radius = baroclinic_wave_parameters.R,
        horizontal_direction_halo = baroclinic_wave_parameters.H,
        non_uniform_conformal_mapping = baroclinic_wave_parameters.non_uniform_conformal_mapping,
        partition = CubedSpherePartition(; R = 1))

    return baroclinic_wave_grid
end
