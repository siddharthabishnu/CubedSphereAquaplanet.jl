function BaroclinicInstabilityGrid(baroclinic_instability_parameters;
                                   arch = CPU(),
                                   FT::DataType = Oceananigans.defaults.FloatType)
    #####
    ##### Grid generation
    #####

    @info "Generating grid..."

    baroclinic_instability_grid = ConformalCubedSphereGrid(
        arch, FT;
        panel_size = (baroclinic_instability_parameters.Nx, baroclinic_instability_parameters.Ny,
                      baroclinic_instability_parameters.Nz),
        z = (-baroclinic_instability_parameters.Lz, 0),
        radius = baroclinic_instability_parameters.R,
        horizontal_direction_halo = baroclinic_instability_parameters.H,
        non_uniform_conformal_mapping = baroclinic_instability_parameters.non_uniform_conformal_mapping,
        partition = CubedSpherePartition(; R = 1))

    return baroclinic_instability_grid
end
