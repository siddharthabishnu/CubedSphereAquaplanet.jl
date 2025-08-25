function BickleyJetGrid(bickley_jet_parameters; arch = CPU(), FT::DataType = Oceananigans.defaults.FloatType)
    bickley_jet_grid = ConformalCubedSphereGrid(
        arch, FT;
        panel_size = (bickley_jet_parameters.Nx, bickley_jet_parameters.Ny, bickley_jet_parameters.Nz),
        z = (-bickley_jet_parameters.Lz, 0),
        radius = bickley_jet_parameters.R,
        horizontal_direction_halo = bickley_jet_parameters.H,
        partition = CubedSpherePartition(; R = 1))

    return bickley_jet_grid
end
