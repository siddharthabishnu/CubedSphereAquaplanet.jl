function TracerAdvectionGrid(tracer_advection_parameters; arch = CPU(), FT::DataType = Oceananigans.defaults.FloatType)
    #####
    ##### Grid generation
    #####

    @info "Generating grid..."

    tracer_advection_grid = ConformalCubedSphereGrid(
        arch, FT;
        panel_size = (tracer_advection_parameters.Nx, tracer_advection_parameters.Ny, tracer_advection_parameters.Nz),
        z = (-tracer_advection_parameters.Lz, 0),
        radius = tracer_advection_parameters.R,
        horizontal_direction_halo = tracer_advection_parameters.H,
        non_uniform_conformal_mapping = tracer_advection_parameters.non_uniform_conformal_mapping,
        partition = CubedSpherePartition(; R = 1))

    return tracer_advection_grid
end
