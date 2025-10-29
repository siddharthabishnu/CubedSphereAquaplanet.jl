function BaroclinicInstabilityInitialConditions!(baroclinic_instability_parameters, baroclinic_instability_model)
    baroclinic_instability_grid = baroclinic_instability_model.grid
    arch = baroclinic_instability_grid.architecture
    Nx, Ny, Nz = size(baroclinic_instability_grid)
    
    @inline Tᵢ(λ, φ, z) =
        baroclinic_instability_parameters.T0 * (1 - tanh((abs(φ) - baroclinic_instability_parameters.φ0)
                                                         / baroclinic_instability_parameters.Δφ)) / 2
        + baroclinic_instability_parameters.ϵT * randn()
    @inline Sᵢ(λ, φ, z) =
        (baroclinic_instability_parameters.S0 - baroclinic_instability_parameters.γS * z
         + baroclinic_instability_parameters.ϵS * randn())

    #####
    ##### Initial condition
    #####

    @info "Setting initial condition..."
    
    set!(baroclinic_instability_model.tracers.T, Tᵢ)
    set!(baroclinic_instability_model.tracers.S, Sᵢ)

    fill_halo_regions!(baroclinic_instability_model.tracers.T)
    fill_halo_regions!(baroclinic_instability_model.tracers.S)
end
