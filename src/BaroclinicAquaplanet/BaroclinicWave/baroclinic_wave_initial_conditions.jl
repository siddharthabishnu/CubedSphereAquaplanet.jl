using KernelAbstractions: @kernel, @index
using Oceananigans.Grids: node
using Oceananigans.Utils: Iterate, launch!

function BaroclinicWaveInitialConditions!(baroclinic_wave_parameters, baroclinic_wave_model)
    baroclinic_wave_grid = baroclinic_wave_model.grid
    arch = baroclinic_wave_grid.architecture
    Nx, Ny, Nz = size(baroclinic_wave_grid)
    
    @inline Tᵢ(λ, φ, z) =
        baroclinic_wave_parameters.T0 * (1 - tanh((abs(φ) - baroclinic_wave_parameters.φ0)
                                                  / baroclinic_wave_parameters.Δφ)) / 2
        + baroclinic_wave_parameters.ϵT * randn()
    @inline Sᵢ(λ, φ, z) =
        baroclinic_wave_parameters.S0 - baroclinic_wave_parameters.γS * z + baroclinic_wave_parameters.ϵS * randn()

    #####
    ##### Initial condition
    #####

    @info "Setting initial condition..."
    
    set!(baroclinic_wave_model.tracers.T, Tᵢ)
    set!(baroclinic_wave_model.tracers.S, Sᵢ)

    fill_halo_regions!(baroclinic_wave_model.tracers.T)
    fill_halo_regions!(baroclinic_wave_model.tracers.S)
end
