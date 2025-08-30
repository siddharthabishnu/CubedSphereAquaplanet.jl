using KernelAbstractions: @kernel, @index
using Oceananigans.Grids: node
using Oceananigans.Utils: Iterate, launch!

function BaroclinicWaveInitialConditions!(baroclinic_wave_parameters, baroclinic_wave_model)
    baroclinic_wave_grid = baroclinic_wave_model.grid
    arch = baroclinic_wave_grid.architecture
    Nx, Ny, Nz = size(baroclinic_wave_grid)
    
    @inline Tᵢ(λ, φ, z) = 30 * (1 - tanh((abs(φ) - 45) / 8)) / 2 + rand()
    @inline Sᵢ(λ, φ, z) = 28 - 5e-3 * z + rand()
    
    #####
    ##### Initial condition
    #####

    @info "Setting initial condition..."
    
    set!(baroclinic_wave_model.tracers.T, Tᵢ)
    set!(baroclinic_wave_model.tracers.S, Sᵢ)

    fill_halo_regions!(baroclinic_wave_model.tracers.T)
    fill_halo_regions!(baroclinic_wave_model.tracers.S)
end
