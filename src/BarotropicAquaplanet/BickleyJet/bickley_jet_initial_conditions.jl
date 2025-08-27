using KernelAbstractions: @kernel, @index
using Oceananigans.Grids: node
using Oceananigans.Utils: Iterate, launch!

function BickleyJetInitialConditions!(bickley_jet_parameters, bickley_jet_model)
    bickley_jet_grid = bickley_jet_model.grid
    arch = bickley_jet_grid.architecture
    Nx, Ny, Nz = size(bickley_jet_grid)
    
    @inline Ψ(y) = -tanh(y)
    @inline U(y) = sech(y)^2

    # A sinusoidal tracer
    @inline C(y, L) = sin(2π * y / L)

    # Slightly off-center vortical perturbations
    @inline ψ̃(x, y, ℓ, k) = exp(-(y + ℓ/10)^2 / 2ℓ^2) * cos(k * x) * cos(k * y)

    # Vortical velocity fields (ũ, ṽ) = (-∂_y, +∂_x) ψ̃
    @inline ũ(x, y, ℓ, k) = + ψ̃(x, y, ℓ, k) * (k * tan(k * y) + (y + ℓ/10) / ℓ^2)
    @inline ṽ(x, y, ℓ, k) = - ψ̃(x, y, ℓ, k) * k * tan(k * x)
    
    dr(x) = deg2rad(x)
    
    ϵ       = bickley_jet_parameters.ϵ
    ℓ       = bickley_jet_parameters.ℓ
    k       = bickley_jet_parameters.k
    a_jet   = bickley_jet_parameters.a_jet
    a_pert  = bickley_jet_parameters.a_pert
    b_pert  = bickley_jet_parameters.b_pert
    ψ₀_jet  = bickley_jet_parameters.ψ₀_jet
    ψ₀_pert = bickley_jet_parameters.ψ₀_pert

    @inline ψᵢ(λ, φ, z) = ψ₀_jet * Ψ(a_jet * dr(φ)) + ϵ * ψ₀_pert * ψ̃(b_pert * dr(λ), a_pert * dr(φ), ℓ, k)
    @inline cᵢ(λ, φ, z) = C(dr(φ), π)
    
    #####
    ##### Initial condition
    #####

    @info "Setting initial condition..."

    ψ = Field{Face, Face, Center}(bickley_jet_grid)
    
    set!(ψ, ψᵢ)
    
    # Note that set! fills only interior points; to compute u and v, we need information in the halo regions.
    fill_halo_regions!(ψ)

    # Note that fill_halo_regions! works for (Face, Face, Center) field, *except* for the two corner points that do not
    # correspond to an interior point! We need to manually fill the Face-Face halo points of the two corners that do not
    # have a corresponding interior point.

    @kernel function _set_ψ_missing_corners!(ψ, grid, region)
        k = @index(Global, Linear)
        i = 1
        j = grid.Ny+1
        if region in (1, 3, 5)
            λ, φ, z = node(i, j, k, grid, Face(), Face(), Center())
            @inbounds ψ[i, j, k] = ψᵢ(λ, φ, z)
        end
        i = grid.Nx+1
        j = 1
        if region in (2, 4, 6)
            λ, φ, z = node(i, j, k, grid, Face(), Face(), Center())
            @inbounds ψ[i, j, k] = ψᵢ(λ, φ, z)
        end
    end
    
    region = Iterate(1:6)
    @apply_regionally launch!(arch, bickley_jet_grid, (Nz,), _set_ψ_missing_corners!, ψ, bickley_jet_grid, region)

    @kernel function _set_initial_velocities!(ψ, u, v)
        i, j, k = @index(Global, NTuple)
        @inbounds u[i, j, k] = - ∂y(ψ)[i, j, k]
        @inbounds v[i, j, k] = + ∂x(ψ)[i, j, k]
    end

    @apply_regionally launch!(arch, bickley_jet_grid, (Nx, Ny, Nz), _set_initial_velocities!, ψ,
                              bickley_jet_model.velocities.u, bickley_jet_model.velocities.v)
    fill_halo_regions!((bickley_jet_model.velocities.u, bickley_jet_model.velocities.v))

    @kernel function _set_initial_tracer_distribution!(c, grid)
        i, j, k = @index(Global, NTuple)
        λ, φ, z = node(i, j, k, grid, Center(), Center(), Center())
        @inbounds c[i, j, k] = cᵢ(λ, φ, z)
    end

    @apply_regionally launch!(bickley_jet_model.architecture, bickley_jet_grid, (Nx, Ny, Nz),
                              _set_initial_tracer_distribution!, bickley_jet_model.tracers.c, bickley_jet_grid)

    fill_halo_regions!(bickley_jet_model.tracers.c)
end
