using KernelAbstractions: @kernel, @index
using Oceananigans.Grids: node
using Oceananigans.Utils: Iterate, launch!

function RossbyHaurwitzWaveInitialConditions!(rossby_haurwitz_wave_parameters, rossby_haurwitz_wave_model)
    rossby_haurwitz_wave_grid = rossby_haurwitz_wave_model.grid
    arch = rossby_haurwitz_wave_grid.architecture
    Nx, Ny, Nz = size(rossby_haurwitz_wave_grid)

    ## Rossby-Haurwitz initial condition from Williamson et al. (§3.6, 1992)
    ## # Here: θ ∈ [-π/2, π/2] is latitude and ϕ ∈ [0, 2π) is longitude.

    R  = rossby_haurwitz_wave_parameters.R
    Lz = rossby_haurwitz_wave_parameters.Lz

    n  = rossby_haurwitz_wave_parameters.n
    K  = rossby_haurwitz_wave_parameters.K
    ω  = rossby_haurwitz_wave_parameters.ω

    g  = rossby_haurwitz_wave_parameters.g
    Ω  = rossby_haurwitz_wave_parameters.Ω

    @inline A(θ) = (0.5ω * (2Ω + ω) * cos(θ)^2 
                    + 0.25K^2 * cos(θ)^(2n) * ((n + 1) * cos(θ)^2 + (2n^2 - n - 2) - 2n^2 * sec(θ)^2))
    @inline B(θ) = 2K * (Ω + ω) * ((n + 1) * (n + 2))^(-1) * cos(θ)^n * (n^2 + 2n + 2 - (n + 1)^2 * cos(θ)^2)
    @inline C(θ) = 0.25K^2 * cos(θ)^(2n) * ((n + 1) * cos(θ)^2 - (n + 2))

    @inline ψ_function(θ, ϕ) = -R^2 * ω * sin(θ) + R^2 * K * cos(θ)^n * sin(θ) * cos(n * ϕ)

    @inline u_function(θ, ϕ) =  R * ω * cos(θ) + R * K * cos(θ)^(n - 1) * (n * sin(θ)^2 - cos(θ)^2) * cos(n * ϕ)
    @inline v_function(θ, ϕ) = -n * K * R * cos(θ)^(n - 1) * sin(θ) * sin(n * ϕ)

    @inline h_function(θ, ϕ) = Lz + R^2/g * (A(θ) + B(θ) * cos(n * ϕ) + C(θ) * cos(2n * ϕ))

    # Initial conditions
    # Previously: θ ∈ [-π/2, π/2] is latitude and ϕ ∈ [0, 2π) is longitude
    # Oceananigans: ϕ ∈ [-90, 90] and λ ∈ [-180, 180]

    @inline rescale¹(λ) = (λ + 180) / 360 * 2π # λ to θ
    @inline rescale²(ϕ) = ϕ / 180 * π # θ to ϕ

    # Arguments were u(θ, ϕ), λ |-> ϕ, θ |-> ϕ
    #=
    @inline u₀(λ, ϕ, z) = u_function(rescale²(ϕ), rescale¹(λ))
    @inline v₀(λ, ϕ, z) = v_function(rescale²(ϕ), rescale¹(λ))
    =#
    @inline η₀(λ, ϕ, z) = h_function(rescale²(ϕ), rescale¹(λ)) - Lz
    #=
    set!(rossby_haurwitz_wave_model, u=u₀, v=v₀, η = η₀)
    =#
    @inline ψ₀(λ, φ, z) = ψ_function(rescale²(φ), rescale¹(λ))
    
    #####
    ##### Initial condition
    #####

    @info "Setting initial condition..."

    ψ = Field{Face, Face, Center}(rossby_haurwitz_wave_grid)

    set!(ψ, ψ₀)

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
            @inbounds ψ[i, j, k] = ψ₀(λ, φ, z)
        end
        i = grid.Nx+1
        j = 1
        if region in (2, 4, 6)
            λ, φ, z = node(i, j, k, grid, Face(), Face(), Center())
            @inbounds ψ[i, j, k] = ψ₀(λ, φ, z)
        end
    end
    
    region = Iterate(1:6)
    @apply_regionally launch!(arch, rossby_haurwitz_wave_grid, (Nz,), _set_ψ_missing_corners!, ψ,
                              rossby_haurwitz_wave_grid, region)

    @kernel function _set_initial_velocities!(ψ, u, v)
        i, j, k = @index(Global, NTuple)
        @inbounds u[i, j, k] = - ∂y(ψ)[i, j, k]
        @inbounds v[i, j, k] = + ∂x(ψ)[i, j, k]
    end

    @apply_regionally launch!(arch, rossby_haurwitz_wave_grid, (Nx, Ny, Nz), _set_initial_velocities!, ψ,
                              rossby_haurwitz_wave_model.velocities.u, rossby_haurwitz_wave_model.velocities.v)
    fill_halo_regions!((rossby_haurwitz_wave_model.velocities.u, rossby_haurwitz_wave_model.velocities.v))

    set!(rossby_haurwitz_wave_model.free_surface.η, η₀)
    fill_halo_regions!(rossby_haurwitz_wave_model.free_surface.η)
end
