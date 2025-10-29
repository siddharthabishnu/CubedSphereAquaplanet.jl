using KernelAbstractions: @kernel, @index
using Oceananigans.Grids: node
using Oceananigans.Operators: ∇²hᶜᶜᶜ
using Oceananigans.Utils: Iterate, launch!

# --- Diffusion step for Gaussian smoothing: cⁿ⁺¹ = cⁿ + (Δτν) ∇²cⁿ ---
# Pass the *effective* step (Δτν = νΔτ) to satisfy FTCS stability without re-multiplying by ν inside the kernel.
@kernel function _apply_gaussian_filter!(c_old, c_new, grid, Δτν)
    i, j, k = @index(Global, NTuple)
    # Horizontal Laplacian on the conformal cubed sphere (metric-aware)
    ∇²c_old = ∇²hᶜᶜᶜ(i, j, k, grid, c_old)
    @inbounds c_new[i, j, k] = c_old[i, j, k] + Δτν * ∇²c_old
end

@kernel function _area_weighted_field!(c, c_area_weighted, Azᶜᶜᵃ, A, k)
    i, j = @index(Global, NTuple)
    @inbounds c_area_weighted[i, j, k] = c[i, j, k] * Azᶜᶜᵃ[i, j] / A
end

# --- Enforce horizontal mean per level (protect against drift) ---
# c̄ is the current area-weighted mean at level k; c̄₀ is the target mean.
@kernel function _enforce_horizontal_mean!(c, c̄, c̄₀, k)
    i, j = @index(Global, NTuple)
    δ = c̄ - c̄₀
    @inbounds c[i, j, k] -= δ
end

"""
    apply_gaussian_filter!(c, params)

Apply Gaussian smoothing to scalar field `c` by integrating the horizontal diffusion
equation up to a total "diffusion time" τ_target, such that

    σ² = 2ντ_target = L²   (with e-folding length L).

# Parameters
- `L`        :: target smoothing length [m] (informational if τ_target is given directly)
- `ν`        :: (artificial) diffusivity, used to define τ_target via L² = 2ντ_target
- `safety`   :: 0 < safety ≤ 0.5, CFL margin for explicit FTCS
- `τ_target` :: total diffusion time to reach desired σ (e.g., L² / (2ν))
"""
function apply_gaussian_filter!(c, params)
    grid   = c.grid
    arch   = grid.architecture

    L        = params.L
    ν        = params.ν
    safety   = params.safety
    τ_target = params.τ_target

    # --- Stable explicit step (FTCS) ---
    # For equal spacings in x and y, FTCS stability requires νΔτ ≤ Δs_min² / 4.
    # We therefore advance with Δτν = νΔτ = safety * Δs_min² / 4.
    # If spacings differ, use the standard bound with 2/(Δx²)+2/(Δy²) cellwise and take the minimum.
    Δs_min = min(minimum_xspacing(grid, Face(), Face(), Center()),
                 minimum_yspacing(grid, Face(), Face(), Center()))
    Δτν = 0.25 * safety * Δs_min^2  # Effective step for ∇² (already includes ν)

    # Number of diffusion steps so that ΣΔτ = τ_target ⇒ ΣΔτν = ντ_target
    Nsteps = ceil(Int, (ν * τ_target) / Δτν)

    # --- Precompute target horizontal means per level (area-weighted) ---
    Nx, Ny, Nz = size(grid)
    c̄₀ = [] # Vector{Float64}(undef, Nz)
    c_area_weighted = CenterField(grid)
    A = 4π * params.R^2 / 6  # Surface area of a cubed sphere panel

    for k in 1:Nz
        @apply_regionally launch!(arch, grid, (Nx, Ny), _area_weighted_field!, c, c_area_weighted, grid.Azᶜᶜᵃ, A, k)
        @apply_regionally c̄₀_k = sum(view(c_area_weighted, :, :, k))  # Target mean
        push!(c̄₀, c̄₀_k)
    end

    c_old = deepcopy(c)
    c_new = CenterField(grid)

    # --- Diffusive Gaussian smoothing loop ---
    for n = 1:Nsteps
        # Halos must be up-to-date for metric Laplacian across panel seams
        fill_halo_regions!(c_old)

        # Advance one explicit diffusion step (panel-parallel)
        @apply_regionally launch!(arch, grid, (Nx, Ny, Nz), _apply_gaussian_filter!, c_old, c_new, grid, Δτν)

        # Re-enforce horizontal mean at each level to avoid drift
        c̄ = []
        for k in 1:Nz
            @apply_regionally launch!(arch, grid, (Nx, Ny), _area_weighted_field!, c_new, c_area_weighted, grid.Azᶜᶜᵃ,
                                      A, k)
            @apply_regionally c̄_k = sum(view(c_area_weighted, :, :, k))
            push!(c̄, c̄_k)
            @apply_regionally launch!(arch, grid, (Nx, Ny), _enforce_horizontal_mean!, c_new, c̄_k, c̄₀[k], k)
        end

        c_old, c_new = c_new, c_old
    end

    c = deepcopy(c_old)
end

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

    apply_gaussian_filter!(baroclinic_instability_model.tracers.T, baroclinic_instability_parameters)
    apply_gaussian_filter!(baroclinic_instability_model.tracers.S, baroclinic_instability_parameters)

    fill_halo_regions!(baroclinic_instability_model.tracers.T)
    fill_halo_regions!(baroclinic_instability_model.tracers.S)
end
