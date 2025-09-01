#=
using Oceananigans.Grids: λnode, φnode

# 4 Gaussians with width δR (degrees) and magnitude θ₀
δR = 2
θ₀ = 1

# Gaussian 1
i₁ = 1
j₁ = 1
panel = 1
λ₁ = λnode(i₁, j₁, grid[panel], Center(), Center())
φ₁ = φnode(i₁, j₁, grid[panel], Center(), Center())

# Gaussian 2
i₂ = Nx÷4 + 1
j₂ = 3Ny÷4 + 1
panel = 4
λ₂ = λnode(i₂, j₂, grid[panel], Center(), Center())
φ₂ = φnode(i₂, j₂, grid[panel], Center(), Center())

# Gaussian 3
i₃ = 3Nx÷4 + 1
j₃ = 3Ny÷4 + 1
panel = 3
λ₃ = λnode(i₃, j₃, grid[panel], Center(), Center())
φ₃ = φnode(i₃, j₃, grid[panel], Center(), Center())

# Gaussian 4
i₄ = 3Nx÷4+1
j₄ = 3Ny÷4+1
panel = 6
λ₄ = λnode(i₄, j₄, grid[panel], Center(), Center())
φ₄ = φnode(i₄, j₄, grid[panel], Center(), Center())

θᵢ(λ, φ, z) = θ₀ * exp(-((λ - λ₁)^2 + (φ - φ₁)^2) / 2δR^2) +
              θ₀ * exp(-((λ - λ₂)^2 + (φ - φ₂)^2) / 2δR^2) +
              θ₀ * exp(-((λ - λ₃)^2 + (φ - φ₃)^2) / 2δR^2) +
              θ₀ * exp(-((λ - λ₄)^2 + (φ - φ₄)^2) / 2δR^2)
=#

function TracerAdvectionInitialConditions!(tracer_advection_parameters, tracer_advection_model)
    tracer_advection_grid = tracer_advection_model.grid
    arch = tracer_advection_grid.architecture
    Nx, Ny, Nz = size(tracer_advection_grid)

    θ₀ = tracer_advection_parameters.θ₀
    Δφ = tracer_advection_parameters.Δφ
    
    @inline θᵢ(λ, φ, z) = θ₀ * cosd(4λ) * exp(-φ^2 / 2Δφ^2)

    #####
    ##### Initial condition
    #####

    @info "Setting initial condition..."
    set!(tracer_advection_model.tracers.θ, θᵢ)
    fill_halo_regions!(tracer_advection_model.tracers.θ)
end