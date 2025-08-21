using KernelAbstractions: @kernel, @index
using Oceananigans.Utils: Iterate, launch!

function TracerAdvectionPrescribedVelocities!(tracer_advection_parameters, tracer_advection_grid)
    arch = tracer_advection_grid.architecture
    Nx, Ny, Nz = size(tracer_advection_grid)
    
    # Tracer parameters
    U = tracer_advection_parameters.U
    R = tracer_advection_parameters.R
    α = tracer_advection_parameters.α

    @inline ψᵣ(λ, φ, z) = - U * R * (sind(φ) * cosd(α) - cosd(λ) * cosd(φ) * sind(α))
    
    ψ = Field{Face, Face, Center}(tracer_advection_grid)
    set!(ψ, ψᵣ)

    # Note that set! fills only interior points; to compute u and v, we need information in the halo regions.
    fill_halo_regions!(ψ)
    
    # Note that fill_halo_regions! works for (Face, Face, Center) field, *except* for the two corner points that do not
    # correspond to an interior point! We need to manually fill the Face-Face halo points of the two corners that do not 
    # have a corresponding interior point.
    
    @kernel function _set_ψ_missing_corners!(ψ, grid, region)
        k = @index(Global, Linear)
        i = 1
        j = Ny+1
        if region in (1, 3, 5)
            λ, φ, z = node(i, j, k, grid, Face(), Face(), Center())
            @inbounds ψ[i, j, k] = ψᵣ(λ, φ, z)
        end
        i = Nx+1
        j = 1
        if region in (2, 4, 6)
            λ, φ, z = node(i, j, k, grid, Face(), Face(), Center())
            @inbounds ψ[i, j, k] = ψᵣ(λ, φ, z)
        end
    end
    
    region = Iterate(1:6)
    @apply_regionally launch!(arch, tracer_advection_grid, (Nz,), _set_ψ_missing_corners!, ψ, tracer_advection_grid,
                              region)

    u = XFaceField(tracer_advection_grid)
    v = YFaceField(tracer_advection_grid)
    
    @kernel function _set_prescribed_velocities!(ψ, u, v)
        i, j, k = @index(Global, NTuple)
        u[i, j, k] = - ∂y(ψ)[i, j, k]
        v[i, j, k] = + ∂x(ψ)[i, j, k]
    end

    @apply_regionally launch!(arch, tracer_advection_grid, (Nx, Ny, Nz), _set_prescribed_velocities!, ψ, u, v)

    return u, v
end