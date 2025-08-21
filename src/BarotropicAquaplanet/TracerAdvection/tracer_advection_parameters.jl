function TracerAdvectionParameters()
    # Domain extents and resolution
    R = 1    # Radius of the sphere (m)
    Lz = 1   # Depth of the domain (m)
    Nx = 32  # Number of grid cells in the local x direction
    Ny = 32  # Number of grid cells in the local y direction
    Nz = 1   # Number of grid cells in the z direction
    H = 6    # Number of halo cells in each horizontal direction
    
    # Solid body rotation
    U = 1         # Velocity scale (m s⁻¹)
    φʳ = 0        # Latitude pierced by the axis of rotation
    α  = 90 - φʳ  # Angle between axis of rotation and north pole (degrees)
    θ₀ = 1        # Tracer amplitude
    Δφ = 20       # Range of the exponential tracer profile

    tracer_advection_parameters = (
        R = R,
        Lz = Lz,
        Nx = Nx,
        Ny = Ny,
        Nz = Nz,
        H = H,
        U = U,
        φʳ = φʳ,
        α = α,
        θ₀ = θ₀,
        Δφ = Δφ
    )

    return tracer_advection_parameters
end
