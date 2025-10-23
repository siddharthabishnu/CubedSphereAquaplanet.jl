using Oceananigans: defaults

function TracerAdvectionParameters()
    non_uniform_conformal_mapping = true
    # If true, applies a stretched conformal map (exponential or geometric spacing) to make the cubed-sphere grid more
    # uniform. This enlarges corner cells, relaxes the CFL constraint, and permits larger time steps.

    # Domain extents and resolution
    R = defaults.planet_radius  # Radius of the sphere (m)
    Lz = 1000                   # Depth of the domain (m)
    Nx = 128                    # Number of grid cells in the local x direction
    Ny = 128                    # Number of grid cells in the local y direction
    Nz = 1                      # Number of grid cells in the z direction
    H = 6                       # Number of halo cells in each horizontal direction

    # Solid body rotation
    U = R * 2π/86400  # Velocity scale (m s⁻¹)
    φʳ = 0            # Latitude pierced by the axis of rotation
    α  = 90 - φʳ      # Angle between axis of rotation and north pole (degrees)
    θ₀ = 1            # Tracer amplitude
    Δφ = 20           # Range of the exponential tracer profile

    g = defaults.gravitational_acceleration  # Earth's gravitational acceleration (m s⁻²)

    tracer_advection_parameters = (
        non_uniform_conformal_mapping = non_uniform_conformal_mapping,
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
        Δφ = Δφ,
        g = g
    )

    return tracer_advection_parameters
end
