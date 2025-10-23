using Oceananigans: defaults

function RossbyHaurwitzWaveParameters()
    non_uniform_conformal_mapping = true
    # If true, applies a stretched conformal map (exponential or geometric spacing) to make the cubed-sphere grid more
    # uniform. This enlarges corner cells, relaxes the CFL constraint, and permits larger time steps.

    # Domain extents and resolution
    R = defaults.planet_radius  # Radius of the sphere (m)
    Lz = 8000                   # Depth of the domain (m)
    Nx = 128                    # Number of grid cells in the local x direction
    Ny = 128                    # Number of grid cells in the local y direction
    Nz = 1                      # Number of grid cells in the z direction
    H = 6                       # Number of halo cells in each horizontal direction
    
    # Rossby-Haurwitz wave initial condition parameters
    n = 4
    K = 7.848e-6
    ω = 0

    g = defaults.gravitational_acceleration  # Earth's gravitational acceleration (m s⁻²)
    Ω = defaults.planet_rotation_rate        # Earth's rotational rate (s⁻¹)

    rossby_haurwitz_wave_parameters = (
        non_uniform_conformal_mapping = non_uniform_conformal_mapping,
        R = R,    
        Lz = Lz,
        Nx = Nx,
        Ny = Ny,
        Nz = Nz,
        H = H,
        n = n,
        K = K,
        ω = ω,
        g = g,
        Ω = Ω
    )

    return rossby_haurwitz_wave_parameters
end
