using Oceananigans.BuoyancyFormulations: g_Earth
using Oceananigans.Coriolis: Ω_Earth
using Oceananigans.Grids: R_Earth

function RossbyHaurwitzWaveParameters()
    # Domain extents and resolution
    R = R_Earth  # Radius of the sphere (m)
    Lz = 8000    # Depth of the domain (m)
    Nx = 128     # Number of grid cells in the local x direction
    Ny = 128     # Number of grid cells in the local y direction
    Nz = 1       # Number of grid cells in the z direction
    H = 6        # Number of halo cells in each horizontal direction
    
    # Rossby-Haurwitz wave initial condition parameters
    n = 4
    K = 7.848e-6
    ω = 0
    
    g = g_Earth  # Earth's gravitational acceleration (m s⁻²)
    Ω = Ω_Earth  # Earth's rotational rate (s⁻¹)

    rossby_haurwitz_wave_parameters = (
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
