using Oceananigans.BuoyancyFormulations: g_Earth
using Oceananigans.Coriolis: Ω_Earth
using Oceananigans.Grids: R_Earth

function BaroclinicWaveParameters(; unit_sphere::Bool = true)
    non_uniform_conformal_mapping = true
    # If true, applies a stretched conformal map (exponential or geometric spacing) to make the cubed-sphere grid more
    # uniform. This enlarges corner cells, relaxes the CFL constraint, and permits larger time steps.

    # Domain extents and resolution
    R = R_Earth  # Radius of the sphere (m)
    Lz = 3000    # Depth of the domain (m)
    Nx = 180     # Number of grid cells in the local x direction corresponding to 0.5° resolution
    Ny = 180     # Number of grid cells in the local y direction corresponding to 0.5° resolution
    Nz = 10      # Number of grid cells in the z direction
    H = 6        # Number of halo cells in each horizontal direction

    g = g_Earth  # Earth's gravitational acceleration (m s⁻²)
    Ω = Ω_Earth  # Earth's rotational rate (s⁻¹)

    Cᴰ = 1e-3    # Drag coefficient

    # Derived quantities
    c = sqrt(g * Lz)  # Gravity wave speed (m s⁻¹)

    # Initial condition parameters
    T0  = 30     # Reference surface conservative temperature (°C)
    φ0  = 45     # Central latitude of temperature front (degrees)
    Δφ  = 8      # Meridional half-width of temperature front (degrees)
    ϵT  = 1e-3   # Amplitude of small temperature perturbations (°C)

    S0  = 35     # Reference surface absolute salinity (g/kg)
    γS  = 3e-4   # Vertical salinity gradient (g/kg per meter)
    ϵS  = 1e-3   # Amplitude of small salinity perturbations (g/kg)

    baroclinic_wave_parameters = (
        non_uniform_conformal_mapping = non_uniform_conformal_mapping,
        R = R,    
        Lz = Lz,
        Nx = Nx,
        Ny = Ny,
        Nz = Nz,
        H = H,
        g = g,
        Ω = Ω,
        Cᴰ = Cᴰ,
        c = c,
        T0 = T0,
        φ0 = φ0,
        Δφ = Δφ,
        ϵT = ϵT,
        S0 = S0,
        γS = γS,
        ϵS = ϵS
    )

    return baroclinic_wave_parameters
end
