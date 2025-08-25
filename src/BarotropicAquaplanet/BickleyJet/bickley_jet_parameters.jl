using Oceananigans.BuoyancyFormulations: g_Earth
using Oceananigans.Coriolis: Ω_Earth
using Oceananigans.Grids: R_Earth

function BickleyJetParameters()
    # Domain extents and resolution
    R = R_Earth  # Radius of the sphere (m)
    Lz = 1000    # Depth of the domain (m)
    Nx = 64      # Number of grid cells in the local x direction
    Ny = 64      # Number of grid cells in the local y direction
    Nz = 1       # Number of grid cells in the z direction
    H = 6        # Number of halo cells in each horizontal direction

    g = g_Earth  # Earth's gravitational acceleration (m s⁻²)
    Ω = Ω_Earth  # Earth's rotational rate (s⁻¹)
    
    # Bickley jet initial condition parameters
    ϵ = 0.1      # Non-dimensional perturbation amplitude
    ℓ = 0.5      # Gaussian envelope width (in angular y-coordinate units)
    k = 0.5      # Sinusoidal wavenumber (in angular x, y units)

    a = 4        # Meridional angular scaling factor (controls jet width in latitude)
    b = 2        # Zonal angular scaling factor (controls number of waves in longitude)

    Fr_target = 0.01       # Target Froude number U / sqrt(gH), used to set velocity scale

    # Derived dimensional quantities
    
    c = sqrt(g * Lz)       # Gravity-wave speed (m s⁻¹)
    U_ref = Fr_target * c  # Reference velocity scale (m s⁻¹)
    ψ₀ = U_ref * R / a     # Streamfunction amplitude (m² s⁻¹), sets peak jet speed

    bickley_jet_parameters = (
        R = R,    
        Lz = Lz,
        Nx = Nx,
        Ny = Ny,
        Nz = Nz,
        H = H,
        g = g,
        Ω = Ω,
        ϵ = ϵ,
        ℓ = ℓ,
        k = k,
        a = a,
        b = b,
        Fr_target = Fr_target,
        c = c,
        U_ref = U_ref,
        ψ₀ = ψ₀
    )

    return bickley_jet_parameters
end
