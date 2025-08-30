using Oceananigans.BuoyancyFormulations: g_Earth
using Oceananigans.Coriolis: Ω_Earth
using Oceananigans.Grids: R_Earth

function BickleyJetParameters(; unit_sphere::Bool = true)
    non_uniform_conformal_mapping = true
    # If true, applies a stretched conformal map (exponential or geometric spacing) to make the cubed-sphere grid more
    # uniform. This enlarges corner cells, relaxes the CFL constraint, and permits larger time steps.

    # Domain extents and resolution
    R = R_Earth  # Radius of the sphere (m)
    Lz = 1000    # Depth of the domain (m)
    Nx = 128     # Number of grid cells in the local x direction
    Ny = 128     # Number of grid cells in the local y direction
    Nz = 1       # Number of grid cells in the z direction
    H = 6        # Number of halo cells in each horizontal direction

    g = g_Earth  # Earth's gravitational acceleration (m s⁻²)
    Ω = Ω_Earth  # Earth's rotational rate (s⁻¹)
    
    # Bickley jet initial condition parameters
    ϵ = 0.1      # Perturbation amplitude (non-dimensional multiplier of U_ref)
    ℓ = 0.5      # Gaussian envelope width (in angular y-coordinate units)
    k = 0.5      # Wavenumber used in cos(kx)cos(ky) factors (in angular units)

    if unit_sphere
        # Non-dimensional test case
        R  = 1                    # Nondimensional sphere radius
        Lz = 1                    # Nondimensional depth
        g  = 1                    # Nondimensional gravitational acceleration
        Ω  = 1                    # Nondimensional rotational rate

        a_jet  = 1                # Jet width parameter (L = R / a_jet)
        b_pert = 1 / k            # Sets zonal wave number: m = k * b_pert

    else
        # Earth-like dimensional case
        a_jet  = 24               # Jet width parameter (L = R / 24 ≈ 265 km), gives β* ~ O(1)
        b_pert = 12 / k           # Sets zonal wavenumber: m = k * b_pert = 24 (within unstable range)
    end

    U_ref   = 1                   # Target reference velocity scale (m s⁻¹)

    # Derived quantities
    a_pert  = 4 / k               # Meridional scaling → 2 cosine periods pole-to-pole
    c       = sqrt(g * Lz)        # Gravity wave speed (m s⁻¹)
    
    ψ₀_jet  = U_ref * R / a_jet   # Base jet streamfunction amplitude → peak jet speed ≈ U_ref
    ψ₀_pert = U_ref * R / a_pert  # Perturbation streamfunction amplitude (scales u_pert ≈ ϵ U_ref)

    bickley_jet_parameters = (
        non_uniform_conformal_mapping = non_uniform_conformal_mapping,
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
        a_jet = a_jet,
        a_pert = a_pert,
        b_pert = b_pert,
        c = c,
        ψ₀_jet = ψ₀_jet,
        ψ₀_pert = ψ₀_pert
    )

    return bickley_jet_parameters
end
