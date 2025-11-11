using Oceananigans.Operators: Δx, Δy
using Oceananigans.TurbulenceClosures: HorizontalVectorInvariantScalarBiharmonicDiffusivity

function BaroclinicInstabilityClosure(baroclinic_instability_parameters)
    νz = 1e-3
    κz = 1e-4

    # Filter width squared, expressed as a harmonic mean of x and y spacings
    @inline Δ²ᶜᶜᶜ(i, j, k, grid, lx, ly, lz) =  2 * (1 / (1 / Δx(i, j, k, grid, lx, ly, lz)^2
                                                          + 1 / Δy(i, j, k, grid, lx, ly, lz)^2))

    # Use a biharmonic diffusivity for momentum. Define the diffusivity function as gridsize^4 divided by the timescale.
    @inline νhb(i, j, k, grid, lx, ly, lz, clock, fields, p) = Δ²ᶜᶜᶜ(i, j, k, grid, lx, ly, lz)^2 / p.λ_rts

    baroclinic_instability_horizontal_viscosity =
        HorizontalVectorInvariantScalarBiharmonicDiffusivity(
            ν = νhb, discrete_form = true, parameters = (; λ_rts = baroclinic_instability_parameters.λ_rts))

    baroclinic_instability_vertical_diffusivity = 
        VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), ν = νz, κ = κz)

    baroclinic_instability_convective_adjustment = 
        ConvectiveAdjustmentVerticalDiffusivity(VerticallyImplicitTimeDiscretization(), convective_κz = 1)

    if baroclinic_instability_parameters.vector_invariant_momentum_advection
        baroclinic_instability_closure = (baroclinic_instability_horizontal_viscosity,
                                          baroclinic_instability_vertical_diffusivity,
                                          baroclinic_instability_convective_adjustment)
    else
        baroclinic_instability_closure = (baroclinic_instability_vertical_diffusivity,
                                          baroclinic_instability_convective_adjustment)
    end

    return baroclinic_instability_closure
end
