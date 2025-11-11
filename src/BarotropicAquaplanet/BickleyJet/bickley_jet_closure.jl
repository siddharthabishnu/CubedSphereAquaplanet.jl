
using Oceananigans.Operators: Δx, Δy
using Oceananigans.TurbulenceClosures: HorizontalVectorInvariantScalarBiharmonicDiffusivity

function BickleyJetClosure(bickley_jet_parameters)
    # Filter width squared, expressed as a harmonic mean of x and y spacings
    @inline Δ²ᶜᶜᶜ(i, j, k, grid, lx, ly, lz) =  2 * (1 / (1 / Δx(i, j, k, grid, lx, ly, lz)^2
                                                          + 1 / Δy(i, j, k, grid, lx, ly, lz)^2))

    # Use a biharmonic diffusivity for momentum. Define the diffusivity function as gridsize^4 divided by the timescale.
    @inline νhb(i, j, k, grid, lx, ly, lz, clock, fields, p) = Δ²ᶜᶜᶜ(i, j, k, grid, lx, ly, lz)^2 / p.λ_rts

    bickley_jet_horizontal_viscosity = 
        HorizontalVectorInvariantScalarBiharmonicDiffusivity(
            ν = νhb, discrete_form = true, parameters = (; λ_rts = bickley_jet_parameters.λ_rts))

    if bickley_jet_parameters.vector_invariant_momentum_advection
        bickley_jet_closure = bickley_jet_horizontal_viscosity
    else
        bickley_jet_closure = nothing
    end

    return bickley_jet_closure
end
