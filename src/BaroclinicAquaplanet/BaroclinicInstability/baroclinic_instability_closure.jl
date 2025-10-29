function BaroclinicInstabilityClosure()
    νz = 1e-3
    κz = 1e-4

    baroclinic_instability_vertical_diffusivity = 
        VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), ν = νz, κ = κz)

    baroclinic_instability_convective_adjustment = 
        ConvectiveAdjustmentVerticalDiffusivity(VerticallyImplicitTimeDiscretization(), convective_κz = 1)

    baroclinic_instability_closure = (baroclinic_instability_vertical_diffusivity,
                                      baroclinic_instability_convective_adjustment)

    return baroclinic_instability_closure
end
