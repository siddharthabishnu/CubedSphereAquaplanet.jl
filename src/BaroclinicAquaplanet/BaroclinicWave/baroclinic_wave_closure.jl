function BaroclinicWaveClosure()
    νz = 1e-3
    κz = 1e-4

    baroclinic_wave_vertical_diffusivity = 
        VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), ν = νz, κ = κz)

    baroclinic_wave_convective_adjustment = 
        ConvectiveAdjustmentVerticalDiffusivity(VerticallyImplicitTimeDiscretization(), convective_κz = 1)

    baroclinic_wave_closure = (baroclinic_wave_vertical_diffusivity, baroclinic_wave_convective_adjustment)

    return baroclinic_wave_closure
end
