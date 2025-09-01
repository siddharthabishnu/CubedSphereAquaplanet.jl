function BaroclinicWaveBoundaryConditions(baroclinic_wave_parameters)
    @inline ϕ²(i, j, k, grid, ϕ) = @inbounds ϕ[i, j, k]^2

    @inline speedᶠᶜᶜ(i, j, k, grid, u, v) = @inbounds sqrt(u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, ϕ², v))
    @inline speedᶜᶠᶜ(i, j, k, grid, u, v) = @inbounds sqrt(ℑxyᶜᶠᵃ(i, j, k, grid, ϕ², u) + v[i, j, k]^2)

    @inline u_drag(i, j, grid, clock, fields, p) = (
    @inbounds - p.Cᴰ * speedᶠᶜᶜ(i, j, 1, grid, fields.u, fields.v) * fields.u[i, j, 1])
    @inline v_drag(i, j, grid, clock, fields, p) = (
    @inbounds - p.Cᴰ * speedᶜᶠᶜ(i, j, 1, grid, fields.u, fields.v) * fields.v[i, j, 1])
    
    u_bot_bc = FluxBoundaryCondition(u_drag, discrete_form = true, parameters = (; Cᴰ = baroclinic_wave_parameters.Cᴰ))
    v_bot_bc = FluxBoundaryCondition(v_drag, discrete_form = true, parameters = (; Cᴰ = baroclinic_wave_parameters.Cᴰ))
    
    u_bcs = FieldBoundaryConditions(bottom = u_bot_bc)
    v_bcs = FieldBoundaryConditions(bottom = v_bot_bc)

    baroclinic_wave_boundary_conditions = (u = u_bcs, v = v_bcs)
    
    return baroclinic_wave_boundary_conditions
end
