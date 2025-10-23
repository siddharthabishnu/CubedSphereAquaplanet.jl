using CubedSphere: spherical_area_triangle
using Oceananigans.Grids: size, lat_lon_to_cartesian
using Oceananigans.Utils: getregion

function BickleyJetGrid(bickley_jet_parameters; arch = CPU(), FT::DataType = Oceananigans.defaults.FloatType)
    #####
    ##### Grid generation
    #####

    @info "Generating grid..."

    bickley_jet_grid = ConformalCubedSphereGrid(
        arch, FT;
        panel_size = (bickley_jet_parameters.Nx, bickley_jet_parameters.Ny, bickley_jet_parameters.Nz),
        z = (-bickley_jet_parameters.Lz, 0),
        radius = bickley_jet_parameters.R,
        horizontal_direction_halo = bickley_jet_parameters.H,
        non_uniform_conformal_mapping = bickley_jet_parameters.non_uniform_conformal_mapping,
        partition = CubedSpherePartition(; R = 1))

    Nx, Ny, Nz = size(bickley_jet_grid)

    a = lat_lon_to_cartesian(bickley_jet_grid[1].φᶜᶜᵃ[1, Ny], bickley_jet_grid[1].λᶜᶜᵃ[1, Ny], 1)
    b = lat_lon_to_cartesian(bickley_jet_grid[3].φᶜᶜᵃ[1, Ny], bickley_jet_grid[3].λᶜᶜᵃ[1, Ny], 1)
    c = lat_lon_to_cartesian(bickley_jet_grid[5].φᶜᶜᵃ[1, Ny], bickley_jet_grid[5].λᶜᶜᵃ[1, Ny], 1)

    Azᶠᶠᵃ = spherical_area_triangle(a, b, c) * bickley_jet_parameters.R^2
    for region in 1:6
        getregion(bickley_jet_grid, region).Azᶠᶠᵃ[1,    1]    = Azᶠᶠᵃ
        getregion(bickley_jet_grid, region).Azᶠᶠᵃ[Nx+1, 1]    = Azᶠᶠᵃ
        getregion(bickley_jet_grid, region).Azᶠᶠᵃ[Nx+1, Ny+1] = Azᶠᶠᵃ
        getregion(bickley_jet_grid, region).Azᶠᶠᵃ[1,    Ny+1] = Azᶠᶠᵃ
    end

    return bickley_jet_grid
end
