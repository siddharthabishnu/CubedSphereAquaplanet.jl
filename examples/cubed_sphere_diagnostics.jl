using KernelAbstractions: @kernel, @index

using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: halo_size, total_size
using Oceananigans.MultiRegion: number_of_regions, ConformalCubedSphereGridOfSomeKind
using Oceananigans.Operators: ζ₃ᶠᶠᶜ
using Oceananigans.Utils
# using Oceananigans.Utils: getregion, launch!

using Oceananigans: CubedSphereField

function compute_size_metrics(grid::ConformalCubedSphereGridOfSomeKind, field::CubedSphereField, ssh::Bool,
                              consider_all_levels::Bool, levels::UnitRange{Int}, read_parent_field_data::Bool)
    Nx, Ny, Nz = size(field)
    Hx, Hy, Hz = halo_size(grid)

    hx = read_parent_field_data ? Hx : 0
    hy = read_parent_field_data ? Hy : 0
    hz = read_parent_field_data ? Hz : 0
    
    if ssh
        consider_all_levels = false
        levels = read_parent_field_data ? (1:1) : (Nz+1:Nz+1)
    else
        if consider_all_levels
            levels = 1:Nz .+ hz
        else
            levels = levels .+ hz
        end
    end

    return Nx, Ny, Nz, hx, hy, hz, consider_all_levels, levels
end

function compute_vertical_index(grid::ConformalCubedSphereGridOfSomeKind, field::CubedSphereField, k::Int, ssh::Bool,
                                read_parent_field_data::Bool)
    Nz = size(field, 3)
    Hz = halo_size(grid)[3]
    hz = read_parent_field_data ? Hz : 0
    if ssh
        k = read_parent_field_data ? 1 : Nz + 1
    else
        k += hz
    end
    return k
end

function extract_field_time_series_array(grid, field_time_series;
                                         with_halos::Bool = false,
                                         ssh::Bool = false,
                                         consider_all_levels::Bool = true,
                                         levels::UnitRange{Int} = 1:1,
                                         read_parent_field_data::Bool = false,
                                         Δ::Int = 1)
    Nx, Ny, Nz, hx, hy, hz, consider_all_levels, levels = (
        compute_size_metrics(grid, field_time_series[1], ssh, consider_all_levels, levels, read_parent_field_data))
    Hx, Hy, Hz = halo_size(grid)
    
    n = length(field_time_series)
    m = floor(Int, (length(field_time_series) - 1)/Δ + 1)

    if with_halos
        field_time_series_array = zeros(Nx+2Hx, Ny+2Hy, length(levels), 6, m)
    else
        field_time_series_array = zeros(Nx, Ny, length(levels), 6, m)
    end

    j = 0
    for i in 1:Δ:n
        j += 1
        for region in 1:number_of_regions(grid)
            if with_halos
                field_time_series_array[:, :, :, region, j] = (
                    getregion(field_time_series[i], region)[1-Hx+hx:Nx+Hx+hx, 1-Hy+hy:Ny+Hy+hy, levels])
            else
                field_time_series_array[:, :, :, region, j] = (
                    getregion(field_time_series[i], region)[1+hx:Nx+hx, 1+hy:Ny+hy, levels])
            end
        end
    end
    
    return field_time_series_array
end

function interpolate_cubed_sphere_field_to_cell_centers(grid, field, field_location;
                                                        ssh::Bool = false,
                                                        consider_all_levels::Bool = true,
                                                        levels::UnitRange{Int} = 1:1,
                                                        read_parent_field_data::Bool = false)
    if field_location == "cc" && !read_parent_field_data
        return field
    end
    Nx, Ny, Nz, hx, hy, hz, consider_all_levels, levels = compute_size_metrics(
        grid, field, ssh, consider_all_levels, levels, read_parent_field_data)

    interpolated_field = Field{Center, Center, Center}(grid, indices = (:, :, levels))

    set!(interpolated_field, 0)

    @inbounds for region in 1:number_of_regions(grid), j in 1:Ny, i in 1:Nx
        dest = getregion(interpolated_field, region)
        src  = getregion(field, region)

        if field_location == "fc"
            dest[i, j, levels] = 0.5(src[i+hx, j+hy, levels] + src[i+1+hx, j+hy, levels])

        elseif field_location == "cf"
            dest[i, j, levels] = 0.5(src[i+hx, j+hy, levels] + src[i+hx, j+1+hy, levels])

        elseif field_location == "ff"
            dest[i, j, levels] = 0.25(src[i+hx, j+hy, levels] + src[i+1+hx, j+hy, levels] + src[i+1+hx, j+1+hy, levels]
                                      + src[i+hx, j+1+hy, levels])

        elseif field_location == "cc"
            dest[i, j, levels] = src[i+hx, j+hy, levels]
        end
    end

    return interpolated_field
end

function compute_vorticity_time_series(grid, u_time_series, v_time_series;
                                       interpolate_to_cell_centers::Bool = true)
    arch = grid.architecture
    ζ = Field{Face, Face, Center}(grid)

    @kernel function _compute_vorticity!(ζ, grid, u, v)
        i, j, k = @index(Global, NTuple)
        @inbounds ζ[i, j, k] = ζ₃ᶠᶠᶜ(i, j, k, grid, u, v)
    end

    offset = -1 .* halo_size(grid)

    ζ_time_series = Field[]
    
    for i in 1:length(u_time_series)
        u = u_time_series[i]
        v = v_time_series[i]
        fill_halo_regions!((u, v))
        @apply_regionally begin
            kernel_parameters = KernelParameters(total_size(ζ[1]), offset)
            launch!(arch, grid, kernel_parameters, _compute_vorticity!, ζ, grid, u, v)
        end
        if interpolate_to_cell_centers
            ζᶜᶜ = interpolate_cubed_sphere_field_to_cell_centers(grid, ζ, "ff")
            push!(ζ_time_series, deepcopy(ζᶜᶜ))
        else
            push!(ζ_time_series, deepcopy(ζ))
        end
    end
    
    return ζ_time_series
end