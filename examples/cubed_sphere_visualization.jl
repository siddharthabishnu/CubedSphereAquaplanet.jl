using Oceananigans
using Oceananigans.MultiRegion: number_of_regions
using Oceananigans.Utils: getregion
#=
Install Imaginocean.jl in your home environment (not the current project):
using Pkg; Pkg.add(url="https://github.com/navidcy/Imaginocean.jl", rev="main")
=#
using CairoMakie, Imaginocean

import Imaginocean: heatlatlon!

using Oceananigans: location, CubedSphereField
using Imaginocean: get_grid, get_lat_lon_nodes_and_vertices

include("cubed_sphere_diagnostics.jl")

function heatlatlon!(ax::Axis, field, k_index::Int=1; suppress_shading_warning::Bool=true, kwargs...)
    LX, LY, LZ = location(field)

    grid = get_grid(field)

    _, (λvertices, φvertices) = get_lat_lon_nodes_and_vertices(grid, LX(), LY(), LZ())

    quad_points = vcat([Point2.(λvertices[:, i, j], φvertices[:, i, j])
                        for i in axes(λvertices, 2), j in axes(λvertices, 3)]...)
    quad_faces = vcat([begin; j = (i-1) * 4 + 1; [j j+1  j+2; j+2 j+3 j]; end for i in 1:length(quad_points)÷4]...)
    colors_per_point = vcat(fill.(vec(interior(field, :, :, k_index)), 4)...)
    shading = suppress_shading_warning ? Makie.NoShading : false

    mesh!(ax, quad_points, quad_faces; color = colors_per_point, shading, kwargs...)
end

function heatlatlon!(ax::Axis, field::CubedSphereField, k_index::Int=1; suppress_shading_warning::Bool=true, kwargs...)
    apply_regionally!(heatlatlon!, ax, field, k_index; suppress_shading_warning, kwargs...)

    xlims!(ax, (-180, 180))
    ylims!(ax, (-90, 90))
end

heatlatlon!(ax::Axis, field::Observable{<:CubedSphereField}, k_index::Int=1;
            suppress_shading_warning::Bool=true, kwargs...) = (
    heatlatlon!(ax, field.val, k_index; suppress_shading_warning, kwargs...))

function specify_colorrange(grid, field;
                            k::Int = 1,
                            ssh::Bool = false,
                            consider_all_levels::Bool = true,
                            levels::UnitRange{Int} = k:k,
                            read_parent_field_data::Bool = false,
                            use_symmetric_colorrange::Bool = true)
    Nx, Ny, Nz, hx, hy, hz, consider_all_levels, levels = compute_size_metrics(
        grid, field, ssh, consider_all_levels, levels, read_parent_field_data)
    
    field_array = zeros(Nx, Ny, length(levels), 6)
    for region in 1:number_of_regions(grid)
        field_array[:, :, :, region] = getregion(field, region)[1+hx:Nx+hx, 1+hy:Ny+hy, levels]
    end
    
    field_maximum = maximum(field_array)
    field_minimum = minimum(field_array)
    
    if use_symmetric_colorrange
        colorrange_limit = max(abs(field_maximum), abs(field_minimum))
        colorrange = [-colorrange_limit, colorrange_limit]
    else
        colorrange = [field_minimum, field_maximum]
    end
    
    return colorrange
end

function specify_colorrange_time_series(grid, field_time_series;
                                        ssh::Bool = false,
                                        consider_all_levels::Bool = true,
                                        levels::UnitRange{Int} = 1:1,
                                        read_parent_field_data::Bool = false,
                                        Δ::Int = 1,
                                        use_symmetric_colorrange::Bool = true)
    field_time_series_array = (
        extract_field_time_series_array(grid, field_time_series;
                                        ssh, consider_all_levels, levels, read_parent_field_data, Δ))

    field_maximum = maximum(field_time_series_array)
    field_minimum = minimum(field_time_series_array)

    if use_symmetric_colorrange
        colorrange_limit = max(abs(field_maximum), abs(field_minimum))
        colorrange = [-colorrange_limit, colorrange_limit]
    else
        colorrange = [field_minimum, field_maximum]
    end
    
    return colorrange
end

function panelwise_visualization(grid, field;
                                 with_halos::Bool = false,
                                 k::Int = 1,
                                 use_symmetric_colorrange::Bool = true,
                                 ssh::Bool = false,
                                 consider_all_levels::Bool = true,
                                 levels::UnitRange{Int} = k:k,
                                 read_parent_field_data::Bool = false,
                                 colorrange::Union{Nothing, Vector} = nothing,
                                 colormap::Union{Nothing, Symbol} = nothing)
    fig = Figure(size = (2450, 1400))

    axis_kwargs = (xlabelsize = 22.5, ylabelsize = 22.5, xlabelpadding = 10, ylabelpadding = 10, xticklabelsize = 17.5,
                   yticklabelsize = 17.5, xticklabelpad = 20, yticklabelpad = 20, aspect = 1.0, titlesize = 27.5,
                   titlegap = 15, titlefont = :bold, xlabel = "Local x direction", ylabel = "Local y direction")

    if isnothing(colorrange)
        colorrange = specify_colorrange(grid, field;
                                        k, ssh, consider_all_levels, levels, read_parent_field_data,
                                        use_symmetric_colorrange)
    end

    colormap = something(colormap, use_symmetric_colorrange ? :balance : :amp)

    Nx, Ny, Nz, hx, hy, hz, consider_all_levels, levels = compute_size_metrics(
        grid, field, ssh, consider_all_levels, levels, read_parent_field_data)
    k = compute_vertical_index(grid, field, k, ssh, read_parent_field_data)

    function slice_panel_data(f, panel)
        if with_halos
            data = getregion(f, panel)[:, :, k]
        else
            data = getregion(f, panel)[1+hx:Nx+hx, 1+hy:Ny+hy, k]
        end
        return parent(data)
    end

    panel_positions = [(3, 1), (3, 3), (2, 3), (2, 5), (1, 5), (1, 7)]

    for (i, pos) in enumerate(panel_positions)
        ax = Axis(fig[pos...]; title = "Panel $i", axis_kwargs...)
        hm = heatmap!(ax, slice_panel_data(field, i); colorrange, colormap)
        Colorbar(fig[pos[1], pos[2] + 1], hm)
    end

    return fig
end

function geo_heatlatlon_visualization(grid, field, field_location, title;
                                      k::Int = 1,
                                      use_symmetric_colorrange::Bool = true,
                                      ssh::Bool = false,
                                      consider_all_levels::Bool = true,
                                      levels::UnitRange{Int} = k:k,
                                      read_parent_field_data::Bool = false,
                                      colorbarlabel::String = "",
                                      colorrange::Union{Nothing, Vector} = nothing,
                                      colormap::Union{Nothing, Symbol} = nothing)
    fig = Figure(size = (2700, 1300))

    axis_kwargs = (xlabelsize = 37.5, ylabelsize = 37.5, xlabelpadding = 25, ylabelpadding = 25, xticklabelsize = 32.5,
                   yticklabelsize = 32.5, xticklabelpad = 20, yticklabelpad = 20, titlesize = 45, titlegap = 30,
                   titlefont = :bold)

    interpolated_field = interpolate_cubed_sphere_field_to_cell_centers(grid, field, field_location;
                                                                        ssh, consider_all_levels, levels,
                                                                        read_parent_field_data)

    if isnothing(colorrange)
        colorrange = specify_colorrange(grid, interpolated_field;
                                        k, ssh, consider_all_levels, levels, read_parent_field_data,
                                        use_symmetric_colorrange)
    end

    colormap = something(colormap, use_symmetric_colorrange ? :balance : :amp)

    ax = Axis(fig[1, 1]; title, axis_kwargs...)
    heatlatlon!(ax, interpolated_field, k; colorrange, colormap)

    Colorbar(fig[1, 2]; limits = colorrange, colormap, label = colorbarlabel, labelsize = 37.5, labelpadding = 25,
             ticklabelsize = 30, ticksize = 22.5, width = 35, height = Relative(1))
    colsize!(fig.layout, 1, Auto(0.8))
    colgap!(fig.layout, 75)

    return fig
end

function panelwise_visualization_animation_Makie(grid, field_time_series;
                                                 with_halos::Bool = false,
                                                 start_index::Int = 1,
                                                 k::Int = 1,
                                                 use_symmetric_colorrange::Bool = true,
                                                 ssh::Bool = false,
                                                 consider_all_levels::Bool = true,
                                                 levels::UnitRange{Int} = k:k,
                                                 read_parent_field_data::Bool = false,
                                                 Δ::Int = 1,
                                                 colorrange::Union{Nothing, Vector} = nothing,
                                                 colormap::Union{Nothing, Symbol} = nothing,
                                                 framerate::Int = 10,
                                                 output_directory::AbstractString = "output_directory",
                                                 filename::AbstractString = "filename",
                                                 format::AbstractString = ".mp4")
    field_time_series_array = extract_field_time_series_array(grid, field_time_series;
                                                              with_halos, ssh, consider_all_levels, levels,
                                                              read_parent_field_data, Δ)
    k = compute_vertical_index(grid, field_time_series[1], k, ssh, read_parent_field_data)

    # Observables
    n = Observable(start_index) # the current index
    field = @lift field_time_series_array[:, :, :, :, $n]

    # Create the initial visualization.
    fig = Figure(size = (2450, 1400))

    axis_kwargs = (xlabelsize = 22.5, ylabelsize = 22.5, xlabelpadding = 10, ylabelpadding = 10, xticklabelsize = 17.5,
                   yticklabelsize = 17.5, xticklabelpad = 20, yticklabelpad = 20, aspect = 1.0, titlesize = 27.5,
                   titlegap = 15, titlefont = :bold, xlabel = "Local x direction", ylabel = "Local y direction")

    if isnothing(colorrange)
        colorrange = specify_colorrange_time_series(grid, field_time_series;
                                                    ssh, consider_all_levels, levels, read_parent_field_data, Δ,
                                                    use_symmetric_colorrange)
    end

    colormap = something(colormap, use_symmetric_colorrange ? :balance : :amp)

    panel_positions = [(3, 1), (3, 3), (2, 3), (2, 5), (1, 5), (1, 7)]

    for (i, pos) in enumerate(panel_positions)
        ax = Axis(fig[pos...]; title = "Panel $i", axis_kwargs...)
        hm = heatmap!(ax, parent(field[][:, :, k, i]); colorrange, colormap)
        Colorbar(fig[pos[1], pos[2] + 1], hm)
    end

    frames = start_index:size(field_time_series_array, 5)
    CairoMakie.record(fig, joinpath(output_directory, filename * format), frames, framerate = framerate) do i
        print("Plotting frame $i of $(frames[end]) \r")
        field[] = field_time_series_array[:, :, :, :, i]
        for (i, pos) in enumerate(panel_positions)
            ax = Axis(fig[pos...]; title = "Panel $i", axis_kwargs...)
            hm = heatmap!(ax, parent(field[][:, :, k, i]); colorrange, colormap)
            Colorbar(fig[pos[1], pos[2] + 1], hm)
        end
    end
end

function panelwise_visualization_animation_frames(grid, field_time_series;
                                                  with_halos::Bool = false,
                                                  start_index::Int = 1,
                                                  k::Int = 1,
                                                  use_symmetric_colorrange::Bool = true,
                                                  ssh::Bool = false,
                                                  consider_all_levels::Bool = true,
                                                  levels::UnitRange{Int} = k:k,
                                                  read_parent_field_data::Bool = false,
                                                  Δ::Int = 1,
                                                  colorrange::Union{Nothing, Vector} = nothing,
                                                  colormap::Union{Nothing, Symbol} = nothing,
                                                  output_directory::AbstractString = "output_directory",
                                                  filename::AbstractString = "filename",
                                                  format::AbstractString = ".png")
    field_time_series_array = extract_field_time_series_array(grid, field_time_series;
                                                              with_halos, ssh, consider_all_levels, levels,
                                                              read_parent_field_data, Δ)
    k = compute_vertical_index(grid, field_time_series[1], k, ssh, read_parent_field_data)

    fig = Figure(size = (2450, 1400))

    axis_kwargs = (xlabelsize = 22.5, ylabelsize = 22.5, xlabelpadding = 10, ylabelpadding = 10, xticklabelsize = 17.5,
                   yticklabelsize = 17.5, xticklabelpad = 20, yticklabelpad = 20, aspect = 1.0, titlesize = 27.5,
                   titlegap = 15, titlefont = :bold, xlabel = "Local x direction", ylabel = "Local y direction")

    if isnothing(colorrange)
        colorrange = specify_colorrange_time_series(grid, field_time_series;
                                                    ssh, consider_all_levels, levels, read_parent_field_data, Δ,
                                                    use_symmetric_colorrange)
    end

    colormap = something(colormap, use_symmetric_colorrange ? :balance : :amp)

    panel_positions = [(3, 1), (3, 3), (2, 3), (2, 5), (1, 5), (1, 7)]

    frames = start_index:size(field_time_series_array, 5)
    for i in frames
        print("Plotting frame $i of $(frames[end]) \r")
        field = field_time_series_array[:, :, :, :, i]
        for (j, pos) in enumerate(panel_positions)
            ax = Axis(fig[pos...]; title = "Panel $j", axis_kwargs...)
            hm = heatmap!(ax, parent(field[:, :, k, j]); colorrange, colormap)
            Colorbar(fig[pos[1], pos[2] + 1], hm)
        end
        save(joinpath(output_directory, filename * "_$i" * format), fig)
    end
end

function panelwise_visualization_animation(grid, field_time_series;
                                           plot_frames::Bool = false,
                                           with_halos::Bool = false,
                                           start_index::Int = 1,
                                           k::Int = 1,
                                           use_symmetric_colorrange::Bool = true,
                                           ssh::Bool = false,
                                           consider_all_levels::Bool = true,
                                           levels::UnitRange{Int} = k:k,
                                           read_parent_field_data::Bool = false,
                                           Δ::Int = 1,
                                           colorrange::Union{Nothing, Vector} = nothing,
                                           colormap::Union{Nothing, Symbol} = nothing,
                                           framerate::Int = 10,
                                           output_directory::AbstractString = "output_directory",
                                           filename::AbstractString = "filename",
                                           figure_format::AbstractString = ".png",
                                           movie_format::AbstractString = ".mp4")
    if plot_frames
        format = figure_format
        path = joinpath(pwd(), output_directory * "/" * filename * "_frames")
        isdir(path) || mkdir(path)
        filename = joinpath(filename * "_frames/", filename)
        panelwise_visualization_animation_frames(grid, field_time_series;
                                                 with_halos, start_index, k, use_symmetric_colorrange, ssh,
                                                 consider_all_levels, levels, read_parent_field_data, Δ, colorrange,
                                                 colormap, output_directory, filename, format)
    else
        format = movie_format
        panelwise_visualization_animation_Makie(grid, field_time_series;
                                                with_halos, start_index, k, use_symmetric_colorrange, ssh,
                                                consider_all_levels, levels, read_parent_field_data, Δ, colorrange,
                                                colormap, framerate, output_directory, filename, format)
    end
end

function geo_heatlatlon_visualization_animation_Makie(grid, field_time_series, field_location, prettytimes,
                                                      title_prefix;
                                                      start_index::Int = 1,
                                                      k::Int = 1,
                                                      use_symmetric_colorrange::Bool = true,
                                                      ssh::Bool = false,
                                                      consider_all_levels::Bool = true,
                                                      levels::UnitRange{Int} = k:k,
                                                      read_parent_field_data::Bool = false,
                                                      Δ::Int = 1,
                                                      colorrange::Union{Nothing, Vector} = nothing,
                                                      colormap::Union{Nothing, Symbol} = nothing,
                                                      colorbarlabel::String = "",
                                                      framerate::Int = 10,
                                                      output_directory::AbstractString = "output_directory",
                                                      filename::AbstractString = "filename",
                                                      format::AbstractString = ".mp4")
    # Observables
    n = Observable(start_index) # the current index
    field = @lift field_time_series[$n]
    prettytime = @lift prettytimes[$n]

    # Create the initial visualization.
    fig = Figure(size = (2700, 1300))
    
    axis_kwargs = (xlabelsize = 37.5, ylabelsize = 37.5, xlabelpadding = 25, ylabelpadding = 25, xticklabelsize = 32.5,
                   yticklabelsize = 32.5, xticklabelpad = 20, yticklabelpad = 20, titlesize = 45, titlegap = 30,
                   titlefont = :bold)

    if isnothing(colorrange)
        colorrange = specify_colorrange_time_series(grid, field_time_series;
                                                    ssh, consider_all_levels, levels, read_parent_field_data, Δ,
                                                    use_symmetric_colorrange)
    end

    colormap = something(colormap, use_symmetric_colorrange ? :balance : :amp)

    ax = Axis(fig[1, 1]; axis_kwargs...)
    ax.title = title_prefix * " after " * prettytime[]
    interpolated_field = interpolate_cubed_sphere_field_to_cell_centers(grid, field[], field_location;
                                                                        ssh, consider_all_levels, levels,
                                                                        read_parent_field_data)
    heatlatlon!(ax, interpolated_field, k; colorrange, colormap)

    Colorbar(fig[1, 2]; limits = colorrange, colormap, label = colorbarlabel, labelsize = 37.5, labelpadding = 25,
             ticklabelsize = 30, ticksize = 22.5, width = 35, height = Relative(1))
    colsize!(fig.layout, 1, Auto(0.8))
    colgap!(fig.layout, 75)

    frames = start_index:length(field_time_series)
    CairoMakie.record(fig, joinpath(output_directory, filename * format), frames, framerate = framerate) do i
        print("Plotting frame $i of $(frames[end]) \r")

        field[] = field_time_series[i]
        prettytime[] = prettytimes[i]

        # Update the title of the plot.
        ax.title = title_prefix * " after " * prettytime[]

        # Update the plot.
        interpolated_field = interpolate_cubed_sphere_field_to_cell_centers(grid, field[], field_location;
                                                                            ssh, consider_all_levels, levels,
                                                                            read_parent_field_data)
        heatlatlon!(ax, interpolated_field, k; colorrange, colormap)

        Colorbar(fig[1, 2]; limits = colorrange, colormap, label = colorbarlabel, labelsize = 37.5, labelpadding = 25,
                 ticklabelsize = 30, ticksize = 22.5, width = 35, height = Relative(1))
        colsize!(fig.layout, 1, Auto(0.8))
        colgap!(fig.layout, 75)
    end
end

function geo_heatlatlon_visualization_animation_frames(grid, field_time_series, field_location, prettytimes, 
                                                       title_prefix;
                                                       start_index::Int = 1,
                                                       k::Int = 1,
                                                       use_symmetric_colorrange::Bool = true,
                                                       ssh::Bool = false,
                                                       consider_all_levels::Bool = true,
                                                       levels::UnitRange{Int} = k:k,
                                                       read_parent_field_data::Bool = false,
                                                       Δ::Int = 1,
                                                       colorrange::Union{Nothing, Vector} = nothing,
                                                       colormap::Union{Nothing, Symbol} = nothing,
                                                       colorbarlabel::String = "",
                                                       output_directory::AbstractString = "output_directory",
                                                       filename::AbstractString = "filename",
                                                       format::AbstractString = ".png")
    # Create the initial visualization.
    fig = Figure(size = (2700, 1300))

    axis_kwargs = (xlabelsize = 37.5, ylabelsize = 37.5, xlabelpadding = 25, ylabelpadding = 25, xticklabelsize = 32.5,
                   yticklabelsize = 32.5, xticklabelpad = 20, yticklabelpad = 20, titlesize = 45, titlegap = 30,
                   titlefont = :bold)

    if isnothing(colorrange)
        colorrange = specify_colorrange_time_series(grid, field_time_series;
                                                    ssh, consider_all_levels, levels, read_parent_field_data, Δ,
                                                    use_symmetric_colorrange)
    end

    colormap = something(colormap, use_symmetric_colorrange ? :balance : :amp)

    ax = Axis(fig[1, 1]; axis_kwargs...)

    Colorbar(fig[1, 2]; limits = colorrange, colormap, label = colorbarlabel, labelsize = 37.5, labelpadding = 25,
             ticklabelsize = 30, ticksize = 22.5, width = 35, height = Relative(1))
    colsize!(fig.layout, 1, Auto(0.8))
    colgap!(fig.layout, 75)

    frames = start_index:length(field_time_series)
    for i in frames
        print("Plotting frame $i of $(frames[end]) \r")
        field = field_time_series[i]
        prettytime = prettytimes[i]
        ax.title = title_prefix * " after " * prettytime
        interpolated_field = interpolate_cubed_sphere_field_to_cell_centers(grid, field, field_location;
                                                                            ssh, consider_all_levels, levels,
                                                                            read_parent_field_data)
        heatlatlon!(ax, interpolated_field, k; colorrange, colormap)
        save(joinpath(output_directory, filename * "_$i" * format), fig)
    end
end

function geo_heatlatlon_visualization_animation(grid, field_time_series, field_location, prettytimes, title_prefix;
                                                plot_frames::Bool = false,
                                                start_index::Int = 1,
                                                k::Int = 1,
                                                use_symmetric_colorrange::Bool = true,
                                                ssh::Bool = false,
                                                consider_all_levels::Bool = true,
                                                levels::UnitRange{Int} = k:k,
                                                read_parent_field_data::Bool = false,
                                                Δ::Int = 1,
                                                colorrange::Union{Nothing, Vector} = nothing,
                                                colormap::Union{Nothing, Symbol} = nothing,
                                                colorbarlabel::String = "",
                                                framerate::Int = 10,
                                                output_directory::AbstractString = "output",
                                                filename::AbstractString = "filename",
                                                figure_format::AbstractString = ".png",
                                                movie_format::AbstractString = ".mp4")
    if plot_frames
        format = figure_format
        path = joinpath(pwd(), output_directory * "/" * filename * "_frames")
        isdir(path) || mkdir(path)
        filename = joinpath(filename * "_frames/", filename)
        geo_heatlatlon_visualization_animation_frames(grid, field_time_series, field_location, prettytimes,
                                                      title_prefix;
                                                      start_index, k, use_symmetric_colorrange, ssh,
                                                      consider_all_levels, levels, read_parent_field_data, Δ,
                                                      colorrange, colormap, colorbarlabel, output_directory, filename,
                                                      format)
    else
        format = movie_format
        geo_heatlatlon_visualization_animation_Makie(grid, field_time_series, field_location, prettytimes,
                                                     title_prefix;
                                                     start_index, k, use_symmetric_colorrange, ssh, consider_all_levels,
                                                     levels, read_parent_field_data, Δ, colorrange, colormap,
                                                     colorbarlabel, framerate, output_directory, filename, format)
    end
end
