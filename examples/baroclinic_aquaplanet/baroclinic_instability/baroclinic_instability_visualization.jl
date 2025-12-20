include("../../cubed_sphere_visualization.jl")

function BaroclinicInstabilityVisualization!(
    baroclinic_instability_grid, Nplots, Δt, plot_iteration_interval,
    prettytimes, framerate;
    output_directory::String = "baroclinic_instability_output",
    prognostic_output_filename::String = "baroclinic_instability_surface_prognostic_fields_output.jld2",
    diagnostic_output_filename::String = "baroclinic_instability_surface_diagnostic_fields_output.jld2",
    iPlot_Start::Int = 0,
    iPlot_Δ::Int = 1,
    plot_states::Dict{Symbol,Bool} = Dict(:u => true,
                                          :v => true,
                                          :w => true,
                                          :η => true,
                                          :T => true,
                                          :S => true,
                                          :ζ => true),
    make_panelwise_visualization_plots_with_halos::Bool = false,
    make_panelwise_visualization_plots::Bool = false,
    geo_heatmap_type::String = "heatlatlon",
    make_geo_heatmap_visualization_plots::Bool = false,
    plot_frames::Bool = false,
    make_panelwise_visualization_animation_with_halos::Bool = false,
    make_panelwise_visualization_animation::Bool = false,
    make_geo_heatmap_visualization_animation::Bool = false)
    (plot_states[:u] || plot_states[:ζ]) &&
        (u_time_series = FieldTimeSeries(joinpath(output_directory, prognostic_output_filename), "u"))
    (plot_states[:v] || plot_states[:ζ]) &&
        (v_time_series = FieldTimeSeries(joinpath(output_directory, prognostic_output_filename), "v"))
    plot_states[:w] && (w_time_series = FieldTimeSeries(joinpath(output_directory, diagnostic_output_filename), "w"))
    plot_states[:η] && (η_time_series = FieldTimeSeries(joinpath(output_directory, diagnostic_output_filename), "η"))
    plot_states[:T] && (T_time_series = FieldTimeSeries(joinpath(output_directory, prognostic_output_filename), "T"))
    plot_states[:S] && (S_time_series = FieldTimeSeries(joinpath(output_directory, prognostic_output_filename), "S"))

    @info "Computing ζ_time_series..."
    plot_states[:ζ] &&
        (ζ_time_series = compute_vorticity_time_series(baroclinic_instability_grid, u_time_series, v_time_series))
    @info "Interpolating to cell centers for u time series..."
    plot_states[:u] &&
        (u_time_series = interpolate_cubed_sphere_field_time_series_to_cell_centers(baroclinic_instability_grid,
                                                                                    u_time_series, "fc"))
    @info "Interpolating to cell centers for v time series..."
    plot_states[:v] &&
        (v_time_series = interpolate_cubed_sphere_field_time_series_to_cell_centers(baroclinic_instability_grid,
                                                                                    v_time_series, "cf"))

    plot_states[:u] && (colorrange_u = specify_colorrange_time_series(baroclinic_instability_grid, u_time_series))
    plot_states[:v] && (colorrange_v = specify_colorrange_time_series(baroclinic_instability_grid, v_time_series))
    plot_states[:w] && (colorrange_w = specify_colorrange_time_series(baroclinic_instability_grid, w_time_series))
    plot_states[:η] && (colorrange_η = specify_colorrange_time_series(baroclinic_instability_grid, η_time_series;
                                                                      ssh = true))
    plot_states[:T] && (colorrange_T = specify_colorrange_time_series(baroclinic_instability_grid, T_time_series))
    plot_states[:S] && (colorrange_S = specify_colorrange_time_series(baroclinic_instability_grid, S_time_series))
    plot_states[:ζ] && (colorrange_ζ = specify_colorrange_time_series(baroclinic_instability_grid, ζ_time_series))
    
    colormap = :balance
    
    u_time_series_plots = []
    v_time_series_plots = []
    w_time_series_plots = []
    η_time_series_plots = []
    T_time_series_plots = []
    S_time_series_plots = []
    ζ_time_series_plots = []

    for iPlot in iPlot_Start:iPlot_Δ:Nplots
        plot_iteration = iPlot * plot_iteration_interval + 1
        plot_states[:u] && push!(u_time_series_plots, deepcopy(u_time_series[plot_iteration]))
        plot_states[:v] && push!(v_time_series_plots, deepcopy(v_time_series[plot_iteration]))
        plot_states[:w] && push!(w_time_series_plots, deepcopy(w_time_series[plot_iteration]))
        plot_states[:η] && push!(η_time_series_plots, deepcopy(η_time_series[plot_iteration]))
        plot_states[:T] && push!(T_time_series_plots, deepcopy(T_time_series[plot_iteration]))
        plot_states[:S] && push!(S_time_series_plots, deepcopy(S_time_series[plot_iteration]))
        plot_states[:ζ] && push!(ζ_time_series_plots, deepcopy(ζ_time_series[plot_iteration]))
    end
    
    plot_states[:u] && (colorrange_plots_u = specify_colorrange_time_series(baroclinic_instability_grid,
                                                                            u_time_series_plots))
    plot_states[:v] && (colorrange_plots_v = specify_colorrange_time_series(baroclinic_instability_grid,
                                                                            v_time_series_plots))
    plot_states[:w] && (colorrange_plots_w = specify_colorrange_time_series(baroclinic_instability_grid,
                                                                            w_time_series_plots))
    plot_states[:η] && (colorrange_plots_η = specify_colorrange_time_series(baroclinic_instability_grid,
                                                                            η_time_series_plots;
                                                                            ssh = true))
    plot_states[:T] && (colorrange_plots_T = specify_colorrange_time_series(baroclinic_instability_grid,
                                                                            T_time_series_plots))
    plot_states[:S] && (colorrange_plots_S = specify_colorrange_time_series(baroclinic_instability_grid,
                                                                            S_time_series_plots))
    plot_states[:ζ] && (colorrange_plots_ζ = specify_colorrange_time_series(baroclinic_instability_grid,
                                                                            ζ_time_series_plots))

    Nx, Ny, Nz = size(baroclinic_instability_grid)

    for iPlot in iPlot_Start:iPlot_Δ:Nplots
        plot_iteration = iPlot * plot_iteration_interval + 1
        title_u = "Zonal velocity after $(prettytimes[plot_iteration])"
        title_v = "Meridional velocity after $(prettytimes[plot_iteration])"
        title_w = "Vertical velocity after $(prettytimes[plot_iteration])"
        title_η = "Surface elevation after $(prettytimes[plot_iteration])"
        title_T = "Temperature profile after $(prettytimes[plot_iteration])"
        title_S = "Salinity profile after $(prettytimes[plot_iteration])"
        title_ζ = "Relative vorticity after $(prettytimes[plot_iteration])"

        if make_panelwise_visualization_plots_with_halos && plot_states[:u]
            filename = joinpath(
                output_directory,
                @sprintf("cubed_sphere_baroclinic_instability_panelwise_visualization_with_halos_u_%d.png",
                         plot_iteration))
            fig = panelwise_visualization(baroclinic_instability_grid, u_time_series[plot_iteration];
                                          with_halos = true, k = Nz, colorrange = colorrange_plots_u, colormap)
            save(filename, fig)
        end
        if make_panelwise_visualization_plots_with_halos && plot_states[:v]
            filename = joinpath(
                output_directory,
                @sprintf("cubed_sphere_baroclinic_instability_panelwise_visualization_with_halos_v_%d.png",
                         plot_iteration))
            fig = panelwise_visualization(baroclinic_instability_grid, v_time_series[plot_iteration];
                                          with_halos = true, k = Nz, colorrange = colorrange_plots_v, colormap)
            save(filename, fig)
        end
        if make_panelwise_visualization_plots_with_halos && plot_states[:w]
            filename = joinpath(
                output_directory,
                @sprintf("cubed_sphere_baroclinic_instability_panelwise_visualization_with_halos_w_%d.png",
                         plot_iteration))
            fig = panelwise_visualization(baroclinic_instability_grid, w_time_series[plot_iteration];
                                          with_halos = true, k = Nz + 1, colorrange = colorrange_plots_w, colormap)
            save(filename, fig)
        end
        if make_panelwise_visualization_plots_with_halos && plot_states[:η]
            filename = joinpath(
                output_directory,
                @sprintf("cubed_sphere_baroclinic_instability_panelwise_visualization_with_halos_η_%d.png",
                         plot_iteration))
            fig = panelwise_visualization(baroclinic_instability_grid, η_time_series[plot_iteration];
                                          with_halos = true, k = Nz + 1, ssh = true, colorrange = colorrange_plots_η,
                                          colormap)
            save(filename, fig)
        end
        if make_panelwise_visualization_plots_with_halos && plot_states[:T]
            filename = joinpath(
                output_directory,
                @sprintf("cubed_sphere_baroclinic_instability_panelwise_visualization_with_halos_T_%d.png",
                         plot_iteration))
            fig = panelwise_visualization(baroclinic_instability_grid, T_time_series[plot_iteration];
                                          with_halos = true, k = Nz, colorrange = colorrange_plots_T, colormap)
            save(filename, fig)
        end
        if make_panelwise_visualization_plots_with_halos && plot_states[:S]
            filename = joinpath(
                output_directory,
                @sprintf("cubed_sphere_baroclinic_instability_panelwise_visualization_with_halos_S_%d.png",
                         plot_iteration))
            fig = panelwise_visualization(baroclinic_instability_grid, S_time_series[plot_iteration];
                                          with_halos = true, k = Nz, colorrange = colorrange_plots_S, colormap)
            save(filename, fig)
        end
        if make_panelwise_visualization_plots_with_halos && plot_states[:ζ]
            filename = joinpath(
                output_directory,
                @sprintf("cubed_sphere_baroclinic_instability_panelwise_visualization_with_halos_ζ_%d.png",
                         plot_iteration))
            fig = panelwise_visualization(baroclinic_instability_grid, ζ_time_series[plot_iteration];
                                          with_halos = true, k = Nz, colorrange = colorrange_plots_ζ, colormap)
            save(filename, fig)
        end
        
        if make_panelwise_visualization_plots && plot_states[:u]
            filename = joinpath(output_directory,
                                @sprintf("cubed_sphere_baroclinic_instability_panelwise_visualization_u_%d.png",
                                         plot_iteration))
            fig = panelwise_visualization(baroclinic_instability_grid, u_time_series[plot_iteration];
                                          k = Nz, colorrange = colorrange_plots_u, colormap)
            save(filename, fig)
        end
        if make_panelwise_visualization_plots && plot_states[:v]
            filename = joinpath(output_directory,
                                @sprintf("cubed_sphere_baroclinic_instability_panelwise_visualization_v_%d.png",
                                         plot_iteration))
            fig = panelwise_visualization(baroclinic_instability_grid, v_time_series[plot_iteration];
                                          k = Nz, colorrange = colorrange_plots_v, colormap)
            save(filename, fig)
        end
        if make_panelwise_visualization_plots && plot_states[:w]
            filename = joinpath(output_directory,
                                @sprintf("cubed_sphere_baroclinic_instability_panelwise_visualization_w_%d.png",
                                         plot_iteration))
            fig = panelwise_visualization(baroclinic_instability_grid, w_time_series[plot_iteration];
                                          k = Nz + 1, colorrange = colorrange_plots_w, colormap)
            save(filename, fig)
        end
        if make_panelwise_visualization_plots && plot_states[:η]
            filename = joinpath(output_directory,
                                @sprintf("cubed_sphere_baroclinic_instability_panelwise_visualization_η_%d.png",
                                         plot_iteration))
            fig = panelwise_visualization(baroclinic_instability_grid, η_time_series[plot_iteration];
                                          k = Nz + 1, ssh = true, colorrange = colorrange_plots_η, colormap)
            save(filename, fig)
        end
        if make_panelwise_visualization_plots && plot_states[:T]
            filename = joinpath(output_directory,
                                @sprintf("cubed_sphere_baroclinic_instability_panelwise_visualization_T_%d.png",
                                         plot_iteration))
            fig = panelwise_visualization(baroclinic_instability_grid, T_time_series[plot_iteration];
                                          k = Nz, colorrange = colorrange_plots_T, colormap)
            save(filename, fig)
        end
        if make_panelwise_visualization_plots && plot_states[:S]
            filename = joinpath(output_directory,
                                @sprintf("cubed_sphere_baroclinic_instability_panelwise_visualization_S_%d.png",
                                         plot_iteration))
            fig = panelwise_visualization(baroclinic_instability_grid, S_time_series[plot_iteration];
                                          k = Nz, colorrange = colorrange_plots_S, colormap)
            save(filename, fig)
        end
        if make_panelwise_visualization_plots && plot_states[:ζ]
            filename = joinpath(
                output_directory,
                @sprintf("cubed_sphere_baroclinic_instability_panelwise_visualization_ζ_%d.png", plot_iteration))
            fig = panelwise_visualization(baroclinic_instability_grid, ζ_time_series[plot_iteration];
                                          k = Nz, colorrange = colorrange_plots_ζ, colormap)
            save(filename, fig)
        end
        
        if make_geo_heatmap_visualization_plots && plot_states[:u]
            filename = joinpath(output_directory,
                                @sprintf("cubed_sphere_baroclinic_instability_geo_%s_visualization_u_%d.png",
                                         geo_heatmap_type, plot_iteration))
            fig = geo_heatmap_visualization(baroclinic_instability_grid, u_time_series[plot_iteration], "cc", title_u;
                                            geo_heatmap_type, k = Nz, colorrange = colorrange_plots_u,
                                            colorbarlabel = "zonal velocity", colormap)
            save(filename, fig)
        end
        if make_geo_heatmap_visualization_plots && plot_states[:v]
            filename = joinpath(output_directory,
                                @sprintf("cubed_sphere_baroclinic_instability_geo_%s_visualization_v_%d.png",
                                         geo_heatmap_type, plot_iteration))
            fig = geo_heatmap_visualization(baroclinic_instability_grid, v_time_series[plot_iteration], "cc", title_v;
                                            geo_heatmap_type, k = Nz, colorrange = colorrange_plots_v,
                                            colorbarlabel = "meridional velocity", colormap)
            save(filename, fig)
        end
        if make_geo_heatmap_visualization_plots && plot_states[:w]
            filename = joinpath(output_directory,
                                @sprintf("cubed_sphere_baroclinic_instability_geo_%s_visualization_w_%d.png",
                                         geo_heatmap_type, plot_iteration))
            fig = geo_heatmap_visualization(baroclinic_instability_grid, w_time_series[plot_iteration], "cc", title_w;
                                            geo_heatmap_type, k = Nz + 1, colorrange = colorrange_plots_w,
                                            colorbarlabel = "vertical velocity", colormap)
            save(filename, fig)
        end
        if make_geo_heatmap_visualization_plots && plot_states[:η]
            filename = joinpath(output_directory,
                                @sprintf("cubed_sphere_baroclinic_instability_geo_%s_visualization_η_%d.png",
                                         geo_heatmap_type, plot_iteration))
            fig = geo_heatmap_visualization(baroclinic_instability_grid, η_time_series[plot_iteration], "cc", title_η;
                                            geo_heatmap_type, ssh = true, colorrange = colorrange_plots_η,
                                            colorbarlabel = "surface elevation", colormap)
            save(filename, fig)
        end
        if make_geo_heatmap_visualization_plots && plot_states[:T]
            filename = joinpath(output_directory,
                                @sprintf("cubed_sphere_baroclinic_instability_geo_%s_visualization_T_%d.png",
                                         geo_heatmap_type, plot_iteration))
            fig = geo_heatmap_visualization(baroclinic_instability_grid, T_time_series[plot_iteration], "cc", title_T;
                                            geo_heatmap_type, k = Nz, colorrange = colorrange_plots_T,
                                            colorbarlabel = "temperature", colormap)
            save(filename, fig)
        end
        if make_geo_heatmap_visualization_plots && plot_states[:S]
            filename = joinpath(output_directory,
                                @sprintf("cubed_sphere_baroclinic_instability_geo_%s_visualization_S_%d.png",
                                         geo_heatmap_type, plot_iteration))
            fig = geo_heatmap_visualization(baroclinic_instability_grid, S_time_series[plot_iteration], "cc", title_S;
                                            geo_heatmap_type, k = Nz, colorrange = colorrange_plots_S,
                                            colorbarlabel = "salinity", colormap)
            save(filename, fig)
        end
        if make_geo_heatmap_visualization_plots && plot_states[:ζ]
            filename = joinpath(output_directory,
                                @sprintf("cubed_sphere_baroclinic_instability_geo_%s_visualization_ζ_%d.png",
                                         geo_heatmap_type, plot_iteration))
            fig = geo_heatmap_visualization(baroclinic_instability_grid, ζ_time_series[plot_iteration], "cc", title_ζ;
                                            geo_heatmap_type, k = Nz, colorrange = colorrange_plots_ζ,
                                            colorbarlabel = "relative vorticity", colormap)
            save(filename, fig)
        end
    end
    
    if make_panelwise_visualization_animation_with_halos && plot_states[:u]
        filename = "cubed_sphere_baroclinic_instability_panelwise_visualization_with_halos_u"
        panelwise_visualization_animation(baroclinic_instability_grid, u_time_series;
                                          plot_frames, with_halos = true, k = Nz, colorrange = colorrange_u, colormap,
                                          framerate, output_directory, filename)
    end
    if make_panelwise_visualization_animation_with_halos && plot_states[:v]
        filename = "cubed_sphere_baroclinic_instability_panelwise_visualization_with_halos_v"
        panelwise_visualization_animation(baroclinic_instability_grid, v_time_series;
                                          plot_frames, with_halos = true, k = Nz, colorrange = colorrange_v, colormap,
                                          framerate, output_directory, filename)
    end
    if make_panelwise_visualization_animation_with_halos && plot_states[:w]
        filename = "cubed_sphere_baroclinic_instability_panelwise_visualization_with_halos_w"
        panelwise_visualization_animation(baroclinic_instability_grid, w_time_series;
                                          plot_frames, with_halos = true, k = Nz + 1, colorrange = colorrange_w,
                                          colormap, framerate, output_directory, filename)
    end
    if make_panelwise_visualization_animation_with_halos && plot_states[:η]
        filename = "cubed_sphere_baroclinic_instability_panelwise_visualization_with_halos_η"
        panelwise_visualization_animation(baroclinic_instability_grid, η_time_series;
                                          plot_frames, with_halos = true, ssh = true, colorrange = colorrange_η,
                                          colormap, framerate, output_directory, filename)
    end
    if make_panelwise_visualization_animation_with_halos && plot_states[:T]
        filename = "cubed_sphere_baroclinic_instability_panelwise_visualization_with_halos_T"
        panelwise_visualization_animation(baroclinic_instability_grid, T_time_series;
                                          plot_frames, with_halos = true, k = Nz, colorrange = colorrange_T, colormap,
                                          framerate, output_directory, filename)
    end
    if make_panelwise_visualization_animation_with_halos && plot_states[:S]
        filename = "cubed_sphere_baroclinic_instability_panelwise_visualization_with_halos_S"
        panelwise_visualization_animation(baroclinic_instability_grid, S_time_series;
                                          plot_frames, with_halos = true, k = Nz, colorrange = colorrange_S, colormap,
                                          framerate, output_directory, filename)
    end
    if make_panelwise_visualization_animation_with_halos && plot_states[:ζ]
        filename = "cubed_sphere_baroclinic_instability_panelwise_visualization_with_halos_ζ"
        panelwise_visualization_animation(baroclinic_instability_grid, ζ_time_series;
                                          plot_frames, with_halos = true, k = Nz, colorrange = colorrange_ζ, colormap,
                                          framerate, output_directory, filename)
    end
    
    if make_panelwise_visualization_animation && plot_states[:u]
        filename = "cubed_sphere_baroclinic_instability_panelwise_visualization_u"
        panelwise_visualization_animation(baroclinic_instability_grid, u_time_series;
                                          plot_frames, k = Nz, colorrange = colorrange_u, colormap, framerate,
                                          output_directory, filename)
    end
    if make_panelwise_visualization_animation && plot_states[:v]
        filename = "cubed_sphere_baroclinic_instability_panelwise_visualization_v"
        panelwise_visualization_animation(baroclinic_instability_grid, v_time_series;
                                          plot_frames, k = Nz, colorrange = colorrange_v, colormap, framerate,
                                          output_directory, filename)
    end
    if make_panelwise_visualization_animation && plot_states[:w]
        filename = "cubed_sphere_baroclinic_instability_panelwise_visualization_w"
        panelwise_visualization_animation(baroclinic_instability_grid, w_time_series;
                                          plot_frames, k = Nz + 1, colorrange = colorrange_w, colormap, framerate,
                                          output_directory, filename)
    end
    if make_panelwise_visualization_animation && plot_states[:η]
        filename = "cubed_sphere_baroclinic_instability_panelwise_visualization_η"
        panelwise_visualization_animation(baroclinic_instability_grid, η_time_series;
                                          plot_frames, ssh = true, colorrange = colorrange_η, colormap, framerate,
                                          output_directory, filename)
    end
    if make_panelwise_visualization_animation && plot_states[:T]
        filename = "cubed_sphere_baroclinic_instability_panelwise_visualization_T"
        panelwise_visualization_animation(baroclinic_instability_grid, T_time_series;
                                          plot_frames, k = Nz, colorrange = colorrange_T, colormap, framerate,
                                          output_directory, filename)
    end
    if make_panelwise_visualization_animation && plot_states[:S]
        filename = "cubed_sphere_baroclinic_instability_panelwise_visualization_S"
        panelwise_visualization_animation(baroclinic_instability_grid, S_time_series;
                                          plot_frames, k = Nz, colorrange = colorrange_S, colormap, framerate,
                                          output_directory, filename)
    end
    if make_panelwise_visualization_animation && plot_states[:ζ]
        filename = "cubed_sphere_baroclinic_instability_panelwise_visualization_ζ"
        panelwise_visualization_animation(baroclinic_instability_grid, ζ_time_series;
                                          plot_frames, k = Nz, colorrange = colorrange_ζ, colormap, framerate,
                                          output_directory, filename)
    end

    if make_geo_heatmap_visualization_animation && plot_states[:u]
        filename = @sprintf("cubed_sphere_baroclinic_instability_geo_%s_visualization_u", geo_heatmap_type)
        geo_heatmap_visualization_animation(baroclinic_instability_grid, u_time_series, "cc", prettytimes,
                                            "Zonal velocity";
                                            plot_frames, geo_heatmap_type, k = Nz, colorrange = colorrange_u,
                                            colormap, framerate, output_directory, filename)
    end
    if make_geo_heatmap_visualization_animation && plot_states[:v]
        filename = @sprintf("cubed_sphere_baroclinic_instability_geo_%s_visualization_v", geo_heatmap_type)
        geo_heatmap_visualization_animation(baroclinic_instability_grid, v_time_series, "cc", prettytimes,
                                            "Meridional velocity";
                                            plot_frames, geo_heatmap_type, k = Nz, colorrange = colorrange_v,
                                            colormap, framerate, output_directory, filename)
    end
    if make_geo_heatmap_visualization_animation && plot_states[:w]
        filename = @sprintf("cubed_sphere_baroclinic_instability_geo_%s_visualization_w", geo_heatmap_type)
        geo_heatmap_visualization_animation(baroclinic_instability_grid, w_time_series, "cc", prettytimes,
                                            "Vertical velocity";
                                            plot_frames, geo_heatmap_type, k = Nz + 1, colorrange = colorrange_w,
                                            colormap, framerate, output_directory, filename)
    end
    if make_geo_heatmap_visualization_animation && plot_states[:η]
        filename = @sprintf("cubed_sphere_baroclinic_instability_geo_%s_visualization_η", geo_heatmap_type)
        geo_heatmap_visualization_animation(baroclinic_instability_grid, η_time_series, "cc", prettytimes,
                                            "Surface elevation";
                                            plot_frames, geo_heatmap_type, ssh = true, colorrange = colorrange_η,
                                            colormap, framerate, output_directory, filename)
    end
    if make_geo_heatmap_visualization_animation && plot_states[:T]
        filename = @sprintf("cubed_sphere_baroclinic_instability_geo_%s_visualization_T", geo_heatmap_type)
        geo_heatmap_visualization_animation(baroclinic_instability_grid, T_time_series, "cc", prettytimes,
                                            "Temperature profile";
                                            plot_frames, geo_heatmap_type, k = Nz, colorrange = colorrange_T,
                                            colormap, framerate, output_directory, filename)
    end
    if make_geo_heatmap_visualization_animation && plot_states[:S]
        filename = @sprintf("cubed_sphere_baroclinic_instability_geo_%s_visualization_S", geo_heatmap_type)
        geo_heatmap_visualization_animation(baroclinic_instability_grid, S_time_series, "cc", prettytimes,
                                            "Salinity profile";
                                            plot_frames, geo_heatmap_type, k = Nz, colorrange = colorrange_S,
                                            colormap, framerate, output_directory, filename)
    end
    if make_geo_heatmap_visualization_animation && plot_states[:ζ]
        filename = @sprintf("cubed_sphere_baroclinic_instability_geo_%s_visualization_ζ", geo_heatmap_type)
        geo_heatmap_visualization_animation(baroclinic_instability_grid, ζ_time_series, "cc", prettytimes,
                                            "Relative vorticity";
                                            plot_frames, geo_heatmap_type, k = Nz, colorrange = colorrange_ζ,
                                            colormap, framerate, output_directory, filename)
    end
end