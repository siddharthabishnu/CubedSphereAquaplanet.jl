include("../../cubed_sphere_visualization.jl")

function RossbyHaurwitzWaveVisualization!(rossby_haurwitz_wave_grid, Nplots, Δt, plot_iteration_interval,
                                          prettytimes, framerate;
                                          output_directory = "rossby_haurwitz_wave_output",
                                          output_filename = "rossby_haurwitz_wave_output.jld2",
                                          make_panelwise_visualization_plots_with_halos = false,
                                          make_panelwise_visualization_plots = false,
                                          make_geo_heatlatlon_visualization_plots = false,
                                          make_panelwise_visualization_animation_with_halos = false,
                                          make_panelwise_visualization_animation = false,
                                          make_geo_heatlatlon_visualization_animation = false)
    u_time_series = FieldTimeSeries(joinpath(output_directory, output_filename), "u")
    v_time_series = FieldTimeSeries(joinpath(output_directory, output_filename), "v")
    η_time_series = FieldTimeSeries(joinpath(output_directory, output_filename), "η")

    ζ_time_series = compute_vorticity_time_series(rossby_haurwitz_wave_grid, u_time_series, v_time_series)

    colorrange_η = specify_colorrange_time_series(rossby_haurwitz_wave_grid, η_time_series;
                                                  ssh = true)
    colorrange_ζ = specify_colorrange_time_series(rossby_haurwitz_wave_grid, ζ_time_series)
    colormap = :balance
    
    for iPlot in 0:Nplots
        plot_iteration = iPlot * plot_iteration_interval + 1
        simulation_time = (plot_iteration - 1) * Δt
        title_η = "Surface elevation after $(prettytime(simulation_time))"
        title_ζ = "Relative vorticity after $(prettytime(simulation_time))"
        
        if make_panelwise_visualization_plots_with_halos
            fig = panelwise_visualization(rossby_haurwitz_wave_grid, η_time_series[plot_iteration];
                                          with_halos = true, ssh = true, colorrange = colorrange_η, colormap)
            filename = joinpath(
                output_directory,
                @sprintf("cubed_sphere_rossby_haurwitz_wave_panelwise_visualization_with_halos_η_%d.png",
                         plot_iteration))
            save(filename, fig)
            fig = panelwise_visualization(rossby_haurwitz_wave_grid, ζ_time_series[plot_iteration];
                                          with_halos = true, colorrange = colorrange_ζ, colormap)
            filename = joinpath(
                output_directory,
                @sprintf("cubed_sphere_rossby_haurwitz_wave_panelwise_visualization_with_halos_ζ_%d.png",
                         plot_iteration))
            save(filename, fig)
        end

        if make_panelwise_visualization_plots
            fig = panelwise_visualization(rossby_haurwitz_wave_grid, η_time_series[plot_iteration];
                                          ssh = true, colorrange = colorrange_η, colormap)
            filename = joinpath(output_directory,
                                @sprintf("cubed_sphere_rossby_haurwitz_wave_panelwise_visualization_η_%d.png",
                                         plot_iteration))
            save(filename, fig)
            fig = panelwise_visualization(rossby_haurwitz_wave_grid, ζ_time_series[plot_iteration];
                                          colorrange = colorrange_ζ, colormap)
            filename = joinpath(
                output_directory,
                @sprintf("cubed_sphere_rossby_haurwitz_wave_panelwise_visualization_ζ_%d.png",
                         plot_iteration))
            save(filename, fig)
        end

        if make_geo_heatlatlon_visualization_plots
            fig = geo_heatlatlon_visualization(rossby_haurwitz_wave_grid, η_time_series[plot_iteration], "cc", title_η;
                                               ssh = true, colorrange = colorrange_η,
                                               colorbarlabel = "surface elevation", colormap)
            filename = joinpath(output_directory,
                                @sprintf("cubed_sphere_rossby_haurwitz_wave_geo_heatlatlon_visualization_η_%d.png",
                                         plot_iteration))
            save(filename, fig)
            fig = geo_heatlatlon_visualization(rossby_haurwitz_wave_grid, ζ_time_series[plot_iteration], "cc", title_ζ;
                                               colorrange = colorrange_ζ, colorbarlabel = "relative vorticity",
                                               colormap)
            filename = joinpath(output_directory,
                                @sprintf("cubed_sphere_rossby_haurwitz_wave_geo_heatlatlon_visualization_ζ_%d.png",
                                         plot_iteration))
            save(filename, fig)
        end
    end

    if make_panelwise_visualization_animation_with_halos
        filename = "cubed_sphere_rossby_haurwitz_wave_panelwise_visualization_with_halos_η"
        panelwise_visualization_animation(rossby_haurwitz_wave_grid, η_time_series;
                                          with_halos = true, ssh = true, colorrange = colorrange_η, colormap, framerate,
                                          output_directory, filename)
        filename = "cubed_sphere_rossby_haurwitz_wave_panelwise_visualization_with_halos_ζ"
        panelwise_visualization_animation(rossby_haurwitz_wave_grid, ζ_time_series;
                                          with_halos = true, colorrange = colorrange_ζ, colormap, framerate,
                                          output_directory, filename)
    end

    if make_panelwise_visualization_animation
        filename = "cubed_sphere_rossby_haurwitz_wave_panelwise_visualization_η"
        panelwise_visualization_animation(rossby_haurwitz_wave_grid, η_time_series;
                                          ssh = true, colorrange = colorrange_η, colormap, framerate, output_directory,
                                          filename)
        filename = "cubed_sphere_rossby_haurwitz_wave_panelwise_visualization_ζ"
        panelwise_visualization_animation(rossby_haurwitz_wave_grid, ζ_time_series;
                                          colorrange = colorrange_ζ, colormap, framerate, output_directory, filename)
    end

    if make_geo_heatlatlon_visualization_animation
        filename = "cubed_sphere_rossby_haurwitz_wave_geo_heatlatlon_visualization_η"
        geo_heatlatlon_visualization_animation(rossby_haurwitz_wave_grid, η_time_series, "cc", prettytimes,
                                               "Surface elevation";
                                               ssh = true, colorrange = colorrange_η, colormap, framerate,
                                               output_directory, filename)
        filename = "cubed_sphere_rossby_haurwitz_wave_geo_heatlatlon_visualization_ζ"
        geo_heatlatlon_visualization_animation(rossby_haurwitz_wave_grid, ζ_time_series, "cc", prettytimes,
                                               "Relative vorticity";
                                               colorrange = colorrange_ζ, colormap, framerate, output_directory,
                                               filename)
    end
end