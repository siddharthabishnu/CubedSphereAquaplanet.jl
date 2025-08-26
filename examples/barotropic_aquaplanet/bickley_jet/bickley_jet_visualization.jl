include("../../cubed_sphere_visualization.jl")

function BickleyJetVisualization!(bickley_jet_grid, Nplots, Δt, plot_iteration_interval,
                                  prettytimes, framerate;
                                  output_directory::String = "bickley_jet_output",
                                  output_filename::String = "bickley_jet_output.jld2",
                                  make_panelwise_visualization_plots_with_halos::Bool = false,
                                  make_panelwise_visualization_plots::Bool = false,
                                  geo_heatmap_type::String = "heatlatlon",
                                  make_geo_heatmap_visualization_plots::Bool = false,
                                  plot_frames::Bool = false,
                                  make_panelwise_visualization_animation_with_halos::Bool = false,
                                  make_panelwise_visualization_animation::Bool = false,
                                  make_geo_heatmap_visualization_animation::Bool = false)
    u_time_series = FieldTimeSeries(joinpath(output_directory, output_filename), "u")
    v_time_series = FieldTimeSeries(joinpath(output_directory, output_filename), "v")
    η_time_series = FieldTimeSeries(joinpath(output_directory, output_filename), "η")
    c_time_series = FieldTimeSeries(joinpath(output_directory, output_filename), "c")

    ζ_time_series = compute_vorticity_time_series(bickley_jet_grid, u_time_series, v_time_series)

    colorrange_η = specify_colorrange_time_series(bickley_jet_grid, η_time_series;
                                                  ssh = true)
    colorrange_c = specify_colorrange_time_series(bickley_jet_grid, c_time_series)
    colorrange_ζ = specify_colorrange_time_series(bickley_jet_grid, ζ_time_series)
    colormap = :balance
    
    for iPlot in 0:Nplots
        plot_iteration = iPlot * plot_iteration_interval + 1
        title_η = "Surface elevation after $(prettytimes[plot_iteration])"
        title_c = "Tracer concentration after $(prettytimes[plot_iteration])"
        title_ζ = "Relative vorticity after $(prettytimes[plot_iteration])"

        if make_panelwise_visualization_plots_with_halos
            fig = panelwise_visualization(bickley_jet_grid, η_time_series[plot_iteration];
                                          with_halos = true, ssh = true, colorrange = colorrange_η, colormap)
            filename = joinpath(
                output_directory,
                @sprintf("cubed_sphere_bickley_jet_panelwise_visualization_with_halos_η_%d.png", plot_iteration))
            save(filename, fig)
            filename = joinpath(
                output_directory,
                @sprintf("cubed_sphere_bickley_jet_panelwise_visualization_with_halos_c_%d.png", plot_iteration))
            fig = panelwise_visualization(bickley_jet_grid, c_time_series[plot_iteration];
                                          with_halos = true, colorrange = colorrange_c, colormap)
            save(filename, fig)
            fig = panelwise_visualization(bickley_jet_grid, ζ_time_series[plot_iteration];
                                          with_halos = true, colorrange = colorrange_ζ, colormap)
            filename = joinpath(
                output_directory,
                @sprintf("cubed_sphere_bickley_jet_panelwise_visualization_with_halos_ζ_%d.png", plot_iteration))
            save(filename, fig)
        end

        if make_panelwise_visualization_plots
            fig = panelwise_visualization(bickley_jet_grid, η_time_series[plot_iteration];
                                          ssh = true, colorrange = colorrange_η, colormap)
            filename = joinpath(output_directory,
                                @sprintf("cubed_sphere_bickley_jet_panelwise_visualization_η_%d.png", plot_iteration))
            save(filename, fig)
            fig = panelwise_visualization(bickley_jet_grid, c_time_series[plot_iteration];
                                          colorrange = colorrange_c, colormap)
            filename = joinpath(output_directory,
                                @sprintf("cubed_sphere_bickley_jet_panelwise_visualization_c_%d.png", plot_iteration))
            save(filename, fig)
            fig = panelwise_visualization(bickley_jet_grid, ζ_time_series[plot_iteration];
                                          colorrange = colorrange_ζ, colormap)
            filename = joinpath(
                output_directory,
                @sprintf("cubed_sphere_bickley_jet_panelwise_visualization_ζ_%d.png", plot_iteration))
            save(filename, fig)
        end

        if make_geo_heatmap_visualization_plots
            fig = geo_heatmap_visualization(bickley_jet_grid, η_time_series[plot_iteration], "cc", title_η;
                                            geo_heatmap_type, ssh = true, colorrange = colorrange_η,
                                            colorbarlabel = "surface elevation", colormap)
            filename = joinpath(output_directory,
                                @sprintf("cubed_sphere_bickley_jet_geo_%s_visualization_η_%d.png",
                                         geo_heatmap_type, plot_iteration))
            save(filename, fig)
            fig = geo_heatmap_visualization(bickley_jet_grid, c_time_series[plot_iteration], "cc", title_c;
                                            geo_heatmap_type, colorrange = colorrange_c,
                                            colorbarlabel = "tracer concentration", colormap)
            filename = joinpath(output_directory,
                                @sprintf("cubed_sphere_bickley_jet_geo_%s_visualization_c_%d.png",
                                         geo_heatmap_type, plot_iteration))
            save(filename, fig)
            fig = geo_heatmap_visualization(bickley_jet_grid, ζ_time_series[plot_iteration], "cc", title_ζ;
                                            geo_heatmap_type, colorrange = colorrange_ζ,
                                            colorbarlabel = "relative vorticity", colormap)
            filename = joinpath(output_directory,
                                @sprintf("cubed_sphere_bickley_jet_geo_%s_visualization_ζ_%d.png",
                                         geo_heatmap_type, plot_iteration))
            save(filename, fig)
        end
    end

    if make_panelwise_visualization_animation_with_halos
        filename = "cubed_sphere_bickley_jet_panelwise_visualization_with_halos_η"
        panelwise_visualization_animation(bickley_jet_grid, η_time_series;
                                          plot_frames, with_halos = true, ssh = true, colorrange = colorrange_η,
                                          colormap, framerate, output_directory, filename)
        filename = "cubed_sphere_bickley_jet_panelwise_visualization_with_halos_c"
        panelwise_visualization_animation(bickley_jet_grid, c_time_series;
                                          plot_frames, with_halos = true, colorrange = colorrange_c, colormap,
                                          framerate, output_directory, filename)
        filename = "cubed_sphere_bickley_jet_panelwise_visualization_with_halos_ζ"
        panelwise_visualization_animation(bickley_jet_grid, ζ_time_series;
                                          plot_frames, with_halos = true, colorrange = colorrange_ζ, colormap,
                                          framerate, output_directory, filename)
    end

    if make_panelwise_visualization_animation
        filename = "cubed_sphere_bickley_jet_panelwise_visualization_η"
        panelwise_visualization_animation(bickley_jet_grid, η_time_series;
                                          plot_frames, ssh = true, colorrange = colorrange_η, colormap, framerate,
                                          output_directory, filename)
        filename = "cubed_sphere_bickley_jet_panelwise_visualization_c"
        panelwise_visualization_animation(bickley_jet_grid, c_time_series;
                                          plot_frames, colorrange = colorrange_c, colormap, framerate, output_directory,
                                          filename)
        filename = "cubed_sphere_bickley_jet_panelwise_visualization_ζ"
        panelwise_visualization_animation(bickley_jet_grid, ζ_time_series;
                                          plot_frames, colorrange = colorrange_ζ, colormap, framerate, output_directory,
                                          filename)
    end

    if make_geo_heatmap_visualization_animation
        filename = @sprintf("cubed_sphere_bickley_jet_geo_%s_visualization_η", geo_heatmap_type)
        geo_heatmap_visualization_animation(bickley_jet_grid, η_time_series, "cc", prettytimes, "Surface elevation";
                                            plot_frames, geo_heatmap_type, ssh = true, colorrange = colorrange_η,
                                            colormap, framerate, output_directory, filename)
        filename = @sprintf("cubed_sphere_bickley_jet_geo_%s_visualization_c", geo_heatmap_type)
        geo_heatmap_visualization_animation(bickley_jet_grid, c_time_series, "cc", prettytimes, "Tracer concentration";
                                            plot_frames, geo_heatmap_type, colorrange = colorrange_c, colormap,
                                            framerate, output_directory, filename)
        filename = @sprintf("cubed_sphere_bickley_jet_geo_%s_visualization_ζ", geo_heatmap_type)
        geo_heatmap_visualization_animation(bickley_jet_grid, ζ_time_series, "cc", prettytimes, "Relative vorticity";
                                            plot_frames, geo_heatmap_type, colorrange = colorrange_ζ, colormap,
                                            framerate, output_directory, filename)
    end
end