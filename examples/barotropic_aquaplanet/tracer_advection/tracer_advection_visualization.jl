include("../../cubed_sphere_visualization.jl")

function TracerAdvectionVisualization!(tracer_advection_grid, θ₀, Nplots, Δt, plot_iteration_interval, prettytimes,
                                       framerate;
                                       output_directory::String = "tracer_advection_output",
                                       output_filename::String = "tracer_advection_output.jld2",
                                       make_panelwise_visualization_plots_with_halos::Bool = false,
                                       make_panelwise_visualization_plots::Bool = false,
                                       geo_heatmap_type::String = "heatlatlon",
                                       make_geo_heatmap_visualization_plots::Bool = false,
                                       plot_frames::Bool = false,
                                       make_panelwise_visualization_animation_with_halos::Bool = false,
                                       make_panelwise_visualization_animation::Bool = false,
                                       make_geo_heatmap_visualization_animation::Bool = false)
    θ_time_series = FieldTimeSeries(joinpath(output_directory, output_filename), "θ")
    
    colorrange = float([-θ₀, θ₀])
    colormap = :balance
    
    for iPlot in 0:Nplots
        plot_iteration = iPlot * plot_iteration_interval + 1
        title = "Tracer concentration after $(prettytimes[plot_iteration])"

        if make_panelwise_visualization_plots_with_halos
            fig = panelwise_visualization(tracer_advection_grid, θ_time_series[plot_iteration];
                                          with_halos = true, colorrange, colormap)
            filename = joinpath(output_directory,
                                @sprintf("cubed_sphere_tracer_advection_panelwise_visualization_with_halos_θ_%d.png",
                                         plot_iteration))
            save(filename, fig)
        end

        if make_panelwise_visualization_plots
            fig = panelwise_visualization(tracer_advection_grid, θ_time_series[plot_iteration];
                                          colorrange, colormap)
            filename = joinpath(output_directory,
                                @sprintf("cubed_sphere_tracer_advection_panelwise_visualization_θ_%d.png",
                                         plot_iteration))
            save(filename, fig)
        end

        if make_geo_heatmap_visualization_plots
            fig = geo_heatmap_visualization(tracer_advection_grid, θ_time_series[plot_iteration], "cc", title;
                                            geo_heatmap_type, colorbarlabel = "tracer concentration", colorrange,
                                            colormap)
            filename = joinpath(output_directory,
                                @sprintf("cubed_sphere_tracer_advection_geo_%s_visualization_θ_%d.png",
                                         geo_heatmap_type, plot_iteration))
            save(filename, fig)
        end
    end

    if make_panelwise_visualization_animation_with_halos
        filename = "cubed_sphere_tracer_advection_panelwise_visualization_with_halos_θ"
        panelwise_visualization_animation(tracer_advection_grid, θ_time_series;
                                          plot_frames, with_halos=true, colorrange, colormap, framerate,
                                          output_directory, filename)
    end

    if make_panelwise_visualization_animation
        filename = "cubed_sphere_tracer_advection_panelwise_visualization_θ"
        panelwise_visualization_animation(tracer_advection_grid, θ_time_series;
                                          plot_frames, colorrange, colormap, framerate, output_directory, filename)
    end

    if make_geo_heatmap_visualization_animation
        filename = @sprintf("cubed_sphere_tracer_advection_geo_%s_visualization_θ", geo_heatmap_type)
        geo_heatmap_visualization_animation(tracer_advection_grid, θ_time_series, "cc", prettytimes,
                                            "Tracer concentration";
                                            plot_frames, geo_heatmap_type, colorrange, colormap, framerate,
                                            output_directory, filename)
    end
end