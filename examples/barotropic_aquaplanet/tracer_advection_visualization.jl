include("../cubed_sphere_visualization.jl")

function TracerAdvectionVisualization!(tracer_advection_grid, θ₀, Nplots, Δt, plot_iteration_interval, prettytimes,
                                       framerate;
                                       output_directory = "tracer_advection",
                                       output_filename = "tracer_advection_output.jld2",
                                       make_panelwise_visualization_plots_with_halos = false,
                                       make_panelwise_visualization_plots = false,
                                       make_geo_heatlatlon_visualization_plots = false,
                                       make_panelwise_visualization_animation_with_halos = false,
                                       make_panelwise_visualization_animation = false,
                                       make_geo_heatlatlon_visualization_animation = false)
    θ_time_series = FieldTimeSeries(joinpath(output_directory, output_filename), "θ")
    
    colorrange = float([-θ₀, θ₀])
    colormap = :balance
    
    for iPlot in 0:Nplots
        plot_iteration = iPlot * plot_iteration_interval + 1
        simulation_time = (plot_iteration - 1) * Δt
        title = "Tracer concentration after $(prettytime(simulation_time))"

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

        if make_geo_heatlatlon_visualization_plots
            fig = geo_heatlatlon_visualization(tracer_advection_grid, θ_time_series[plot_iteration], "cc", title;
                                               colorbarlabel = "tracer concentration", colorrange, colormap)
            filename = joinpath(output_directory,
                                @sprintf("cubed_sphere_tracer_advection_geo_heatlatlon_visualization_θ_%d.png",
                                         plot_iteration))
            save(filename, fig)
        end
    end

    if make_panelwise_visualization_animation_with_halos
        filename = "cubed_sphere_tracer_advection_panelwise_visualization_with_halos_θ"
        panelwise_visualization_animation(tracer_advection_grid, θ_time_series;
                                          with_halos=true, colorrange, colormap, framerate, output_directory, filename)
    end

    if make_panelwise_visualization_animation
        filename = "cubed_sphere_tracer_advection_panelwise_visualization_θ"
        panelwise_visualization_animation(tracer_advection_grid, θ_time_series;
                                          colorrange, colormap, framerate, output_directory, filename)
    end

    if make_geo_heatlatlon_visualization_animation
        filename = "cubed_sphere_tracer_advection_geo_heatlatlon_visualization_θ"
        geo_heatlatlon_visualization_animation(tracer_advection_grid, θ_time_series, "cc", prettytimes,
                                               "Tracer concentration";
                                               colorrange, colormap, framerate, output_directory, filename)
    end
end