"""
Plot simulation statistics module
"""
module PlotSimStats

import Base.CoreLogging: @info, @warn, @error, @debug
using DataFrames
using CSV
using Plots
import Statistics
using ArgParse
using Printf

const PLOT_3D = "3d"
const PLOT_HMAP = "hmap"

const CONTACT_RADIUS_PALETTE = :tab10

# Map matplotlib colormap names to Plots.jl colormap symbols
function get_plots_colormap(cmap::String)::Symbol
    colormap_map = Dict(
        "hot_r" => :heat,
        "hot" => :heat,
        "gist_heat_r" => :heat,
        "afmhot_r" => :heat,
        "YlOrRd" => :YlOrRd,
        "Greys" => :grays,
        "gist_yarg" => :grays,
        "cool" => :cool,
        "viridis" => :viridis,
        "plasma" => :plasma,
        "inferno" => :inferno,
        "magma" => :magma,
        "YlGnBu_r" => :YlGnBu,
        "YlGnBu" => :YlGnBu,
    )
    
    if haskey(colormap_map, cmap)
        return colormap_map[cmap]
    end
    
    try
        return Symbol(cmap)
    catch
        @warn "Unknown colormap '$cmap', using 'heat' as fallback"
        return :heat
    end
end

function plot_for_radius(subplot, r::String, z_col::String, color, sim_stats_all_pd::DataFrame, args::Dict,
                         list_nlef::Union{Vector{Int}, Nothing} = nothing)
    # Filter by radius_contact, handling missing values
    # Convert missing values to empty string for comparison
    sim_stats_pd = filter(row -> begin
        rc = coalesce(row.radius_contact, "")
        rc == r
    end, sim_stats_all_pd)
    
    if list_nlef !== nothing
        sim_stats_pd = filter(row -> row.nlef in list_nlef, sim_stats_pd)
    end
    
    sim_stats_pd = sort(sim_stats_pd, :sim_hic_file)
    
    x = sim_stats_pd.nlef
    y = sim_stats_pd.km
    z = sim_stats_pd[!, Symbol(z_col)]
    
    if args["plot_mode"] == PLOT_HMAP
        # Check if we have data after filtering
        if nrow(sim_stats_pd) == 0
            @warn "No data found for radius_contact='$r'. Skipping plot."
            return nothing
        end
        
        # Create pivot table similar to pandas pivot_table
        # Handle duplicates by taking the mean (like pandas pivot_table default)
        # First, group by nlef and km, then take mean of z_col
        grouped = groupby(sim_stats_pd, [:nlef, :km])
        aggregated = combine(grouped, Symbol(z_col) => Statistics.mean => Symbol(z_col))
        
        # Now unstack should work without duplicates
        heatmap_data = unstack(aggregated, :nlef, :km, Symbol(z_col))
        
        # Check if unstack produced valid data
        if ncol(heatmap_data) < 2 || nrow(heatmap_data) == 0
            @warn "Insufficient data for heatmap after unstack. Skipping plot."
            return nothing
        end
        
        # Convert DataFrame to matrix for Plots.jl heatmap
        # First column is nlef (index), remaining columns are km values
        nlef_values = heatmap_data[!, 1]  # First column is nlef
        km_values = names(heatmap_data)[2:end]  # Column names are km values
        
        # Extract the data matrix (excluding the first column which is nlef)
        data_matrix = Matrix(heatmap_data[:, 2:end])
        
        # Handle missing values in the matrix
        data_matrix = coalesce.(data_matrix, 0.0)
        
        # Convert column names (km values) to numeric
        km_numeric = [parse(Float64, string(k)) for k in km_values]
        
        # Create heatmap with proper axis labels
        return heatmap(data_matrix, 
                      xticks=(1:length(km_numeric), [@sprintf("%.2e", k) for k in km_numeric]),
                      yticks=(1:length(nlef_values), string.(nlef_values)),
                      xlabel="km", ylabel="Nlef", 
                      title=z_col,
                      colormap=get_plots_colormap(args["cmap"]))
    else
        return scatter3d!(subplot, x, y, z, color=color, marker=:circle)
    end
end

function plot_3d_stats(z_col::String, list_nlef::Union{Vector{Int}, Nothing}, subplot, args::Dict, sim_stats_all_pd::DataFrame)
    legend = String[]
    color_shift = 0
    len_radii = length(args["list_contact_radii"])
    loop_radii(args["list_contact_radii"], list_nlef, color_shift, legend, subplot, z_col, args, sim_stats_all_pd)
    
    if args["plot_mode"] == PLOT_3D
        xlabel!(subplot, "Nlef")
        ylabel!(subplot, "km")
        zlabel!(subplot, z_col)
        plot!(subplot, legend=legend)
    end
end

function loop_radii(radii::Vector{String}, list_nlef::Union{Vector{Int}, Nothing}, color_shift::Int,
                    legend::Vector{String}, subplot, z_col::String, args::Dict, sim_stats_all_pd::DataFrame)
    for (i, r) in enumerate(radii)
        @info "plot radius:$r with color:$(color_shift + i)"
        if args["plot_mode"] == PLOT_HMAP
            fig = plot()
        end
        subplot = plot_for_radius(subplot, r, z_col, CONTACT_RADIUS_PALETTE, sim_stats_all_pd, args, list_nlef)
        
        if args["plot_mode"] == PLOT_HMAP && subplot !== nothing
            basename_file = basename(args["stats_file"])
            title!(fig, "Heatmap of '$z_col' with contact radius:$r for: \n$basename_file")
            xlabel!(fig, "km")
            ylabel!(fig, "Nlef")
            
            file_name = "$(basename_file)_$(args["plot_mode"])_$(z_col)_r$(r).$(args["file_extension"])"
            savefig(fig, joinpath(args["output_folder"], file_name))
            @info "Plot saved in $file_name"
        end
    end
end

function main(args::Dict)
    sim_stats_all_pd = CSV.read(args["stats_file"], DataFrame, types=Dict(:radius_contact => String))
    
    # Extract radius_contact from folder path if it's missing in CSV
    # Pattern: r2.84, r3.55, etc. in the path
    for i in 1:nrow(sim_stats_all_pd)
        if ismissing(sim_stats_all_pd.radius_contact[i]) || isempty(coalesce(sim_stats_all_pd.radius_contact[i], ""))
            # Try to extract from sim_hic_file path (e.g., ./out/demo_run/r2.84/hic_003.hdf5)
            hic_file = sim_stats_all_pd.sim_hic_file[i]
            # Match pattern r<number> or r<number>p in the path
            radius_match = match(r"r(\d+\.?\d*)(p?)", hic_file)
            if radius_match !== nothing
                radius_val = radius_match.captures[1]
                prob_suffix = radius_match.captures[2]
                sim_stats_all_pd.radius_contact[i] = "$radius_val$prob_suffix"
            end
        end
    end
    
    if args["plot_mode"] == PLOT_3D
        fig = plot()
        plt_1 = fig
    else
        plt_1 = nothing
    end
    
    plot_3d_stats(args["z_column"], args["list_nlef"], plt_1, args, sim_stats_all_pd)
    
    if args["plot_mode"] == PLOT_3D
        file_name = "$(basename(args["stats_file"]))_$(args["plot_mode"])_$(args["z_column"])_r$(join(args["list_contact_radii"], "_r")).$(args["file_extension"])"
        savefig(plt_1, joinpath(args["output_folder"], file_name))
        @info "Plot saved in $file_name"
    end
end

function main()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--stats_file", "-f"
            help = "Simulation statistics' repository file (sim_stats.csv)"
            default = "./sim_stats.csv"
        "--output_folder", "-o"
            help = "Output folder to save plots"
            default = "."
        "--file_extension", "-e"
            help = "Image file format extension: png, tif, svg"
            default = "png"
        "--z_column", "-z"
            help = "Column from a --stats_file to be plotted on the z-axis"
            default = "chi2_log"
        "--plot_mode", "-p"
            help = "Plot method. 3d - 3-dimensional scatter plot; hmap - 2D heatmap plot"
            default = PLOT_HMAP
        "--cmap", "-c"
            help = "Color map"
            default = "YlGnBu_r"
        "--list_nlef"
            help = "Show data only for the list of LEFs occupancy (Nlefs)"
            nargs = '*'
            arg_type = Int
            default = Int[]
        "--list_contact_radii", "-k"
            help = "Show data only for the list of contact radii"
            nargs = '+'
            default = ["2.84", "3.55"]
    end
    
    args = parse_args(s)
    # Convert vector arguments from Vector{Any} to proper types
    args["list_nlef"] = length(args["list_nlef"]) > 0 ? [Int(x) for x in args["list_nlef"]] : nothing
    args["list_contact_radii"] = [String(x) for x in args["list_contact_radii"]]
    main(args)
end

end # module

