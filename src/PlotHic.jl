"""
Plot Hi-C module
"""
module PlotHic

import Base.CoreLogging: @info, @warn, @error, @debug
using HDF5
using Plots
using Glob
using PyCall
using ArgParse
using LinearAlgebra
using Base: log10

# Use Julia cooler implementation
include(joinpath(@__DIR__, "cooler", "CoolerModule.jl"))
using .CoolerModule
@info "Using Julia cooler implementation"

const DEFAULT_CMAP = "hot_r"
const DEFAULT_RESOLUTION = 2000
const DEFAULT_HIC_WILDCARD = "hic*.hdf5"
const DEFAULT_CLIM = [-2.75, 0]
const DEFAULT_PLOT_FORMAT = "png"
const CHR_X_SYNONYMS = ["6", "chrX", "X"]
const DEFAULT_TITLE = "Hi-C for measurement {hic_file}"
const DEFAULT_DPI = 150
const SIM_RESOLUTION = 2000

# Helper function to fill diagonal with a value (equivalent to np.fill_diagonal)
function fill_diagonal!(matrix::Matrix, value::Number)
    n = min(size(matrix, 1), size(matrix, 2))
    for i in 1:n
        matrix[i, i] = value
    end
end

# Map matplotlib colormap names to Plots.jl colormap symbols
# Plots.jl uses different names and doesn't support "_r" suffix for reversed colormaps
function get_plots_colormap(cmap::String)::Symbol
    # Map common matplotlib colormaps to Plots.jl equivalents
    colormap_map = Dict(
        "hot_r" => :heat,           # reversed hot -> heat (similar appearance)
        "hot" => :heat,
        "gist_heat_r" => :heat,
        "afmhot_r" => :heat,
        "YlOrRd" => :YlOrRd,        # Available in Plots.jl
        "Greys" => :grays,
        "gist_yarg" => :grays,
        "cool" => :cool,
        "viridis" => :viridis,
        "plasma" => :plasma,
        "inferno" => :inferno,
        "magma" => :magma,
        "YlGnBu_r" => :YlGnBu,      # reversed -> normal
        "YlGnBu" => :YlGnBu,
    )
    
    # Check if we have a direct mapping
    if haskey(colormap_map, cmap)
        return colormap_map[cmap]
    end
    
    # Try to convert to Symbol directly (might work for some colormaps)
    try
        return Symbol(cmap)
    catch
        # Fallback to heat if unknown
        @warn "Unknown colormap '$cmap', using 'heat' as fallback"
        return :heat
    end
end

function get_hic(hic_file::String, resolution::Int, balanced::Bool, hic_chrs::Vector{String})::Matrix
    if endswith(hic_file, ".hdf5")
        h5open(hic_file, "r") do f
            a_group_key = first(keys(f))
            @info "get_hic for $hic_file HDF5 Keys $(keys(f)) and dimensions $a_group_key"
            data = read(f[a_group_key])
            hic = Array(data)
            fill_diagonal!(hic, 0)
            if resolution > SIM_RESOLUTION
                hic = hic_coarsen(hic, hic_file, resolution)
            end
            return hic
        end
    elseif endswith(hic_file, ".cool") || endswith(hic_file, ".mcool")
        cooler_ref = endswith(hic_file, ".cool") ? "$hic_file::/" : "$hic_file::/resolutions/$resolution"
        hic_cooler = CoolerFile(cooler_ref)
        hic_chr_names = chromnames(hic_cooler)
        hic_chr = first(intersect(hic_chr_names, hic_chrs))
        balance = balanced && (bins(hic_cooler)["weights"] !== nothing)
        mat_obj = matrix(hic_cooler, balance=balance)
        hic = fetch(mat_obj, hic_chr)
        if balanced
            hic = replace(hic, NaN => 0.0)
        end
        fill_diagonal!(hic, 0)
        return hic
    else
        error("Unsupported file format: $hic_file")
    end
end

function hic_coarsen(hic::Matrix, hic_h5::String, resolution::Int, from_resolution::Int = SIM_RESOLUTION)::Matrix
    factor = resolution ÷ from_resolution
    if factor == 1
        return hic
    end
    if resolution % from_resolution != 0
        @error "Cannot do binning from resolution $from_resolution to $resolution for hic file $hic_h5"
    end
    
    N = size(hic, 1)
    new_N = N ÷ factor
    hic_binned = zeros(new_N, new_N)
    
    for i in 1:new_N
        for j in 1:new_N
            i_start = (i-1)*factor + 1
            i_end = i*factor
            j_start = (j-1)*factor + 1
            j_end = j*factor
            hic_binned[i, j] = sum(hic[i_start:i_end, j_start:j_end])
        end
    end
    
    return hic_binned
end

function run(output_folder::String; hic_wildcard::String = DEFAULT_HIC_WILDCARD,
             resolution::Int = DEFAULT_RESOLUTION, balanced::Bool = false,
             hic_chrs::Vector{String} = CHR_X_SYNONYMS, cmap::String = DEFAULT_CMAP,
             clim::Vector{Float64} = DEFAULT_CLIM, title::String = DEFAULT_TITLE,
             plot_format::String = DEFAULT_PLOT_FORMAT)
    @info "Plotting HiC for $output_folder with color map: $cmap in file format: $plot_format ..."
    
    if isempty(cmap)
        cmap = DEFAULT_CMAP
    end
    if isempty(plot_format)
        plot_format = DEFAULT_PLOT_FORMAT
    end
    
    pattern = joinpath(output_folder, hic_wildcard)
    hic_files = sort(glob(pattern))
    
    for hic_file in hic_files
        hic = get_hic(hic_file, resolution, balanced, hic_chrs)
        @info "data.shape: $(size(hic))"
        
        hic_log = log10.(replace(hic, 0.0 => 1e-10))
        
        res_kb = resolution / 1000
        n = size(hic_log, 1)
        # Axis formatters: show bin indices as genomic position in kb (like Python FuncFormatter)
        tick_step = max(1, n ÷ 10)
        tick_pos = 1:tick_step:n
        tick_labels = round.(Int, (tick_pos .* res_kb))
        
        # Title: clean hic_file like Python (replace '_hic_003.hdf5', '' and './')
        hic_file_clean = replace(hic_file, "_hic_003.hdf5" => "")
        hic_file_clean = replace(hic_file_clean, "./" => "")
        title_str = replace(title, "{hic_file}" => hic_file_clean)
        title_str = replace(title_str, "{output_folder}" => output_folder)
        title_str = replace(title_str, "{resolution}" => string(resolution))
        title_str = replace(title_str, "{balanced}" => string(balanced))
        
        p = heatmap(hic_log, colormap=get_plots_colormap(cmap), title=title_str,
                    xticks=(tick_pos, string.(tick_labels)),
                    yticks=(tick_pos, string.(tick_labels)))
        
        if length(clim) > 1
            plot!(p, clims=(clim[1], clim[2]))
        end
        
        basename_file = basename(hic_file)
        extension = split(basename_file, ".")[end]
        plot_file = joinpath(output_folder, replace(basename_file, ".$extension" => "_$(cmap).$(plot_format)"))
        
        savefig(p, plot_file; dpi=DEFAULT_DPI)
        @info "Hic plot saved in file $plot_file"
    end
end

function main()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--output_folder", "-o"
            help = "'Analysis' step output folder containing raw hic_*.hdf5 files"
            default = "."
        "--hic_wildcard", "-w"
            help = "Hi-C files wildcard"
            default = DEFAULT_HIC_WILDCARD
        "--hic_chrs"
            help = "Synonyms of the chromosome from Hi-C matrixes to be plotted"
            nargs = '*'
            default = CHR_X_SYNONYMS
        "--resolution", "-r"
            help = "Hi-C data resolution in bp"
            arg_type = Int
            default = DEFAULT_RESOLUTION
        "--balanced", "-b"
            help = "Get balanced if available, for .cool and .mcool"
            action = :store_true
        "--cmap", "-c"
            help = "Color map"
            default = DEFAULT_CMAP
        "--clim"
            help = "Set the color limits of the current image"
            nargs = '*'
            arg_type = Float64
            default = DEFAULT_CLIM
        "--plot_format", "-f"
            help = "Image file format extension: png, tif, svg"
            default = DEFAULT_PLOT_FORMAT
        "--title"
            help = "Plot title template"
            default = DEFAULT_TITLE
    end
    
    args = parse_args(s)
    
    # Convert vectors from Vector{Any} to proper types
    clim = length(args["clim"]) > 0 ? [Float64(x) for x in args["clim"]] : DEFAULT_CLIM
    hic_chrs = [String(x) for x in args["hic_chrs"]]
    
    run(args["output_folder"],
        hic_wildcard=args["hic_wildcard"],
        resolution=args["resolution"],
        balanced=args["balanced"],
        hic_chrs=hic_chrs,
        cmap=args["cmap"],
        clim=clim,
        title=args["title"],
        plot_format=args["plot_format"])
end

end # module

