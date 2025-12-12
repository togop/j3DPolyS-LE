"""
Hi-C Analysis module
"""
module HicAnalysis

import Base.CoreLogging: @info, @warn, @error, @debug
import Base: splitext, basename, dirname, intersect
using HDF5
using DataFrames
using CSV
using Statistics
using Distributions
using LinearAlgebra
using Glob
using PyCall
using Plots
using Base: searchsorted, log10

# Try to import cooler via PyCall
try
    global cooler = PyNULL()
    # Add bioconda channel to conda config if not already present
    try
        # Check if conda is available and add bioconda channel silently
        result = run(`conda config --add channels bioconda`, wait=true, stdout=devnull, stderr=devnull)
    catch
        # Channel might already be added or conda not available, continue anyway
    end
    pyimport_conda("cooler", "cooler", "bioconda")
    @info "Using cooler via PyCall"
catch
    @warn "Could not import cooler. Some functions may not work."
end

# Constants
const SIM_RESOLUTION = 2000
const EXP_RESOLUTION = 10000
const RESOLUTION = 10000
const RES_FACTORS = [1, 2, 4]

const CMAP = "YlOrRd"
const DECAY_CMAP = "tab20"
const DECAY_PALETTE_SHIFT = 2
const DECAY_COLOR_MAX = 20
const DECAY_PLOT_ALPHA = 0.5
const DECAY_USE_HDF5 = true
const DECAY_USE_HDF5_COOL = false

const SAVE_DECAY_PROBABILITY = true

const CHI2_MODE_LOG = "log"
const CHI2_MODE_LINEAR = "linear"
const CHI2_MODE_LOG_ZOOM = false

const CHI2_RANGE_START = 20
const CHI2_RANGE_END = 2000
const CHI2_RANGE_NUM = 100

const CHI2_USE_SEM = true
const CHI2_USE_BALANCED = false

const PLOT_FORMAT = "png"

# Output file names
const DR_OUT = "dr.out"
const CONTACT_OUT = "contact.out"
const CONFIG_OUT = "config.out"
const NLEF_OUT = "Nlef.out"
const PROCESS_OUT = "process.out"
const CHIP_OUT = "chip_lef.out"
const CHIP_BED_GRAPH = "chip_lef.bedGraph"
const XYZCONFIG_OUT = "xyzconfig_001.out"

const ALPHA = "α"

const SIM_CHR_SYNONYMS = ["6", "chrX", "X", "1", "chrI", "I"]
const SIM_CHR = "chrX"

const CHR_X_SYNONYMS = ["6", "chrX", "X"]
const chr_x_size = 17718942
const CHR_A_SYNONYMS = ["1", "chrI", "I"]
const chr_a_size = 15072434

global CHR_SYNONYMS = CHR_X_SYNONYMS
global chr_size = chr_x_size

const DEFAULT_CHIP_CORRELATION = "spearmanr"

const PLOTS_FOLDER = "plots"
const DEMO_CLIM = [-1, 3]
const DEFAULT_CLIM = DEMO_CLIM
const PLOT_COMP_TADS = false

# Logger
const logger = Base.CoreLogging.current_logger()

function get_last_hic(output_folder::String)::String
    if !isdir(output_folder)
        @error "Directory does not exist: $output_folder"
        @error "Please ensure the analysis folder exists or run the analysis step first"
        exit(1)
    end
    files = filter(f -> occursin(r"^hic_.*\.hdf5$", f), readdir(output_folder))
    if isempty(files)
        @error "Not found hic_*.hdf5 files in $output_folder found: $(readdir(output_folder))"
        exit(1)
    end
    sim_hic_snapshot = sort(files)[end]
    return joinpath(output_folder, sim_hic_snapshot)
end

function remove_duplicates(list::Vector, max_val::Union{Int, Nothing} = nothing)::Vector
    final_list = Int[]
    for el in list
        if (max_val === nothing || el <= max_val) && !(el in final_list)
            push!(final_list, el)
        end
    end
    return final_list
end

function hic_to_cooler(hic_file::String, chr::String = SIM_CHR, resolution::Int = SIM_RESOLUTION)::String
    cool_file = "$(hic_file).$(resolution).cool"
    h5open(hic_file, "r") do f
        @info "Converting hic file $hic_file to $cool_file with Keys: $(keys(f)) ..."
        a_group_key = first(keys(f))
        @info "a_group_key0: $a_group_key"
        
        data = read(f[a_group_key])
        hic = Array(data)
        
        @info "data.shape: $(size(hic))"
        
        N = size(hic, 1)
        
        bins_index = [[chr, i * resolution, i * resolution + resolution] for i in 0:(N-1)]
        # 'end' is a reserved keyword, so we use Dict and convert
        bins_dict = Dict(
            :chrom => [b[1] for b in bins_index],
            :start => [b[2] for b in bins_index],
            :end => [b[3] for b in bins_index]
        )
        bins = DataFrame(bins_dict)
        
        pixels_bin1_id = Int64[]
        pixels_bin2_id = Int64[]
        pixels_count = Float64[]
        
        tot_iter = (N - 1) * N / 2
        iter = 0
        for bin1_id in 0:(N-2)
            for bin2_id in (bin1_id+1):(N-1)
                iter += 1
                progress = (iter / tot_iter) * 100
                if progress % 10 == 0
                    @info "pixels progress: $progress%"
                end
                count = hic[bin1_id+1, bin2_id+1]
                if count != 0
                    push!(pixels_bin1_id, bin1_id)
                    push!(pixels_bin2_id, bin2_id)
                    push!(pixels_count, count)
                end
            end
        end
        
        pixels_dic = Dict(
            "bin1_id" => pixels_bin1_id,
            "bin2_id" => pixels_bin2_id,
            "count" => pixels_count
        )
        
        metadata = Dict(
            "format" => "HDF5::Cooler",
            "format-version" => "0.8.6",
            "bin-type" => "fixed",
            "bin-size" => resolution,
            "storage-mode" => "symmetric-upper",
            "assembly" => "ce11",
            "generated-by" => "j3DPolySLE"
        )
        
        if !isfile(cool_file)
            # Use PyCall to create cooler file
            py_cooler = pyimport("cooler")
            py_cooler.create_cooler(cool_file, bins=bins, pixels=pixels_dic, dtypes=Dict("count" => "float64"), ordered=true, metadata=metadata)
        else
            @warn "Cooler file $cool_file already exists and will not be replaced!"
        end
    end
    return cool_file
end

function average_contact_prob(prob_mat::Matrix, dist::Int, plot::Bool = false)
    count = 0
    probs = Float64[]
    x = Int[]
    y = Int[]
    
    for i in (dist+1):size(prob_mat, 1)
        j = i - dist
        count += 1
        push!(probs, prob_mat[i, j])
        if plot
            push!(x, i)
            push!(y, j)
        end
    end
    
    avrg_prob = count > 0 ? mean(probs) : 0.0
    if CHI2_USE_SEM
        # Standard error of the mean: std / sqrt(n)
        sd_sem = count > 1 ? std(probs) / sqrt(count) : 0.0
    else
        sd_sem = count > 1 ? std(probs) : 0.0
    end
    
    return Float64(avrg_prob), Float64(sd_sem)
end

function get_chi2_dist_range(chi2_mode::String, res::Int, max_val::Union{Int, Nothing} = nothing)::Vector{Int}
    if chi2_mode == CHI2_MODE_LINEAR
        linear_end = max_val !== nothing ? min(CHI2_RANGE_END, Int(round(max_val * res / SIM_RESOLUTION))) : CHI2_RANGE_END
        dist_range = range(CHI2_RANGE_START, linear_end, length=CHI2_RANGE_NUM)
    else  # chi2_mode == CHI2_MODE_LOG
        log_range = range(0, log10(CHI2_RANGE_END/CHI2_RANGE_START), length=CHI2_RANGE_NUM)
        dist_range = CHI2_RANGE_START .* 10 .^ log_range
    end
    dist_range = [Int(round(x)) for x in (dist_range .* SIM_RESOLUTION ./ res)]
    dist_range = remove_duplicates(dist_range, max_val)
    return dist_range
end

function get_hic(hic_h5::String, resolution::Int = SIM_RESOLUTION, balance::Bool = CHI2_USE_BALANCED)::Matrix
    if !DECAY_USE_HDF5_COOL
        h5open(hic_h5, "r") do f
            a_group_key = first(keys(f))
            @info "get_hic for $hic_h5 HDF5 Keys $(keys(f)) and dimensions $a_group_key"
            
            data = read(f[a_group_key])
            hic = Array(data)
            fill_diagonal!(hic, 0)
            
            if resolution > SIM_RESOLUTION
                hic = hic_coarsen(hic, hic_h5, resolution)
            end
            return hic
        end
    else
        # Use cooler approach
        cool_file = "$(hic_h5).$(SIM_RESOLUTION).cool"
        if !isfile(cool_file)
            @info "Creating hic $cool_file..."
            cool_file = hic_to_cooler(hic_h5, SIM_CHR, SIM_RESOLUTION)
        end
        @info "Extracting hic from $cool_file..."
        py_cooler = pyimport("cooler")
        hic_cooler = py_cooler.Cooler("$cool_file::/")
        hic_chr_names = [String(x) for x in hic_cooler.chromnames]
        hic_chr = first(intersect(hic_chr_names, SIM_CHR_SYNONYMS))
        balanced = balance && (hic_cooler.bins()["weights"] !== nothing)
        hic = hic_cooler.matrix(balance=balanced).fetch(hic_chr)
        if balanced
            hic = replace(hic, NaN => 0.0)
        end
        hic = hic_coarsen(hic, hic_h5, resolution)
        fill_diagonal!(hic, 0)
        return hic
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
    
    # Coarsen by summing over blocks
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

function get_decay_distribution(hic::Matrix, confidence::Float64 = 0.0, min_dist::Int = 1, max_dist::Union{Int, Nothing} = nothing)::Tuple{Vector{Int}, Vector{Float64}, Vector{Float64}, Vector{Float64}}
    @info "get_decay_distribution: data.shape: $(size(hic)), confidence:$confidence, max_dist:$max_dist"
    
    max_dist = max_dist !== nothing ? min(size(hic, 1), max_dist) : size(hic, 1)
    dists = collect((min_dist+1):max_dist)
    probs = Float64[]
    probs_confi_u = Float64[]
    probs_confi_l = Float64[]
    
    for d in dists
        avrg_prob, sd_sem_prob = average_contact_prob(hic, d)
        push!(probs, avrg_prob)
        
        if confidence != 0.0
            n = max_dist - d
            t_dist = TDist(n - 1)
            test_stat = quantile(t_dist, (confidence + 1) / 2)
            if CHI2_USE_SEM
                sem_prob = sd_sem_prob
            else
                sem_prob = sd_sem_prob / sqrt(n)
            end
            h = test_stat * sem_prob
            push!(probs_confi_u, avrg_prob - h)
            push!(probs_confi_l, avrg_prob + h)
        end
    end
    
    return dists, probs, probs_confi_l, probs_confi_u
end

function get_exp_sim_mcool(hic::String, chrs::Vector{String}, res::Int = RESOLUTION)::String
    exp_sim_cool = endswith(hic, ".cool") ? hic : "$(hic).cool"
    if isfile(exp_sim_cool) && !isfile("$(hic).hdf5")
        # Experimental cooler
        py_cooler = pyimport("cooler")
        hic_cooler = py_cooler.Cooler("$exp_sim_cool::/")
        # Convert PyObject to Int - PyCall should handle this automatically
        root_res = convert(Int, hic_cooler.binsize)
        factors = copy(RES_FACTORS)
        factor_val = res ÷ root_res
        idx = searchsortedfirst(factors, factor_val)
        if idx > length(factors) || factors[idx] != factor_val
            insert!(factors, idx, factor_val)
        end
        resolutions = [Int(i * root_res) for i in factors]
        exp_sim_mcool = replace(exp_sim_cool, r"\.cool$" => ".$(resolutions[1]).mcool")
        
        res_cool = nothing
        try
            res_cool = py_cooler.Cooler("$exp_sim_mcool::/resolutions/$res")
        catch
            @info "Missing resolution $res in $exp_sim_mcool so it will be generated again"
        end
        
        if !isfile(exp_sim_mcool) || res_cool === nothing
            exp_sim_mcool = balance_mcool(exp_sim_cool, resolutions, exp_sim_mcool)
        end
    else
        # Simulation
        root_res = SIM_RESOLUTION
        factors = RES_FACTORS
        insert!(factors, searchsorted(factors, res ÷ root_res), res ÷ root_res)
        resolutions = [Int(i * root_res) for i in factors]
        exp_sim_mcool = "$(hic).$(resolutions[1]).mcool"
        if !isfile(exp_sim_mcool)
            @info "Generating cooler files for $hic"
            exp_sim_mcool = hic_to_mcool(hic, chrs[1], root_res, factors)
        end
    end
    return exp_sim_mcool
end

function hic_to_mcool(hic_file::String, chr::String, resolution::Int, factors::Vector{Int})::String
    cool_file = "$(hic_file).$(resolution).cool"
    if !isfile(cool_file)
        cool_file = hic_to_cooler(hic_file, chr, resolution)
    end
    resolutions = [Int(i * resolution) for i in factors]
    mcool_file = "$(hic_file).$(resolutions[1]).mcool"
    try
        py_cooler = pyimport("cooler")
        py_cooler.zoomify_cooler(cool_file, mcool_file, resolutions=resolutions, chunksize=Int(10e6))
    catch
        @warn "Problem to zoomify file $cool_file so will regenerate it!"
        rm(cool_file, force=true)
        cool_file = hic_to_cooler(hic_file, chr, resolution)
        py_cooler = pyimport("cooler")
        py_cooler.zoomify_cooler(cool_file, mcool_file, resolutions=resolutions, chunksize=Int(10e6))
    end
    return mcool_file
end

function balance_mcool(cool_file::String, resolutions::Vector{Int}, mcool_file::String)::String
    res_str = join(string.(resolutions), ", ")
    cmd = `cooler zoomify --balance --balance-args '--convergence-policy store_nan' -o $mcool_file -c 10000000 -r '$res_str' $cool_file`
    run(cmd)
    
    @info "check balance for $mcool_file"
    py_cooler = pyimport("cooler")
    for res in resolutions
        hic_cooler = py_cooler.Cooler("$mcool_file::/resolutions/$res")
        chr_names = [String(x) for x in hic_cooler.chromnames]
        chr = first(intersect(chr_names, CHR_SYNONYMS))
        hic_mat = hic_cooler.matrix(balance=true).fetch(chr)
        @info "resolution: $res hic_mat.shape: $(size(hic_mat))"
    end
    return mcool_file
end

# Placeholder functions for remaining functionality
# These would need to be fully implemented based on the Python code

function read_tads_bed(tads_bed_file::Union{String, Nothing})::Union{Matrix{Int}, Nothing}
    # Implementation needed
    return nothing
end

function compare_hic_chromosome(hic_file::String, cmp_hic::String;
                                hic_chrs::Union{Vector{String}, Nothing} = nothing,
                                chrs::Vector{String} = CHR_SYNONYMS, res::Int = RESOLUTION,
                                tads_boundary::Union{String, Nothing} = nothing,
                                plots_folder::Union{String, Nothing} = nothing, norm::Bool = true,
                                chi2_mode::String = CHI2_MODE_LOG, hic_balance::Bool = CHI2_USE_BALANCED,
                                plot_cmap::String = CMAP, plot_format::String = PLOT_FORMAT)::Tuple{Float64, Float64}
    # Implementation needed - this is a complex function
    return (0.0, 1.0)
end

# Helper function to get decay distribution from HDF5
function get_decay_distribution_hdf5(hic_h5::String, resolution::Int, confidence::Float64 = 0.0, min_dist::Int = 1, max_dist::Union{Int, Nothing} = nothing)::Tuple{Tuple{Vector{Int}, Vector{Float64}, Vector{Float64}, Vector{Float64}}, Int}
    hic = get_hic(hic_h5, resolution)
    dists, probs, probs_confi_l, probs_confi_u = get_decay_distribution(hic, confidence, min_dist, max_dist)
    return (dists, probs, probs_confi_l, probs_confi_u), size(hic, 1)
end

# Helper function to get decay distribution from cooler
function get_decay_distribution_cool(hic_cooler, chr::String, resolution::Int, confidence::Float64 = 0.0, min_dist::Int = 1, max_dist::Union{Int, Nothing} = nothing)::Tuple{Vector{Int}, Vector{Float64}, Vector{Float64}, Vector{Float64}}
    hic = hic_cooler.matrix(balance=false).fetch(chr)
    fill_diagonal!(hic, 0)
    
    # Convert PyObject binsize to Int - access as Python value and convert
    binsize_py = hic_cooler.binsize
    # Try to convert directly - PyCall should handle Python int to Julia Int
    binsize_val = try
        Int(binsize_py)
    catch
        # Fallback: use Python int() function
        py_int = pyimport("builtins").int
        Int(py_int(binsize_py))
    end
    resolution_factor = Float64(binsize_val) / Float64(resolution)
    min_dist_corrected = min_dist !== nothing ? Int(round(Float64(min_dist) / resolution_factor)) : min_dist
    max_dist_corrected = max_dist !== nothing ? Int(round(Float64(max_dist) / resolution_factor)) : max_dist
    
    dists, probs, probs_confi_l, probs_confi_u = get_decay_distribution(hic, confidence, min_dist_corrected, max_dist_corrected)
    dists_corrected = resolution_factor != 1.0 ? [Int(round(Float64(d) * resolution_factor)) for d in dists] : dists
    return dists_corrected, probs, probs_confi_l, probs_confi_u
end

# Helper function to get cooler and chromosome names
function get_hic_cooler_res(hic_file::String, hic_chrs::Vector{String}, res::Int)
    py_cooler = pyimport("cooler")
    hic_cooler = nothing
    # Try to open as mcool first (multi-resolution cooler)
    try
        hic_cooler = py_cooler.Cooler("$hic_file::/resolutions/$res")
    catch
        # If that fails, try as single-resolution cooler
        try
            hic_cooler = py_cooler.Cooler("$hic_file::/")
        catch e
            @error "Failed to open cooler file $hic_file: $e"
            rethrow(e)
        end
    end
    
    if hic_cooler === nothing
        error("Failed to open cooler file $hic_file")
    end
    
    hic_chr_names = [String(x) for x in hic_cooler.chromnames]
    hic_chrs_select = intersect(hic_chr_names, hic_chrs)
    return hic_cooler, collect(hic_chrs_select)
end

# Helper function to fill diagonal with zeros
function fill_diagonal!(matrix::Matrix, value::Number = 0)
    n = min(size(matrix, 1), size(matrix, 2))
    for i in 1:n
        matrix[i, i] = value
    end
end

function plot_distance_contact_prob_decay(hic_list::Vector{String};
                                           hic_chrs::Vector{String} = CHR_SYNONYMS,
                                           tads::Union{String, Nothing} = nothing,
                                           exp_cool::Union{String, Nothing} = nothing,
                                           output_folder::Union{String, Nothing} = nothing,
                                           res::Int = RESOLUTION,
                                           confidence::Float64 = 0.0,
                                           replace::Bool = true,
                                           chi2_mode::String = CHI2_MODE_LOG,
                                           format::String = PLOT_FORMAT)::String
    if isempty(hic_chrs)
        # If hic_chrs is empty, use CHR_SYNONYMS, but if that's also empty, use default CHR_X_SYNONYMS
        hic_chrs = !isempty(CHR_SYNONYMS) ? CHR_SYNONYMS : CHR_X_SYNONYMS
    end
    
    # Determine output folder
    if output_folder === nothing
        output_folder = dirname(hic_list[1])
    end
    
    # Build filename
    hic_names = [splitext(basename(hic))[1] for hic in hic_list]
    exp_base_name = exp_cool !== nothing ? splitext(basename(exp_cool))[1] : ""
    
    hics = ""
    for (i, hic_r) in enumerate(hic_names)
        hic = hic_list[i]
        hics_suffix = (isfile(hic) && endswith(hic, ".hdf5") && DECAY_USE_HDF5) ? ".h5" : ""
        hics *= "_$(hic_r)$(hics_suffix)"
    end
    
    tads_pref = tads !== nothing ? "_t$(splitext(basename(tads))[1])" : ""
    conf_suffix = confidence > 0.0 ? "_c$(round(confidence, digits=2))" : ""
    chi2_suffix = "_chi2_$(chi2_mode[1:min(3, length(chi2_mode))])"
    sem_suffix = !CHI2_USE_SEM ? "_sd" : ""
    zoom_suffix = CHI2_MODE_LOG_ZOOM ? "_zoom" : ""
    
    # Build filename - we'll determine actual chromosome names from the data
    # Initialize with defaults, will be updated when we process the data
    # Python uses CHR_SYNONYMS[-1] for experimental and hic_chrs[-1] for simulation
    # Use fallback to CHR_X_SYNONYMS if CHR_SYNONYMS is empty
    actual_chr_synonyms = !isempty(CHR_SYNONYMS) ? CHR_SYNONYMS : CHR_X_SYNONYMS
    exp_chr_name = actual_chr_synonyms[end]
    # For simulation, use hic_chrs[-1] (which should be set from parameter or CHR_SYNONYMS/CHR_X_SYNONYMS)
    sim_chr_name = hic_chrs[end]
    
    # Build initial filename - will be updated with correct chromosome names after processing
    hic_decay_plot = joinpath(output_folder,
                              "contact-decay_$(exp_base_name)_$(exp_chr_name)_vs_$(sim_chr_name)$(hics)$(tads_pref)$(conf_suffix)$(chi2_suffix)$(sem_suffix)$(zoom_suffix).$(res).$(format)")
    
    if isfile(hic_decay_plot)
        if !replace
            @info "Distance-contact decay plot already existing: $hic_decay_plot, so skip it"
            return hic_decay_plot
        else
            @warn "Replacing existing Distance-contact decay plot: $hic_decay_plot"
        end
    end
    
    @info "Distance-contact decay plot: $hic_decay_plot ..."
    @info "hic_list: $hic_list"
    @info "exp_cool: $exp_cool"
    @info "output_folder: $output_folder"
    
    max_tad_size = chr_size ÷ SIM_RESOLUTION
    
    min_dist = CHI2_RANGE_START * SIM_RESOLUTION ÷ res
    max_dist = CHI2_RANGE_END * SIM_RESOLUTION ÷ res
    
    # Determine scale attributes - use :log10 for xscale when in log mode (not :log)
    xscale_attr = chi2_mode == CHI2_MODE_LOG ? :log10 : :identity
    if chi2_mode == CHI2_MODE_LOG && !CHI2_MODE_LOG_ZOOM
        min_dist = 1
        max_dist = nothing
    end
    
    # Create empty plot with scales - use :log10 for yscale (not :log)
    # Don't create empty plot - we'll create it with the first data point
    fig = nothing
    
    lines = []
    legend_entries = String[]
    cmp_hic = exp_cool !== nothing ? Base.replace(exp_cool, r".cool$" => "") : ""
    chr_i = 0
    first_plot = true  # Track if this is the first plot to create the figure
    
    for (i, hic) in enumerate(hic_list)
        if !isfile(hic)
            @error "Missing file $hic, skip it!"
            continue
        end
        
        hic_cooler = nothing
        if isfile(hic) && endswith(hic, ".hdf5") && DECAY_USE_HDF5
            # For HDF5 files, use the last element of hic_chrs (matching Python's hic_chrs[-1])
            # If hic_chrs is empty or was set to empty, use CHR_SYNONYMS instead
            actual_hic_chrs = !isempty(hic_chrs) ? hic_chrs : CHR_SYNONYMS
            hic_chrs_select = !isempty(actual_hic_chrs) ? [actual_hic_chrs[end]] : ["chrX"]
            # Update sim_chr_name from actual data - use last element to match Python's hic_chrs[-1]
            sim_chr_name = String(hic_chrs_select[1])
        else
            try
                # hic_chrs should already be set (either from parameter or CHR_SYNONYMS)
                hic_cooler, hic_chrs_select = get_hic_cooler_res(hic, hic_chrs, res)
                # Update sim_chr_name from actual data (use last element, matching Python's hic_chrs[-1])
                if !isempty(hic_chrs_select)
                    sim_chr_name = String(hic_chrs_select[end])
                end
            catch e
                @error "Failed to get cooler for $hic: $e"
                continue
            end
        end
        
        for hic_chr in hic_chrs_select
            if isfile(hic) && endswith(hic, ".hdf5") && DECAY_USE_HDF5
                max_dist_val = max_dist !== nothing ? max_dist : nothing
                (dists, probs, probs_confi_l, probs_confi_u), hic_chr_size = get_decay_distribution_hdf5(hic, res, confidence, min_dist, max_dist_val)
            else
                if hic_cooler === nothing
                    @error "hic_cooler is not defined for $hic"
                    continue
                end
                dists, probs, probs_confi_l, probs_confi_u = get_decay_distribution_cool(hic_cooler, hic_chr, res, confidence, min_dist, max_dist)
                # Convert PyObject chromsize to Int
                try
                    hic_chr_size = Int(hic_cooler.chromsizes[hic_chr])
                catch
                    py_int = pyimport("builtins").int
                    hic_chr_size = Int(py_int(hic_cooler.chromsizes[hic_chr]))
                end
            end
            
            max_tad_size = min(max_tad_size, hic_chr_size ÷ res)
            
            # Calculate chi2 and alpha
            plots_folder = joinpath(output_folder, PLOTS_FOLDER)
            if !isempty(cmp_hic)
                chi2, alpha = compare_hic_chromosome(hic, cmp_hic, hic_chrs=[hic_chr], chrs=CHR_SYNONYMS, res=res,
                                                      tads_boundary=tads, plots_folder=plots_folder, norm=true,
                                                      chi2_mode=chi2_mode)
            else
                chi2, alpha = 0.0, 1.0
            end
            
            probs_adjust = probs .* alpha
            color_idx = chr_i * DECAY_PALETTE_SHIFT
            # Use color palette - tab20 has 20 colors, cycle through them
            color_num = (color_idx % DECAY_COLOR_MAX) + 1
            color = palette(:tab20)[color_num]
            
            legend_label = "$(hic_names[i]).$hic_chr chi2_min: $(round(chi2, digits=3))"
            @info "Plotting $(length(dists)) points: dists range [$(minimum(dists)), $(maximum(dists))], probs range [$(minimum(probs_adjust)), $(maximum(probs_adjust))]"
            # Filter out zero or negative values for log scale - keep only positive values
            valid_idx = probs_adjust .> 0
            num_valid = sum(valid_idx)
            @info "Valid (non-zero) points: $num_valid out of $(length(dists))"
            
            if any(valid_idx)
                dists_plot = dists[valid_idx]
                probs_plot = probs_adjust[valid_idx]
                
                # Create plot with first data point if this is the first plot
                if first_plot
                    @info "Creating plot with $(length(dists_plot)) points: dists [$(minimum(dists_plot)), $(maximum(dists_plot))], probs [$(minimum(probs_plot)), $(maximum(probs_plot))]"
                    fig = plot(dists_plot, probs_plot, 
                              xscale=xscale_attr, yscale=:log10,
                              xlabel="genomic distance in $(res ÷ 1000)kb",
                              ylabel="average #contacts ~ contact probability",
                              title="Distance-contact decay chi2:x$chi2_mode\n$output_folder",
                              color=color, alpha=DECAY_PLOT_ALPHA, label=legend_label,
                              showaxis=true, grid=true, legend=true,
                              framestyle=:box, minorgrid=false)
                    first_plot = false
                else
                    plot!(fig, dists_plot, probs_plot, color=color, alpha=DECAY_PLOT_ALPHA, label=legend_label)
                end
            else
                @warn "All probabilities are zero or negative for $hic, skipping plot"
            end
            push!(lines, (dists, probs_adjust))
            
            if confidence != 0.0 && fig !== nothing
                conf_color_num = ((color_idx + 1) % DECAY_COLOR_MAX) + 1
                conf_color = palette(:tab20)[conf_color_num]
                # Replace zeros in confidence intervals
                probs_confi_u_plot = copy(probs_confi_u)
                probs_confi_u_plot[probs_confi_u_plot .<= 0] .= 1e-10
                probs_confi_l_plot = copy(probs_confi_l)
                probs_confi_l_plot[probs_confi_l_plot .<= 0] .= 1e-10
                plot!(fig, dists, probs_confi_u_plot, color=conf_color, alpha=DECAY_PLOT_ALPHA, label="", linestyle=:dash)
                plot!(fig, dists, probs_confi_l_plot, color=conf_color, alpha=DECAY_PLOT_ALPHA, label="", linestyle=:dash)
            end
            
            push!(legend_entries, legend_label)
            
            if SAVE_DECAY_PROBABILITY
                f_name = "$(hic).$(res)_decay_probs_$(chi2_mode[1:min(3, length(chi2_mode))]).txt"
                @info "Save simulation data from Distance-contact decay plot: $hic to $f_name"
                open(f_name, "w") do f
                    write(f, join([string(elem) for elem in probs], "\n"))
                end
            end
            
            chr_i += 1
        end
    end
    
    # Plot experimental decay
    if exp_cool !== nothing && isfile(exp_cool)
        # Use the original exp_cool path, not cmp_hic (which has .cool removed)
        exp_cooler, exp_chr_select = get_hic_cooler_res(exp_cool, CHR_SYNONYMS, res)
        exp_chr = !isempty(exp_chr_select) ? exp_chr_select[1] : nothing
        
        if exp_chr !== nothing
            # Update exp_chr_name from actual experimental data
            exp_chr_name = String(exp_chr)
            
            dists, probs, probs_confi_l, probs_confi_u = get_decay_distribution_cool(exp_cooler, exp_chr, res, confidence, min_dist, max_dist)
            # Convert PyObject chromsize to Int
            exp_chr_size = try
                Int(exp_cooler.chromsizes[exp_chr])
            catch
                py_int = pyimport("builtins").int
                Int(py_int(exp_cooler.chromsizes[exp_chr]))
            end
            max_tad_size = min(max_tad_size, Int(round(exp_chr_size / res)))
            
            if SAVE_DECAY_PROBABILITY
                f_name = "$(cmp_hic).$(res)_decay_probs_$(chi2_mode[1:min(3, length(chi2_mode))]).txt"
                @info "Save experimental data from Distance-contact decay plot: $exp_cool to $f_name"
                open(f_name, "w") do f
                    write(f, join([string(elem) for elem in probs], "\n"))
                end
            end
            
            exp_legend_label = "$(basename(exp_cool)).$exp_chr"
            @info "Plotting experimental $(length(dists)) points: dists range [$(minimum(dists)), $(maximum(dists))], probs range [$(minimum(probs)), $(maximum(probs))]"
            # Filter out zero or negative values in experimental data
            valid_exp_idx = probs .> 0
            if any(valid_exp_idx)
                dists_exp = dists[valid_exp_idx]
                probs_exp = probs[valid_exp_idx]
                
                # Create plot if it doesn't exist yet (shouldn't happen, but just in case)
                if fig === nothing
                    fig = plot(dists_exp, probs_exp,
                              xscale=xscale_attr, yscale=:log10,
                              xlabel="genomic distance in $(res ÷ 1000)kb",
                              ylabel="average #contacts ~ contact probability",
                              title="Distance-contact decay chi2:x$chi2_mode\n$output_folder",
                              color=:navy, alpha=DECAY_PLOT_ALPHA, label=exp_legend_label,
                              showaxis=true, grid=true, legend=true,
                              framestyle=:box, minorgrid=false)
                else
                    plot!(fig, dists_exp, probs_exp, color=:navy, alpha=DECAY_PLOT_ALPHA, label=exp_legend_label)
                end
            else
                @warn "All experimental probabilities are zero or negative, skipping plot"
            end
            push!(legend_entries, exp_legend_label)
            
            if confidence != 0.0
                # Replace zeros in experimental confidence intervals too
                probs_confi_u_exp = copy(probs_confi_u)
                probs_confi_u_exp[probs_confi_u_exp .<= 0] .= 1e-10
                probs_confi_l_exp = copy(probs_confi_l)
                probs_confi_l_exp[probs_confi_l_exp .<= 0] .= 1e-10
                plot!(fig, dists, probs_confi_u_exp, color=:blue, alpha=DECAY_PLOT_ALPHA, label="", linestyle=:dash)
                plot!(fig, dists, probs_confi_l_exp, color=:blue, alpha=DECAY_PLOT_ALPHA, label="", linestyle=:dash)
            end
        end
    end
    
    # Plot chi2 dist range
    @info "call get_chi2_dist_range($chi2_mode, $res) # calculated max_tad_size:$max_tad_size"
    dist_range = get_chi2_dist_range(chi2_mode, res, max_tad_size)
    if fig !== nothing
        # Add vertical lines for chi2 range boundaries
        plot!(fig, [CHI2_RANGE_START * SIM_RESOLUTION / res], seriestype=:vline, color=:gray, linewidth=1, label="", legend=false)
        plot!(fig, [CHI2_RANGE_END * SIM_RESOLUTION / res], seriestype=:vline, color=:gray, linewidth=1, label="", legend=false)
        
        for d in dist_range
            plot!(fig, [d], seriestype=:vline, color=:gray, linewidth=1, label="", linestyle=:dash, alpha=0.5, legend=false)
        end
    end
    
    # Rebuild filename with correct chromosome names now that we've processed the data
    hic_decay_plot = joinpath(output_folder,
                              "contact-decay_$(exp_base_name)_$(exp_chr_name)_vs_$(sim_chr_name)$(hics)$(tads_pref)$(conf_suffix)$(chi2_suffix)$(sem_suffix)$(zoom_suffix).$(res).$(format)")
    
    # Ensure plot exists and has data
    if fig === nothing
        @error "No plot was created! No data to plot."
        return ""
    end
    
    # Update title if needed (labels are already set when creating the plot)
    plot!(fig, title="Distance-contact decay chi2:x$chi2_mode\n$output_folder")
    
    # Ensure plot has proper axis settings - explicitly set to show axes
    plot!(fig, showaxis=true, grid=true, framestyle=:box, minorgrid=false, 
          xguide="genomic distance in $(res ÷ 1000)kb",
          yguide="average #contacts ~ contact probability")
    
    savefig(fig, hic_decay_plot)
    @info " save plot: $hic_decay_plot"
    return hic_decay_plot
end

function chip_out_to_bedgraph(chip_out_file::String, bed_graph_file::Union{String, Nothing} = nothing,
                               chrom::String = SIM_CHR, resolution::Int = SIM_RESOLUTION)::String
    # Implementation needed
    return ""
end

function plot_chip_seq(chip_out_file::String, exp_chip::Union{String, Nothing}, boundary::Union{String, Nothing},
                      resolution::Int = SIM_RESOLUTION, correlation::String = DEFAULT_CHIP_CORRELATION,
                      plot::Bool = true, replace::Bool = true)::Tuple{Float64, Float64, Float64}
    # Implementation needed
    return (0.0, 0.0, 0.0)
end

end # module

