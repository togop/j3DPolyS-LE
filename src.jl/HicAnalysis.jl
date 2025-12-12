"""
Hi-C Analysis module
"""
module HicAnalysis

import Base.CoreLogging: @info, @warn, @error, @debug
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

function plot_distance_contact_prob_decay(hic_list::Vector{String}, hic_chrs::Vector{String} = CHR_SYNONYMS,
                                           tads::Union{String, Nothing} = nothing, exp_cool::Union{String, Nothing} = nothing,
                                           output_folder::Union{String, Nothing} = nothing, res::Int = RESOLUTION,
                                           confidence::Float64 = 0.0, replace::Bool = true,
                                           chi2_mode::String = CHI2_MODE_LOG, format::String = PLOT_FORMAT)::String
    # Implementation needed
    return ""
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

