"""
Hi-C converters module
"""
module HicConverters

import Base.CoreLogging: @info, @warn, @error, @debug
using HDF5
using DataFrames
using CSV
using PyCall
using ArgParse

# Use Julia cooler implementation
include(joinpath(@__DIR__, "cooler", "CoolerModule.jl"))
using .CoolerModule
@info "Using Julia cooler implementation"

function read_hic_hdf5(hic_hdf5::String)::Matrix
    h5open(hic_hdf5, "r") do f
        @info "Read hic file $hic_hdf5 with Keys: $(keys(f)) ..."
        a_group_key = first(keys(f))
        data = read(f[a_group_key])
        hic = Array(data)
        return hic
    end
end

function hic_to_cool(hic::Matrix, chr::String, resolution::Int, cool_file::String)::String
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
        "genome-assembly" => "ce11",
        "generated-by" => "j3DPolySLE-2026.1"
    )
    
    create_cooler(cool_file, bins=bins, pixels=pixels_dic, dtypes=Dict("count" => "float64"), ordered=true, metadata=metadata)
    return cool_file
end

function hic_to_mcool(hic::Matrix, chr::String, resolutions::Vector{Int}, mcool_file::String)::String
    cool_file = replace(mcool_file, ".mcool" => ".cool")
    hic_to_cool(hic, chr, resolutions[1], cool_file)
    
    res_str = join(string.(resolutions), ", ")
    cmd = `cooler zoomify -o $mcool_file -c 10000000 -r '$res_str' $cool_file`
    @info "call: $cmd"
    run(cmd)
    return mcool_file
end

function main(input_file::String, output_file::String, chr::String = "chrS", resolutions::Vector{Int} = [2000])
    if endswith(input_file, ".hdf5")
        hic = read_hic_hdf5(input_file)
    elseif endswith(input_file, ".mat")
        # Note: SciPy.jl may not be available, would need PyCall to scipy.io
        @error "MATLAB .mat file support requires SciPy via PyCall - not implemented"
        exit(1)
    else
        @error "Unsupported input format $input_file. Supported file types: .hdf5 (3DPolyS_LE format), .mat (MATLAB 2D matrix)"
        exit(1)
    end
    
    if endswith(output_file, ".cool")
        hic_to_cool(hic, chr, resolutions[1], output_file)
    elseif endswith(output_file, ".mcool")
        hic_to_mcool(hic, chr, resolutions, output_file)
    elseif endswith(output_file, ".mat")
        # Note: SciPy.jl may not be available, would need PyCall to scipy.io
        @error "MATLAB .mat file support requires SciPy via PyCall - not implemented"
        exit(1)
    else
        @error "Unsupported output format $output_file. Supported file types: .cool, .mcool, .mat"
        exit(1)
    end
end

function main()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--input_file", "-i"
            help = "Input file name. Supported file types: .hdf5 (3DPolyS_LE format), .mat (MATLAB 2D matrix)"
            required = true
        "--output_file", "-o"
            help = "Output file name. Allowed file types: .cool, .mcool"
            required = true
        "--chr"
            help = "Chromosome name to be used for"
            default = "chrS"
        "--resolutions", "-r"
            help = "List of resolutions for .mcool output file"
            nargs = '+'
            arg_type = Int
            default = [2000]
    end
    
    args = parse_args(s)
    main(args["input_file"], args["output_file"], args["chr"], args["resolutions"])
end

function hdf5_to_cooler()
    main()
end

end # module

