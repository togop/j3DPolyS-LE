"""
Zoomify.jl - Multi-resolution cooler file generation

Creates .mcool files with multiple resolutions by coarsening Hi-C data.
"""
module Zoomify

import Base.CoreLogging: @info, @warn, @error
using HDF5
using DataFrames
using ..Cooler: CoolerFile, chromnames, chromsizes, binsize, bins, create_cooler

export zoomify_cooler

function coarsen_pixels(pixels::Dict{String, Vector}, factor::Int, n_bins_original::Int)::Dict{String, Vector}
    """Coarsen pixels by aggregating into larger bins"""
    n_bins_new = div(n_bins_original, factor)

    # Map old bins to new bins
    bin1_new = div.(pixels["bin1_id"], factor)
    bin2_new = div.(pixels["bin2_id"], factor)

    # Aggregate counts
    aggregated = Dict{Tuple{Int, Int}, Float64}()
    for i in 1:length(bin1_new)
        key = (bin1_new[i], bin2_new[i])
        aggregated[key] = get(aggregated, key, 0.0) + pixels["count"][i]
    end

    # Convert back to vectors
    new_pixels = Dict{String, Vector}()
    new_pixels["bin1_id"] = Int[]
    new_pixels["bin2_id"] = Int[]
    new_pixels["count"] = Float64[]

    for ((b1, b2), count) in sort(collect(aggregated))
        push!(new_pixels["bin1_id"], b1)
        push!(new_pixels["bin2_id"], b2)
        push!(new_pixels["count"], count)
    end

    return new_pixels
end

function zoomify_cooler(input_uri::String, output_uri::String;
                        resolutions::Vector{Int},
                        chunksize::Int=10000000,
                        balance::Bool=false,
                        balance_args::Union{Dict, Nothing}=nothing)

    @info "Zoomifying cooler: $input_uri -> $output_uri with resolutions: $resolutions"

    # Open input cooler
    input_cooler = CoolerFile(input_uri)
    base_resolution = binsize(input_cooler)

    # Check that resolutions are multiples of base resolution
    for res in resolutions
        if res % base_resolution != 0
            error("Resolution $res is not a multiple of base resolution $base_resolution")
        end
    end

    # Sort resolutions
    resolutions = sort(resolutions)

    # Parse output URI
    parts = split(output_uri, "::")
    output_filepath = parts[1]

    # Create output mcool file
    mode = isfile(output_filepath) ? "r+" : "w"
    output_h5 = h5open(output_filepath, mode)

    try
        # Create resolutions group
        if !haskey(output_h5, "resolutions")
            create_group(output_h5, "resolutions")
        end

        # Read original data
        input_grp = input_cooler.h5file[input_cooler.root]
        input_bins = bins(input_cooler)
        input_bins_group = input_grp["bins"]

        # Read original bins as DataFrame
        # Note: 'end' is a reserved keyword, so we use Symbol("end") for the column name
        bins_df = DataFrame()
        bins_df[!, :chrom] = input_bins["chrom"]
        bins_df[!, :start] = input_bins["start"]
        bins_df[!, Symbol("end")] = input_bins["end"]

        # Map chrom IDs back to names
        chrom_names_list = chromnames(input_cooler)
        bins_df[!, :chrom] = [chrom_names_list[id + 1] for id in bins_df[!, :chrom]]

        # Read original pixels
        input_pixels_group = input_grp["pixels"]
        original_pixels = Dict{String, Vector}()
        original_pixels["bin1_id"] = read(input_pixels_group["bin1_id"])
        original_pixels["bin2_id"] = read(input_pixels_group["bin2_id"])
        original_pixels["count"] = read(input_pixels_group["count"])

        n_bins_original = length(input_bins["start"])

        # Process each resolution
        for res in resolutions
            @info "Creating resolution: $res"

            factor = div(res, base_resolution)
            resolution_path = "/resolutions/$res"

            # Skip if resolution already exists
            if haskey(output_h5, resolution_path)
                @warn "Resolution $res already exists, skipping"
                continue
            end

            if factor == 1
                # Copy original data
                new_bins_df = bins_df
                new_pixels = original_pixels
            else
                # Coarsen bins
                n_bins_new = div(n_bins_original, factor)
                new_starts = [bins_df[i*factor + 1, :start] for i in 0:(n_bins_new-1)]
                new_ends = [bins_df[min((i+1)*factor, n_bins_original), Symbol("end")] for i in 0:(n_bins_new-1)]
                new_chroms = [bins_df[i*factor + 1, :chrom] for i in 0:(n_bins_new-1)]

                new_bins_df = DataFrame()
                new_bins_df[!, :chrom] = new_chroms
                new_bins_df[!, :start] = new_starts
                new_bins_df[!, Symbol("end")] = new_ends

                # Coarsen pixels
                new_pixels = coarsen_pixels(original_pixels, factor, n_bins_original)
            end

            # Create metadata
            metadata = Dict{String, Any}(
                "bin-size" => res,
                "bin-type" => "fixed",
                "format" => "HDF5::Cooler",
                "format-version" => 3,
                "storage-mode" => "symmetric-upper"
            )

            # Read original metadata if available
            if haskey(attributes(input_grp), "genome-assembly")
                metadata["genome-assembly"] = read(attributes(input_grp)["genome-assembly"])
            end

            # Write cooler at this resolution
            close(output_h5)  # Close to allow create_cooler to open
            create_cooler("$output_filepath::$resolution_path",
                         bins=new_bins_df,
                         pixels=new_pixels,
                         metadata=metadata,
                         ordered=true)
            output_h5 = h5open(output_filepath, "r+")  # Reopen
        end

    finally
        close(output_h5)
        close(input_cooler.h5file)
    end

    @info "Zoomify complete: $output_uri"
end

end # module
