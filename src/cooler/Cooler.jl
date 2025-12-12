"""
Cooler.jl - Pure Julia implementation of cooler file format

This module provides functionality to read and write cooler/mcool files,
which are HDF5-based formats for storing genomic interaction matrices (Hi-C data).

Based on the cooler library: https://github.com/open2c/cooler
"""
module Cooler

import Base.CoreLogging: @info, @warn, @error, @debug
using HDF5
using DataFrames
using Statistics

export CoolerFile, CoolerMatrix
export create_cooler, zoomify_cooler
export chromnames, chromsizes, binsize, bins, matrix, fetch

# Type definitions
struct CoolerFile
    uri::String
    h5file::HDF5.File
    root::String

    function CoolerFile(uri::String)
        # Parse URI: file.cool::/path or file.mcool::/resolutions/1000
        parts = split(uri, "::")
        if length(parts) != 2
            error("Invalid cooler URI format. Expected 'file.cool::/' or 'file.mcool::/resolutions/resolution'")
        end

        filepath = parts[1]
        root = parts[2]

        if !isfile(filepath)
            error("Cooler file not found: $filepath")
        end

        h5file = h5open(filepath, "r")

        # Verify this is a valid cooler file
        if !haskey(h5file, root)
            close(h5file)
            error("Root path $root not found in $filepath")
        end

        # Check for required groups
        grp = h5file[root]
        required = ["chroms", "bins", "pixels"]
        for req in required
            if !haskey(grp, req)
                close(h5file)
                error("Required group '$req' not found in cooler at $root")
            end
        end

        new(uri, h5file, root)
    end
end

struct CoolerMatrix
    cooler::CoolerFile
    balance::Bool

    function CoolerMatrix(cooler::CoolerFile, balance::Bool=false)
        new(cooler, balance)
    end
end

# Property accessors
function chromnames(cooler::CoolerFile)::Vector{String}
    grp = cooler.h5file[cooler.root]
    chrom_names = read(grp["chroms/name"])
    return String.(chrom_names)
end

function chromsizes(cooler::CoolerFile)::Dict{String, Int}
    grp = cooler.h5file[cooler.root]
    chrom_names = read(grp["chroms/name"])
    chrom_lengths = read(grp["chroms/length"])

    sizes = Dict{String, Int}()
    for (name, length) in zip(chrom_names, chrom_lengths)
        sizes[String(name)] = Int(length)
    end
    return sizes
end

function binsize(cooler::CoolerFile)::Int
    grp = cooler.h5file[cooler.root]

    # Try to read bin-size from metadata
    if haskey(attributes(grp), "bin-size")
        return Int(read(attributes(grp)["bin-size"]))
    end

    # Calculate from bins
    bins_group = grp["bins"]
    starts = read(bins_group["start"])
    ends = read(bins_group["end"])

    if length(starts) > 0
        return Int(ends[1] - starts[1])
    end

    error("Could not determine bin size")
end

function bins(cooler::CoolerFile)::Dict{String, Any}
    grp = cooler.h5file[cooler.root]
    bins_group = grp["bins"]

    result = Dict{String, Any}()
    result["chrom"] = read(bins_group["chrom"])
    result["start"] = read(bins_group["start"])
    result["end"] = read(bins_group["end"])

    # Check for weights (ICE balancing)
    if haskey(bins_group, "weight")
        result["weights"] = read(bins_group["weight"])
    else
        result["weights"] = nothing
    end

    return result
end

function matrix(cooler::CoolerFile; balance::Bool=false)::CoolerMatrix
    return CoolerMatrix(cooler, balance)
end

function fetch(mat::CoolerMatrix, region::String)::Matrix{Float64}
    return fetch(mat, (region,))
end

function fetch(mat::CoolerMatrix, region::Tuple)::Matrix{Float64}
    cooler = mat.cooler
    grp = cooler.h5file[cooler.root]

    # Parse region
    if length(region) == 1
        # Whole chromosome: "chr1"
        chrom = String(region[1])
        start_bp = nothing
        end_bp = nothing
    elseif length(region) == 3
        # Region with coordinates: ("chr1", 100000, 200000)
        chrom = String(region[1])
        start_bp = Int(region[2])
        end_bp = Int(region[3])
    else
        error("Invalid region format. Use \"chr\" or (\"chr\", start, end)")
    end

    # Read bins
    bins_group = grp["bins"]
    bin_chroms = read(bins_group["chrom"])
    bin_starts = read(bins_group["start"])
    bin_ends = read(bins_group["end"])

    # Read chromosome mapping
    chroms_group = grp["chroms"]
    chrom_names = String.(read(chroms_group["name"]))

    # Find chromosome index
    chrom_idx = findfirst(==(chrom), chrom_names)
    if chrom_idx === nothing
        error("Chromosome $chrom not found in cooler file")
    end
    chrom_id = chrom_idx - 1  # 0-based indexing

    # Find bin range
    chrom_mask = bin_chroms .== chrom_id
    chrom_bin_indices = findall(chrom_mask)

    if isempty(chrom_bin_indices)
        error("No bins found for chromosome $chrom")
    end

    # Filter by coordinates if provided
    if start_bp !== nothing && end_bp !== nothing
        coord_mask = (bin_starts .>= start_bp) .& (bin_ends .<= end_bp)
        bin_indices = findall(chrom_mask .& coord_mask)
    else
        bin_indices = chrom_bin_indices
    end

    if isempty(bin_indices)
        error("No bins found in region")
    end

    bin_start = minimum(bin_indices)
    bin_end = maximum(bin_indices)
    n_bins = bin_end - bin_start + 1

    # Read pixels
    pixels_group = grp["pixels"]
    pixel_bin1 = read(pixels_group["bin1_id"])
    pixel_bin2 = read(pixels_group["bin2_id"])
    pixel_count = read(pixels_group["count"])

    # Initialize matrix
    matrix_data = zeros(Float64, n_bins, n_bins)

    # Fill matrix with pixel data
    for i in 1:length(pixel_bin1)
        b1 = pixel_bin1[i] + 1  # Convert to 1-based
        b2 = pixel_bin2[i] + 1
        count = pixel_count[i]

        # Check if pixels are in our range
        if b1 >= bin_start && b1 <= bin_end && b2 >= bin_start && b2 <= bin_end
            row = b1 - bin_start + 1
            col = b2 - bin_start + 1

            matrix_data[row, col] = count
            # Cooler stores upper triangle, make symmetric
            if row != col
                matrix_data[col, row] = count
            end
        end
    end

    # Apply balancing if requested
    if mat.balance
        bins_data = bins(cooler)
        if bins_data["weights"] !== nothing
            weights = bins_data["weights"]
            bin_weights = weights[bin_start:bin_end]

            # Apply ICE balancing: M_ij = count_ij * weight_i * weight_j
            for i in 1:n_bins
                for j in 1:n_bins
                    w_i = bin_weights[i]
                    w_j = bin_weights[j]

                    # Handle NaN weights
                    if isnan(w_i) || isnan(w_j) || isinf(w_i) || isinf(w_j)
                        matrix_data[i, j] = NaN
                    else
                        matrix_data[i, j] *= w_i * w_j
                    end
                end
            end
        end
    end

    return matrix_data
end

function create_cooler(output_uri::String;
                       bins::DataFrame,
                       pixels::Dict{String, Vector},
                       dtypes::Dict{String, String}=Dict{String, String}(),
                       ordered::Bool=true,
                       metadata::Union{Dict{String, Any}, Nothing}=nothing,
                       assembly::Union{String, Nothing}=nothing,
                       h5opts::Union{Dict, Nothing}=nothing)

    # Parse URI
    parts = split(output_uri, "::")
    if length(parts) == 2
        filepath = parts[1]
        root = parts[2]
    else
        filepath = output_uri
        root = "/"
    end

    # Validate inputs
    if !hasproperty(bins, :chrom) || !hasproperty(bins, :start) || !hasproperty(bins, Symbol("end"))
        error("bins DataFrame must have columns: chrom, start, end")
    end

    if !haskey(pixels, "bin1_id") || !haskey(pixels, "bin2_id") || !haskey(pixels, "count")
        error("pixels Dict must have keys: bin1_id, bin2_id, count")
    end

    # Create/open HDF5 file
    # For cooler files, always use "w" mode to create fresh or "r+" for existing with new resolution
    # Check if this is adding a new resolution to existing mcool
    is_new_resolution = isfile(filepath) && root != "/"
    mode = is_new_resolution ? "r+" : "w"
    h5file = h5open(filepath, mode)

    try
        # Create root group if it doesn't exist
        if root != "/" && !haskey(h5file, root)
            create_group(h5file, root)
        end

        grp = root == "/" ? h5file : h5file[root]

        # Store metadata (only if attribute doesn't already exist)
        if metadata !== nothing
            for (key, value) in metadata
                if !haskey(attributes(grp), key)
                    attributes(grp)[key] = value
                end
            end
        end

        # Create chroms group
        if !haskey(grp, "chroms")
            create_group(grp, "chroms")
        end
        chroms_group = grp["chroms"]

        # Get unique chromosomes from bins
        unique_chroms = unique(bins.chrom)
        chrom_lengths = Int[]
        for chrom in unique_chroms
            chrom_bins = bins[bins.chrom .== chrom, :]
            push!(chrom_lengths, maximum(chrom_bins[!, Symbol("end")]))
        end

        # Write chromosome data
        chroms_group["name"] = collect(String.(unique_chroms))
        chroms_group["length"] = chrom_lengths

        # Create bins group
        if !haskey(grp, "bins")
            create_group(grp, "bins")
        end
        bins_group = grp["bins"]

        # Map chromosome names to IDs
        chrom_to_id = Dict(chrom => i-1 for (i, chrom) in enumerate(unique_chroms))
        bin_chrom_ids = [chrom_to_id[c] for c in bins.chrom]

        bins_group["chrom"] = bin_chrom_ids
        bins_group["start"] = Int.(bins.start)
        bins_group["end"] = Int.(bins[!, Symbol("end")])

        # Create pixels group
        if !haskey(grp, "pixels")
            create_group(grp, "pixels")
        end
        pixels_group = grp["pixels"]

        # Sort pixels if ordered
        if ordered
            # Create sorted indices
            sort_order = sortperm(collect(zip(pixels["bin1_id"], pixels["bin2_id"])))
            pixels_group["bin1_id"] = Int.(pixels["bin1_id"][sort_order])
            pixels_group["bin2_id"] = Int.(pixels["bin2_id"][sort_order])
            pixels_group["count"] = Float64.(pixels["count"][sort_order])
        else
            pixels_group["bin1_id"] = Int.(pixels["bin1_id"])
            pixels_group["bin2_id"] = Int.(pixels["bin2_id"])
            pixels_group["count"] = Float64.(pixels["count"])
        end

        # Add indexes for efficient querying
        bin1_offset = zeros(Int, length(bin_chrom_ids) + 1)
        current_bin1 = -1
        pixel_idx = 0

        for i in 1:length(pixels_group["bin1_id"])
            bin1 = read(pixels_group["bin1_id"])[i]
            if bin1 != current_bin1
                for b in (current_bin1+2):(bin1+1)
                    bin1_offset[b] = pixel_idx
                end
                current_bin1 = bin1
            end
            pixel_idx += 1
        end
        # Fill remaining
        for b in (current_bin1+2):length(bin1_offset)
            bin1_offset[b] = pixel_idx
        end

        indexes_group = create_group(grp, "indexes")
        indexes_group["bin1_offset"] = bin1_offset
        indexes_group["chrom_offset"] = zeros(Int, length(unique_chroms) + 1)

    finally
        close(h5file)
    end

    @info "Created cooler file: $output_uri"
end

# Cleanup
function Base.close(cooler::CoolerFile)
    close(cooler.h5file)
end

end # module
