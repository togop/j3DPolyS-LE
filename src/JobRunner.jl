"""
JobRunner module for managing job execution
"""
module JobRunner

# Note: doc! patch is defined in the main j3DPolySLE module to avoid duplicate definitions

import Base.CoreLogging: @info, @warn, @error, @debug
using TOML
using Base: searchsorted

export CfgJobRunner, CFG_SECTION_3DPOLYS_LE, convert_dat_to_cfg, _normalize_config_file

const CFG_SECTION_3DPOLYS_LE = "3dpolys_le"

# Abstract type - renamed to avoid conflict with module name and documentation issues
abstract type AbstractJobRunner end

# Alias for backward compatibility (not exported)
const JobRunner = AbstractJobRunner

function run_cmd(jr::AbstractJobRunner, cmd::String, dep_jobid::String, profile::String = "")::String
    start_cmd = _get_start_cmd(jr, dep_jobid, profile)
    full_cmd = "$start_cmd $cmd"
    
    jobid = ""
    if _cmd_run_shell(jr, profile)
        # Set MPI environment variables to avoid network interface errors
        # Completely disable OFI and use only shared memory transport
        env_vars = Dict(
            "OMPI_MCA_btl" => "^openib,usnic,ofi,tcp",
            "OMPI_MCA_btl_vader_single_copy_mechanism" => "none",
            "OMPI_MCA_pml" => "ob1",
            "OMPI_MCA_ofi_interface" => "",
            "OMPI_MCA_btl_tcp_if_exclude" => "bridge101,lo,docker0"
        )
        result = run(setenv(`sh -c $full_cmd`, env_vars), wait = true)
        @info "call: $full_cmd"
        if result.exitcode == 0
            jobout = read(result, String)
            jobid = _get_jobid(jr, jobout, profile)
            @info "JobID is: $jobid"
        else
            @error "Error submitting Job: $full_cmd"
        end
    elseif _cmd_run_stdout(jr, profile)
        println(full_cmd)
    else
        cmd_run_file = _cmd_run_file(jr, profile)
        if !isempty(cmd_run_file)
            if !isfile(cmd_run_file)
                open(cmd_run_file, "w") do f
                    println(f, "#! /bin/bash")
                end
            end
            open(cmd_run_file, "a") do f
                println(f, full_cmd)
            end
        end
    end
    return jobid
end

mutable struct CfgJobRunner <: AbstractJobRunner
    input_cfg::String
    cmd_run_file::Union{String, Nothing}
    _config::Dict{String, Any}
    
    function CfgJobRunner(input_cfg::String, cmd_run_file::Union{String, Nothing} = nothing)
        config = Dict{String, Any}()
        if isfile(input_cfg)
            # Pre-process config file to fix TOML-incompatible floating point numbers
            # TOML requires digits after decimal point (0. -> 0.0)
            # This maintains compatibility with existing config files
            normalized_cfg = _normalize_config_file(input_cfg)
            try
                config = TOML.parsefile(normalized_cfg)
            finally
                # Clean up temporary file if we created one
                if normalized_cfg != input_cfg
                    rm(normalized_cfg, force=true)
                end
            end
        end
        new(input_cfg, cmd_run_file, config)
    end
end

# Helper function to normalize config file for TOML parsing
# TOML is stricter than Python's configparser, so we need to:
# 1. Fix floating point numbers: 0. -> 0.0
# 2. Quote unquoted string values
function _normalize_config_file(cfg_path::String)::String
    # Read the file line by line
    lines = readlines(cfg_path)
    normalized_lines = String[]
    needs_normalization = false
    
    for line in lines
        original_line = line
        # Skip comments, blank lines, and section headers
        stripped = strip(line)
        if isempty(stripped) || startswith(stripped, "#") || startswith(stripped, "[")
            push!(normalized_lines, line)
            continue
        end
        
        # Check for key = value pattern
        if occursin(r"^\s*[^=]+\s*=\s*", line)
            # Split on = (but only the first one)
            parts = split(line, "=", limit=2)
            if length(parts) == 2
                key = parts[1]
                value = strip(parts[2])
                
                # Remove inline comments from value
                comment_pos = findfirst("#", value)
                if comment_pos !== nothing
                    value = strip(value[1:prevind(value, first(comment_pos))])
                    comment = value[first(comment_pos):end]
                else
                    comment = ""
                end
                
                # Normalize floating point: X. -> X.0
                if match(r"^\d+\.\s*$", value) !== nothing
                    value = replace(value, r"^(\d+)\.\s*$" => s"\1.0")
                    needs_normalization = true
                end
                
                # Quote unquoted string values (not already quoted, not numbers, not booleans)
                if !isempty(value) && 
                   !startswith(value, "\"") && !startswith(value, "'") &&
                   !occursin(r"^[\d\.\+\-eE]+$", value) &&  # Not a number
                   !occursin(r"^(true|false)$", lowercase(value)) &&  # Not a boolean
                   !startswith(value, "[") && !startswith(value, "{")  # Not array/table
                    value = "\"$value\""
                    needs_normalization = true
                end
                
                line = key * " = " * value * (isempty(comment) ? "" : " " * comment)
            end
        end
        
        push!(normalized_lines, line)
    end
    
    # If no normalization needed, return original path
    if !needs_normalization
        return cfg_path
    end
    
    # Write to temporary file
    temp_file = tempname() * ".cfg"
    write(temp_file, join(normalized_lines, "\n"))
    return temp_file
end

function _get_start_cmd(jr::CfgJobRunner, dep_jobid::String, profile::String)::String
    cmd_job_dependency = get_property(jr, profile, "cmd_job_dependency", "")
    cmd_dep = isempty(dep_jobid) ? "" : replace(cmd_job_dependency, "{jobid}" => dep_jobid)
    cmd_prefix = get_property(jr, profile, "cmd_prefix", "")
    return replace(cmd_prefix, "{cmd_job_dependency}" => cmd_dep)
end

function _get_jobid(jr::CfgJobRunner, jobout::String, profile::String)::String
    jobid_re = get_property(jr, profile, "jobid_re", "")
    m = match(Regex(jobid_re), jobout)
    return m === nothing ? "" : m.match
end

function get_sim_property(jr::CfgJobRunner, name::String)::String
    section = get(jr._config, CFG_SECTION_3DPOLYS_LE, Dict{String, Any}())
    return get(section, name, "")
end

function get_property(jr::CfgJobRunner, profile::String, name::String, default::String = "")::String
    value = default
    if !isempty(profile)
        section_name = "job_runner_$profile"
        section = get(jr._config, section_name, Dict{String, Any}())
        value = get(section, name, default)
    end
    if isempty(value)
        section = get(jr._config, "job_runner", Dict{String, Any}())
        value = get(section, name, default)
    end
    return value
end

function _cmd_run_shell(jr::CfgJobRunner, profile::String)::Bool
    cmd_run = get_property(jr, profile, "cmd_run", "")
    return cmd_run == "shell" && (jr.cmd_run_file === nothing || isempty(jr.cmd_run_file))
end

function _cmd_run_stdout(jr::CfgJobRunner, profile::String)::Bool
    cmd_run = get_property(jr, profile, "cmd_run", "")
    return cmd_run == "stdout" && (jr.cmd_run_file === nothing || isempty(jr.cmd_run_file))
end

function _cmd_run_file(jr::CfgJobRunner, profile::String)::String
    if jr.cmd_run_file !== nothing && !isempty(jr.cmd_run_file)
        return jr.cmd_run_file
    end
    cmd_run = get_property(jr, profile, "cmd_run", "")
    if startswith(cmd_run, "file:")
        return split(cmd_run, ":")[2]
    end
    return ""
end

function convert_dat_to_cfg(input_dat::String)::String
    config = Dict{String, Any}()
    config[CFG_SECTION_3DPOLYS_LE] = Dict{String, Any}()
    
    open(input_dat, "r") do fp
        for line in eachline(fp)
            if isempty(strip(line)) || startswith(strip(line), "#")
                continue
            end
            parts = split(line, "::")
            if length(parts) == 2
                val = strip(parts[1])
                name = strip(parts[2])
                config[CFG_SECTION_3DPOLYS_LE][name] = val
            end
        end
    end
    
    input_cfg = replace(input_dat, ".dat" => ".cfg")
    open(input_cfg, "w") do cf
        TOML.print(cf, config)
    end
    return input_cfg
end

end # module

