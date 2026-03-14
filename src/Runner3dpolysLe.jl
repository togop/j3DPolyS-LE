"""
3DPolyS-LE Runner module - Main entry point
"""
module Runner3dpolysLe

import Base.CoreLogging: @info, @warn, @error, @debug
using ArgParse
using DataFrames
using CSV
using Glob
using TOML
using ..HicAnalysis
using ..JobRunner
import ..JobRunner: CfgJobRunner, get_property, run_cmd, CFG_SECTION_3DPOLYS_LE, _normalize_config_file  # Import directly to avoid documentation issues
using ..Stats3dpolysLe
using ..PlotHic
using ..PlotSimStats

export main

const EXP_COOL_AS_STATS = "."

const logger = Base.CoreLogging.current_logger()

mutable struct DccExtrusionArgs
    input_cfg::String
    output_folder::String
    analyse::String
    boundary::String
    lef_loading_sites::String
    tads_boundary::String
    stats_file::String
    exp_cool::String
    exp_chip::String
    nlef::Int
    km::Float64
    radius_contact::Float64
    contact_probability::Union{Bool, Nothing}
    boundary_direction::Union{Int, Nothing}
    z_loop::Union{Bool, Nothing}
    unidirectional::Union{Bool, Nothing}
    init_mode::String
    stats::Bool
    all_stats::Bool
    cmp_chrs::Union{Vector{String}, Nothing}
    resolution::Int
    _config::Dict{String, Any}
    
    function DccExtrusionArgs(;
        stats::Bool = false,
        all_stats::Bool = false,
        boundary::String = "",
        lef_loading_sites::String = "",
        input_cfg::String = "./input.cfg",
        tads_boundary::String = "",
        stats_file::String = "./py3dpolys_le_stats.csv",
        exp_cool::String = "",
        exp_chip::String = "",
        nlef::Int = 0,
        km::Float64 = 0.0,
        radius_contact::Float64 = 0.0,
        contact_probability::Union{Bool, Nothing} = nothing,
        boundary_direction::Union{Int, Nothing} = nothing,
        z_loop::Union{Bool, Nothing} = nothing,
        unidirectional::Union{Bool, Nothing} = nothing,
        init_mode::String = "",
        output_folder::String = "",
        analyse::String = "",
        cmp_chrs::Union{Vector{String}, Nothing} = nothing,
        resolution::Int = HicAnalysis.RESOLUTION
    )
        config = Dict{String, Any}()
        if isfile(input_cfg)
            # Use the same normalization function as CfgJobRunner
            normalized_cfg = _normalize_config_file(input_cfg)
            try
                config = TOML.parsefile(normalized_cfg)
            finally
                # Clean up temporary file if we created one
                if normalized_cfg != input_cfg
                    rm(normalized_cfg, force=true)
                end
            end
        else
            @warn "Input configuration file $input_cfg does not exist or is not accessible, all parameters will be read from the CLI!"
        end
        
        output_folder_val = output_folder
        analyse_val = analyse
        
        boundary_val = !isempty(boundary) || !isfile(input_cfg) ? boundary : get_property_from_config(config, "boundary", "")
        if !isempty(boundary_val) && !isfile(boundary_val)
            @error "Boundary file $boundary_val not found!"
            exit(1)
        end
        
        lef_loading_sites_val = !isempty(lef_loading_sites) || !isfile(input_cfg) ? lef_loading_sites : get_property_from_config(config, "lef_loading_sites", "")
        if !isempty(lef_loading_sites_val) && !isfile(lef_loading_sites_val)
            @error "LEF loading-sites file $lef_loading_sites_val not found!"
            exit(1)
        end
        
        tads_boundary_val = !isempty(tads_boundary) || !isfile(input_cfg) ? tads_boundary : get_property_from_config(config, "tads_boundary", "")
        if !isempty(tads_boundary_val) && !isfile(tads_boundary_val)
            @error "TADs-boundary file $tads_boundary_val not found!"
            exit(1)
        end
        
        exp_cool_val = !isempty(exp_cool) || !isfile(input_cfg) ? exp_cool : get_property_from_config(config, "exp_cool", "")
        exp_chip_val = exp_chip
        
        nlef_val = nlef != 0 || !isfile(input_cfg) ? nlef : parse(Int, get_property_from_config(config, "Nlef", "0"))
        km_val = km != 0.0 || !isfile(input_cfg) ? km : parse(Float64, get_property_from_config(config, "km", "0.0"))
        radius_contact_val = radius_contact != 0.0 || !isfile(input_cfg) ? radius_contact : parse(Float64, get_property_from_config(config, "radius_contact", "0.0"))
        
        contact_probability_val = contact_probability !== nothing || !isfile(input_cfg) ? contact_probability : str2bool(get_property_from_config(config, "contact_probability", "false"))
        boundary_direction_val = boundary_direction !== nothing || !isfile(input_cfg) ? boundary_direction : parse(Int, get_property_from_config(config, "boundary_direction", "0"))
        
        z_loop_val = z_loop !== nothing || !isfile(input_cfg) ? z_loop : str2bool(get_property_from_config(config, "z_loop", "false"))
        unidirectional_val = unidirectional !== nothing || !isfile(input_cfg) ? unidirectional : str2bool(get_property_from_config(config, "unidirectional", "false"))
        init_mode_val = !isempty(init_mode) || !isfile(input_cfg) ? init_mode : get_property_from_config(config, "init_mode", "")
        
        if cmp_chrs !== nothing
            cmp_chrs_val = cmp_chrs
        elseif isfile(input_cfg)
            cmp_chrs_str = get_property_from_config(config, "cmp_chrs", "")
            cmp_chrs_val = !isempty(cmp_chrs_str) ? [String(x) for x in split(cmp_chrs_str, r"[\s;,]+")] : nothing
        else
            cmp_chrs_val = nothing
        end
        
        if cmp_chrs_val === nothing || isempty(cmp_chrs_val)
            cmp_chrs_val = HicAnalysis.CHR_SYNONYMS
        else
            global HicAnalysis.CHR_SYNONYMS = cmp_chrs_val
        end
        
        if !isempty(exp_cool_val)
            if !isfile(exp_cool_val)
                @error "Experimental cool file $exp_cool_val not found!"
                exit(1)
            else
                @info "Ensure .mcool file of $exp_cool_val is created..."
                cmp_hic = replace(exp_cool_val, r"\.cool$" => "")
                HicAnalysis.get_exp_sim_mcool(cmp_hic, cmp_chrs_val, resolution)
            end
        end
        
        new(input_cfg, output_folder_val, analyse_val, boundary_val, lef_loading_sites_val,
            tads_boundary_val, stats_file, exp_cool_val, exp_chip_val, nlef_val, km_val,
            radius_contact_val, contact_probability_val, boundary_direction_val, z_loop_val,
            unidirectional_val, init_mode_val, stats, all_stats, cmp_chrs_val, resolution, config)
    end
end

function get_property_from_config(config::Dict{String, Any}, name::String, default::String = "")::String
    section = get(config, CFG_SECTION_3DPOLYS_LE, Dict{String, Any}())
    value = get(section, name, default)
    # Convert to String if it's not already a String
    if value isa String
        return value
    elseif value === nothing
        return default
    else
        return string(value)
    end
end

function str2bool(b::String)::Bool
    return lowercase(b) in ["true", "t", "yes", "y", "1"]
end

function default_analysis_subfolder(args::DccExtrusionArgs)::String
    cp = args.contact_probability === true ? "p" : ""
    return args.radius_contact > 0 ? "r$(round(args.radius_contact, digits=2))$cp" : ""
end

function default_analysis_folder(args::DccExtrusionArgs)::String
    return joinpath(args.output_folder, default_analysis_subfolder(args))
end

function default_output_folder(args::DccExtrusionArgs)::String
    z_opt = args.z_loop === true ? "_z-loop" : ""
    u_opt = args.unidirectional === true ? "_unidir" : ""
    im_opt = !isempty(args.init_mode) ? "_im-$(args.init_mode[1])" : ""
    return "out-Nlef$(args.nlef)-km$(args.km)-bd$(args.boundary_direction)$(im_opt)$(z_opt)$(u_opt)"
end

function is_analysis_output_folder(folder::String)::Bool
    chip_out_file = joinpath(folder, HicAnalysis.CHIP_OUT)
    xyzconfig_out = joinpath(folder, HicAnalysis.XYZCONFIG_OUT)
    return isdir(folder) && isfile(chip_out_file) && isfile(xyzconfig_out)
end

function is_simulation_output_folder(folder::String)::Bool
    return isdir(folder) && !isempty(readdir(folder)) &&
           isfile(joinpath(folder, HicAnalysis.CONFIG_OUT)) &&
           isfile(joinpath(folder, HicAnalysis.CONTACT_OUT)) &&
           isfile(joinpath(folder, HicAnalysis.NLEF_OUT)) &&
           isfile(joinpath(folder, HicAnalysis.PROCESS_OUT)) &&
           isfile(joinpath(folder, HicAnalysis.DR_OUT))
end

mutable struct DccExtrusionRunner
    _running_jobids::Vector{String}
    _job_runner::CfgJobRunner  # Use imported name directly
    
    function DccExtrusionRunner(job_runner::CfgJobRunner)
        new(String[], job_runner)
    end
end

function read_stats_file(stats_file::String)::DataFrame
    return CSV.read(stats_file, DataFrame, types=Dict(
        :boundary => String,
        :boundary_direction => String,
        :km => String,
        :radius_contact => String,
        :chi2_log => String,
        :alpha_log => String,
        :chi2_lin => String,
        :alpha_lin => String
    ), missingstring="")
end

function split_radius_contact_probability(rcp::String)::Tuple{Bool, Float64}
    if isdigit(last(rcp))
        r = parse(Float64, rcp)
        cp = false
    else
        r = parse(Float64, rcp[1:end-1])
        cp = true
    end
    return cp, r
end

function run(runner::DccExtrusionRunner, dcc_args::DccExtrusionArgs;
             stats_only::Bool = false, radii::Vector{String} = String[],
             dep_jobid::Union{String, Nothing} = nothing, replace::Bool = false)::String
    prev_jobid = dep_jobid !== nothing ? dep_jobid : (isempty(runner._running_jobids) ? "" : runner._running_jobids[end])
    
    cp = dcc_args.contact_probability === true ? "-cp" : ""
    z_loop = dcc_args.z_loop === true ? "-z" : ""
    u_opt = dcc_args.unidirectional === true ? "-u" : ""
    init_mode = !isempty(dcc_args.init_mode) ? "-im:$(dcc_args.init_mode)" : ""
    r = dcc_args.radius_contact
    r_opt = r > 0 ? "-r:$r $cp" : ""
    r_opt_py = r > 0 ? "-r $r $cp" : ""
    
    boundary = dcc_args.boundary
    input_cfg = dcc_args.input_cfg
    stats_file = dcc_args.stats_file
    exp_cool = dcc_args.exp_cool
    nlef = dcc_args.nlef
    km = dcc_args.km
    a_opt = !isempty(dcc_args.analyse) ? "-a:$(dcc_args.analyse)" : ""
    
    if isempty(dcc_args.output_folder)
        dcc_args.output_folder = default_output_folder(dcc_args)
    end
    
    jobid = ""
    if !dcc_args.stats && !dcc_args.all_stats && !stats_only
        run_sim_or_analysis = true
        if isempty(dcc_args.analyse)
            if is_simulation_output_folder(dcc_args.output_folder)
                if replace
                    @warn "Simulation output folder $(dcc_args.output_folder) already exists and data files will be replaced!"
                else
                    @warn "Simulation output folder $(dcc_args.output_folder) already exists so simulation will be skipped!"
                    run_sim_or_analysis = false
                end
            end
        end
        
        if !isempty(dcc_args.analyse)
            if is_analysis_output_folder(dcc_args.analyse)
                if replace
                    @warn "Analysis output folder $(dcc_args.analyse) exists and data files will be replaced!"
                else
                    @warn "Analysis output folder $(dcc_args.analyse) exists so analysis will be skipped!"
                    run_sim_or_analysis = false
                end
            end
        end
        
        if run_sim_or_analysis
            cmd_prefix = get_cmd_prefix(runner, dcc_args, mpirun=true)
            boundary_opt_f = !isempty(boundary) ? "-b:$boundary" : ""
            lef_loading_sites_opt_f = !isempty(dcc_args.lef_loading_sites) ? "-lls:$(dcc_args.lef_loading_sites)" : ""
            boundary_direction_opt_f = dcc_args.boundary_direction !== nothing ? "-bd:$(dcc_args.boundary_direction)" : ""
            cmd = "$cmd_prefix 3dpolys_le -o:$(dcc_args.output_folder) --km:$km --nlef:$nlef $boundary_opt_f $lef_loading_sites_opt_f $boundary_direction_opt_f $z_loop $u_opt $init_mode $a_opt $r_opt $input_cfg"
            
            if !isempty(dcc_args.analyse)
                jobid = run_cmd(runner._job_runner, cmd, prev_jobid, "analysis")
            else
                jobid = run_cmd(runner._job_runner, cmd, prev_jobid, "sim")
            end
        end
    end
    
    if dep_jobid === nothing && !isempty(jobid)
        push!(runner._running_jobids, jobid)
    end
    
    # Contact radius analysis
    stats_dep_jobid = jobid
    if isempty(dcc_args.analyse) && !stats_only && !isempty(radii)
        # Filter out radii that match the simulation radius (no need to re-analyze)
        sim_radius_str = string(dcc_args.radius_contact)
        additional_radii = filter(r -> !startswith(r, sim_radius_str), radii)

        # Only run additional radius analysis if batch mode is configured
        # In shell mode, the Fortran program has issues with analysis-only mode
        cmd_run = get_property(runner._job_runner, "", "cmd_run", "")
        is_batch_mode = !isempty(cmd_run) && cmd_run != "shell" && cmd_run != "stdout"

        if !isempty(additional_radii) && is_batch_mode
            dcc_args_analysis = deepcopy(dcc_args)
            stats_dep_jobid = sim_contact_radius_analysis(runner, dcc_args_analysis, radii=additional_radii, dep_jobid=jobid)
        elseif !isempty(additional_radii) && !is_batch_mode
            @warn "Skipping additional radius analysis ($(join(additional_radii, ", "))) in shell mode. " *
                  "Set cmd_run to batch mode (e.g., 'sbatch') in config to enable multi-radius analysis."
        end
    end
    
    # Statistics
    analyse_folder = dcc_args.analyse
    hic_mcool = filter(f -> occursin(r"^hic_.*\.mcool$", f), isdir(analyse_folder) ? readdir(analyse_folder) : String[])
    if !isempty(analyse_folder) || (isfile(joinpath(analyse_folder, HicAnalysis.CHIP_OUT)) && (isempty(hic_mcool) || dcc_args.all_stats || stats_only))
        s_cmp_chrs = dcc_args.cmp_chrs !== nothing ? "--cmp_chrs $(join(dcc_args.cmp_chrs, " "))" : ""
        t_opt = !isempty(dcc_args.tads_boundary) ? "-t $(dcc_args.tads_boundary)" : ""
        boundary_opt_py = !isempty(boundary) ? "-b $boundary" : ""
        boundary_direction_opt_py = dcc_args.boundary_direction !== nothing ? "-bd $(dcc_args.boundary_direction)" : ""
        cmd_prefix = get_cmd_prefix(runner, dcc_args, mpirun=false)
        cmd = "$cmd_prefix 3dpolys_le_stats -o $(dcc_args.output_folder) -a $analyse_folder --km $km --nlef $nlef -e $exp_cool $boundary_opt_py $boundary_direction_opt_py $t_opt $r_opt_py -i $input_cfg -f $stats_file $s_cmp_chrs"
        jobid = run_cmd(runner._job_runner, cmd, stats_dep_jobid, "stats")
    end
    
    return jobid
end

# Check if mpirun is available and get its full path
function find_mpirun()::String
    try
        # Try to find mpirun in PATH
        result = read(`which mpirun`, String)
        mpirun_path = strip(result)
        if !isempty(mpirun_path) && isfile(mpirun_path)
            return mpirun_path
        end
    catch
        # If which fails, try common locations
        common_paths = [
            "/usr/bin/mpirun",
            "/usr/local/bin/mpirun",
            "/opt/homebrew/bin/mpirun",
            "/opt/local/bin/mpirun"
        ]
        for path in common_paths
            if isfile(path)
                return path
            end
        end
    end
    return ""  # MPI not available
end

# Cache the mpirun path to avoid repeated lookups
const MPIRUN_PATH = find_mpirun()

function get_cmd_prefix(runner::DccExtrusionRunner, dcc_args::DccExtrusionArgs; mpirun::Bool = false)::String
    container_prefix = get_property(runner._job_runner, "", "container_prefix", "")

    # Only create cmd.sh wrapper if using a container
    if !isempty(container_prefix)
        cmd_sh = joinpath(dcc_args.output_folder, "cmd.sh")
        if !isfile(cmd_sh)
            if !isdir(dcc_args.output_folder)
                mkdir(dcc_args.output_folder)
            end
            open(cmd_sh, "w") do f
                println(f, "#! /bin/bash")
                # Set MPI environment variables to avoid network interface errors
                println(f, "export OMPI_MCA_btl_vader_single_copy_mechanism=none")
                println(f, "export OMPI_MCA_btl=^openib,usnic,ofi")
                println(f, "export OMPI_MCA_btl_tcp_if_exclude=bridge101,lo,docker0")
                println(f, "export OMPI_MCA_pml=ob1")
                println(f, "\"\$@\"")  # Escape $ to prevent string interpolation
            end
        end

        # Build command prefix with MPI flags if using mpirun
        if mpirun
            # In container mode, mpirun should be in the container's PATH
            mpi_flags = "-mca btl ^openib,usnic,ofi,tcp -mca btl_vader_single_copy_mechanism none -mca pml ob1 -mca ofi_interface \"\""
            cmd_prefix = "$cmd_sh mpirun $mpi_flags $container_prefix"
        else
            cmd_prefix = "$cmd_sh $container_prefix"
        end
    else
        # No container - execute directly
        if mpirun
            # Check if mpirun is available
            if isempty(MPIRUN_PATH)
                @warn "mpirun requested but not found in PATH. Running without MPI..."
                cmd_prefix = ""
            else
                mpi_flags = "-mca btl ^openib,usnic,ofi,tcp -mca btl_vader_single_copy_mechanism none -mca pml ob1 -mca ofi_interface \"\""
                cmd_prefix = "$MPIRUN_PATH $mpi_flags"
            end
        else
            cmd_prefix = ""
        end
    end

    return cmd_prefix
end

function sim_contact_radius_analysis(runner::DccExtrusionRunner, dcc_args::DccExtrusionArgs;
                                      radii::Vector{String} = String[], dep_jobid::Union{String, Nothing} = nothing)::String
    last_jobid = ""
    for rcp in radii
        dcc_args_r = deepcopy(dcc_args)
        cp, r = split_radius_contact_probability(rcp)
        dcc_args_r.radius_contact = r
        dcc_args_r.contact_probability = cp
        dcc_args_r.analyse = default_analysis_folder(dcc_args_r)

        # Create analysis folder if it doesn't exist
        # The Fortran program requires the folder to exist before running
        if !isdir(dcc_args_r.analyse)
            mkpath(dcc_args_r.analyse)
        end

        last_jobid = run(runner, dcc_args_r, stats_only=false, radii=radii, dep_jobid=dep_jobid)
        last_jobid = run_multi_decay_plot(runner, dcc_args, last_jobid)
    end
    return last_jobid
end

function run_multi_decay_plot(runner::DccExtrusionRunner, dcc_args::DccExtrusionArgs, dep_jobid::String)::String
    s_cmp_chrs = dcc_args.cmp_chrs !== nothing ? "--cmp_chrs $(join(dcc_args.cmp_chrs, " "))" : ""
    cmd_prefix = get_cmd_prefix(runner, dcc_args, mpirun=false)
    cmd = "$cmd_prefix 3dpolys_le_runner multi_decay_plot -o $(dcc_args.output_folder) -i $(dcc_args.input_cfg) -e $(dcc_args.exp_cool) $s_cmp_chrs"
    return run_cmd(runner._job_runner, cmd, dep_jobid, "")
end

function multi_decay_plot(dcc_args::DccExtrusionArgs; replace::Bool = false)
    @info "multi_decay_plot called with output_folder=$(dcc_args.output_folder), analyse=$(dcc_args.analyse), exp_cool=$(dcc_args.exp_cool)"
    sub_folders = !isempty(dcc_args.analyse) ? [dcc_args.analyse] : sort([f.path for f in readdir(dcc_args.output_folder, join=true) if isdir(f)])
    plots_folder = joinpath(dcc_args.output_folder, HicAnalysis.PLOTS_FOLDER)
    if plots_folder in sub_folders
        filter!(x -> x != plots_folder, sub_folders)
    end
    @info "sub_folders: $sub_folders"
    hic_h5_list = String[]
    for analysis_folder in sub_folders
        last_hic = HicAnalysis.get_last_hic(analysis_folder)
        @info "Found HIC file: $last_hic"
        push!(hic_h5_list, last_hic)
    end
    @info "hic_h5_list: $hic_h5_list"
    
    # Call plot_distance_contact_prob_decay for both log and linear modes
    @info "Calling plot_distance_contact_prob_decay for log mode..."
    hic_multi_decay_plot_log = HicAnalysis.plot_distance_contact_prob_decay(hic_h5_list,
                                                                              exp_cool=dcc_args.exp_cool,
                                                                              output_folder=dcc_args.output_folder,
                                                                              res=dcc_args.resolution,
                                                                              confidence=0.0,
                                                                              replace=replace,
                                                                              chi2_mode=HicAnalysis.CHI2_MODE_LOG)
    @info "Log mode plot saved to: $hic_multi_decay_plot_log"
    
    @info "Calling plot_distance_contact_prob_decay for linear mode..."
    hic_multi_decay_plot_lin = HicAnalysis.plot_distance_contact_prob_decay(hic_h5_list,
                                                                             exp_cool=dcc_args.exp_cool,
                                                                             output_folder=dcc_args.output_folder,
                                                                             res=dcc_args.resolution,
                                                                             confidence=0.0,
                                                                             replace=replace,
                                                                             chi2_mode=HicAnalysis.CHI2_MODE_LINEAR)
    @info "Linear mode plot saved to: $hic_multi_decay_plot_lin"
end

function grid_nlef_km(runner::DccExtrusionRunner, dcc_args::DccExtrusionArgs,
                      nlef_list::Vector{Int}, km_list::Vector{Float64},
                      radii::Vector{String}, replace::Bool)
    for nlef in nlef_list
        for km in km_list
            dcc_args_grid = deepcopy(dcc_args)
            dcc_args_grid.nlef = nlef
            dcc_args_grid.km = km
            run(runner, dcc_args_grid, radii=radii, replace=replace)
        end
    end
end

function cli_parser()
    s = ArgParseSettings(description="Running 3dpolys_le_runner.", epilog="DISCLAIMER: As almost not sufficient tests prove the correctness of all possible parameter combinations (no comprehensive test coverage), please check your output data and log files, and make sure all went as you have expected.")
    
    @add_arg_table! s begin
        "run_command"
            help = "Run batch command"
            arg_type = String
            required = true
        "--boundary_direction", "-d"
            help = "Impermeability direction applied to all boundaries"
            arg_type = Int
            default = 0
        "--boundary", "-b"
            help = "Boundary sites file in CSV format"
            default = ""
        "--lef_loading_sites", "-s"
            help = "LEFs loading sites file in a csv format"
            default = ""
        "--hic_chrs"
            help = "Synonyms of the chromosome from Hi-C matrixes to be compared"
            nargs = '*'
        "--cmp_chrs"
            help = "Synonyms of the Hi-C chromosome with which simulation to be compared"
            nargs = '*'
        "--z_loop", "-z"
            help = "Allow z_loop for LEFs move in a simulation"
            action = :store_true
        "--unidirectional", "-u"
            help = "Unidirectional mode for LEFs move otherwise bidirectional"
            action = :store_true
        "--init_mode", "-n"
            help = "Initial folding mode"
            default = ""
        "--tads_boundary", "-t"
            help = "TADs boundary file"
            default = ""
        "--exp_cool", "-e"
            help = "Experimental cooler (.cool) file"
            default = ""
        "--exp_cools", "-x"
            help = "List of experimental cooler files"
            nargs = '+'
        "--exp_chip", "-p"
            help = "Experimental ChIP-seq file"
        "--resolution", "-q"
            help = "Resolution for distance-contact-decay plots"
            arg_type = Int
            default = HicAnalysis.RESOLUTION
        "--input_cfg", "-i"
            help = "Input.cfg file used from a simulation"
            default = "./input.cfg"
        "--nlef", "-l"
            help = "Nlef value to use in a simulation"
            arg_type = Int
        "--km", "-m"
            help = "km value to use in a simulation"
            arg_type = Float64
        "--stats_file", "-f"
            help = "Simulation statistics' repository file"
            default = "./sim_stats.csv"
        "--stats"
            help = "Run statistical analysis only for entries missing statistics plots"
            action = :store_true
        "--all_stats"
            help = "Run statistical analysis for all entries"
            action = :store_true
        "--radius_contact", "-r"
            help = "Contact radius in lattice units"
            arg_type = Float64
            default = 0.0
        "--list_contact_radii", "-lr"
            help = "List of contact radii"
            nargs = '+'
            default = ["2.84"]
        "--replace"
            help = "Replace of any existing output data files"
            action = :store_true
        "--threading"
            help = "Use Python multi-threading"
            action = :store_true
        "--contact_probability", "-y"
            help = "Use contact radius probability"
            action = :store_true
        "--correlation"
            help = "Correlation method to use"
            default = "spearmanr"
        "--output_folder", "-o"
            help = "Simulation output folder"
            default = ""
        "--cmd_run_file"
            help = "File where to save all commands"
            default = ""
        "--analysis_folder", "-a"
            help = "Analysis output folder"
            default = ""
        "--nlef_list"
            help = "List of Nlef values"
            nargs = '+'
            arg_type = Int
            default = [50, 100, 150, 200, 300, 400, 500, 600, 800, 1000, 1200]
        "--km_list"
            help = "List of km values"
            nargs = '+'
            arg_type = Float64
            default = [5.4e-4, 3*5.4e-4, 5*5.4e-4, 7*5.4e-4, 9*5.4e-4, 12*5.4e-4]
    end
    return s
end

function main(cmd_args::Union{Vector{String}, Nothing} = nothing)
    # Allow passing arguments directly (useful for REPL usage)
    # If cmd_args is provided, use it; otherwise use global ARGS
    settings = cli_parser()
    if cmd_args !== nothing
        args = parse_args(cmd_args, settings)
    else
        args = parse_args(settings)
    end
    
    if haskey(args, "cmp_chrs") && args["cmp_chrs"] !== nothing
        cmp_chrs = [String(x) for x in args["cmp_chrs"]]
        @info "Overwriting the default hic_analysis.CHR_SYNONYMS: $(join(cmp_chrs, ",")) with which HiCs will be compared!"
        global HicAnalysis.CHR_SYNONYMS = cmp_chrs
    end
    
    job_runner = CfgJobRunner(args["input_cfg"], args["cmd_run_file"])
    
    dcc_args = DccExtrusionArgs(
        boundary=args["boundary"],
        lef_loading_sites=args["lef_loading_sites"],
        input_cfg=args["input_cfg"],
        tads_boundary=args["tads_boundary"],
        stats_file=args["stats_file"],
        exp_cool=args["exp_cool"],
        nlef=(haskey(args, "nlef") && args["nlef"] !== nothing) ? Int(args["nlef"]) : 0,
        km=(haskey(args, "km") && args["km"] !== nothing) ? Float64(args["km"]) : 0.0,
        radius_contact=args["radius_contact"],
        contact_probability=get(args, "contact_probability", false) ? true : nothing,
        boundary_direction=args["boundary_direction"],
        z_loop=get(args, "z_loop", false) ? true : nothing,
        unidirectional=get(args, "unidirectional", false) ? true : nothing,
        init_mode=args["init_mode"],
        stats=get(args, "stats", false),
        all_stats=get(args, "all_stats", false),
        cmp_chrs=(haskey(args, "cmp_chrs") && args["cmp_chrs"] !== nothing && !isempty(args["cmp_chrs"])) ? Vector{String}([String(x) for x in args["cmp_chrs"]]) : nothing,
        resolution=args["resolution"],
        output_folder=args["output_folder"],
        analyse=args["analysis_folder"]
    )
    
    dcc_run = DccExtrusionRunner(job_runner)
    
    plot_format = get_property(job_runner, "stats", "plot_format", HicAnalysis.PLOT_FORMAT)
    if isempty(plot_format)
        plot_format = HicAnalysis.PLOT_FORMAT
    end
    
    run_command = args["run_command"]
    # Convert vector arguments from Vector{Any} to proper types
    nlef_list = [Int(x) for x in args["nlef_list"]]
    km_list = [Float64(x) for x in args["km_list"]]
    list_contact_radii = [String(x) for x in args["list_contact_radii"]]
    
    if run_command == "grid_nlef_km"
        grid_nlef_km(dcc_run, dcc_args, nlef_list, km_list, list_contact_radii, args["replace"])
    elseif run_command == "run"
        run(dcc_run, dcc_args, radii=list_contact_radii, replace=args["replace"])
    elseif run_command == "multi_decay_plot"
        multi_decay_plot(dcc_args, replace=args["replace"])
    # Additional commands would be implemented here
    end
end

end # module

