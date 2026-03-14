"""
3DPolyS-LE Statistics module
"""
module Stats3dpolysLe

import Base.CoreLogging: @info, @warn, @error, @debug
using CSV
using DataFrames
using ArgParse
using ..HicAnalysis
using ..PlotHic
using ..JobRunner
import ..JobRunner: CfgJobRunner, get_sim_property, get_property, convert_dat_to_cfg  # Import directly to avoid documentation issues

const logger = Base.CoreLogging.current_logger()

function cli_parser()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--output_folder", "-o"
            help = "Simulation's output folder containing result files"
            default = "."
        "--analyse", "-a"
            help = "'Analyse' step output folder containing hic_*.hdf5 files"
            default = "./out/analyse"
        "--boundary", "-b"
            help = "Boundary file in TSV format"
            arg_type = String
        "--boundary_direction", "-bd"
            help = "Impermeability direction applied to all boundaries: -1: opposite, 0: both, 1: same direction"
            arg_type = Int
            default = 0
        "--tads_boundary", "-t"
            help = "TADs boundary file in CSV format"
            arg_type = String
        "--exp_cool", "-e"
            help = "Experimental cool file with which simulation data to be compared"
            default = ""
        "--nlef", "-l"
            help = "Nlef value used in a simulation"
            default = "200"
        "--km", "-m"
            help = "km value used in a simulation"
            default = "2.7e-3"
        "--input_cfg", "-i"
            help = "input.cfg file used in a simulation"
            default = "./input.cfg"
        "--radius_contact", "-r"
            help = "Contact radius in lattice units"
            arg_type = Float64
            default = 0.0
        "--contact_probability", "-cp"
            help = "Use contact radius probability: (1 - r^2 / max_r^2)"
            action = :store_true
        "--stats_file", "-f"
            help = "Simulation statistics' repository file"
            default = "./sim_stats.csv"
        "--resolution", "-res"
            help = "Resolution to downscale Chip-seq output data in .bedGraph format (default 2kb)"
            arg_type = Int
            default = HicAnalysis.SIM_RESOLUTION
        "--replace"
            help = "Whether to replace existing files"
            action = :store_true
        "--chi2_mode"
            help = "Chi2-min mode for sampling contact distances"
            default = HicAnalysis.CHI2_MODE_LINEAR
        "--cmp_chrs"
            help = "Synonyms of the Hi-C chromosome with which simulation data to be compared"
            nargs = '*'
    end
    return s
end

function main()
    logging_config = Dict(
        "version" => 1,
        "disable_existing_loggers" => false,
        "formatters" => Dict(
            "default" => Dict(
                "format" => "%(levelname)s: %(message)s"
            )
        ),
        "handlers" => Dict(
            "console" => Dict(
                "class" => "logging.StreamHandler",
                "level" => "INFO",
                "formatter" => "default"
            )
        ),
        "root" => Dict(
            "level" => "INFO",
            "handlers" => ["console"]
        )
    )
    
    args = parse_args(cli_parser())
    
    input_cfg = args["input_cfg"]
    if endswith(input_cfg, ".dat")
        input_cfg = convert_dat_to_cfg(input_cfg)
    end
    
    cfg_job_runner = CfgJobRunner(input_cfg)
    
    @info "start with parameters: $args"
    
    if haskey(args, "cmp_chrs") && args["cmp_chrs"] !== nothing
        cmp_chrs = [String(x) for x in args["cmp_chrs"]]
        @info "Overwriting the default hic_analysis.CHR_SYNONYMS: $(join(cmp_chrs, ",")) with which HiCs will be compared!"
        global HicAnalysis.CHR_SYNONYMS = cmp_chrs  # Already Vector{String} from line 109
    else
        cmp_chrs_str = get_sim_property(cfg_job_runner, "cmp_chrs")
        if !isempty(cmp_chrs_str)
            global HicAnalysis.CHR_SYNONYMS = [String(x) for x in split(cmp_chrs_str, r"[\s;,]+")]
        end
    end
    
    exp_cool = haskey(args, "exp_cool") && !isempty(args["exp_cool"]) ? args["exp_cool"] : get_sim_property(cfg_job_runner, "exp_cool")
    interaction_sites = get_sim_property(cfg_job_runner, "interaction_sites")
    lef_loading_sites = get_sim_property(cfg_job_runner, "lef_loading_sites")
    basal_loading_factor = get_sim_property(cfg_job_runner, "basal_loading_factor")
    boundary = haskey(args, "boundary") && args["boundary"] !== nothing ? args["boundary"] : get_sim_property(cfg_job_runner, "boundary")
    tads_boundary = haskey(args, "tads_boundary") && args["tads_boundary"] !== nothing ? args["tads_boundary"] : get_sim_property(cfg_job_runner, "tads_boundary")
    
    cmp_hic_file = isempty(exp_cool) ? nothing : replace(exp_cool, r"\.(cool|mcool)" => "")
    
    # Find last HIC.hdf5 do analysis and store
    sim_hic_file = HicAnalysis.get_last_hic(args["analyse"])
    
    hic_folder = dirname(sim_hic_file)
    plots_folder = joinpath(hic_folder, HicAnalysis.PLOTS_FOLDER)
    
    plot_cmap = get_property(cfg_job_runner, "stats", "plot_cmap", HicAnalysis.CMAP)
    plot_format = get_property(cfg_job_runner, "stats", "plot_format", HicAnalysis.PLOT_FORMAT)
    
    if !isempty(exp_cool)
        # Ensure CHR_SYNONYMS is Vector{String} (convert from Vector{Any} if needed)
        chrs_vec = Vector{String}([String(x) for x in HicAnalysis.CHR_SYNONYMS])
        (chi2_lin, alpha_lin) = HicAnalysis.compare_hic_chromosome(sim_hic_file, cmp_hic_file;
                                                                   chrs=chrs_vec,
                                                                   res=HicAnalysis.RESOLUTION, tads_boundary=tads_boundary,
                                                                   norm=true, chi2_mode=HicAnalysis.CHI2_MODE_LINEAR,
                                                                   plot_cmap=plot_cmap, plot_format=plot_format)
        (chi2_log, alpha_log) = HicAnalysis.compare_hic_chromosome(sim_hic_file, cmp_hic_file;
                                                                   chrs=chrs_vec,
                                                                   res=HicAnalysis.RESOLUTION, tads_boundary=tads_boundary,
                                                                   plots_folder=plots_folder, norm=true,
                                                                   chi2_mode=HicAnalysis.CHI2_MODE_LOG,
                                                                   plot_cmap=plot_cmap, plot_format=plot_format)
    else
        (chi2_lin, alpha_lin) = (0.0, 1.0)
        (chi2_log, alpha_log) = (0.0, 1.0)
    end
    
    @info "comparing $sim_hic_file with $exp_cool result chi2: $chi2_log"
    
    # Save stats in repository file
    if !isfile(args["stats_file"])
        # Create header - create empty DataFrame with column names
        header = ["sim_hic_file", "sim_out_folder", "exp_cool", "resolution",
                  "boundary", "boundary_direction", "tads_boundary", "input.cfg",
                  "nlef", "km", "radius_contact", "chi2_log", "alpha_log",
                  "chi2_lin", "alpha_lin", "chr"]
        # Create empty DataFrame with column names
        df = DataFrame([[] for _ in header], Symbol.(header))
        CSV.write(args["stats_file"], df, append=false)
    end
    
    # Simple file-based lock for concurrent access protection
    lock_file = "$(args["stats_file"]).lock"
    max_wait = 60  # Maximum seconds to wait for lock
    wait_time = 0.1  # Seconds between lock checks
    elapsed = 0.0
    while isfile(lock_file) && elapsed < max_wait
        sleep(wait_time)
        elapsed += wait_time
    end
    if isfile(lock_file)
        @warn "Lock file still exists after $max_wait seconds, proceeding anyway"
    end
    # Create lock file
    touch(lock_file)
    try
        cp = args["contact_probability"] ? "p" : ""
        radius_p = args["radius_contact"] > 0 ? "$(args["radius_contact"])$cp" : ""
        
        row = DataFrame(
            sim_hic_file = [sim_hic_file],
            sim_out_folder = [args["output_folder"]],
            exp_cool = [exp_cool],
            resolution = [HicAnalysis.RESOLUTION],
            boundary = [boundary],
            boundary_direction = [args["boundary_direction"]],
            tads_boundary = [tads_boundary],
            input_cfg = [args["input_cfg"]],
            nlef = [args["nlef"]],
            km = [args["km"]],
            radius_contact = [radius_p],
            chi2_log = [chi2_log],
            alpha_log = [alpha_log],
            chi2_lin = [chi2_lin],
            alpha_lin = [alpha_lin],
            chr = [isempty(HicAnalysis.CHR_SYNONYMS) ? "" : String(HicAnalysis.CHR_SYNONYMS[end])]
        )
        CSV.write(args["stats_file"], row, append=true)
    finally
        # Remove lock file
        if isfile(lock_file)
            rm(lock_file, force=true)
        end
    end
    
    # Generate hic_*_hot_r.png plot files if not already done
    hic_plot_file = replace(sim_hic_file, ".hdf5" => "_$(plot_cmap).$(plot_format)")
    if !isfile(hic_plot_file) || args["replace"]
        PlotHic.run(args["analyse"]; resolution=HicAnalysis.RESOLUTION, cmap=plot_cmap, plot_format=plot_format)
    else
        @info "Hic plot file already created $hic_plot_file so skip it"
    end
end

end # module

