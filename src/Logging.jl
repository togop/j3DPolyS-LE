"""
Logging utilities for 3DPolyS-LE
"""
module LoggingUtils

# Import from Base.CoreLogging to avoid conflicts
import Base.CoreLogging: AbstractLogger, @info, @warn, @error, @debug, SimpleLogger, Info

const _loggers = Dict{String, AbstractLogger}()

function get_logger(name::String = "j3DPolySLE")
    global _loggers
    
    if !haskey(_loggers, name)
        _loggers[name] = SimpleLogger(stdout, Info)
    end
    
    return _loggers[name]
end

end # module

