function buildsmgstats(v::Vector{String})
    return Tuple{ (getfield(SMGStats, Symbol(s)) for s in v)... }
end

function build_namedtuple(data)
    if data isa Dict
        # Dict → NamedTuple
        return NamedTuple(Symbol(k) => v for (k, v) in data)
    elseif data isa Vector
        # Vector of Dicts → Vector of NamedTuples
        return [build_namedtuple(d) for d in data]
    else
        error("Unsupported config element: $(typeof(data))")
    end
end

function warncompatible(f, a)

    if compatible(f, a)
        return true
    else
        smglog("Warning stat: ", f, " removed due to incompatibility with file's ", a)
        return false
    end
end
function loadstatconfig(file, auxmap=nothing)
    yamlconfig = YAML.load_file(file)

    stattypes = [getfield(SMGStats, Symbol(s)) for s in yamlconfig["stats"]]
    
    if !isnothing(auxmap)
        filter!(f -> warncompatible(f, typeof(auxmap)), stattypes)
    end
    stats = Tuple{stattypes...}
    
    config = NamedTuple(Symbol(k) => build_namedtuple(c) for (k, c) in yamlconfig["configs"])

    stats, config
end