### online stats for streaming reads


abstract type RecordStat end

struct RecordUpdates end
struct ModUpdates end
struct PostModUpdates end

recordupdates(::Type) = nothing
modupdates(::Type) = nothing
postmodupdates(::Type) = nothing


plotstat(stat::String, dict::Dict{T, V}) where {T, V} = plotstat(stat, dict[stat])
plotstat(stat::String, data::DataFrame) = plotstat(getfield(SMGStats, Symbol(stat)), data)
statname(stat::S) where {S <: AbstractString} = statname(getfield(SMGStats, Symbol(stat)))
statname(stat) = string("Unk: ", stat)

### Each stat should define
# 1. A subtype of RecordStat: struct MyStat <: RecordStat
# 2. An instantiate(::Type{MyStat}, config) function
# 3. A function that returns the name: statname(::Type{MyStat}) = 
# 4. Functions to determine when the stat should be updated
# 4a. isrecordstat(::MyStat) returns true if stat should be updated with each record with `updaterecord!`
# 4b. ismodstat(::MyStat) returns true if stat should be updated with each mod with `updatemod!`
# 4c. postmodupdate(::MyStat) returns true if it need updating after iterating over mods
# Update functions
#   updaterecord!(stat::MyStat, record, recorddata)
#   updatemod!(stat::MyStat, mod, record, recorddata)
#   updatepost!(stat::MyStat, record, recorddata)
# 5. A output function with default filename: writestats(stat::MyStat, path, file="mystat.tsv.gz")
#     returns written file
#     writes a dataframe
# 6. A plotting function (using CairoMakie) defined on the dataframe written in 4:
#    plotstat(::Type{MyStat}, data::DataFrame)

include("statconfig.jl")

include("chromcounts.jl")
# include("readlengthhist.jl")
include("modhist.jl")
include("modrate.jl")
include("nucmsprate.jl")
include("nucmsplength.jl")
include("modmeta.jl")
include("readlengthcounter.jl")
include("modcrosscor.jl")
include("metaplot.jl")


# isrecordstat(x) = true

####################################################################
####################################################################
####################################################################

# statnames = (ChromStat, ReadLengthCounter, ModHist, ModRate, NucMSPRate, NucMSPLenHist, ModMetaHist)

instantiate(x, reader, mods) = error(string("unrecognised: ", x))
      
macro smgstats(types...)
    return Expr(:curly, :Tuple, types...)
end


@inline updatestat!(::Nothing, args...) = nothing
@inline updatestat!(::Type{RecordUpdates}, stat, record, recorddata)    = updaterecord!(stat, record, recorddata)
@inline updatestat!(::Type{ModUpdates}, stat, mod, record, recorddata)  = updatemod!(stat, mod, record, recorddata)
@inline updatestat!(::Type{PostModUpdates}, stat, record, recorddata)   = updatepostmod!(stat, record, recorddata)


@generated function instantiate(::Type{T}, config) where {T <: Tuple}
    exprs = [:(instantiate($(T.parameters[i]), config)) for i in 1:length(T.parameters)]
    return :(($(exprs...),))
end 

@generated function instantiate(::Type{T}, reader, mods) where {T <: Tuple}
    exprs = [:(instantiate($(T.parameters[i]), reader, mods)) for i in 1:length(T.parameters)]
    return :(($(exprs...),))
end 

@generated function update_record_stats!(stats::TT, record, recorddata) where {TT <: Tuple}
    exprs = [:(updatestat!(recordupdates($(TT.parameters[i])), stats[$i], record, recorddata)) for i in 1:length(TT.parameters)]
    return Expr(:block, exprs...)
end


@generated function update_mod_stats!(stats::TT, mod, record, recorddata) where {TT <: Tuple}
    exprs = [:(updatestat!(modupdates($(TT.parameters[i])), stats[$i], mod, record, recorddata)) for i in 1:length(TT.parameters)]
    return Expr(:block, exprs...)
end


@generated function update_postmod_stats!(stats::TT, record, recorddata) where {TT <: Tuple}
    exprs = [:(updatestat!(postmodupdates($(TT.parameters[i])), stats[$i], record, recorddata)) for i in 1:length(TT.parameters)]
    return Expr(:block, exprs...)
end


function calculatestats(file,  stattypes::Tuple{<:Tuple}, mods, recordata; nr=100_000, config=(;))

    reader = open(HTSFileReader, file)
end


function firestats(file, stattypes::Type{<:Tuple}; nr=100_000, config=(;))
    reader = open(HTSFileReader, file)
    recorddata = StencillingData(AuxMapModFire())
    mods = (mod_5mC, mod_5hmC, mod_6mA)

    fullconfig = (; reader, mods, config...)

    stats = instantiate(stattypes, fullconfig)

    r = 0
    for record in eachrecord(reader)
        processread!(record, recorddata)

        update_record_stats!(stats, record, recorddata)
        for mi in ModIterator(record, recorddata)
            update_mod_stats!(stats, mi, record, recorddata)
        end
        update_postmod_stats!(stats, record, recorddata)

        r += 1
        (r == nr) && break
    end

    close(reader)
    stats
end



function writeallstats(stats, path::String="", file="meta.tsv"; bamfile="", statdir="stats")
    if isempty(path) && !isempty(bamfile)
        path = joinpath(dirname(bamfile), statdir)
    end
    mkpath(path)
    filepath = joinpath(path, file)
    statnames = String[]
    files = String[]
    for s in stats
        f = writestats(s, path)
        println("Written: ", f)
        push!(statnames, string(nameof(typeof(s))))
        push!(files, f)
    end
    f = CSV.write(filepath, DataFrame(Stat=statnames, File=files), delim='\t')
    println("Written: ", f)
end

function readstats(path::String="", file="meta.tsv" ; bamfile="", statdir="stats")
    if isempty(path) && !isempty(bamfile)
        path = joinpath(dirname(bamfile), statdir)
    end
    meta = CSV.read(joinpath(path, file), DataFrame)
    stats = Dict(statlabel => CSV.read(file, DataFrame) for (statlabel, file) in zip(meta.Stat, meta.File))

    (;path, meta, stats)
end


function samplesummary(data, thresh=0.9)

    summary = DataFrame(Stat=String[], value=Float64[], Units=String[])
    modstats = DataFrame(Mod=String[], Percent=Float64[], Total_Gb=Float64[])
    ### total reads
    if haskey(data, "ChromStat")
        totalreads = sum(data["ChromStat"].count)/1e+6
        push!(summary, ("Total Reads", totalreads, "M"))
    end

    ### total bases and N50 mean +/- std read length
    if haskey(data, "ReadLengthCounter")
        df = data["ReadLengthCounter"]

        totalseq = df.readlength'*df.count
        push!(summary, ("Total Bases", totalseq/1e+6, "Gb"))
        cumlat = 0
        n50_i = 0
        for i = size(df, 1):-1:1
            cumlat += df.readlength[i]*df.count[i]
            if cumlat > totalseq*0.5
                n50_i = i - 1
                break
            end
        end
    
        n50 = df.readlength[n50_i]/1000
        push!(summary, ("N50*", n50, "Kb"))

        w = Weights(df.count)
        mean_read_length = mean(df.readlength, w)
        std_read_length = std(df.readlength, w)
        
        push!(summary, ("Mean Read Length", mean_read_length/1000, "Kb"))
        push!(summary, ("Std Read Length", std_read_length/1000, "Kb"))
    

    end

    ### total modified bases
    if haskey(data, "ModHist")
        df = data["ModHist"]
        sdf = stack(df, names(df, r"mod"))
        sdf.ismod = sdf.prob .> thresh
        transform!(groupby(sdf, :variable), :value => (v -> v/sum(v)) => :proportion)

        mdf = combine(groupby(sdf, [:variable, :ismod]), :proportion => sum => :proportion, :value => sum => :Total_Gb)
        mdf = mdf[mdf.ismod, :]
        mdf.Total_Gb ./= 1e+6
        mdf.proportion .*= 100
        display(mdf)

        rename!(mdf,  :variable => :Mod, :proportion => :Percent)
        modstats = mdf[!, [:Mod, :Percent, :Total_Gb]]
        
        modstats.Mod = replace.(modstats.Mod, "mod_" => "")
        sort!(modstats, :Mod)

    end

    ### modification rates
    summary, modstats
end