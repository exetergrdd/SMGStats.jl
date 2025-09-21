### online stats for streaming reads


abstract type RecordStat end

plotstat(stat::String, dict::Dict{T, V}) where {T, V} = plotstat(stat, dict[stat])
plotstat(stat::String, data::DataFrame) = plotstat(getfield(SMGStats, Symbol(stat)), data)
statname(stat::S) where {S <: AbstractString} = statname(getfield(SMGStats, Symbol(stat)))
statname(stat) = string("Unk: ", stat)

### Each stat should define
# 1. A subtype of RecordStat: struct MyStat <: RecordStat
# 2. A constructor: MyStat()
# 3. A function that returns the name: statname(::Type{MyStat}) = 
# 3. An update function: update!(stat::MyStat, reader, record, recorddata)
# 4. A output function with default filename: writestats(stat::MyStat, path, file="mystat.tsv.gz")
#     returns written file
#     writes a dataframe
# 5. A plotting function (using CairoMakie) defined on the dataframe written in 4:
#    plotstat(::Type{MyStat}, data::DataFrame)


include("chromcounts.jl")
include("readlengthhist.jl")
include("modhist.jl")
include("modrate.jl")
include("nucmsprate.jl")
include("nucmsplength.jl")
include("modmeta.jl")





####################################################################
####################################################################
####################################################################

function firestats(file, nr = 100_000)
    
    reader = open(HTSFileReader, file)
    @show nrecords(reader)
    recorddata = StencillingData(AuxMapModFire())

    bamstats = (ChromStat(), ReadLengthHist(), ModHist(), ModRate(), NucMSPRate(), NucMSPLenHist(), ModMetaHist())
    r = 0
    it = eachrecord(reader)
    nt = length(it)
    @show nt
    p = Progress(nr == -1 ? nt : min(nt, nr))

    for record in it
        processread!(record, recorddata)

        for bamstat in bamstats
            update!(bamstat, reader, record, recorddata)
        end
        r += 1
        next!(p)
        (nr == r) && break
    end
    close(reader)
    bamstats
end

function writeallstats(stats, path::String, file="meta.tsv")
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

function readstats(path::String, file="meta.tsv")
    meta = CSV.read(joinpath(path, file), DataFrame)
    Dict(statlabel => CSV.read(file, DataFrame) for (statlabel, file) in zip(meta.Stat, meta.File))
end


function samplesummary(data, thresh=0.9)

    summary = DataFrame(Stat=String[], value=Float64[], Units=String[])
    ### total reads
    if haskey(data, "ChromStat")
        totalreads = sum(data["ChromStat"].count)/1e+6
        push!(summary, ("Total Reads", totalreads, "M"))
    end

    ### total bases and N50 mean +/- std read length
    if haskey(data, "ReadLengthHist")
        df = data["ReadLengthHist"]
        df.rl = 10.0.^(df.readlength)
        df.totalbases = df.rl.*df.count

        sort!(df.readlength, rev=true)
        df.cumlative = cumsum(df.totalbases)
        push!(summary, ("Total Bases", df.cumlative[1]/1e+6, "Gb"))
    
        n50 = df.rl[findfirst(df.cumlative .> 0.5*df.cumlative[1])]/1000
        
        push!(summary, ("N50", n50, "Kb"))


        mean_read_length = df.rl'*(df.count/sum(df.count))
        std_read_length  = sqrt((df.rl.*df.rl)'*(df.count/sum(df.count)) - mean_read_length*mean_read_length)
        
        push!(summary, ("Mean Read Length", mean_read_length/1000, "Kb"))
        push!(summary, ("Std Read Length", std_read_length/1000, "Kb"))

    end

    ### total modified bases
    if haskey(data, "ModHist")
        df = data["ModHist"]
        sdf = stack(df, names(df, r"mod"))
        sdf.ismod = sdf.prob .> thresh
        transform!(groupby(sdf, :variable), :value => (v -> v/sum(v)) => :proportion)

        mdf = combine(groupby(sdf, [:variable, :ismod]), :proportion => sum => :proportion, :value => sum => :value)
        mdf = mdf[mdf.ismod, :]
        mdf.value ./= 1e+6
        rename!(mdf, :value => :total, :variable => :mod)
        
        ss = stack(mdf, [:proportion, :total])
        ss.Stat = string.(ss.mod, " ", ss.variable)
        ss.Units = ifelse.(ss.variable .== "total", "Gb", "proportion")

        append!(summary, ss[!, [:Stat, :value, :Units]])

    end

    ### modification rates
    summary
end