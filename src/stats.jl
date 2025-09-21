

### Calculate stats

# 1. Distribution of reads over chroms
# 2. Read length histogram
# 3. histogram of all mod counts for each mod
# 4. Mean rate of each mod per read
# 5. msp/nucs per reat
# 6. mean msp size/mean nuc size/mean fire size



abstract type RecordStat end

plotstat(stat::String, dict::Dict{T, V}) where {T, V} = plotstat(stat, dict[stat])
plotstat(stat::String, data::DataFrame) = plotstat(getfield(SMGStats, Symbol(stat)), data)
statname(stat::S) where {S <: AbstractString} = statname(getfield(SMGStats, Symbol(stat)))
statname(stat) = string("Unk: ", stat)
############################# chromosome counts #############################
struct ChromStat{S} <: RecordStat
    stat::S
end

ChromStat() = ChromStat(CountMap(String))
update!(stat::ChromStat, reader::HTSFileReader, record::BamRecord, recorddata) = fit!(stat.stat, refname(reader, record))
function chromlt(chrA, chrB)

    if occursin("_", chrA) && occursin("_", chrB)
        chrA < chrB
    elseif occursin("_", chrB)
        true
    elseif occursin("_", chrA)
        false
    else
        if isnumeric(chrA[4]) && isnumeric(chrB[4])
            parse(Int, chrA[4:end]) < parse(Int, chrB[4:end])
        elseif isnumeric(chrA[4])
            true
        else
            false
        end
    end
end
function unicodeplot(stat::ChromStat)
    chroms = sort(collect(keys(stat.stat)), lt=chromlt)
    counts = [stat.stat.value[c] for c in chroms]
    [barplot(chroms, counts, title="Chromosome counts")]
end
function writestats(stat::ChromStat, path::String, file="chromstats.tsv")
    filepath = joinpath(path, file)
    chroms = sort(collect(keys(stat.stat)), lt=chromlt)
    counts = [stat.stat.value[c] for c in chroms]
    df = DataFrame(chrom=chroms, count=counts)
    CSV.write(filepath, df, delim='\t', compress=endswith(file, ".gz"))
end
#  TODO:: define this if a specialised stat reader is necessary
# readstats(::Type{ChromStat}, file) =
function plotstat(::Type{ChromStat}, data::DataFrame, scale=1000)
    n = size(data, 1)
    f = Figure(size=(800, 800))
    ax = Axis(f[1, 1], xgridvisible=false, ygridvisible=false, xlabel="Counts (k)", yticks=(1:n, data.chrom), xticklabelrotation=π/4, yreversed=true, title="Total reads per chromosome (thousands)")
    barplot!(1:n, data.count/scale, direction=:x, bar_labels=:y, flip_labels_at=0.8*maximum(data.count)/scale, color_over_bar=:white)
    # text!(data.count/scale, 1:n, text=string.(round.(data.count/scale, digits=2)), align=(:left, :middle))
    f



end
statname(::Type{ChromStat}) = "Chromosome Counts"

############################# read length hist #############################
struct ReadLengthHist{S} <: RecordStat
    stat::S
end
ReadLengthHist() = ReadLengthHist(KHist(250))
update!(stat::ReadLengthHist, reader::HTSFileReader, record::BamRecord, recorddata) = fit!(stat.stat, log10(Float64(querylength(record))))
unicodeplot(stat::ReadLengthHist) = [lineplot(first.(stat.stat.bins), last.(stat.stat.bins), title="Read Length Histogram")]
function writestats(stat::ReadLengthHist, path::String, file="readlengthhist.tsv.gz")
    filepath = joinpath(path, file)
    readlengths = first.(stat.stat.bins)
    counts = last.(stat.stat.bins)
    df = DataFrame(readlength=readlengths, count=counts)
    CSV.write(filepath, df, delim='\t', compress=endswith(file, ".gz"))
end
function plotstat(::Type{ReadLengthHist}, data::DataFrame, scale=1, xtrans=x -> x)
    f = Figure()
    ax = Axis(f[1, 1], xgridvisible=false, ygridvisible=false, ylabel="Counts (k)")
    lines!(xtrans.(data.readlength), data.count, color=:black)
    mean_read_length = (xtrans.(data.readlength)'*(data.count/sum(data.count)))
    vlines!([mean_read_length/scale], color=:red, label=string("Mean = ", round(mean_read_length/scale, digits=2), " kb"))
    # Legend(f[1, 2], ax)
    axislegend(ax)
    if xtrans == identity
        xmax = 1.1*maximum(data.readlength)
    else
        xmax = min(1.1*xtrans(maximum(data.readlength)), 50_000)
    end
    # xlims!(0, xmax)
    
    # xlims!(0, 100)
    f
end
statname(::Type{ReadLengthHist}) = "Read length histogram"
############################# modifcation hist #############################
struct ModHist <: RecordStat
    data::Dict{Modification, Vector{Int}}
end
ModHist() = ModHist(Dict{Modification, Vector{Int}}())
function update!(stat::ModHist, reader::HTSFileReader, record::BamRecord, recorddata)
    for mi in ModIterator(record, recorddata)
        haskey(stat.data, mi.mod) || (stat.data[mi.mod] = zeros(Int, 256))
        stat.data[mi.mod][mi.prob + 1] += 1
    end
end
filtunibarplot(labels, values, thresh = 0, ind=values .> thresh; kwargs...) = barplot(labels[ind], values[ind]; kwargs...)
function unicodeplot(stat::ModHist)
    labels = string.(0:255)
    bps = [filtunibarplot(labels, vec, title=string(mod)) for (mod, vec) in stat.data]
    bps
end

function writestats(stat::ModHist, path::String, file="modhist.tsv.gz")
    filepath = joinpath(path, file)
    df = DataFrame(Dict(string(m) => v for (m, v) in stat.data))
    modlabels = names(df)
    df.ml = 0:255
    df.prob = df.ml./255
    df = df[!, ["ml" ; "prob" ; modlabels]]
    # df = mapreduce((mod, vec) -> DataFrame(mod=mod, ml=ml, prob=prob, count=count), stat.data)
    CSV.write(filepath, df, delim='\t', compress=endswith(file, ".gz"))
end
function plotstat(::Type{ModHist}, df::DataFrame, scale=1e+6, thresh=0.9)
    sdf = stack(df, names(df, r"mod"))
    sdf.ismod = sdf.prob .> thresh
    transform!(groupby(sdf, :variable), :value => (v -> v/sum(v)) => :proportion)


    mdf = combine(groupby(sdf, [:variable, :ismod]), :proportion => sum => :proportion)
    mdf = mdf[mdf.ismod, :]

   

    f = Figure(; size=(1000, 400))
    prate = data(mdf) * mapping(:variable, :proportion, color=:variable) * visual(BarPlot, bar_labels=:y)
    draw!(f[1, 1], prate, axis= (; xgridvisible=false, ygridvisible=false, ylabel="Proportion modified", xlabel="Modification"))
    ylims!(-0.005, 1.2*maximum(mdf.proportion))

    plt = data(sdf) * mapping(:prob, :proportion , color=:variable) * mapping(row=:variable) * visual(BarPlot)
    ag = draw!(f[1, 2], plt, axis= (; xgridvisible=false, ygridvisible=false, ylabel="Mods / million", xlabel="Mod Probability"))
    legend!(f[1, 3], ag)
    f
end
statname(::Type{ModHist}) = "Modification histogram"



############################# mod rate #############################
struct ModRate <: RecordStat
    n::Int
    thresh::Float64
    stat::Dict{Modification, KHist{Float64}}
end
ModRate(n=500, thresh=0.9*255) = ModRate(n, thresh, Dict{Modification, KHist{Float64}}())
function update!(stat::ModRate, reader::HTSFileReader, record::BamRecord, recorddata)
    totalmods = Dict{Modification, Int}()
    for mi in ModIterator(record, recorddata)
        (mi.prob > stat.thresh)   || continue
        haskey(totalmods, mi.mod) || (totalmods[mi.mod] = 0)
        totalmods[mi.mod] = totalmods[mi.mod] + 1
    end
    for (mod, total) in totalmods
        haskey(stat.stat, mod) || (stat.stat[mod] = KHist(stat.n))
        fit!(stat.stat[mod], total/querylength(record))
    end
end
unicodeplot(stat::ModRate) = [lineplot(first.(kh.bins), last.(kh.bins), title=string(mod)) for (mod, kh) in stat.stat]

function writestats(stat::ModRate, path::String, file="modrate.tsv.gz")
    filepath = joinpath(path, file)
    df = mapreduce(((k, v), ) -> DataFrame(mod=k, rate=first.(v.bins), count=last.(v.bins)), vcat, collect(stat.stat))
    CSV.write(filepath, df, delim='\t', compress=endswith(file, ".gz"))
end

function plotstat(::Type{ModRate}, df::DataFrame, scale=1)
    transform!(groupby(df, :mod), :count => (c -> c/sum(c)) => :proportion)

    ### quantiles for boxplot
    qdf = combine(groupby(df, :mod), df -> quantile_from_freq(df.rate, df.proportion))
    sort!(qdf, :mod)
    # display(qdf)
    n = length(unique(df.mod))
    f  = Figure(size=(1000, 400))
    ### draw boxplots from qdf
    ax = Axis(f[1:2, 1:n], xgridvisible=false, ygridvisible=false, ylabel="mod/bp", xticks=(1:n, qdf.mod))
    quantileboxplot!(qdf)

    plt = data(df) * mapping(:rate, :proportion, color=:mod) * mapping(col=:mod) * visual(Lines)
    draw!(f[1, (1:n) .+ n], plt, axis= (; xgridvisible=false, ygridvisible=false))
    plt_overlay = data(df) * mapping(:rate, :proportion, color=:mod) * visual(Lines)
    ag = draw!(f[2, (1:n) .+ n], plt_overlay, axis=(; xlabel="mod/bp", xgridvisible=false, ygridvisible=false))
    legend!(f[1:2, 2*n + 1], ag)

    f
end
statname(::Type{ModRate}) = "Modification Rate Histogram"


############################# nuc msp rate #############################
struct NucMSPRate <: RecordStat
    nuc::KHist{Float64}
    msp::KHist{Float64}
end
NucMSPRate(n=500) = NucMSPRate(KHist(n), KHist(n))
function update!(stat::NucMSPRate, reader::HTSFileReader, record::BamRecord, recorddata)
    nn = length(firenucs(record, recorddata))
    nm = length(firemsps(record, recorddata))

    fit!(stat.nuc, nn/querylength(record))
    fit!(stat.msp, nm/querylength(record))
end
unicodeplot(stat::NucMSPRate) = [lineplot(first.(stat.nuc.bins), last.(stat.nuc.bins), title="Nucleosome Rate"),
                                 lineplot(first.(stat.msp.bins), last.(stat.msp.bins), title="MSP Rate")]

function writestats(stat::NucMSPRate, path::String, file="nuc_msp_rate.tsv.gz")
    filepath = joinpath(path, file)

    features = [fill("nuc", length(stat.nuc.bins)) ; fill("msp", length(stat.msp.bins))]
    rates  = [first.(stat.nuc.bins) ; first.(stat.msp.bins)]
    counts = [last.(stat.nuc.bins) ; last.(stat.msp.bins)]

    df = DataFrame(Feature=features, Rate=rates, Count=counts)
    CSV.write(filepath, df, delim='\t', compress=endswith(file, ".gz"))
end


function plotstat(::Type{NucMSPRate}, df::DataFrame, scale=1)
    # plt = data(df) * mapping(:Rate, :Count, color=:Feature) * mapping(col=:Feature) * visual(Lines)
    # f = draw(plt, axis= (; xgridvisible=false, ygridvisible=false))
    # plt_overlay = data(df) * mapping(:Rate, :Count, color=:Feature) * visual(Lines)
    # n = length(unique(df.Feature))
    
    # draw!(f.figure[2, 1:n], plt_overlay, axis=(; xlabel="feature/bp", xgridvisible=false, ygridvisible=false))
    # f


    #####
    transform!(groupby(df, :Feature), :Count => (c -> c/sum(c)) => :proportion)

    ### quantiles for boxplot
    qdf = combine(groupby(df, :Feature, sort=true), df -> quantile_from_freq(df.Rate, df.proportion))
    
    
    n = length(unique(df.Feature))
    f  = Figure(size=(1000, 400))
    ### draw boxplots from qdf
    ax = Axis(f[1:2, 1:n], xgridvisible=false, ygridvisible=false, ylabel="Feature/bp", xticks=(1:n, qdf.Feature))
    quantileboxplot!(qdf)

    plt = data(df) * mapping(:Rate, :proportion, color=:Feature) * mapping(col=:Feature) * visual(Lines)
    draw!(f[1, (1:n) .+ n], plt, axis= (; xgridvisible=false, ygridvisible=false))
    plt_overlay = data(df) * mapping(:Rate, :proportion, color=:Feature) * visual(Lines)
    ag = draw!(f[2, (1:n) .+ n], plt_overlay, axis=(; xlabel="Feature/bp", xgridvisible=false, ygridvisible=false))
    legend!(f[1:2, 2*n + 1], ag)

    f
end
statname(::Type{NucMSPRate}) = "Nucleosome and MSP rate historgram"


############################# nuc msp size #############################
struct NucMSPLenHist <: RecordStat
    nuc::Vector{Int}
    msp::Vector{Int}
end
NucMSPLenHist(n=1000) = NucMSPLenHist(zeros(n), zeros(n))
function update!(stat::NucMSPLenHist, reader::HTSFileReader, record::BamRecord, recorddata)
    for (ns, nl) in firenucs(record, recorddata)
        if 1 ≤ nl ≤ length(stat.nuc)
            stat.nuc[nl] += 1
        end
    end
    for (as, al, aq) in firemsps(record, recorddata)
        if 1 ≤ al ≤ length(stat.nuc)
            stat.msp[al] += 1
        end
    end
   
end
unicodeplot(stat::NucMSPLenHist) = [barplot(1:length(stat.nuc), stat.nuc, title="Nucleosome Length Hist"),
                                    barplot(1:length(stat.msp), stat.msp, title="MSP Length HIst")]

function writestats(stat::NucMSPLenHist, path::String, file="nuc_msp_length_hist.tsv.gz")
    filepath = joinpath(path, file)

    features = [fill("nuc", length(stat.nuc)) ; fill("msp", length(stat.msp))]
    lengths  = [1:length(stat.nuc) ; 1:length(stat.msp)]
    counts   = [stat.nuc ; stat.msp]

    df = DataFrame(Feature=features, Length=lengths, Count=counts)
    CSV.write(filepath, df, delim='\t', compress=endswith(file, ".gz"))
end

function plotstat(::Type{NucMSPLenHist}, df::DataFrame, scale=1)
    df = df[df.Length .> 1, :]
    plt = data(df) * mapping(:Length, :Count, color=:Feature) * mapping(col=:Feature) * visual(Lines)
    f = draw(plt, axis= (; xgridvisible=false, ygridvisible=false))
    plt_overlay = data(df) * mapping(:Length, :Count, color=:Feature) * visual(Lines)
    n = length(unique(df.Feature))
    
    draw!(f.figure[2, 1:n], plt_overlay, axis=(; xgridvisible=false, ygridvisible=false))
    f
end
statname(::Type{NucMSPLenHist}) = "Nucleosome and MSP length histogram"


############################# mod meta read size #############################
struct ModMetaHist <: RecordStat
    n::Int
    thresh::Float64
    data::Dict{Modification, Vector{Int}}
end
ModMetaHist(n=1000) = ModMetaHist(n, 0.9*255, Dict{Modification, Vector{Int}}())
function update!(stat::ModMetaHist, reader::HTSFileReader, record::BamRecord, recorddata)
    for mi in ModIterator(record, recorddata)
        (mi.prob > stat.thresh)   || continue
        haskey(stat.data, mi.mod) || (stat.data[mi.mod] = zeros(Int, stat.n))

        i = Int(floor(stat.n*mi.pos/querylength(record)))
        i = max(i, 1)
        stat.data[mi.mod][i] += 1
    end
end
# unicodeplot(stat::ModMetaHist) = ...


function writestats(stat::ModMetaHist, path::String, file="mod_meta_hist.tsv.gz")
    filepath = joinpath(path, file)
    df = mapreduce(((mod, counts), ) -> DataFrame(mod=mod, bin=range(0, 1, length=stat.n), count=counts), vcat, collect(stat.data))
    CSV.write(filepath, df, delim='\t', compress=endswith(file, ".gz"))
end

function plotstat(::Type{ModMetaHist}, df::DataFrame, scale=1)

    plt = data(df) * mapping(:bin, :count, color=:mod) * mapping(col=:mod) * visual(Lines)
    f = draw(plt, axis= (; xgridvisible=false, ygridvisible=false))
    plt_overlay = data(df) * mapping(:bin, :count, color=:mod) * visual(Lines)
    n = length(unique(df.mod))
    
    draw!(f.figure[2, 1:n], plt_overlay, axis=(; xlabel="mod/bp", xgridvisible=false, ygridvisible=false))
    f
end
statname(::Type{ModMetaHist}) = "Modification Distribution along read"


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