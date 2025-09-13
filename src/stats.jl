

### Calculate stats

# 1. Distribution of reads over chroms
# 2. histogram of all mod counts for each mod
# 3. Mean rate of each mod per read
# 4. nucs per read
# 5. msp per re
# 6. mean msp size
# 7. mean nuc size

# struct ModStat{T, F, O}
#     stat::T
#     name::String
#     fn::F
#     onlinestat::O
# end

# struct ChromCounter end
# chromcounter() = ModStat(ChromCounter(), "Chromosome Counts", (reader, record, recorddata) -> refname(reader, record), CountMap(String))

abstract type RecordStat end
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

############################# read length hist #############################
struct ReadLengthHist{S} <: RecordStat
    stat::S
end
ReadLengthHist() = ReadLengthHist(KHist(100))
update!(stat::ReadLengthHist, reader::HTSFileReader, record::BamRecord, recorddata) = fit!(stat.stat, Float64(querylength(record)))
unicodeplot(stat::ReadLengthHist) = [lineplot(first.(stat.stat.bins), last.(stat.stat.bins), title="Read Length Histogram")]


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

############################# mod rate #############################
struct ModRate <: RecordStat
    n::Int
    thresh::Float64
    stat::Dict{Modification, KHist{Float64}}
end
ModRate(n=100, thresh=0.9*255) = ModRate(n, thresh, Dict{Modification, KHist{Float64}}())
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


############################# nuc msp rate #############################
struct NucMSPRate <: RecordStat
    nuc::KHist{Float64}
    msp::KHist{Float64}
end
NucMSPRate(n=100) = NucMSPRate(KHist(n), KHist(n))
function update!(stat::NucMSPRate, reader::HTSFileReader, record::BamRecord, recorddata)
    nn = length(firenucs(record, recorddata))
    nm = length(firemsps(record, recorddata))

    fit!(stat.nuc, nn/querylength(record))
    fit!(stat.msp, nm/querylength(record))
end
unicodeplot(stat::NucMSPRate) = [lineplot(first.(stat.nuc.bins), last.(stat.nuc.bins), title="Nucleosome Rate"),
                              lineplot(first.(stat.msp.bins), last.(stat.msp.bins), title="MSP Rate")]

function firestats(file)
    reader = open(HTSFileReader, file)
    recorddata = StencillingData(AuxMapModFire())

    bamstats = [ChromStat(), ReadLengthHist(), ModHist(), ModRate(), NucMSPRate()]
    for record in eachrecord(reader)
        processread!(record, recorddata)

        for bamstat in bamstats
            update!(bamstat, reader, record, recorddata)
        end
    end

    close(reader)
    bamstats
end

############################# nuc msp size #############################
struct NucMSPLenHist <: RecordStat
    nuc::Vector{Int}
    msp::Vector{Int}
end
NucMSPLenHist(n=500) = NucMSPLenHist(zeros(n), zeros(n))
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

function firestats(file)
    reader = open(HTSFileReader, file)
    recorddata = StencillingData(AuxMapModFire())

    bamstats = [ChromStat(), ReadLengthHist(), ModHist(), ModRate(), NucMSPRate(), NucMSPLenHist()]
    for record in eachrecord(reader)
        processread!(record, recorddata)

        for bamstat in bamstats
            update!(bamstat, reader, record, recorddata)
        end
    end

    close(reader)
    bamstats
end


bamfile = "/Users/ndlo201/julia/dev/SMGReader/test/data/test.bam"
stats = firestats(bamfile)
ups = unicodeplot.(stats);
displayplots(ups)