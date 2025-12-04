############################# chromosome counts #############################
struct ChromStat{S} <: RecordStat
    stat::S
    chroms::Dict{Int32, String}
end

instantiate(::Type{ChromStat}, config) = ChromStat(CountMap(Int32), referencedict(config.reader))
# instantiate(::Type{ChromStat}, reader, mods) = ChromStat(CountMap(Int32), referencedict(reader))
ChromStat(reader) = ChromStat(CountMap(Int32), referencedict(reader))

recordupdates(::Type{<:ChromStat})  = RecordUpdates
modupdates(::Type{<:ChromStat})     = nothing
postmodupdates(::Type{<:ChromStat}) = nothing

@inline updaterecord!(stat::ChromStat, record::BamRecord, recorddata) = fit!(stat.stat, record.core.tid)




statname(::Type{ChromStat}) = "Chromosome Counts"

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

function writestats(stat::ChromStat, path::String, file="chromstats.tsv")
    filepath = joinpath(path, file)
    
    chromids = sort(collect(keys(stat.stat)))

    counts = [stat.stat.value[c] for c in chromids]
    chroms = [stat.chroms[c] for c in chromids]

    df = DataFrame(chrom=chroms, count=counts)
    sort!(df, order(:chrom, by=chromlt))
    CSV.write(filepath, df, delim='\t', compress=endswith(file, ".gz"))
end
#  TODO:: define this if a specialised stat reader is necessary
# readstats(::Type{ChromStat}, file) =
function plotstat(::Type{ChromStat}, data::DataFrame, scale=1000)
    n = size(data, 1)
    f = Figure(size=(800, 800))
    Axis(f[1, 1], xgridvisible=false, ygridvisible=false, xlabel="Counts (k)", yticks=(1:n, data.chrom), xticklabelrotation=π/4, yreversed=true, title="Total reads per chromosome (thousands)")
    barplot!(1:n, data.count/scale, direction=:x, bar_labels=:y, flip_labels_at=0.8*maximum(data.count)/scale, color_over_bar=:white)
    f
end
