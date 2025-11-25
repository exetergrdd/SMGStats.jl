### code for intervals on the genome

### Use a simplified structure rather than overhead of GenomicFeatures to reduce memory
### For use with a sorted sweep of reads and genomic intervals


struct ChromInterval
    start::Int
    stop::Int
    strand::Bool
    regionindex::Int
    group::Int
end

struct ChromIntervalCollection{S}
    intervals::Dict{S, Vector{ChromInterval}}
    numgroups::Int
    grouplabels::Vector{String}
    dw::Vector{Int}
end



function chromintervals_to_tid(cic, reader)
    tid_to_chrom = referencedict(reader)
    chrom_to_tid = Dict(c => tid for (tid, c) in tid_to_chrom)

    intervals = Dict(chrom_to_tid[c] => v for (c, v) in cic.intervals)
    ChromIntervalCollection(intervals, cic.numgroups, cic.grouplabels, cic.dw)
end




function Base.isless(a::ChromInterval, b::ChromInterval)
    (a.start < b.start) && return true
    if a.start == b.start
        (a.stop < b.stop) && return true
        return a.strand < b.strand
    else
        return false
    end

end

Base.:(==)(a::ChromInterval, b::ChromInterval) = a.start == b.start && a.stop == b.stop && a.strand == b.strand


function loadchromintervals(file; group=1, dw=1500, regions = Dict{String, Vector{ChromInterval}}())
    df = CSV.read(file, DataFrame) 
    strands = "strand" ∈ names(df) ? df.strand : fill("+", size(df, 1))
    for (k, (c, s, e, str)) in enumerate(zip(df.chrom, df.start, df.stop, strands))
        cs = String(c)
        haskey(regions, cs) || (regions[cs] = ChromInterval[])
        mp = div(s + 1 + e, 2) ### bed file is zero based coords
        push!(regions[cs], ChromInterval(mp - dw, mp + dw, str == "+", k, group))
    end

    for v in values(regions)
        sort!(v)
    end
    regions
end



function loadmultichromintervals(config)
    regions = Dict{String, Vector{ChromInterval}}()
    for (k, bedconfig) in enumerate(config)
        loadchromintervals(bedconfig.file, group=k, dw=bedconfig.width, regions=regions)
    end

    ChromIntervalCollection(regions, length(config), [c.label for c in config], [c.width for c in config])
end

