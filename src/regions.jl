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
    intervals::Dict{S,Vector{ChromInterval}}
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


function loadchromintervals(file; group=1, dw=1500, regions=Dict{String,Vector{ChromInterval}}())
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
    regions = Dict{String,Vector{ChromInterval}}()
    for (k, bedconfig) in enumerate(config)
        loadchromintervals(bedconfig.file, group=k, dw=bedconfig.width, regions=regions)
    end

    ChromIntervalCollection(regions, length(config), [c.label for c in config], [c.width for c in config])
end



function bedhasheader(bedfile)

    if endswith(bedfile, ".gz")
        io = open(bedfile) |> GzipDecompressorStream
    else
        io = open(bedfile)
    end
    hasheader = false
    while !eof(io)
        line = readline(io)

        if startswith(line, "chrom") || startswith(line, "#chrom")
            hasheader = true
            break
        elseif startswith(line, "#")
            continue
        else
            break
        end

    end
    close(io)
    hasheader
end

"""
    loadbed(bedfile)

Loads a BED file into a DataFrame. Handles standard BED files and those with headers.
If no header is present, columns are named `:chrom`, `:start`, `:stop`, and optionally `:name`.

# Arguments
- `bedfile`: Path to the BED file.

# Returns
A `DataFrame` containing the BED file data.
"""
function loadbed(bedfile)
    if bedhasheader(bedfile)
        bdf = CSV.read(bedfile, DataFrame)
        # normalise typical chrom header variants
        rename!(bdf, replace.(names(bdf), "#chrom" => "chrom"))
    else
        bdf = CSV.read(bedfile, DataFrame, header=false)
        rename!(bdf, 1 => :chrom, 2 => :start, 3 => :stop)
        if size(bdf, 2) >= 4
            rename!(bdf, 4 => :name)
        else
            bdf[!, :name] = string.(replace(basename(bedfile), ".bed" => ""), "_", 1:size(bdf, 1))
        end
    end
    bdf
end

"""
    bedintervals(bdf; chromlabel=:chrom, startlabel=:start, stoplabel=:stop, namelabel=:name)

Converts a DataFrame of BED data (from `loadbed`) into a `GenomicIntervalCollection`.
Supports flexible column names.

# Arguments
- `bdf`: DataFrame containing BED data.
- `chromlabel`: Column name for chromosome (default `:chrom`).
- `startlabel`: Column name for start position (default `:start`).
- `stoplabel`: Column name for stop position (default `:stop`).
- `namelabel`: Column name for interval name/metadata (default `:name`).
"""
function bedintervals(bdf; chromlabel=:chrom, startlabel=:start, stoplabel=:stop, namelabel=:name)
    # Validate columns exist
    for col in [chromlabel, startlabel, stoplabel]
        if col ∉ propertynames(bdf)
            error("Column :$col not found in DataFrame. Available columns: $(propertynames(bdf))")
        end
    end

    sort!(bdf, [chromlabel, startlabel, stoplabel])

    metadata = if namelabel ∈ propertynames(bdf)
        bdf[!, namelabel]
    else
        # If namelabel not found, create simple index-based metadata
        [i for i in 1:size(bdf, 1)]
    end

    giv = GenomicInterval.(bdf[!, chromlabel], bdf[!, startlabel], bdf[!, stoplabel], '.', metadata)
    GenomicIntervalCollection(giv)
end

"""
    buildintervalset(intervalfile, label="", chromlabel=:chrom, startlabel=:start, stoplabel=:stop)

Wrapper for `loadbed` and `bedintervals` for compatibility.
"""
function buildintervalset(intervalfile, label="", chromlabel=:chrom, startlabel=:start, stoplabel=:stop)
    df = loadbed(intervalfile)
    # For buildintervalset, if label is specific, we might handle metadata differently,
    # but strictly following previous behavior:

    # Check if we need to rename or if the default loadbed names match user expectation.
    # If user passes generic :chrom, :start, :stop, loadbed usually handles it.
    # If file has header with specific names, they used to be passed here.

    bedintervals(df, chromlabel=chromlabel, startlabel=startlabel, stoplabel=stoplabel)
end


"""
    bedchromintervals(bdf; chromlabel=:chrom, startlabel=:start, stoplabel=:stop, strandlabel=:strand)

Converts a DataFrame of BED data (from `loadbed`) into a `ChromIntervalCollection`.
Uses `ChromInterval` structure.
"""
function bedchromintervals(bdf; chromlabel=:chrom, startlabel=:start, stoplabel=:stop, strandlabel=:strand)
    for col in [chromlabel, startlabel, stoplabel]
        if col ∉ propertynames(bdf)
            error("Column :$col not found in DataFrame. Available columns: $(propertynames(bdf))")
        end
    end

    sort!(bdf, [chromlabel, startlabel, stoplabel])

    regions = Dict{String,Vector{ChromInterval}}()
    
    # Handle strand
    strands = if strandlabel ∈ propertynames(bdf)
        bdf[!, strandlabel]
    else
        fill("+", size(bdf, 1))
    end
    
    for (i, (c, s, e, str)) in enumerate(zip(bdf[!, chromlabel], bdf[!, startlabel], bdf[!, stoplabel], strands))
        cs = String(c)
        if !haskey(regions, cs)
            regions[cs] = ChromInterval[]
        end
        # Use simple 1-based indexing for regionindex (i)
        # strand: "+" -> true, else false
        push!(regions[cs], ChromInterval(s, e, str == "+" || str == true, i, 1))
    end

    ChromIntervalCollection(regions, 1, ["regions"], Int[])
end


"""
    buildchromintervalset(intervalfile; chromlabel=:chrom, startlabel=:start, stoplabel=:stop)

Wrapper for `loadbed` and `bedchromintervals`.
"""
function buildchromintervalset(intervalfile; chromlabel=:chrom, startlabel=:start, stoplabel=:stop)
    df = loadbed(intervalfile)
    bedchromintervals(df, chromlabel=chromlabel, startlabel=startlabel, stoplabel=stoplabel)
end
