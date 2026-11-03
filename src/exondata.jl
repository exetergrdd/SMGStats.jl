#=
    exondata.jl: Fast annotation and gene model plotting using exondata (BED-like GTF table)
=#

struct TranscriptModel
    chrom::String
    strand::Bool            # true for '+', false for '-'
    genetype::String
    geneid::String
    genename::String
    transcriptid::String
    start::Int
    stop::Int
    exon_starts::Vector{Int}
    exon_stops::Vector{Int}
end

struct ExonData
    by_chrom::Dict{String, Vector{TranscriptModel}}
end

"""
    parse_coord_list(s)

Parses a comma-separated string of integers into a Vector{Int}.
Returns an empty vector if missing or empty.
"""
function parse_coord_list(s)::Vector{Int}
    (ismissing(s) || isnothing(s)) && return Int[]
    str = strip(string(s))
    isempty(str) && return Int[]
    res = Int[]
    for part in split(str, ',')
        p = strip(part)
        isempty(p) && continue
        val = tryparse(Int, p)
        !isnothing(val) && push!(res, val)
    end
    return res
end

clean_string(s)::String = ismissing(s) ? "" : strip(string(s), ['"', ' ', '\t'])

"""
    load_exondata(file::AbstractString)::ExonData

Loads a gzipped or uncompressed exondata file (e.g. gencode.v49.annotation.exondata.gtf.gz)
into an indexed `ExonData` structure organized by chromosome.
"""
function load_exondata(file::AbstractString)::ExonData
    !isfile(file) && error("ExonData file not found: $file")
    
    df = CSV.read(file, DataFrame; stringtype=String)
    by_chrom = Dict{String, Vector{TranscriptModel}}()

    chroms        = df.chrom
    strands       = df.strand
    gene_types    = hasproperty(df, :GeneType) ? df.GeneType : fill("", nrow(df))
    gene_ids      = hasproperty(df, :GeneID) ? df.GeneID : fill("", nrow(df))
    gene_names    = hasproperty(df, :GeneName) ? df.GeneName : fill("", nrow(df))
    transcript_ids = hasproperty(df, :TranscriptID) ? df.TranscriptID : fill("", nrow(df))
    starts        = df.start
    stops         = df.stop
    exon_starts   = hasproperty(df, :exon_starts) ? df.exon_starts : fill("", nrow(df))
    exon_stops    = hasproperty(df, :exon_stops) ? df.exon_stops : fill("", nrow(df))

    for i in 1:nrow(df)
        c = clean_string(chroms[i])
        str_val = clean_string(strands[i])
        strand = (str_val == "+" || str_val == "1")
        
        gt = clean_string(gene_types[i])
        gid = clean_string(gene_ids[i])
        gn = clean_string(gene_names[i])
        tid = clean_string(transcript_ids[i])
        
        s = Int(starts[i])
        e = Int(stops[i])
        es = parse_coord_list(exon_starts[i])
        ee = parse_coord_list(exon_stops[i])

        tm = TranscriptModel(c, strand, gt, gid, gn, tid, s, e, es, ee)
        if !haskey(by_chrom, c)
            by_chrom[c] = TranscriptModel[]
        end
        push!(by_chrom[c], tm)
    end

    # Sort each chromosome vector by start coordinate for fast searching
    for (c, list) in by_chrom
        sort!(list, by = m -> (m.start, m.stop))
    end

    return ExonData(by_chrom)
end

"""
    annotate_locus(exondata::ExonData, chrom::AbstractString, loc::UnitRange{Int})

Annotates a genomic locus by finding:
1. Overlapping transcripts.
2. The nearest gene TSS and signed distance in bp (negative = upstream promoter, positive = downstream/inside gene).
3. Functional context classification (Promoter, Gene Body, or Distal Intergenic).
"""
function annotate_locus(exondata::ExonData, chrom::AbstractString, loc::UnitRange{Int})
    chr_key = String(chrom)
    if !haskey(exondata.by_chrom, chr_key)
        # Try prepending or stripping 'chr'
        alt_chr = startswith(chr_key, "chr") ? chr_key[4:end] : ("chr" * chr_key)
        if haskey(exondata.by_chrom, alt_chr)
            chr_key = alt_chr
        else
            return (; nearest_gene="None", nearest_gene_id="", gene_type="", tss_distance=0, context="Unknown", transcripts=TranscriptModel[])
        end
    end

    models = exondata.by_chrom[chr_key]
    isempty(models) && return (; nearest_gene="None", nearest_gene_id="", gene_type="", tss_distance=0, context="Unknown", transcripts=TranscriptModel[])

    loc_start = first(loc)
    loc_stop  = last(loc)
    mid_pos   = div(loc_start + loc_stop, 2)

    overlapping = TranscriptModel[]
    nearest_idx = 1
    min_dist_abs = typemax(Int)
    signed_dist = 0

    for (k, m) in enumerate(models)
        # Check overlap
        if m.start <= loc_stop && m.stop >= loc_start
            push!(overlapping, m)
        end

        # Calculate TSS distance
        tss = m.strand ? m.start : m.stop
        # Signed distance: negative = upstream of TSS (promoter), positive = downstream of TSS
        cur_dist = m.strand ? (mid_pos - tss) : (tss - mid_pos)
        cur_dist_abs = abs(cur_dist)
        if cur_dist_abs < min_dist_abs
            min_dist_abs = cur_dist_abs
            signed_dist = cur_dist
            nearest_idx = k
        end
    end

    best_m = models[nearest_idx]
    
    # Classify context
    context = if !isempty(overlapping)
        # Check if overlaps any exon
        has_exon = false
        for tm in overlapping
            for (es, ee) in zip(tm.exon_starts, tm.exon_stops)
                if es <= loc_stop && ee >= loc_start
                    has_exon = true
                    break
                end
            end
            has_exon && break
        end
        if abs(signed_dist) <= 1000
            "Promoter"
        elseif has_exon
            "Exonic"
        else
            "Intronic"
        end
    else
        if abs(signed_dist) <= 1500
            "Promoter"
        else
            "Distal Intergenic"
        end
    end

    return (;
        nearest_gene = best_m.genename,
        nearest_gene_id = best_m.geneid,
        gene_type = best_m.genetype,
        tss_distance = signed_dist,
        context = context,
        transcripts = overlapping
    )
end

"""
    filter_gene_models(models::Vector{TranscriptModel}; longest_only::Bool=true)

Filters transcripts to select only the representative longest transcript per gene
(or all transcripts if `longest_only=false`).
"""
function filter_gene_models(models::Vector{TranscriptModel}; longest_only::Bool=true)::Vector{TranscriptModel}
    !longest_only && return models
    isempty(models) && return TranscriptModel[]

    by_gene = Dict{String, TranscriptModel}()
    for m in models
        key = isempty(m.genename) ? m.geneid : m.genename
        if !haskey(by_gene, key)
            by_gene[key] = m
        else
            cur_span = by_gene[key].stop - by_gene[key].start
            new_span = m.stop - m.start
            if new_span > cur_span
                by_gene[key] = m
            end
        end
    end
    res = collect(values(by_gene))
    sort!(res, by = m -> (m.start, m.stop))
    return res
end

"""
    pack_gene_rows(models::Vector{TranscriptModel}; buffer::Int=1000)

Greedy interval row packing: assigns each transcript model to a vertical row (1, 2, ...)
so non-overlapping genes share the same row. Returns a Vector of row indices matching `models`.
"""
function pack_gene_rows(models::Vector{TranscriptModel}; buffer::Int=1000)::Vector{Int}
    n = length(models)
    n == 0 && return Int[]
    rows = zeros(Int, n)
    row_ends = Int[]

    for i in 1:n
        m = models[i]
        placed = false
        for r in 1:length(row_ends)
            if row_ends[r] + buffer < m.start
                rows[i] = r
                row_ends[r] = m.stop
                placed = true
                break
            end
        end
        if !placed
            push!(row_ends, m.stop)
            rows[i] = length(row_ends)
        end
    end
    return rows
end

"""
    plot_genemodels!(ax, models::Vector{TranscriptModel}, loc::UnitRange{Int}; 
                     longest_only::Bool=true, exon_height=0.4, buffer=1000, color=:darkslateblue)

Plots gene models (intron lines, exon rectangles, and gene name labels) onto the given Axis `ax`.
"""
function plot_genemodels!(ax, models::Vector{TranscriptModel}, loc::UnitRange{Int}; 
                          longest_only::Bool=true, exon_height=0.4, buffer=1000, color=:darkslateblue, rows=nothing)
    filtered = longest_only ? filter_gene_models(models; longest_only=true) : models
    if isempty(filtered)
        # Empty placeholder axis
        hidespines!(ax, :l, :t, :r, :b)
        ax.yticks = ([], [])
        xlims!(ax, first(loc), last(loc))
        return ax
    end

    model_rows = isnothing(rows) ? pack_gene_rows(filtered; buffer=buffer) : rows
    max_row = isempty(model_rows) ? 1 : maximum(model_rows)
    margin = max(50, div(length(loc), 40))

    for (m, r) in zip(filtered, model_rows)
        # 1. Intron backbone line across transcript span (clipped to loc window ± margin)
        intron_start = max(m.start, first(loc) - 5000)
        intron_stop  = min(m.stop, last(loc) + 5000)
        lines!(ax, [intron_start, intron_stop], [r, r], color=:gray50, linewidth=1.5)

        # 2. Exons as solid polygons
        ex_polys = Rect2f[]
        for (es, ee) in zip(m.exon_starts, m.exon_stops)
            # Only draw exons that touch or are near the loc window
            if ee >= (first(loc) - 1000) && es <= (last(loc) + 1000)
                w = max(ee - es, 5) # at least 5 bp width for visual clarity
                push!(ex_polys, Rect2f(es, r - exon_height / 2, w, exon_height))
            end
        end
        if !isempty(ex_polys)
            poly!(ax, ex_polys, color=color)
        end

        # 3. Gene Name label centered over the visible span of the transcript in this window
        vis_start = max(m.start, first(loc))
        vis_stop  = min(m.stop, last(loc))
        mid_label = clamp(div(vis_start + vis_stop, 2), first(loc) + margin, last(loc) - margin)
        strand_sym = m.strand ? " →" : " ←"
        label_text = longest_only ? string(m.genename, strand_sym) : string(m.genename, " (", first(split(m.transcriptid, '.')), ")", strand_sym)
        text!(ax, mid_label, r + exon_height / 2 + 0.12; 
              text=label_text, align=(:center, :bottom), fontsize=9, font=:bold, color=:black)
    end

    hidespines!(ax, :l, :t, :r, :b)
    ax.yticks = ([], [])
    ylims!(ax, 0.4, max_row + 0.7)
    xlims!(ax, first(loc), last(loc))
    return ax
end
