#=
    report_catalog.jl: Multi-region clustering catalog dashboard and batch report builder
=#

"""
    Lightweight, type-stable JSON printer (no extra external dependencies)
"""
function to_json(val)
    io = IOBuffer()
    print_json(io, val)
    return String(take!(io))
end

print_json(io::IO, x::Nothing) = print(io, "null")
print_json(io::IO, x::Missing) = print(io, "null")
print_json(io::IO, x::Bool) = print(io, x ? "true" : "false")
print_json(io::IO, x::Number) = (isnan(x) || isinf(x)) ? print(io, "null") : print(io, x)

function print_json(io::IO, s::AbstractString)
    print(io, "\"")
    for c in s
        if c == '"'
            print(io, "\\\"")
        elseif c == '\\'
            print(io, "\\\\")
        elseif c == '\n'
            print(io, "\\n")
        elseif c == '\r'
            print(io, "\\r")
        elseif c == '\t'
            print(io, "\\t")
        elseif c == '/'
            print(io, "\\/")
        elseif c < ' '
            print(io, "\\u", string(Int(c), base=16, pad=4))
        else
            print(io, c)
        end
    end
    print(io, "\"")
end

print_json(io::IO, s::Symbol) = print_json(io, string(s))

function print_json(io::IO, v::AbstractVector)
    print(io, "[")
    for (i, el) in enumerate(v)
        i > 1 && print(io, ",")
        print_json(io, el)
    end
    print(io, "]")
end

function print_json(io::IO, nt::NamedTuple)
    print(io, "{")
    first = true
    for (k, v) in pairs(nt)
        !first && print(io, ",")
        first = false
        print_json(io, string(k))
        print(io, ":")
        print_json(io, v)
    end
    print(io, "}")
end

function print_json(io::IO, d::AbstractDict)
    print(io, "{")
    first = true
    for (k, v) in d
        !first && print(io, ",")
        first = false
        print_json(io, string(k))
        print(io, ":")
        print_json(io, v)
    end
    print(io, "}")
end

Base.@kwdef struct LocusSummary
    region_id::String
    slug::String
    chrom::String
    peak_start::Int
    peak_stop::Int
    peak_width::Int
    loc_start::Int
    loc_stop::Int
    peak_name::String
    category::String
    nearest_gene::String
    nearest_gene_id::String
    gene_type::String
    tss_distance::Int
    context::String
    optimal_k::Int
    total_reads::Int
    prop_open::Float64
    min_nrl::Union{Nothing, Float64}
    max_nrl::Union{Nothing, Float64}
    has_open_cluster::Bool
    has_closed_cluster::Bool
    has_strand_skew::Bool
    has_hap_skew::Bool
    has_sample_skew::Bool
    sig_cpg_count::Int
    report_url::String
end

function print_json(io::IO, x::LocusSummary)
    print(io, "{")
    for (i, k) in enumerate(fieldnames(LocusSummary))
        i > 1 && print(io, ",")
        print_json(io, string(k))
        print(io, ":")
        print_json(io, getfield(x, k))
    end
    print(io, "}")
end


"""
    update_clustering_nrl!(clustering; phase_threshold::Float64=0.035, nrl_feature::Symbol=:fire)

Ensures that `clustering.clustertable` has accurate, dynamically re-evaluated spectral NRLs
using relative autocorrelation lags (starting at 0) and the specified phasing threshold.
Mutates `clustering.clustertable` in place and returns `clustering`.
"""
function update_clustering_nrl!(clustering; phase_threshold::Float64=0.035, nrl_feature::Symbol=:fire)
    if haskey(clustering, :clustertable) && haskey(clustering, :CM) && haskey(clustering, :optimal_clustering)
        try
            ct = clustering.clustertable
            !issorted(ct.cluster) && sort!(ct, :cluster)
            centers = clustering.optimal_clustering.optimal_model.centers
            k = size(centers, 2)
            cids = ct.cluster

            nrls_fire = _spectral_nrl_columns(clustering.CM, nrl_feature, k; M=centers, phase_threshold=phase_threshold)
            ct.nrl                = [nrls_fire.nrl[c] for c in cids]
            ct.phasing_score      = [nrls_fire.phasing_score[c] for c in cids]
            ct.is_closed          = [nrls_fire.is_closed[c] for c in cids]
            ct.peak_type          = [nrls_fire.peak_type[c] for c in cids]
            ct.nrl_fire           = ct.nrl
            ct.phasing_score_fire = ct.phasing_score
            ct.is_closed_fire     = ct.is_closed
            ct.peak_type_fire     = ct.peak_type

            if :mod_6mA in clustering.CM.features
                nrls_6mA = _spectral_nrl_columns(clustering.CM, :mod_6mA, k; M=centers, phase_threshold=phase_threshold)
                ct.nrl_6mA           = [nrls_6mA.nrl[c] for c in cids]
                ct.phasing_score_6mA = [nrls_6mA.phasing_score[c] for c in cids]
                ct.is_closed_6mA     = [nrls_6mA.is_closed[c] for c in cids]
                ct.peak_type_6mA     = [nrls_6mA.peak_type[c] for c in cids]
            end

            # Backfill read-level 6mA NRLs on pre-existing JLD2 files if missing
            if haskey(clustering, :readtable) && !isnothing(clustering.readtable)
                rt = clustering.readtable
                n_reads = nrow(rt)
                if :mod_6mA in clustering.CM.features && !hasproperty(rt, :nrl_6mA)
                    read_nrls_6mA = _spectral_nrl_columns(clustering.CM, :mod_6mA, n_reads; phase_threshold=phase_threshold)
                    rt.nrl_6mA                    = read_nrls_6mA.nrl
                    rt.phasing_score_6mA          = read_nrls_6mA.phasing_score
                    rt.phasing_implies_closed_6mA = read_nrls_6mA.is_closed
                    rt.phasing_peak_type_6mA      = read_nrls_6mA.peak_type
                end
                if hasproperty(rt, :nrl_6mA) && !hasproperty(ct, :mean_read_nrl_6mA)
                    ct.total_valid_nrls_6mA        = [count(!ismissing, view(rt.nrl_6mA, rt.cluster .== c)) for c in cids]
                    ct.mean_read_nrl_6mA           = [filtmean(view(rt.nrl_6mA, rt.cluster .== c)) for c in cids]
                    ct.mean_read_phasing_score_6mA = [filtmean(view(rt.phasing_score_6mA, rt.cluster .== c)) for c in cids]
                end
            end
        catch err
            # Silently fallback if CM or model cannot be re-evaluated
        end
    end
    return clustering
end

"""
    load_clustering(fpath::AbstractString; phase_threshold::Float64=0.035)

Safely loads a clustering NamedTuple from a `.jld2` file, automatically resolving
either the `"clustering"` key, `"single_stored_object"`, or the primary root group,
and updates spectral NRL properties using relative autocorrelation lags.
"""
function load_clustering(fpath::AbstractString; phase_threshold::Float64=0.035)
    obj = jldopen(fpath, "r") do file
        k = keys(file)
        if "clustering" in k
            return file["clustering"]
        elseif "single_stored_object" in k
            return file["single_stored_object"]
        elseif !isempty(k)
            return file[first(k)]
        else
            error("Empty JLD2 file: $fpath")
        end
    end
    update_clustering_nrl!(obj; phase_threshold=phase_threshold)
    return obj
end

function _build_peak_lookup(peakfile)
    isnothing(peakfile) && return nothing
    bdf = peakfile isa DataFrame ? peakfile : loadbed(peakfile)
    lookup = Dict{Tuple{String, Int, Int}, Tuple{String, String}}()
    for (i, r) in enumerate(eachrow(bdf))
        pname = hasproperty(r, :PeakName) ? String(string(r.PeakName)) : "peak_$i"
        pcat  = hasproperty(r, :Category) ? String(string(r.Category)) : "ClusterPeak"
        lookup[(String(r.chrom), Int(r.start), Int(r.stop))] = (pname, pcat)
    end
    return lookup
end

function _resolve_peak_meta(clustering; idx::Int=1, peak_lookup=nothing,
                            peakname::Union{Nothing, AbstractString}=nothing,
                            category::Union{Nothing, AbstractString}=nothing)
    chrom = String(clustering.chrom)
    pstart = Int(first(clustering.peakloc))
    pstop  = Int(last(clustering.peakloc))

    pname = if !isnothing(peakname)
        String(peakname)
    elseif !isnothing(peak_lookup) && haskey(peak_lookup, (chrom, pstart, pstop))
        peak_lookup[(chrom, pstart, pstop)][1]
    elseif hasproperty(clustering, :PeakName)
        String(clustering.PeakName)
    elseif hasproperty(clustering.clustertable, :PeakName) && nrow(clustering.clustertable) > 0
        String(clustering.clustertable.PeakName[1])
    else
        "peak_$idx"
    end

    pcat = if !isnothing(category)
        String(category)
    elseif !isnothing(peak_lookup) && haskey(peak_lookup, (chrom, pstart, pstop))
        peak_lookup[(chrom, pstart, pstop)][2]
    elseif hasproperty(clustering, :Category)
        String(clustering.Category)
    elseif hasproperty(clustering.clustertable, :Category) && nrow(clustering.clustertable) > 0
        String(clustering.clustertable.Category[1])
    else
        "ClusterPeak"
    end

    # Backfill in-memory tables if loaded from older JLD2 files
    if hasproperty(clustering, :clustertable) && !isnothing(clustering.clustertable)
        ct = clustering.clustertable
        !hasproperty(ct, :PeakName) && (ct.PeakName = fill(pname, nrow(ct)))
        !hasproperty(ct, :Category) && (ct.Category = fill(pcat, nrow(ct)))
    end
    if hasproperty(clustering, :readtable) && !isnothing(clustering.readtable)
        rt = clustering.readtable
        !hasproperty(rt, :PeakName) && (rt.PeakName = fill(pname, nrow(rt)))
        !hasproperty(rt, :Category) && (rt.Category = fill(pcat, nrow(rt)))
    end

    return pname, pcat
end

"""
    extract_locus_summary(clustering; exondata=nothing, slug=nothing, idx=1, peak_lookup=nothing, peakname=nothing, category=nothing)

Extracts key diagnostic metrics, genomic context, and cluster summaries from a `clustering` object
for catalog indexing and client-side filtering.
"""
function extract_locus_summary(clustering; exondata=nothing, slug::Union{Nothing, AbstractString}=nothing,
                               idx::Int=1, peak_lookup=nothing,
                               peakname::Union{Nothing, AbstractString}=nothing,
                               category::Union{Nothing, AbstractString}=nothing)
    chrom = String(clustering.chrom)
    peak_start = first(clustering.peakloc)
    peak_stop  = last(clustering.peakloc)
    peak_width = peak_stop - peak_start + 1
    loc_start  = first(clustering.loc)
    loc_stop   = last(clustering.loc)
    opt_k      = clustering.optimal_clustering.optimal_k
    total_reads = length(clustering.reads)
    ct         = clustering.clustertable
    rt         = clustering.readtable
    
    resolved_slug = isnothing(slug) ? "$(chrom)_$(peak_start)_$(peak_stop)" : String(slug)
    region_id = "$(chrom):$(peak_start)-$(peak_stop)"
    peak_name, peak_cat = _resolve_peak_meta(clustering; idx=idx, peak_lookup=peak_lookup, peakname=peakname, category=category)
    
    # Genomic annotations on the target regulatory peak
    annot = if !isnothing(exondata)
        annotate_locus(exondata, chrom, clustering.peakloc)
    else
        (; nearest_gene="Unannotated", nearest_gene_id="", gene_type="", tss_distance=0, context="Unknown", transcripts=TranscriptModel[])
    end
    
    # Accessibility
    prop_open = hasproperty(rt, :has_fire) ? 
        round(count(rt.has_fire) / max(1, total_reads), digits=3) : 0.0
    
    # NRLs
    valid_nrls = Float64[]
    if hasproperty(ct, :nrl)
        for val in ct.nrl
            if !ismissing(val) && !isnan(val)
                push!(valid_nrls, Float64(val))
            end
        end
    end
    min_nrl = isempty(valid_nrls) ? nothing : round(minimum(valid_nrls), digits=1)
    max_nrl = isempty(valid_nrls) ? nothing : round(maximum(valid_nrls), digits=1)
    
    # Biological flags
    has_open_cluster = hasproperty(ct, :proportion_has_fire) && any(ct.proportion_has_fire .> 0.25)
    has_closed_cluster = hasproperty(ct, :is_closed) && any(ct.is_closed)
    
    # Strand / Haplotype / Sample Fisher skew flags (prefer precomputed clustertable p-values)
    has_strand_skew = if opt_k <= 1
        false
    elseif hasproperty(ct, :strand_pvalue)
        any(p -> !isnan(p) && p < 0.01, ct.strand_pvalue)
    elseif hasproperty(rt, :strand) && length(unique(rt.strand)) > 1
        any(k -> (pv = fisher(rt.cluster .== k, rt.strand).fisher_bpv; !isnan(pv) && pv < 0.01), 1:opt_k)
    else
        false
    end
    
    has_hap_skew = if opt_k <= 1
        false
    elseif hasproperty(ct, :hap_pvalue) && any(!isnan, ct.hap_pvalue)
        any(p -> !isnan(p) && p < 0.01, ct.hap_pvalue)
    elseif hasproperty(rt, :hap) && any(rt.hap .== 1) && any(rt.hap .== 2)
        hap_mask = (rt.hap .== 1) .| (rt.hap .== 2)
        if count(hap_mask) >= 10
            sub_clust = rt.cluster[hap_mask]
            sub_hap = rt.hap[hap_mask] .== 1
            any(k -> (pv = fisher(sub_clust .== k, sub_hap).fisher_bpv; !isnan(pv) && pv < 0.01), 1:opt_k)
        else
            false
        end
    else
        false
    end

    has_sample_skew = if opt_k <= 1
        false
    elseif hasproperty(ct, :sample_pvalue) && any(!isnan, ct.sample_pvalue)
        any(p -> !isnan(p) && p < 0.01, ct.sample_pvalue)
    else
        false
    end
    
    sig_cpg_count = (hasproperty(clustering, :modclustertable) && !isnothing(clustering.modclustertable) && nrow(clustering.modclustertable) > 0) ?
        count(r -> r.pvalue < 0.05, eachrow(clustering.modclustertable)) : 0
        
    return LocusSummary(
        region_id = region_id,
        slug = resolved_slug,
        chrom = chrom,
        peak_start = peak_start,
        peak_stop = peak_stop,
        peak_width = peak_width,
        loc_start = loc_start,
        loc_stop = loc_stop,
        peak_name = peak_name,
        category = peak_cat,
        nearest_gene = String(annot.nearest_gene),
        nearest_gene_id = String(annot.nearest_gene_id),
        gene_type = String(annot.gene_type),
        tss_distance = Int(annot.tss_distance),
        context = String(annot.context),
        optimal_k = Int(opt_k),
        total_reads = Int(total_reads),
        prop_open = Float64(prop_open),
        min_nrl = min_nrl,
        max_nrl = max_nrl,
        has_open_cluster = has_open_cluster,
        has_closed_cluster = has_closed_cluster,
        has_strand_skew = has_strand_skew,
        has_hap_skew = has_hap_skew,
        has_sample_skew = has_sample_skew,
        sig_cpg_count = Int(sig_cpg_count),
        report_url = "regions/$(resolved_slug)/report.html"
    )
end

"""
    build_cluster_catalog(inputs, outdir; exondata=nothing, longest_only=true, generate_reports=true, skip_existing=true, limit=nothing, title="SMGStats Chromatin Cluster Catalog")

Scans clustering results (directory of `.jld2` files or list of file paths), computes summary metrics,
generates individual interactive HTML reports, and builds a unified multi-region catalog dashboard (`index.html` + `catalog.json`).
"""
function build_cluster_catalog(inputs, outdir; 
    exondata=nothing, 
    peakfile=nothing,
    longest_only=true, 
    gene_flank=15_000,
    phase_threshold=0.035,
    generate_reports=true, 
    skip_existing=true, 
    limit=nothing,
    title="SMGStats Chromatin Cluster Catalog"
)
    mkpath(outdir)
    
    # 1. Normalize input files
    files = String[]
    if inputs isa AbstractString
        if isdir(inputs)
            for f in readdir(inputs, join=true)
                if endswith(f, ".jld2")
                    push!(files, f)
                end
            end
        elseif isfile(inputs)
            push!(files, inputs)
        end
    elseif inputs isa AbstractVector
        for f in inputs
            if isfile(f) && endswith(f, ".jld2")
                push!(files, f)
            end
        end
    end
    
    sort!(files)
    if !isnothing(limit) && limit < length(files)
        files = files[1:limit]
    end
    
    n_total = length(files)
    if n_total == 0
        error("No .jld2 clustering files found in inputs: $inputs")
    end
    
    smglog("Catalog: Found $n_total clustering files to index into $outdir")
    
    # 2. Load exondata annotations and optional peakfile lookup
    loaded_exondata = if exondata isa AbstractString
        smglog("Catalog: Loading ExonData annotations from $exondata...")
        load_exondata(exondata)
    else
        exondata
    end
    peak_lookup = _build_peak_lookup(peakfile)
    
    # 3. Process each region
    catalog_entries = LocusSummary[]
    seen_slugs = Dict{String, Int}()
    p = Progress(n_total, desc="Generating Cluster Catalog: ", barlen=40)
    
    for (idx, fpath) in enumerate(files)
        try
            clustering = load_clustering(fpath; phase_threshold=phase_threshold)
            reorder_clusters!(clustering)
            base_slug = "$(clustering.chrom)_$(first(clustering.peakloc))_$(last(clustering.peakloc))"
            slug_count = get(seen_slugs, base_slug, 0) + 1
            seen_slugs[base_slug] = slug_count
            unique_slug = slug_count == 1 ? base_slug : "$(base_slug)_$(slug_count)"

            sum_entry = extract_locus_summary(clustering; exondata=loaded_exondata, slug=unique_slug, idx=idx, peak_lookup=peak_lookup)
            push!(catalog_entries, sum_entry)
            
            if generate_reports
                region_dir = joinpath(outdir, "regions", sum_entry.slug)
                report_html_path = joinpath(region_dir, "report.html")
                if !(skip_existing && isfile(report_html_path))
                    render_cluster_report(clustering, region_dir; 
                        exondata=loaded_exondata, 
                        longest_only=longest_only, 
                        gene_flank=gene_flank,
                        phase_threshold=phase_threshold,
                        peakname=sum_entry.peak_name,
                        category=sum_entry.category,
                        catalog_path="../../index.html",
                        guide_path="../../guide.html",
                        update_nrl=false
                    )
                end
            end
        catch err
            smglog("Warning: Failed processing $fpath: $err")
        end
        next!(p)
    end
    
    # 4. Save catalog.json
    json_path = joinpath(outdir, "catalog.json")
    open(json_path, "w") do io
        print_json(io, catalog_entries)
    end
    smglog("Catalog: Saved JSON summary to $json_path")
    
    # 5. Render Methodology & Guide HTML
    guide_path = joinpath(outdir, "guide.html")
    render_cluster_guide_html(guide_path; back_url="index.html", back_label="Back to Catalog Dashboard")

    # 6. Generate and save master index.html
    html_content = generate_catalog_html(catalog_entries, title)
    index_path = joinpath(outdir, "index.html")
    write(index_path, html_content)
    smglog("Catalog: Rendered Master Dashboard to $index_path")
    
    return index_path
end


"""
    generate_catalog_html(catalog_entries, title)

Constructs the full, standalone, high-performance web dashboard for exploring multi-region clustering results.
Embeds the catalog dataset as JSON for zero-CORS instant offline viewing.
"""
function generate_catalog_html(catalog_entries, title)
    json_data = to_json(catalog_entries)
    n_regions = length(catalog_entries)
    
    # Compute summary KPIs
    multi_state_count = count(r -> r.optimal_k > 1, catalog_entries)
    pct_multi = n_regions > 0 ? round(multi_state_count / n_regions * 100, digits=1) : 0.0
    
    promoter_count = count(r -> r.context == "Promoter", catalog_entries)
    pct_promoter = n_regions > 0 ? round(promoter_count / n_regions * 100, digits=1) : 0.0
    
    hap_skew_count = count(r -> r.has_hap_skew, catalog_entries)
    pct_hap = n_regions > 0 ? round(hap_skew_count / n_regions * 100, digits=1) : 0.0
    
    return """
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>$title</title>
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500;600&display=swap" rel="stylesheet">
  <style>
    :root {
      --bg-main: #f8fafc;
      --bg-surface: #ffffff;
      --bg-card: #f1f5f9;
      --bg-card-hover: #e2e8f0;
      --border-subtle: #e2e8f0;
      --border-glow: rgba(2, 132, 199, 0.3);
      --text-main: #0f172a;
      --text-muted: #475569;
      --text-dim: #94a3b8;
      --accent-blue: #0284c7;
      --accent-indigo: #4f46e5;
      --accent-purple: #7c3aed;
      --accent-emerald: #059669;
      --accent-amber: #d97706;
      --accent-rose: #e11d48;
      --radius-sm: 6px;
      --radius-md: 10px;
      --radius-lg: 16px;
      --font-sans: 'Inter', -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
      --font-mono: 'JetBrains Mono', SFMono-Regular, Menlo, monospace;
    }

    * { box-sizing: border-box; margin: 0; padding: 0; }
    body {
      font-family: var(--font-sans);
      background-color: var(--bg-main);
      color: var(--text-main);
      line-height: 1.5;
      padding-bottom: 5rem;
      min-height: 100vh;
      -webkit-font-smoothing: antialiased;
    }

    .container {
      max-width: 1840px;
      width: 98vw;
      margin: 0 auto;
      padding: 0 0.5rem;
    }

    /* Top Brand Bar */
    .top-nav {
      border-bottom: 1px solid var(--border-subtle);
      background: rgba(255, 255, 255, 0.9);
      backdrop-filter: blur(12px);
      position: sticky;
      top: 0;
      z-index: 50;
    }
    .top-nav-inner {
      display: flex;
      align-items: center;
      justify-content: space-between;
      height: 64px;
    }
    .brand {
      display: flex;
      align-items: center;
      gap: 0.75rem;
      font-weight: 700;
      font-size: 1.15rem;
      color: var(--text-main);
      text-decoration: none;
    }
    .brand-tag {
      background: linear-gradient(135deg, var(--accent-blue), var(--accent-indigo));
      -webkit-background-clip: text;
      -webkit-text-fill-color: transparent;
      font-weight: 800;
    }
    .brand-badge {
      font-size: 0.72rem;
      font-weight: 600;
      text-transform: uppercase;
      letter-spacing: 0.05em;
      background: #e0f2fe;
      border: 1px solid #bae6fd;
      color: #0369a1;
      padding: 0.2rem 0.5rem;
      border-radius: var(--radius-sm);
    }

    /* Hero Banner */
    .hero {
      background: #ffffff;
      border: 1px solid var(--border-subtle);
      border-radius: var(--radius-md);
      padding: 2rem 2.25rem;
      margin: 1.75rem 0 1.25rem 0;
      box-shadow: 0 1px 3px rgba(0,0,0,0.05);
    }
    .hero-title {
      font-size: 2rem;
      font-weight: 800;
      letter-spacing: -0.025em;
      margin-bottom: 0.4rem;
      color: var(--text-main);
    }
    .hero-desc {
      color: var(--text-muted);
      font-size: 0.98rem;
      max-width: 950px;
    }

    /* KPI Grid */
    .kpi-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
      gap: 1.25rem;
      margin: 1.5rem 0;
    }
    .kpi-card {
      background: var(--bg-surface);
      border: 1px solid var(--border-subtle);
      border-radius: var(--radius-md);
      padding: 1.25rem 1.5rem;
      position: relative;
      overflow: hidden;
      box-shadow: 0 1px 3px rgba(0,0,0,0.04);
      transition: transform 0.15s ease, border-color 0.15s ease;
    }
    .kpi-card:hover {
      transform: translateY(-2px);
      border-color: var(--border-glow);
    }
    .kpi-card::before {
      content: '';
      position: absolute;
      top: 0; left: 0; right: 0;
      height: 3px;
    }
    .kpi-card.blue::before { background: linear-gradient(90deg, #0284c7, #4f46e5); }
    .kpi-card.purple::before { background: linear-gradient(90deg, #7c3aed, #e11d48); }
    .kpi-card.emerald::before { background: linear-gradient(90deg, #059669, #0284c7); }
    .kpi-card.amber::before { background: linear-gradient(90deg, #d97706, #e11d48); }

    .kpi-label {
      font-size: 0.8rem;
      font-weight: 600;
      text-transform: uppercase;
      letter-spacing: 0.05em;
      color: var(--text-dim);
      margin-bottom: 0.35rem;
    }
    .kpi-val {
      font-size: 1.85rem;
      font-weight: 700;
      color: var(--text-main);
      display: flex;
      align-items: baseline;
      gap: 0.5rem;
    }
    .kpi-sub {
      font-size: 0.88rem;
      font-weight: 400;
      color: var(--text-muted);
    }

    /* Controls & Filter Toolbar */
    .filter-panel {
      background: var(--bg-surface);
      border: 1px solid var(--border-subtle);
      border-radius: var(--radius-md);
      padding: 1.5rem;
      margin-bottom: 1.75rem;
      box-shadow: 0 1px 3px rgba(0,0,0,0.04);
    }
    .search-row {
      display: flex;
      gap: 1rem;
      margin-bottom: 1.25rem;
      flex-wrap: wrap;
    }
    .search-input-wrapper {
      position: relative;
      flex: 1;
      min-width: 280px;
    }
    .search-input {
      width: 100%;
      background: #ffffff;
      border: 1px solid #cbd5e1;
      border-radius: var(--radius-md);
      padding: 0.65rem 1rem 0.65rem 2.5rem;
      color: var(--text-main);
      font-family: var(--font-sans);
      font-size: 0.92rem;
      outline: none;
      transition: all 0.15s ease;
    }
    .search-input:focus {
      border-color: var(--accent-blue);
      box-shadow: 0 0 0 3px rgba(2, 132, 199, 0.15);
    }
    .search-icon {
      position: absolute;
      left: 0.9rem;
      top: 50%;
      transform: translateY(-50%);
      color: var(--text-dim);
      pointer-events: none;
    }

    .filter-group-row {
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 1.5rem;
      flex-wrap: wrap;
      padding-top: 1rem;
      border-top: 1px solid var(--border-subtle);
    }
    .filter-section {
      display: flex;
      align-items: center;
      gap: 0.75rem;
      flex-wrap: wrap;
    }
    .filter-title {
      font-size: 0.8rem;
      font-weight: 600;
      text-transform: uppercase;
      letter-spacing: 0.05em;
      color: var(--text-dim);
      margin-right: 0.25rem;
    }

    .pill-group {
      display: inline-flex;
      background: #f1f5f9;
      border-radius: var(--radius-md);
      padding: 3px;
      border: 1px solid #e2e8f0;
    }
    .pill-btn {
      background: transparent;
      border: none;
      color: var(--text-muted);
      padding: 0.35rem 0.8rem;
      font-size: 0.82rem;
      font-weight: 500;
      border-radius: var(--radius-sm);
      cursor: pointer;
      transition: all 0.15s ease;
    }
    .pill-btn:hover {
      color: var(--text-main);
    }
    .pill-btn.active {
      background: #ffffff;
      color: var(--accent-blue);
      font-weight: 700;
      box-shadow: 0 1px 2px rgba(0,0,0,0.08);
    }

    .chip-toggle {
      background: #ffffff;
      border: 1px solid #cbd5e1;
      border-radius: var(--radius-sm);
      color: var(--text-muted);
      padding: 0.35rem 0.75rem;
      font-size: 0.82rem;
      font-weight: 500;
      cursor: pointer;
      display: inline-flex;
      align-items: center;
      gap: 0.4rem;
      transition: all 0.15s ease;
    }
    .chip-toggle:hover {
      border-color: #94a3b8;
      color: var(--text-main);
    }
    .chip-toggle.active {
      background: #e0f2fe;
      border-color: var(--accent-blue);
      color: #0369a1;
      font-weight: 600;
    }

    .select-dropdown {
      background: #ffffff;
      border: 1px solid #cbd5e1;
      border-radius: var(--radius-md);
      color: var(--text-main);
      padding: 0.45rem 0.75rem;
      font-size: 0.85rem;
      outline: none;
      cursor: pointer;
    }

    /* Table Container */
    .table-container {
      background: var(--bg-surface);
      border: 1px solid var(--border-subtle);
      border-radius: var(--radius-md);
      overflow: hidden;
      box-shadow: 0 1px 3px rgba(0, 0, 0, 0.04);
    }
    .results-bar {
      display: flex;
      align-items: center;
      justify-content: space-between;
      padding: 0.85rem 1.25rem;
      border-bottom: 1px solid var(--border-subtle);
      font-size: 0.88rem;
      color: var(--text-muted);
      background: #fafafa;
    }
    .results-count strong {
      color: var(--accent-blue);
    }

    .catalog-table {
      width: 100%;
      border-collapse: collapse;
      text-align: left;
      font-size: 0.85rem;
    }
    .catalog-table th {
      background: #f8fafc;
      color: var(--text-muted);
      font-weight: 600;
      text-transform: uppercase;
      font-size: 0.72rem;
      letter-spacing: 0.05em;
      padding: 0.75rem 1rem;
      border-bottom: 1px solid var(--border-subtle);
      white-space: nowrap;
      user-select: none;
      cursor: pointer;
    }
    .catalog-table th:hover {
      color: var(--text-main);
    }
    .catalog-table td {
      padding: 0.75rem 1rem;
      border-bottom: 1px solid var(--border-subtle);
      color: var(--text-main);
      vertical-align: middle;
      transition: background 0.1s ease;
    }
    .catalog-table tbody tr {
      cursor: pointer;
      transition: background 0.15s ease;
    }
    .catalog-table tbody tr:hover td {
      background: #f0f9ff;
    }

    /* Cell Elements */
    .locus-cell {
      font-family: var(--font-mono);
      font-size: 0.85rem;
      font-weight: 600;
      color: var(--accent-blue);
      display: flex;
      align-items: center;
      gap: 0.5rem;
    }
    .copy-btn {
      background: transparent;
      border: none;
      color: var(--text-dim);
      cursor: pointer;
      padding: 2px 4px;
      border-radius: 4px;
      font-size: 0.75rem;
      transition: color 0.15s;
    }
    .copy-btn:hover { color: var(--accent-blue); }

    .gene-title {
      font-weight: 700;
      color: var(--text-main);
      font-size: 0.92rem;
    }
    .gene-dist {
      font-size: 0.75rem;
      font-weight: 500;
      padding: 0.15rem 0.45rem;
      border-radius: var(--radius-sm);
      display: inline-block;
      margin-top: 2px;
    }
    .dist-promoter { background: #dcfce7; color: #15803d; }
    .dist-proximal { background: #e0f2fe; color: #0369a1; }
    .dist-distal   { background: #f1f5f9; color: #475569; }

    .context-badge {
      display: inline-block;
      padding: 0.18rem 0.5rem;
      border-radius: 9999px;
      font-size: 0.74rem;
      font-weight: 600;
    }
    .ctx-promoter  { background: #dcfce7; color: #15803d; border: 1px solid #86efac; }
    .ctx-genebody  { background: #ede9fe; color: #6d28d9; border: 1px solid #ddd6fe; }
    .ctx-intergenic{ background: #f1f5f9; color: #475569; border: 1px solid #e2e8f0; }

    .k-badge {
      font-weight: 700;
      font-size: 0.82rem;
      padding: 0.18rem 0.55rem;
      border-radius: var(--radius-sm);
      display: inline-block;
    }
    .k-1 { background: #f1f5f9; color: #475569; border: 1px solid #e2e8f0; }
    .k-multi { background: #ede9fe; color: #6d28d9; border: 1px solid #ddd6fe; }

    .reads-val {
      font-weight: 600;
      color: var(--text-main);
    }
    .reads-bar {
      height: 4px;
      background: #e2e8f0;
      border-radius: 2px;
      overflow: hidden;
      margin-top: 4px;
      width: 70px;
    }
    .reads-fill {
      height: 100%;
      background: var(--accent-blue);
      border-radius: 2px;
    }

    .openness-pct {
      font-weight: 600;
    }
    .open-high { color: #059669; }
    .open-mid  { color: #0284c7; }
    .open-low  { color: #64748b; }

    .tag-list {
      display: flex;
      flex-wrap: wrap;
      gap: 0.35rem;
    }
    .tag-badge {
      font-size: 0.72rem;
      font-weight: 600;
      padding: 0.15rem 0.45rem;
      border-radius: var(--radius-sm);
    }
    .tag-hap    { background: #ffe4e6; color: #be123c; border: 1px solid #fecdd3; }
    .tag-strand { background: #fef3c7; color: #b45309; border: 1px solid #fde68a; }
    .tag-cpg    { background: #e0f2fe; color: #0369a1; border: 1px solid #bae6fd; }

    .action-btn {
      display: inline-flex;
      align-items: center;
      gap: 0.35rem;
      background: #f0f9ff;
      border: 1px solid #bae6fd;
      color: #0284c7;
      padding: 0.3rem 0.75rem;
      border-radius: var(--radius-sm);
      font-size: 0.8rem;
      font-weight: 600;
      text-decoration: none;
      transition: all 0.15s ease;
      white-space: nowrap;
    }
    .action-btn:hover {
      background: #0284c7;
      color: #ffffff;
      transform: translateX(2px);
    }

    /* Pagination */
    .pagination-bar {
      display: flex;
      align-items: center;
      justify-content: space-between;
      padding: 0.85rem 1.25rem;
      border-top: 1px solid var(--border-subtle);
      background: #fafafa;
    }
    .page-nav {
      display: flex;
      gap: 0.35rem;
    }
    .page-btn {
      background: #ffffff;
      border: 1px solid #cbd5e1;
      color: var(--text-muted);
      min-width: 32px;
      height: 32px;
      display: inline-flex;
      align-items: center;
      justify-content: center;
      border-radius: var(--radius-sm);
      font-size: 0.82rem;
      cursor: pointer;
      transition: all 0.15s ease;
    }
    .page-btn:hover:not(:disabled) {
      color: var(--text-main);
      border-color: #94a3b8;
    }
    .page-btn.active {
      background: var(--accent-blue);
      color: #ffffff;
      font-weight: 700;
      border-color: var(--accent-blue);
    }
    .page-btn:disabled {
      opacity: 0.4;
      cursor: not-allowed;
    }

    /* Empty state */
    .empty-state {
      padding: 4rem 2rem;
      text-align: center;
      color: var(--text-dim);
    }
    .empty-state h3 {
      color: var(--text-muted);
      font-size: 1.15rem;
      margin-bottom: 0.5rem;
    }
  </style>
</head>
<body>

  <!-- Brand Navbar -->
  <nav class="top-nav">
    <div class="container top-nav-inner">
      <div class="brand">
        <span class="brand-tag">SMGStats</span>
        <span>Chromatin Clustering</span>
        <span class="brand-badge">Catalog</span>
      </div>
      <div style="display: flex; align-items: center; gap: 1rem;">
        <a href="guide.html" style="display: inline-flex; align-items: center; gap: 0.4rem; color: #4f46e5; text-decoration: none; font-size: 0.84rem; font-weight: 600; padding: 0.35rem 0.8rem; background: #eef2ff; border: 1px solid #c7d2fe; border-radius: 6px; transition: all 0.15s ease;">📖 Methodology &amp; Guide</a>
        <div style="font-size: 0.85rem; color: var(--text-muted);">
          Single-Molecule Fiber-seq Diagnostics
        </div>
      </div>
    </div>
  </nav>

  <main class="container">
    
    <!-- Hero Header -->
    <header class="hero">
      <h1 class="hero-title">$title</h1>
      <p class="hero-desc">
        Interactive multi-region catalog across <strong>$n_regions</strong> regulatory elements.
        Explore sub-state chromatin heterogeneity, nucleosome repeat lengths (NRL), 5mC differential methylation,
        and allele-specific fiber distributions.
      </p>
    </header>

    <!-- KPI Metric Cards -->
    <section class="kpi-grid">
      <div class="kpi-card blue">
        <div class="kpi-label">Regulatory Elements</div>
        <div class="kpi-val">$n_regions <span class="kpi-sub">analyzed loci</span></div>
      </div>
      <div class="kpi-card purple">
        <div class="kpi-label">Heterogeneous Loci (k* ≥ 2)</div>
        <div class="kpi-val">$pct_multi% <span class="kpi-sub">($multi_state_count loci)</span></div>
      </div>
      <div class="kpi-card emerald">
        <div class="kpi-label">Promoter Proximal (≤1kb)</div>
        <div class="kpi-val">$pct_promoter% <span class="kpi-sub">($promoter_count loci)</span></div>
      </div>
      <div class="kpi-card amber">
        <div class="kpi-label">Allelic / Haplotype Skew</div>
        <div class="kpi-val">$pct_hap% <span class="kpi-sub">($hap_skew_count loci)</span></div>
      </div>
    </section>

    <!-- Search & Filters Panel -->
    <section class="filter-panel">
      <!-- Search Input -->
      <div class="search-row">
        <div class="search-input-wrapper">
          <span class="search-icon">🔍</span>
          <input type="text" id="searchInput" class="search-input" placeholder="Search by PeakName (e.g. ctcf_1), Category (e.g. AccessibleTSS), Gene, or Coordinates...">
        </div>
        <div>
          <select id="sortSelect" class="select-dropdown">
            <option value="reads_desc">Sort: Read Depth (High → Low)</option>
            <option value="reads_asc">Sort: Read Depth (Low → High)</option>
            <option value="k_desc">Sort: Optimal k* (High → Low)</option>
            <option value="tss_asc">Sort: Nearest TSS Distance</option>
            <option value="open_desc">Sort: % Open Molecules</option>
            <option value="coords">Sort: Genomic Position</option>
            <option value="peak_name">Sort: Peak Name</option>
            <option value="category">Sort: Peak Category</option>
          </select>
        </div>
      </div>

      <!-- Filter Controls -->
      <div class="filter-group-row">
        <!-- Context Pills -->
        <div class="filter-section">
          <span class="filter-title">Context:</span>
          <div class="pill-group" id="contextPills">
            <button class="pill-btn active" data-context="all">All</button>
            <button class="pill-btn" data-context="Promoter">Promoter</button>
            <button class="pill-btn" data-context="Gene Body">Gene Body</button>
            <button class="pill-btn" data-context="Intergenic">Intergenic</button>
          </div>
        </div>

        <!-- Category Pills (populated dynamically from dataset categories) -->
        <div class="filter-section">
          <span class="filter-title">Category:</span>
          <div class="pill-group" id="categoryPills">
            <button class="pill-btn active" data-category="all">All</button>
          </div>
        </div>

        <!-- Heterogeneity State -->
        <div class="filter-section">
          <span class="filter-title">Clustering State:</span>
          <div class="pill-group" id="kPills">
            <button class="pill-btn active" data-k="all">All k</button>
            <button class="pill-btn" data-k="single">Single (k=1)</button>
            <button class="pill-btn" data-k="multi">Heterogeneous (k≥2)</button>
            <button class="pill-btn" data-k="complex">Complex (k≥3)</button>
          </div>
        </div>

        <!-- Signal Toggles -->
        <div class="filter-section">
          <span class="filter-title">Signals:</span>
          <button class="chip-toggle" id="toggleOpen">🔥 Open State</button>
          <button class="chip-toggle" id="toggleHap">🧬 Allelic Skew</button>
          <button class="chip-toggle" id="toggleStrand">⮂ Strand Skew</button>
          <button class="chip-toggle" id="toggleSample">🧪 Sample Skew</button>
          <button class="chip-toggle" id="toggleCpG">● Diff 5mC</button>
        </div>
      </div>
    </section>

    <!-- Results Table Container -->
    <section class="table-container">
      <div class="results-bar">
        <div class="results-count">
          Showing <strong id="showingStart">0</strong>–<strong id="showingEnd">0</strong> of <strong id="totalFiltered">0</strong> regions
        </div>
        <div style="display: flex; align-items: center; gap: 0.5rem;">
          <span style="font-size: 0.8rem; color: var(--text-dim);">Per Page:</span>
          <select id="pageSizeSelect" class="select-dropdown" style="padding: 0.25rem 0.5rem; font-size: 0.8rem;">
            <option value="25">25</option>
            <option value="50" selected>50</option>
            <option value="100">100</option>
            <option value="250">250</option>
            <option value="100000">All</option>
          </select>
        </div>
      </div>

      <div style="overflow-x: auto;">
        <table class="catalog-table" id="regionsTable">
          <thead>
            <tr>
              <th data-col="coords">Genomic Locus</th>
              <th data-col="peak_name">Peak Name</th>
              <th data-col="category">Category</th>
              <th data-col="gene">Nearest Gene</th>
              <th data-col="context">Context</th>
              <th data-col="reads">Read Depth</th>
              <th data-col="k">Optimal k*</th>
              <th data-col="open">Openness</th>
              <th data-col="nrl">NRL (bp)</th>
              <th>Biological Flags</th>
              <th style="text-align: right;">Report</th>
            </tr>
          </thead>
          <tbody id="tableBody">
            <!-- Populated dynamically via JS -->
          </tbody>
        </table>
      </div>

      <!-- Empty State -->
      <div id="emptyState" class="empty-state" style="display: none;">
        <h3>No matching genomic regions found</h3>
        <p>Try broadening your search term or clearing active signal filters.</p>
      </div>

      <!-- Pagination Controls -->
      <div class="pagination-bar" id="paginationBar">
        <div style="font-size: 0.82rem; color: var(--text-dim);" id="pageSummary">
          Page 1 of 1
        </div>
        <div class="page-nav" id="pageButtons">
          <!-- Buttons populated by JS -->
        </div>
      </div>
    </section>

  </main>

  <!-- Embedded Catalog Dataset (Zero CORS, 100% offline compatible) -->
  <script id="catalog-data" type="application/json">
  $json_data
  </script>

  <script>
    // --- Application State ---
    const RAW_DATA = JSON.parse(document.getElementById('catalog-data').textContent);
    let maxReads = 1;
    RAW_DATA.forEach(d => { if (d.total_reads > maxReads) maxReads = d.total_reads; });

    let state = {
      searchQuery: '',
      contextFilter: 'all',
      categoryFilter: 'all',
      kFilter: 'all',
      toggleOpen: false,
      toggleHap: false,
      toggleStrand: false,
      toggleSample: false,
      toggleCpG: false,
      sortBy: 'reads_desc',
      pageSize: 50,
      currentPage: 1
    };

    // Populate Category Pills dynamically
    const categoryPills = document.getElementById('categoryPills');
    const uniqueCategories = [...new Set(RAW_DATA.map(d => d.category).filter(Boolean))].sort();
    for (const cat of uniqueCategories) {
      const btn = document.createElement('button');
      btn.className = 'pill-btn';
      btn.dataset.category = cat;
      btn.textContent = cat;
      categoryPills.appendChild(btn);
    }

    // --- DOM Elements ---
    const searchInput = document.getElementById('searchInput');
    const sortSelect = document.getElementById('sortSelect');
    const pageSizeSelect = document.getElementById('pageSizeSelect');
    const tableBody = document.getElementById('tableBody');
    const emptyState = document.getElementById('emptyState');
    const showingStart = document.getElementById('showingStart');
    const showingEnd = document.getElementById('showingEnd');
    const totalFiltered = document.getElementById('totalFiltered');
    const pageButtons = document.getElementById('pageButtons');
    const pageSummary = document.getElementById('pageSummary');

    // --- Filter & Sort Logic ---
    function filterAndSortData() {
      let filtered = RAW_DATA.filter(row => {
        // Search text
        if (state.searchQuery) {
          const q = state.searchQuery.toLowerCase();
          const matchGene = row.nearest_gene && row.nearest_gene.toLowerCase().includes(q);
          const matchGeneId = row.nearest_gene_id && row.nearest_gene_id.toLowerCase().includes(q);
          const matchId = row.region_id && row.region_id.toLowerCase().includes(q);
          const matchChr = row.chrom && row.chrom.toLowerCase().includes(q);
          const matchPeak = row.peak_name && row.peak_name.toLowerCase().includes(q);
          const matchCat = row.category && row.category.toLowerCase().includes(q);
          if (!matchGene && !matchGeneId && !matchId && !matchChr && !matchPeak && !matchCat) return false;
        }

        // Category
        if (state.categoryFilter !== 'all' && row.category !== state.categoryFilter) return false;

        // Context
        if (state.contextFilter !== 'all') {
          if (state.contextFilter === 'Promoter' && row.context !== 'Promoter') return false;
          if (state.contextFilter === 'Gene Body' && row.context !== 'Exonic' && row.context !== 'Intronic' && row.context !== 'Gene Body') return false;
          if (state.contextFilter === 'Intergenic' && !row.context.includes('Intergenic')) return false;
        }

        // K state
        if (state.kFilter === 'single' && row.optimal_k !== 1) return false;
        if (state.kFilter === 'multi' && row.optimal_k < 2) return false;
        if (state.kFilter === 'complex' && row.optimal_k < 3) return false;

        // Toggles
        if (state.toggleOpen && !row.has_open_cluster) return false;
        if (state.toggleHap && !row.has_hap_skew) return false;
        if (state.toggleStrand && !row.has_strand_skew) return false;
        if (state.toggleSample && !row.has_sample_skew) return false;
        if (state.toggleCpG && row.sig_cpg_count === 0) return false;

        return true;
      });

      // Sort
      filtered.sort((a, b) => {
        switch (state.sortBy) {
          case 'reads_desc': return b.total_reads - a.total_reads;
          case 'reads_asc':  return a.total_reads - b.total_reads;
          case 'k_desc':     return b.optimal_k - a.optimal_k || b.total_reads - a.total_reads;
          case 'tss_asc':    return Math.abs(a.tss_distance) - Math.abs(b.tss_distance);
          case 'open_desc':  return b.prop_open - a.prop_open;
          case 'nrl_asc':    return (a.min_nrl ?? 9999) - (b.min_nrl ?? 9999);
          case 'coords':     return a.chrom.localeCompare(b.chrom, undefined, {numeric: true}) || a.peak_start - b.peak_start;
          case 'peak_name':  return (a.peak_name || '').localeCompare(b.peak_name || '', undefined, {numeric: true});
          case 'category':   return (a.category || '').localeCompare(b.category || '') || b.total_reads - a.total_reads;
          default: return 0;
        }
      });

      return filtered;
    }

    // --- Render View ---
    function render() {
      const filtered = filterAndSortData();
      const total = filtered.length;

      // Adjust page bounds
      const maxPages = Math.max(1, Math.ceil(total / state.pageSize));
      if (state.currentPage > maxPages) state.currentPage = maxPages;

      const startIndex = (state.currentPage - 1) * state.pageSize;
      const endIndex = Math.min(startIndex + state.pageSize, total);
      const pageData = filtered.slice(startIndex, endIndex);

      // Summary indicators
      totalFiltered.textContent = total;
      showingStart.textContent = total === 0 ? 0 : startIndex + 1;
      showingEnd.textContent = endIndex;

      if (total === 0) {
        tableBody.innerHTML = '';
        emptyState.style.display = 'block';
        document.getElementById('paginationBar').style.display = 'none';
        return;
      }

      emptyState.style.display = 'none';
      document.getElementById('paginationBar').style.display = 'flex';

      // Table rows
      let html = '';
      for (const row of pageData) {
        // Locus format
        const locusStr = `\${row.chrom}:\${row.peak_start}-\${row.peak_stop}`;
        
        // Gene & TSS distance
        let distBadge = '';
        const absDist = Math.abs(row.tss_distance);
        if (absDist <= 1000) {
          distBadge = `<span class="gene-dist dist-promoter">\${row.tss_distance < 0 ? row.tss_distance + ' bp' : '+' + row.tss_distance + ' bp'}</span>`;
        } else if (absDist <= 10000) {
          distBadge = `<span class="gene-dist dist-proximal">\${(row.tss_distance/1000).toFixed(1)} kb</span>`;
        } else {
          distBadge = `<span class="gene-dist dist-distal">\${(row.tss_distance/1000).toFixed(1)} kb</span>`;
        }

        // Context class
        let ctxClass = 'ctx-intergenic';
        if (row.context === 'Promoter') ctxClass = 'ctx-promoter';
        else if (row.context === 'Exonic' || row.context === 'Intronic' || row.context === 'Gene Body') ctxClass = 'ctx-genebody';

        // Optimal k badge
        const kBadge = row.optimal_k > 1 ? 
          `<span class="k-badge k-multi">k* = \${row.optimal_k}</span>` :
          `<span class="k-badge k-1">k* = 1</span>`;

        // Openness
        const pctOpenVal = (row.prop_open * 100).toFixed(1);
        let openClass = 'open-low';
        if (row.prop_open >= 0.4) openClass = 'open-high';
        else if (row.prop_open >= 0.15) openClass = 'open-mid';

        // Reads bar
        const barPct = Math.min(100, Math.round(row.total_reads / maxReads * 100));

        // NRL display
        let nrlStr = '—';
        if (row.min_nrl !== null && row.max_nrl !== null) {
          nrlStr = row.min_nrl === row.max_nrl ? `\${row.min_nrl}` : `\${row.min_nrl}–\${row.max_nrl}`;
        }

        // Biological flags
        let tagsHtml = '';
        if (row.has_hap_skew) tagsHtml += `<span class="tag-badge tag-hap" title="Allelic / Haplotype Skew (p < 0.01)">Hap Skew</span>`;
        if (row.has_strand_skew) tagsHtml += `<span class="tag-badge tag-strand" title="Strand Skew (p < 0.01)">Strand Skew</span>`;
        if (row.has_sample_skew) tagsHtml += `<span class="tag-badge tag-strand" title="Sample Skew (p < 0.01)">Sample Skew</span>`;
        if (row.sig_cpg_count > 0) tagsHtml += `<span class="tag-badge tag-cpg" title="\${row.sig_cpg_count} differentially methylated CpGs">\${row.sig_cpg_count} Diff 5mC</span>`;
        if (!tagsHtml) tagsHtml = `<span style="color: var(--text-dim); font-size: 0.75rem;">None</span>`;

        html += `
          <tr onclick="window.location.href='\${row.report_url}'">
            <td>
              <div class="locus-cell">
                <span>\${locusStr}</span>
                <button class="copy-btn" onclick="event.stopPropagation(); navigator.clipboard.writeText('\${locusStr}')" title="Copy coordinates">📋</button>
              </div>
              <div style="font-size: 0.75rem; color: var(--text-dim); margin-top: 2px;">\${row.peak_width} bp element</div>
            </td>
            <td>
              <span style="font-family: var(--font-mono); font-size: 0.82rem; font-weight: 600; color: var(--text-main);">\${row.peak_name || '—'}</span>
            </td>
            <td>
              <span class="context-badge \${row.category === 'AccessibleTSS' ? 'ctx-promoter' : (row.category === 'FirePeak' ? 'ctx-genebody' : 'ctx-intergenic')}">\${row.category || 'ClusterPeak'}</span>
            </td>
            <td>
              <div class="gene-title">\${row.nearest_gene || '—'}</div>
              <div>\${distBadge}</div>
            </td>
            <td>
              <span class="context-badge \${ctxClass}">\${row.context}</span>
            </td>
            <td>
              <div class="reads-val">\${row.total_reads}</div>
              <div class="reads-bar"><div class="reads-fill" style="width: \${barPct}%;"></div></div>
            </td>
            <td>\${kBadge}</td>
            <td>
              <span class="openness-pct \${openClass}">\${pctOpenVal}%</span>
            </td>
            <td>
              <span style="font-family: var(--font-mono); font-size: 0.85rem;">\${nrlStr}</span>
            </td>
            <td>
              <div class="tag-list">\${tagsHtml}</div>
            </td>
            <td style="text-align: right;">
              <a href="\${row.report_url}" class="action-btn" onclick="event.stopPropagation();">
                Report &rarr;
              </a>
            </td>
          </tr>
        `;
      }
      tableBody.innerHTML = html;

      // Render pagination
      pageSummary.textContent = `Page \${state.currentPage} of \${maxPages}`;
      let pageHtml = `
        <button class="page-btn" \${state.currentPage === 1 ? 'disabled' : ''} onclick="changePage(1)">&laquo;</button>
        <button class="page-btn" \${state.currentPage === 1 ? 'disabled' : ''} onclick="changePage(\${state.currentPage - 1})">&lsaquo;</button>
      `;

      // Page numbers (at most 5 buttons around current page)
      let pStart = Math.max(1, state.currentPage - 2);
      let pEnd = Math.min(maxPages, pStart + 4);
      if (pEnd - pStart < 4) pStart = Math.max(1, pEnd - 4);

      for (let p = pStart; p <= pEnd; p++) {
        pageHtml += `
          <button class="page-btn \${p === state.currentPage ? 'active' : ''}" onclick="changePage(\${p})">\${p}</button>
        `;
      }

      pageHtml += `
        <button class="page-btn" \${state.currentPage === maxPages ? 'disabled' : ''} onclick="changePage(\${state.currentPage + 1})">&rsaquo;</button>
        <button class="page-btn" \${state.currentPage === maxPages ? 'disabled' : ''} onclick="changePage(\${maxPages})">&raquo;</button>
      `;
      pageButtons.innerHTML = pageHtml;
    }

    function changePage(p) {
      state.currentPage = p;
      render();
      window.scrollTo({ top: 380, behavior: 'smooth' });
    }

    // --- Event Listeners ---
    searchInput.addEventListener('input', (e) => {
      state.searchQuery = e.target.value.trim();
      state.currentPage = 1;
      render();
    });

    sortSelect.addEventListener('change', (e) => {
      state.sortBy = e.target.value;
      render();
    });

    pageSizeSelect.addEventListener('change', (e) => {
      state.pageSize = parseInt(e.target.value, 10);
      state.currentPage = 1;
      render();
    });

    // Context Pills
    document.getElementById('contextPills').addEventListener('click', (e) => {
      if (!e.target.classList.contains('pill-btn')) return;
      document.querySelectorAll('#contextPills .pill-btn').forEach(b => b.classList.remove('active'));
      e.target.classList.add('active');
      state.contextFilter = e.target.dataset.context;
      state.currentPage = 1;
      render();
    });

    // Category Pills
    categoryPills.addEventListener('click', (e) => {
      if (!e.target.classList.contains('pill-btn')) return;
      document.querySelectorAll('#categoryPills .pill-btn').forEach(b => b.classList.remove('active'));
      e.target.classList.add('active');
      state.categoryFilter = e.target.dataset.category;
      state.currentPage = 1;
      render();
    });

    // K State Pills
    document.getElementById('kPills').addEventListener('click', (e) => {
      if (!e.target.classList.contains('pill-btn')) return;
      document.querySelectorAll('#kPills .pill-btn').forEach(b => b.classList.remove('active'));
      e.target.classList.add('active');
      state.kFilter = e.target.dataset.k;
      state.currentPage = 1;
      render();
    });

    // Toggle Chips
    function setupToggle(id, prop) {
      const btn = document.getElementById(id);
      btn.addEventListener('click', () => {
        state[prop] = !state[prop];
        btn.classList.toggle('active', state[prop]);
        state.currentPage = 1;
        render();
      });
    }
    setupToggle('toggleOpen', 'toggleOpen');
    setupToggle('toggleHap', 'toggleHap');
    setupToggle('toggleStrand', 'toggleStrand');
    setupToggle('toggleSample', 'toggleSample');
    setupToggle('toggleCpG', 'toggleCpG');

    // Column Header Click Sorting
    const colSortMap = {
      coords: 'coords',
      peak_name: 'peak_name',
      category: 'category',
      gene: 'tss_asc',
      context: 'tss_asc',
      reads: 'reads_desc',
      k: 'k_desc',
      open: 'open_desc',
      nrl: 'nrl_asc'
    };
    document.querySelectorAll('#regionsTable th[data-col]').forEach(th => {
      th.addEventListener('click', () => {
        const col = th.dataset.col;
        if (col === 'reads' && state.sortBy === 'reads_desc') {
          state.sortBy = 'reads_asc';
        } else if (colSortMap[col]) {
          state.sortBy = colSortMap[col];
        }
        if ([...sortSelect.options].some(o => o.value === state.sortBy)) {
          sortSelect.value = state.sortBy;
        }
        state.currentPage = 1;
        render();
      });
    });

    // Initial render
    render();
  </script>
</body>
</html>
    """
end


function _html_escape(s::AbstractString)
    s1 = replace(s, "&" => "&amp;")
    s2 = replace(s1, "<" => "&lt;")
    return replace(s2, ">" => "&gt;")
end

function _github_heading_slug(s::AbstractString)
    plain = replace(s, r"<[^>]+>" => "")
    plain = replace(plain, "&amp;" => "", "&lt;" => "", "&gt;" => "", "&quot;" => "", "&#39;" => "")
    slug = lowercase(strip(plain))
    slug = replace(slug, r"[^\w\s-]" => "")
    slug = replace(slug, r"\s" => "-")
    return slug
end

function _add_markdown_heading_ids(html_str::AbstractString)
    return replace(html_str, r"<h([1-4])>(.*?)</h\1>"s => function(m)
        mat = match(r"^<h([1-4])>(.*?)</h\1>$"s, m)
        isnothing(mat) && return m
        lvl = mat.captures[1]
        inner = mat.captures[2]
        slug = _github_heading_slug(inner)
        return "<h$(lvl) id=\"$(slug)\">$(inner)</h$(lvl)>"
    end)
end

function _markdown_to_html_with_math(md_str::AbstractString)
    math_blocks = Pair{String, String}[]
    idx = 0

    # 1. Extract display math $$...$$ BEFORE Markdown.parse so backslashes (\\), <, >, _ are untouched
    s1 = replace(md_str, r"\$\$(.+?)\$\$"s => function(m)
        mat = match(r"^\$\$(.+?)\$\$$"s, m)
        isnothing(mat) && return m
        idx += 1
        token = "SMGMATHDISPLAYTOKEN$(idx)END"
        eq = _html_escape(strip(mat.captures[1]))
        push!(math_blocks, token => "<div class=\"math-block\">\\[\n$(eq)\n\\]</div>")
        return token
    end)

    # 2. Extract inline math $...$ BEFORE Markdown.parse (only outside fenced/inline code blocks)
    s2 = replace(s1, r"(```[\s\S]*?```|`[^`\n]+`)|\$([^\$\n]+?)\$" => function(m)
        startswith(m, "`") && return m
        mat = match(r"^\$([^\$\n]+?)\$$", m)
        isnothing(mat) && return m
        idx += 1
        token = "SMGMATHINLINETOKEN$(idx)END"
        eq = _html_escape(strip(mat.captures[1]))
        push!(math_blocks, token => "<span class=\"math-inline\">\\($(eq)\\)</span>")
        return token
    end)

    html_str = sprint(Markdown.html, Markdown.parse(s2))
    for (token, replacement) in math_blocks
        html_str = replace(html_str, token => replacement)
    end
    return html_str
end

const _CLUSTER_GUIDE_MD_DEFAULT_PATH = joinpath(@__DIR__, "..", "docs", "cluster_guide.md")
include_dependency(_CLUSTER_GUIDE_MD_DEFAULT_PATH)
const _CLUSTER_GUIDE_MD_EMBEDDED = isfile(_CLUSTER_GUIDE_MD_DEFAULT_PATH) ? read(_CLUSTER_GUIDE_MD_DEFAULT_PATH, String) : ""

"""
    render_cluster_guide_html(outpath; guide_md_path=pkgdir(SMGStats, "docs", "cluster_guide.md"), back_url="index.html", back_label="Back to Catalog Dashboard")

Converts `docs/cluster_guide.md` into a standalone, styled HTML page (`guide.html`) using
Julia's built-in `Markdown` standard library and lightweight KaTeX CDN scripts to render
inline (`\\(...\\)`) and centered display (`\\[...\\]`) mathematical equations, tables, and diagrams.
"""
function render_cluster_guide_html(outpath::AbstractString;
    guide_md_path::AbstractString = _CLUSTER_GUIDE_MD_DEFAULT_PATH,
    back_url::Union{Nothing, AbstractString} = "index.html",
    back_label::AbstractString = "Back to Catalog Dashboard"
)
    md_text = if isfile(guide_md_path)
        read(guide_md_path, String)
    elseif !isempty(_CLUSTER_GUIDE_MD_EMBEDDED)
        _CLUSTER_GUIDE_MD_EMBEDDED
    else
        "# Comprehensive Guide to `smgstats cluster`\n\nDocumentation file not found at `$(guide_md_path)`."
    end
    raw_body_html = _markdown_to_html_with_math(md_text)
    body_html = _add_markdown_heading_ids(raw_body_html)

    html_doc = """
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>SMGStats Clustering Methodology &amp; Guide</title>
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&family=JetBrains+Mono:wght@400;500;600&display=swap" rel="stylesheet">
  <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/katex.min.css" crossorigin="anonymous">
  <script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/katex.min.js" crossorigin="anonymous"></script>
  <script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/contrib/auto-render.min.js" crossorigin="anonymous"
    onload="renderMathInElement(document.body, {delimiters: [{left: '\\\\[', right: '\\\\]', display: true}, {left: '\$\$', right: '\$\$', display: true}, {left: '\\\\(', right: '\\\\)', display: false}, {left: '\$', right: '\$', display: false}], throwOnError: false});"></script>
  <style>
    :root {
      --bg-main: #f8fafc;
      --bg-surface: #ffffff;
      --border-subtle: #e2e8f0;
      --text-main: #0f172a;
      --text-muted: #475569;
      --text-dim: #64748b;
      --accent-blue: #0284c7;
      --accent-indigo: #4f46e5;
      --radius-md: 12px;
      --font-sans: 'Inter', -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
      --font-mono: 'JetBrains Mono', SFMono-Regular, Menlo, monospace;
    }
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body {
      font-family: var(--font-sans);
      background-color: var(--bg-main);
      color: var(--text-main);
      line-height: 1.65;
      padding-bottom: 5rem;
      -webkit-font-smoothing: antialiased;
    }
    .top-nav {
      border-bottom: 1px solid var(--border-subtle);
      background: rgba(255, 255, 255, 0.92);
      backdrop-filter: blur(12px);
      position: sticky;
      top: 0;
      z-index: 50;
    }
    .top-nav-inner {
      max-width: 1140px;
      margin: 0 auto;
      padding: 0 1.5rem;
      display: flex;
      align-items: center;
      justify-content: space-between;
      height: 60px;
    }
    .brand {
      display: flex;
      align-items: center;
      gap: 0.65rem;
      font-weight: 700;
      font-size: 1.05rem;
      color: var(--text-main);
      text-decoration: none;
    }
    .brand-tag {
      background: linear-gradient(135deg, var(--accent-blue), var(--accent-indigo));
      -webkit-background-clip: text;
      -webkit-text-fill-color: transparent;
      font-weight: 800;
    }
    .brand-badge {
      font-size: 0.7rem;
      font-weight: 600;
      text-transform: uppercase;
      letter-spacing: 0.05em;
      background: #eef2ff;
      border: 1px solid #c7d2fe;
      color: #4f46e5;
      padding: 0.18rem 0.5rem;
      border-radius: 6px;
    }
    .back-btn {
      display: inline-flex;
      align-items: center;
      gap: 0.4rem;
      color: var(--accent-blue);
      text-decoration: none;
      font-size: 0.85rem;
      font-weight: 600;
      padding: 0.35rem 0.8rem;
      background: #f0f9ff;
      border: 1px solid #bae6fd;
      border-radius: 6px;
      transition: all 0.15s ease;
    }
    .back-btn:hover {
      background: var(--accent-blue);
      color: #ffffff;
    }
    .guide-container {
      max-width: 1140px;
      margin: 2rem auto;
      padding: 0 1.5rem;
    }
    .guide-card {
      background: var(--bg-surface);
      border: 1px solid var(--border-subtle);
      border-radius: var(--radius-md);
      padding: 2.5rem 3rem;
      box-shadow: 0 1px 3px rgba(0, 0, 0, 0.04);
    }
    .guide-card h1 {
      font-size: 2rem;
      font-weight: 800;
      letter-spacing: -0.025em;
      margin-bottom: 1rem;
      color: var(--text-main);
    }
    .guide-card h2 {
      font-size: 1.45rem;
      font-weight: 700;
      letter-spacing: -0.015em;
      margin-top: 2.25rem;
      margin-bottom: 0.85rem;
      padding-bottom: 0.4rem;
      border-bottom: 1px solid var(--border-subtle);
      color: var(--text-main);
    }
    .guide-card h3 {
      font-size: 1.15rem;
      font-weight: 700;
      margin-top: 1.75rem;
      margin-bottom: 0.6rem;
      color: var(--text-main);
    }
    .guide-card h4 {
      font-size: 0.98rem;
      font-weight: 600;
      margin-top: 1.25rem;
      margin-bottom: 0.45rem;
      color: var(--text-main);
    }
    .guide-card p {
      margin-bottom: 0.95rem;
      color: #1e293b;
      font-size: 0.95rem;
    }
    .guide-card ul, .guide-card ol {
      margin: 0.6rem 0 1.1rem 1.5rem;
      color: #1e293b;
      font-size: 0.95rem;
    }
    .guide-card li {
      margin-bottom: 0.35rem;
    }
    .guide-card hr {
      border: none;
      border-top: 1px solid var(--border-subtle);
      margin: 2rem 0;
    }
    .guide-card a {
      color: var(--accent-blue);
      text-decoration: none;
      font-weight: 500;
    }
    .guide-card a:hover {
      text-decoration: underline;
    }
    .guide-card code {
      font-family: var(--font-mono);
      font-size: 0.84em;
      background: #f1f5f9;
      color: #0f172a;
      padding: 0.15rem 0.38rem;
      border-radius: 4px;
      border: 1px solid #e2e8f0;
    }
    .guide-card pre {
      background: #0f172a;
      color: #f8fafc;
      padding: 1.1rem 1.25rem;
      border-radius: 8px;
      overflow-x: auto;
      margin: 1rem 0 1.25rem 0;
      line-height: 1.45;
    }
    .guide-card pre code {
      background: transparent;
      color: inherit;
      padding: 0;
      border: none;
      font-size: 0.82rem;
    }
    .guide-card table {
      width: 100%;
      border-collapse: collapse;
      margin: 1rem 0 1.5rem 0;
      font-size: 0.85rem;
      display: block;
      overflow-x: auto;
    }
    .guide-card th {
      background: #f8fafc;
      color: var(--text-muted);
      font-weight: 600;
      text-transform: uppercase;
      font-size: 0.73rem;
      letter-spacing: 0.04em;
      padding: 0.6rem 0.85rem;
      border-bottom: 2px solid var(--border-subtle);
      text-align: left;
      white-space: nowrap;
    }
    .guide-card td {
      padding: 0.55rem 0.85rem;
      border-bottom: 1px solid #f1f5f9;
      vertical-align: top;
      color: #1e293b;
    }
    .guide-card tr:hover td {
      background: #f8fafc;
    }
    .math-block {
      display: block;
      text-align: center;
      margin: 1.15rem 0;
      padding: 0.5rem 0.75rem;
      overflow-x: auto;
      overflow-y: hidden;
    }
    .katex-display {
      margin: 1.1rem 0;
      overflow-x: auto;
      overflow-y: hidden;
      padding: 0.25rem 0;
    }
  </style>
</head>
<body>
  <nav class="top-nav">
    <div class="top-nav-inner">
      <div class="brand">
        <span class="brand-tag">SMGStats</span>
        <span>Chromatin Clustering</span>
        <span class="brand-badge">Methodology &amp; Guide</span>
      </div>
      $(!isnothing(back_url) ? """<a href="$(back_url)" class="back-btn">&larr; $(back_label)</a>""" : "")
    </div>
  </nav>
  <main class="guide-container">
    <article class="guide-card">
      $(body_html)
    </article>
  </main>
</body>
</html>
"""
    mkpath(dirname(abspath(outpath)))
    write(outpath, html_doc)
    return outpath
end
