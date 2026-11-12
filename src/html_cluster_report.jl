#=
    html_cluster_report.jl: Standalone HTML report generator for single-molecule chromatin clustering
=#

"""
    save_cluster_report_assets(clustering, asset_dir; exondata=nothing, longest_only=true, px_per_unit=2)

Generates and saves the hybrid visual suite (PNGs for dense single-molecule tracks, SVGs for statistical panels)
into `asset_dir`. Returns a NamedTuple of relative asset filenames.
"""
function save_cluster_report_assets(clustering, asset_dir; 
    exondata=nothing, 
    longest_only=true, 
    gene_flank=15_000, 
    phase_threshold=0.035, 
    px_per_unit=2,
    only_sample::Bool=false
)
    mkpath(asset_dir)
    
    has_genes = false
    if only_sample && isfile(joinpath(asset_dir, "tracks.png"))
        has_genes = isfile(joinpath(asset_dir, "genemodels.png"))
    else
        # 1. Summary: Cluster sizes, ΔBIC, Schwarz weights (SVG @ 1600x350)
        f_sum = plotclustersummary(clustering; figure=(; size=(1600, 350)))
        summary_path = joinpath(asset_dir, "summary.svg")
        save(summary_path, f_sum)
        
        # 2. Gene Model track (separate figure, spanning ±15kb flanking context with dual loc/peakloc shading)
        if !isnothing(exondata)
            f_gene = plot_genemodel_track(clustering, exondata; gene_flank=gene_flank, longest_only=longest_only)
            if !isnothing(f_gene)
                has_genes = true
                save(joinpath(asset_dir, "genemodels.png"), f_gene, px_per_unit=px_per_unit)
            end
        end

        # 3. Dense Single-Molecule Tracks: Feature Heatmap + FIRE/Nuc + 6mA + 5mC (PNG @ 2x Retina)
        f_tracks = Figure(size=(1600, 650))
        plot_cluster_tracks(clustering; f=f_tracks)
        tracks_path = joinpath(asset_dir, "tracks.png")
        save(tracks_path, f_tracks, px_per_unit=px_per_unit)
        
        # 4. Centroid Profiles: FIRE vs 6mA, Nuc vs 6mA, 5mC vs FIRE (PNG @ 2x)
        f_prof = Figure(size=(1600, 750))
        plot_cluster_profiles(clustering; f=f_prof)
        profiles_path = joinpath(asset_dir, "profiles.png")
        save(profiles_path, f_prof, px_per_unit=px_per_unit)
        
        # 5. NRL: Two separate full-width figures (Row 1: side-by-side FIRE & 6mA read NRL boxplots; Row 2: full-width AC grid)
        k = clustering.optimal_clustering.optimal_k
        f_nrl_box = Figure(size=(1500, 320))
        boxplot_nrl(clustering; f=f_nrl_box, axx=1, axy=1, direction=:vertical, dual=true)
        save(joinpath(asset_dir, "nrl_boxplot.svg"), f_nrl_box)

        f_nrl_spec = Figure(size=(1500, max(280, 230 * cld(k, 4) + 45)))
        plot_spectral_nrl(clustering; f=f_nrl_spec, nrl_feature=:fire, phase_threshold=phase_threshold)
        save(joinpath(asset_dir, "nrl_spectral.svg"), f_nrl_spec)
        
        # 6. Uncertainty / Entropy (SVG)
        f_ent = Figure(size=(1300, 320))
        ax_ent1 = Axis(f_ent[1, 1], title="Population Uncertainty vs k")
        plot_mean_assignment_uncertainty(clustering; f=f_ent, ax=ax_ent1)
        boxplot_read_uncertainty(clustering; f=f_ent, axx=2:3, axy=1, direction=:vertical)
        entropy_path = joinpath(asset_dir, "entropy.svg")
        save(entropy_path, f_ent)
        
        # 7. 5mC Epigenetics: Volcano + Differential Heatmap + Single-Molecule 5mC at Diff Sites (SVG)
        f_5mc = Figure(size=(1600, 360))
        plot_5mC_volcano(clustering; f=f_5mc, axx=1, axy=1, pv_thresh=0.05)
        plot_5mC_diff_heat(clustering; f=f_5mc, axx=3, axy=1, pv_thresh=0.05)
        
        sig_df = filter(r -> r.pvalue < 0.05, clustering.modclustertable)
        sig_pos = unique(sig_df.pos)
        levels = read_plot_levels(clustering)
        ax_mod = Axis(f_5mc[1, 5], title=isempty(sig_pos) ? "Single-Molecule 5mC (No sig sites)" : "Single-Molecule 5mC (Diff Sites in Red)", xgridvisible=false, ygridvisible=false)
        modplot(clustering.chrom, clustering.loc, clustering.reads, levels; 
                mod=:mod_5mC, f=f_5mc, ax=ax_mod, markersize=3, plotline=false, color=:grey75)
        if !isempty(sig_pos)
            modplot(clustering.chrom, clustering.loc, clustering.reads, levels; 
                    mod=:mod_5mC, f=f_5mc, ax=ax_mod, validpos=(Set(sig_pos), Set{Int}()), markersize=4.5, plotline=false, color=:crimson)
        end
        clusterlines!(clustering, ax=ax_mod, plotlines=true, setticks=true, color=:black)
        cpg_path = joinpath(asset_dir, "5mC.svg")
        save(cpg_path, f_5mc)
        
        # 8. Strand Diagnostics (SVG)
        f_strand = strand_diagnostics(clustering)
        strand_path = joinpath(asset_dir, "strand.svg")
        save(strand_path, f_strand)
        
        # 9. Haplotype Diagnostics (SVG)
        f_hap = hap_diagnostics(clustering)
        hap_path = joinpath(asset_dir, "haplotype.svg")
        save(hap_path, f_hap)
    end
    
    # 10. Sample Representation Diagnostics (SVG)
    f_sample = sample_diagnostics(clustering)
    sample_path = joinpath(asset_dir, "sample.svg")
    save(sample_path, f_sample)
    
    return (;
        summary = "assets/summary.svg",
        genemodels = has_genes ? "assets/genemodels.png" : nothing,
        tracks = "assets/tracks.png",
        profiles = "assets/profiles.png",
        nrl_box = "assets/nrl_boxplot.svg",
        nrl_spec = "assets/nrl_spectral.svg",
        entropy = "assets/entropy.svg",
        cpg = "assets/5mC.svg",
        strand = "assets/strand.svg",
        hap = "assets/haplotype.svg",
        sample = "assets/sample.svg"
    )
end

function _format_or_ci_pv(ft; alpha::Float64=0.01)
    pv = ft.fisher_bpv
    is_skew = !isnan(pv) && pv < alpha
    or_str = isinf(ft.or) ? "Inf" : (isnan(ft.or) ? "—" : string(round(ft.or, digits=2)))
    ci_str = if !isfinite(ft.logor_se) || !isfinite(ft.or)
        "—"
    else
        log_or = log(max(1e-6, ft.or))
        ci_low  = round(exp(log_or - 1.96 * ft.logor_se), digits=2)
        ci_high = round(exp(log_or + 1.96 * ft.logor_se), digits=2)
        "$ci_low – $ci_high"
    end
    pv_str = isnan(pv) ? "—" : string(round(pv, sigdigits=3))
    return (; or=or_str, ci=ci_str, pv=pv_str, is_skew=is_skew)
end

@inline _pval_badge(pv_str::AbstractString, is_skew::Bool) =
    is_skew ? "<span class='status-pill status-alert'>p=$pv_str</span>" : "<span class='status-dim'>p=$pv_str</span>"

@inline function _pval_badge(pv::Real; alpha::Float64=0.01)
    isnan(pv) && return "—"
    return _pval_badge(string(round(pv, sigdigits=3)), pv < alpha)
end

function _html_compact_table(headers::AbstractVector{<:AbstractString}, rows_html::AbstractVector{<:AbstractString}; 
                             empty_msg::AbstractString="", margin_top::Bool=false, footer_html::AbstractString="")
    isempty(rows_html) && !isempty(empty_msg) && return "<p class='text-dim'><em>$empty_msg</em></p>"
    io = IOBuffer()
    style = margin_top ? "overflow-x: auto; margin-top: 1rem;" : "overflow-x: auto;"
    println(io, "<div style=\"$style\"><table class=\"compact-table\"><thead><tr>")
    for h in headers
        println(io, "  <th>$h</th>")
    end
    println(io, "</tr></thead><tbody>")
    for r in rows_html
        println(io, r)
    end
    println(io, "</tbody></table></div>")
    !isempty(footer_html) && println(io, footer_html)
    return String(take!(io))
end

"""
    render_cluster_report(clustering, outdir; exondata=nothing, longest_only=true, label="", filename="report.html", catalog_path=nothing)

Generates a complete, self-contained interactive HTML report for a clustered locus.
Includes metadata badges (gene context, TSS distance, open chromatin, allelic skew),
pre-rendered high-resolution figures with zoom lightbox, and sortable compact tables.
"""
function render_cluster_report(clustering, outdir; 
    exondata=nothing, 
    longest_only=true, 
    gene_flank=15_000,
    phase_threshold=0.035,
    peakname::Union{Nothing, AbstractString}=nothing,
    category::Union{Nothing, AbstractString}=nothing,
    label="", 
    filename="report.html", 
    catalog_path=nothing,
    guide_path=nothing,
    update_nrl::Bool=true,
    samples=nothing,
    only_sample::Bool=false
)
    update_nrl && update_clustering_nrl!(clustering; phase_threshold=phase_threshold, samples=samples)
    resolved_peakname, resolved_category = _resolve_peak_meta(clustering; peakname=peakname, category=category)
    mkpath(outdir)
    resolved_guide_path = if !isnothing(guide_path)
        String(guide_path)
    elseif !isnothing(catalog_path)
        "../../guide.html"
    else
        render_cluster_guide_html(joinpath(outdir, "guide.html"); back_url=filename, back_label="Back to Locus Report")
        "guide.html"
    end
    asset_dir = joinpath(outdir, "assets")
    assets = save_cluster_report_assets(clustering, asset_dir; 
        exondata=exondata, 
        longest_only=longest_only, 
        gene_flank=gene_flank,
        phase_threshold=phase_threshold,
        only_sample=only_sample
    )
    
    chrom = clustering.chrom
    peak_start = first(clustering.peakloc)
    peak_stop  = last(clustering.peakloc)
    loc_start  = first(clustering.loc)
    loc_stop   = last(clustering.loc)
    flank_bp   = peak_start - loc_start
    total_reads = length(clustering.reads)
    opt_k      = clustering.optimal_clustering.optimal_k
    rt         = clustering.readtable
    
    # Annotate target peak context if exondata provided
    annot = if !isnothing(exondata)
        annotate_locus(exondata, chrom, clustering.peakloc)
    else
        (; nearest_gene="Unannotated", nearest_gene_id="", gene_type="", tss_distance=0, context="Genomic Region", transcripts=TranscriptModel[])
    end
    
    # Total reads per cluster
    cluster_counts = [count(rt.cluster .== k) for k in 1:opt_k]
    
    # Check for biologically notable flags
    has_open_cluster = any(clustering.clustertable.proportion_has_fire .> 0.25)
    has_closed_cluster = any(clustering.clustertable.is_closed)
    
    # Strand Fisher tests per cluster
    has_strand_skew = false
    strand_by_cluster = Dict{Int, NamedTuple}()
    if opt_k > 1 && hasproperty(rt, :strand) && length(unique(rt.strand)) > 1
        is_pos = rt.strand .== true
        for k in 1:opt_k
            is_k = rt.cluster .== k
            fmt = _format_or_ci_pv(fisher(is_k, is_pos))
            n_pos = count(is_k .& is_pos)
            n_neg = count(is_k .& .!is_pos)
            tot_k = cluster_counts[k]
            pct_pos = tot_k > 0 ? round(n_pos / tot_k * 100, digits=1) : 0.0
            has_strand_skew |= fmt.is_skew
            strand_by_cluster[k] = (; cluster=k, tot=tot_k, pos=n_pos, neg=n_neg, pct_pos=pct_pos, fmt...)
        end
    end
    
    # Haplotype Fisher tests per cluster
    has_hap_skew = false
    hap_by_cluster = Dict{Int, NamedTuple}()
    if opt_k > 1 && hasproperty(rt, :hap) && any(rt.hap .== 1) && any(rt.hap .== 2)
        hap_mask = (rt.hap .== 1) .| (rt.hap .== 2)
        if count(hap_mask) >= 10
            sub_clust = rt.cluster[hap_mask]
            sub_hap = rt.hap[hap_mask] .== 1
            for k in 1:opt_k
                fmt = _format_or_ci_pv(fisher(sub_clust .== k, sub_hap))
                tot_k = cluster_counts[k]
                n_h1 = count((rt.cluster .== k) .& (rt.hap .== 1))
                n_h2 = count((rt.cluster .== k) .& (rt.hap .== 2))
                n_h0 = count((rt.cluster .== k) .& (rt.hap .== 0))
                phased_tot = n_h1 + n_h2
                pct_h1 = phased_tot > 0 ? round(n_h1 / phased_tot * 100, digits=1) : 0.0
                has_hap_skew |= fmt.is_skew
                hap_by_cluster[k] = (; cluster=k, tot=tot_k, h1=n_h1, h2=n_h2, h0=n_h0, pct_h1=pct_h1, fmt...)
            end
        end
    end
    
    # Sample Fisher tests per cluster
    has_sample_skew = false
    sample_by_cluster = Dict{Int, NamedTuple}()
    samples = hasproperty(rt, :sample) ? unique(rt.sample) : String[]
    if opt_k > 1 && length(samples) >= 2
        for k in 1:opt_k
            is_k = rt.cluster .== k
            best_pv = 1.0
            best_s = first(samples)
            best_ft = nothing
            best_n_s = 0
            
            for s in samples
                is_s = rt.sample .== s
                ft = fisher(is_k, is_s)
                if !isnan(ft.fisher_bpv) && ft.fisher_bpv <= best_pv
                    best_pv = ft.fisher_bpv
                    best_s = s
                    best_ft = ft
                    best_n_s = count(is_k .& is_s)
                end
            end
            
            if !isnothing(best_ft)
                fmt = _format_or_ci_pv(best_ft)
                tot_k = cluster_counts[k]
                pct_s = tot_k > 0 ? round(best_n_s / tot_k * 100, digits=1) : 0.0
                has_sample_skew |= fmt.is_skew
                sample_by_cluster[k] = (; cluster=k, tot=tot_k, top_sample=string(best_s), count_sample=best_n_s, pct_sample=pct_s, fmt...)
            end
        end
    end
    
    # 1. Format compact cluster properties table HTML & NRL phasing summary table HTML
    ct_rows = String[]
    nrl_rows_html = String[]
    for row in eachrow(clustering.clustertable)
        cid = row.cluster
        fire_nrl_str = ismissing(row.nrl) ? "—" : string(row.nrl)
        nrl_6mA_val = hasproperty(row, :nrl_6mA) ? row.nrl_6mA : missing
        nrl_6mA_str = ismissing(nrl_6mA_val) ? "—" : string(nrl_6mA_val)
        
        is_closed_6mA = hasproperty(row, :is_closed_6mA) && row.is_closed_6mA === true
        status_pill = row.is_closed ? 
            "<span class='status-pill status-closed'>Periodic</span>" : 
            "<span class='status-pill status-open'>Open</span>"
        status_pill_6mA = is_closed_6mA ?
            "<span class='status-pill status-closed'>Periodic</span>" :
            "<span class='status-pill status-open'>Open</span>"
            
        strand_badge = haskey(strand_by_cluster, cid) ? 
            _pval_badge(strand_by_cluster[cid].pv, strand_by_cluster[cid].is_skew) :
            (hasproperty(row, :strand_pvalue) ? _pval_badge(row.strand_pvalue) : "—")

        hap_badge = haskey(hap_by_cluster, cid) ? 
            _pval_badge(hap_by_cluster[cid].pv, hap_by_cluster[cid].is_skew) :
            (hasproperty(row, :hap_pvalue) ? _pval_badge(row.hap_pvalue) : "—")

        sample_badge = haskey(sample_by_cluster, cid) ? 
            _pval_badge(sample_by_cluster[cid].pv, sample_by_cluster[cid].is_skew) :
            (hasproperty(row, :sample_pvalue) ? _pval_badge(row.sample_pvalue) : "—")

        c_reads = hasproperty(row, :totalreads) ? row.totalreads : (hasproperty(row, :count) ? row.count : (1 <= cid <= length(cluster_counts) ? cluster_counts[cid] : row.counts_sample))

        phasing_fire_str = isfinite(row.phasing_score) ? string(round(row.phasing_score, digits=3)) : "—"
        phasing_6mA_str  = (hasproperty(row, :phasing_score_6mA) && isfinite(row.phasing_score_6mA)) ? string(round(row.phasing_score_6mA, digits=3)) : "—"
        phasing_str = phasing_6mA_str == "—" ? phasing_fire_str : "$phasing_fire_str / $phasing_6mA_str"

        peak_fire_str = (hasproperty(row, :peak_type) && row.peak_type != :none) ? string(row.peak_type) : "—"
        peak_6mA_str  = (hasproperty(row, :peak_type_6mA) && row.peak_type_6mA != :none) ? string(row.peak_type_6mA) : "—"

        mean_r_fire = (hasproperty(row, :mean_read_nrl) && !ismissing(row.mean_read_nrl) && isfinite(row.mean_read_nrl)) ? "$(round(row.mean_read_nrl, digits=1)) bp" : "—"
        n_r_fire    = hasproperty(row, :total_valid_nrls) ? row.total_valid_nrls : 0
        mean_r_6mA  = (hasproperty(row, :mean_read_nrl_6mA) && !ismissing(row.mean_read_nrl_6mA) && isfinite(row.mean_read_nrl_6mA)) ? "$(round(row.mean_read_nrl_6mA, digits=1)) bp" : "—"
        n_r_6mA     = hasproperty(row, :total_valid_nrls_6mA) ? row.total_valid_nrls_6mA : 0

        push!(ct_rows, """
        <tr>
          <td><strong>$cid</strong></td>
          <td><strong>$c_reads</strong></td>
          <td>$(round(row.meancost, digits=2))</td>
          <td>$(round(Int, row.mean_readlength))</td>
          <td>$(round(row.mean_rate_fire_bp * 100, digits=1))%</td>
          <td>$(round(row.mean_rate_mod_6mA * 100, digits=1))%</td>
          <td>$(round(row.mean_rate_mod_5mC * 100, digits=1))%</td>
          <td><strong>$fire_nrl_str</strong></td>
          <td><strong>$nrl_6mA_str</strong></td>
          <td>$phasing_str</td>
          <td>$status_pill</td>
          <td>$strand_badge</td>
          <td>$hap_badge</td>
          <td>$sample_badge</td>
        </tr>""")

        push!(nrl_rows_html, """
        <tr>
          <td><strong>$cid</strong></td>
          <td><strong>$c_reads</strong></td>
          <td><strong>$fire_nrl_str</strong></td>
          <td>$phasing_fire_str</td>
          <td>$peak_fire_str</td>
          <td>$status_pill</td>
          <td><strong>$nrl_6mA_str</strong></td>
          <td>$phasing_6mA_str</td>
          <td>$peak_6mA_str</td>
          <td>$status_pill_6mA</td>
          <td>$mean_r_fire <span class="text-dim">(n=$n_r_fire)</span></td>
          <td>$mean_r_6mA <span class="text-dim">(n=$n_r_6mA)</span></td>
        </tr>""")
    end
    clustertable_html = _html_compact_table(
        ["Cluster", "Reads", "Mean Cost", "Length (bp)", "FIRE (%)", "6mA (%)", "5mC (%)", "FIRE NRL (bp)", "6mA NRL (bp)", "Phasing (FIRE / 6mA)", "Status", "Strand Bias", "Allelic Bias", "Sample Bias"],
        ct_rows
    )
    nrl_table_html = _html_compact_table(
        ["Cluster", "Reads", "FIRE NRL (bp)", "FIRE Phasing", "FIRE Harmonic", "FIRE Call", "6mA NRL (bp)", "6mA Phasing", "6mA Harmonic", "6mA Call", "Read Mean NRL (FIRE)", "Read Mean NRL (6mA)"],
        nrl_rows_html;
        margin_top=true
    )
    
    # 2. Format Strand Fisher table HTML
    strand_rows_html = [
        begin
            r = strand_by_cluster[k]
            eval_pill = r.is_skew ? "<span class='status-pill status-alert'>Strand Skew (p < 0.01)</span>" : "<span class='status-pill status-open'>Balanced</span>"
            "<tr><td><strong>$(r.cluster)</strong></td><td>$(r.tot)</td><td>$(r.pos)</td><td>$(r.neg)</td><td>$(r.pct_pos)%</td><td>$(r.or)</td><td>$(r.ci)</td><td>$(r.pv)</td><td>$eval_pill</td></tr>"
        end for k in sort!(collect(keys(strand_by_cluster)))
    ]
    strand_table_html = _html_compact_table(
        ["Cluster", "Total Reads", "(+) Strand", "(-) Strand", "% (+) Strand", "Odds Ratio", "95% CI", "Fisher p-value", "Assessment"],
        strand_rows_html;
        empty_msg="Strand metadata not available or locus has only 1 cluster.",
        margin_top=true
    )

    # 3. Format Haplotype Fisher table HTML
    hap_rows_html = [
        begin
            r = hap_by_cluster[k]
            eval_pill = r.is_skew ? "<span class='status-pill status-alert'>Allelic Skew (p < 0.01)</span>" : "<span class='status-pill status-open'>Biallelic</span>"
            "<tr><td><strong>$(r.cluster)</strong></td><td>$(r.tot)</td><td>$(r.h1)</td><td>$(r.h2)</td><td>$(r.h0)</td><td>$(r.pct_h1)%</td><td>$(r.or)</td><td>$(r.ci)</td><td>$(r.pv)</td><td>$eval_pill</td></tr>"
        end for k in sort!(collect(keys(hap_by_cluster)))
    ]
    hap_table_html = _html_compact_table(
        ["Cluster", "Total Reads", "Hap 1", "Hap 2", "Unphased (0)", "Hap 1 / (Hap 1 + 2)", "Odds Ratio", "95% CI", "Fisher p-value", "Assessment"],
        hap_rows_html;
        empty_msg="Haplotype phase information not available or insufficient heterozygous SNVs.",
        margin_top=true
    )

    # 4. Format Sample Representation Fisher table HTML
    sample_rows_html = [
        begin
            r = sample_by_cluster[k]
            eval_pill = r.is_skew ? "<span class='status-pill status-alert'>Sample Bias (p < 0.01)</span>" : "<span class='status-pill status-open'>Balanced</span>"
            "<tr><td><strong>$(r.cluster)</strong></td><td>$(r.tot)</td><td><code>$(r.top_sample)</code></td><td>$(r.count_sample)</td><td>$(r.pct_sample)%</td><td>$(r.or)</td><td>$(r.ci)</td><td>$(r.pv)</td><td>$eval_pill</td></tr>"
        end for k in sort!(collect(keys(sample_by_cluster)))
    ]
    sample_table_html = _html_compact_table(
        ["Cluster", "Total Reads", "Top Enriched Sample", "Sample Reads", "Sample Share (%)", "Odds Ratio", "95% CI", "Fisher p-value", "Assessment"],
        sample_rows_html;
        empty_msg="Sample representation metadata not available or locus has only 1 cluster / 1 sample.",
        margin_top=true
    )

    # 5. Format top differentially methylated CpGs table (sorted by most significant p-value)
    sig_cpgs = sort(filter(r -> r.pvalue < 0.05, clustering.modclustertable), :pvalue)
    cpg_rows_html = [
        begin
            cov = hasproperty(row, :total) ? row.total : (hasproperty(row, :depth) ? row.depth : (row.mod + row.unmod))
            "<tr><td><code>$(row.pos)</code></td><td>$(row.cluster)</td><td>$(round(row.rate * 100, digits=1))%</td><td>$(round(row.total_rate * 100, digits=1))%</td><td>$cov</td><td>$(string(round(row.pvalue, sigdigits=3)))</td></tr>"
        end for row in eachrow(sig_cpgs[1:min(15, nrow(sig_cpgs)), :])
    ]
    cpg_footer = nrow(sig_cpgs) > 15 ? "<small style='color: var(--text-dim);'>Showing top 15 of $(nrow(sig_cpgs)) significant CpG sites (p < 0.05)</small>" : ""
    cpg_table_html = _html_compact_table(
        ["Pos", "Cluster", "Cluster Rate", "Background Rate", "Coverage", "p-value (Binomial)"],
        cpg_rows_html;
        empty_msg="No single-CpG sites reached p < 0.05 differential methylation.",
        footer_html=cpg_footer
    )
    
    # Distance formatting
    tss_dist_str = if annot.tss_distance < 0
        "$(-annot.tss_distance) bp Upstream of TSS"
    elseif annot.tss_distance > 0
        "$(annot.tss_distance) bp Downstream of TSS"
    else
        "At TSS"
    end

    html_content = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
      <meta charset="utf-8">
      <meta name="viewport" content="width=device-width, initial-scale=1">
      <title>SMGStats Cluster Report: $chrom:$peak_start-$peak_stop</title>
      <link rel="preconnect" href="https://fonts.googleapis.com">
      <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
      <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500;600&display=swap" rel="stylesheet">
      <style>
        :root {
          --bg-body: #f8fafc;
          --bg-card: #ffffff;
          --border-color: #e2e8f0;
          --border-subtle: #f1f5f9;
          --text-main: #0f172a;
          --text-muted: #475569;
          --text-dim: #94a3b8;
          --accent-blue: #0284c7;
          --accent-indigo: #4f46e5;
          --accent-emerald: #059669;
          --accent-amber: #d97706;
          --accent-rose: #e11d48;
          --radius-md: 10px;
          --font-sans: 'Inter', -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
          --font-mono: 'JetBrains Mono', SFMono-Regular, Menlo, monospace;
        }
        * { box-sizing: border-box; margin: 0; padding: 0; }
        body {
          font-family: var(--font-sans);
          background-color: var(--bg-body);
          color: var(--text-main);
          line-height: 1.5;
          padding-bottom: 5rem;
          -webkit-font-smoothing: antialiased;
        }
        .container {
          max-width: 1840px;
          width: 98vw;
          margin: 0 auto;
          padding: 1.25rem 0.5rem;
        }
        .hero-banner {
          background: #ffffff;
          border: 1px solid var(--border-color);
          border-radius: var(--radius-md);
          padding: 1.25rem 1.5rem;
          margin-bottom: 1.25rem;
          box-shadow: 0 1px 3px rgba(0,0,0,0.05);
        }
        .hero-banner h1 { margin-bottom: 0.35rem; color: var(--text-main); font-size: 1.75rem; font-weight: 700; letter-spacing: -0.02em; }
        .hero-banner p { color: var(--text-muted); margin-bottom: 1rem; font-size: 0.92rem; }
        .badge-bar { display: flex; flex-wrap: wrap; gap: 0.45rem; align-items: center; }
        .badge {
          display: inline-flex;
          align-items: center;
          padding: 0.2rem 0.6rem;
          border-radius: 9999px;
          font-size: 0.76rem;
          font-weight: 600;
          line-height: 1;
        }
        .badge-k     { background: #ede9fe; color: #6d28d9; border: 1px solid #ddd6fe; }
        .badge-reads { background: #e0f2fe; color: #0369a1; border: 1px solid #bae6fd; }
        .badge-gene  { background: #ecfdf5; color: #047857; border: 1px solid #a7f3d0; }
        .badge-ctx   { background: #f1f5f9; color: #334155; border: 1px solid #e2e8f0; }
        .badge-tss   { background: #fef3c7; color: #92400e; border: 1px solid #fde68a; }
        .badge-open  { background: #dcfce7; color: #15803d; border: 1px solid #86efac; }
        .badge-alert { background: #fee2e2; color: #b91c1c; border: 1px solid #fca5a5; }

        .report-section {
          background: var(--bg-card);
          border: 1px solid var(--border-color);
          border-radius: var(--radius-md);
          padding: 1.25rem 1.5rem;
          margin-bottom: 1.5rem;
          box-shadow: 0 1px 3px rgba(0,0,0,0.04);
        }
        .report-section h3 {
          font-size: 1.2rem;
          font-weight: 700;
          color: var(--text-main);
          margin-bottom: 0.4rem;
          display: flex;
          align-items: center;
          gap: 0.5rem;
        }
        .report-section h4 {
          font-size: 0.95rem;
          font-weight: 600;
          color: var(--text-main);
          margin-top: 1.25rem;
          margin-bottom: 0.65rem;
        }
        .section-desc {
          color: var(--text-muted);
          font-size: 0.85rem;
          margin-bottom: 1rem;
        }
        .figure-wrapper {
          overflow-x: auto;
          background: #ffffff;
          border: 1px solid var(--border-color);
          border-radius: 8px;
          padding: 0.5rem;
          margin-bottom: 1rem;
        }
        .figure-wrapper img {
          max-width: 100%;
          height: auto;
          border-radius: 4px;
          cursor: zoom-in;
          display: block;
          margin: 0 auto;
        }

        /* Compact Tables */
        .compact-table {
          width: 100%;
          border-collapse: collapse;
          font-size: 0.78rem;
          text-align: left;
        }
        .compact-table th {
          background: #f8fafc;
          color: var(--text-muted);
          font-weight: 600;
          text-transform: uppercase;
          font-size: 0.70rem;
          letter-spacing: 0.04em;
          padding: 0.35rem 0.55rem;
          border-bottom: 1px solid var(--border-color);
          white-space: nowrap;
        }
        .compact-table td {
          padding: 0.35rem 0.55rem;
          border-bottom: 1px solid var(--border-subtle);
          color: var(--text-main);
          white-space: nowrap;
        }
        .compact-table tbody tr:hover td {
          background: #f1f5f9;
        }
        
        .status-pill {
          display: inline-block;
          padding: 0.1rem 0.35rem;
          border-radius: 3px;
          font-size: 0.70rem;
          font-weight: 600;
        }
        .status-closed { background: #ede9fe; color: #6d28d9; }
        .status-open   { background: #dcfce7; color: #15803d; }
        .status-alert  { background: #fee2e2; color: #b91c1c; font-weight: 700; }
        .status-dim    { color: var(--text-dim); }
        .text-dim      { color: var(--text-dim); }

        .two-col-grid {
          display: grid;
          grid-template-columns: minmax(320px, 1fr) minmax(450px, 2fr);
          gap: 1.5rem;
          align-items: start;
        }

        /* Lightbox Modal */
        .lightbox-modal {
          display: none;
          position: fixed;
          z-index: 9999;
          left: 0; top: 0;
          width: 100%; height: 100%;
          background: rgba(15, 23, 42, 0.85);
          backdrop-filter: blur(4px);
          justify-content: center;
          align-items: center;
        }
        .lightbox-modal img {
          max-width: 95%;
          max-height: 95%;
          border-radius: 8px;
          box-shadow: 0 10px 30px rgba(0,0,0,0.5);
        }
      </style>
    </head>
    <body>
      <main class="container">
        
        <!-- Header Banner -->
        <header class="hero-banner">
          <div style="display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap; gap: 0.5rem; margin-bottom: 0.85rem;">
            <div>
              $(!isnothing(catalog_path) ? """<a href="$(catalog_path)" style="display: inline-flex; align-items: center; gap: 0.4rem; color: #0284c7; text-decoration: none; font-size: 0.85rem; font-weight: 600; padding: 0.3rem 0.75rem; background: #f0f9ff; border: 1px solid #bae6fd; border-radius: 6px; transition: all 0.15s ease;">&larr; Back to Catalog Dashboard</a>""" : "")
            </div>
            <div>
              <a href="$(resolved_guide_path)" style="display: inline-flex; align-items: center; gap: 0.4rem; color: #4f46e5; text-decoration: none; font-size: 0.85rem; font-weight: 600; padding: 0.3rem 0.75rem; background: #eef2ff; border: 1px solid #c7d2fe; border-radius: 6px; transition: all 0.15s ease;">📖 Methodology &amp; Guide</a>
            </div>
          </div>
          <h1>Single-Molecule Cluster Report: <code>$resolved_peakname</code></h1>
          <p>Locus: <strong>$chrom:$loc_start-$loc_stop</strong> &bull; Target Peak: <code>$chrom:$peak_start-$peak_stop</code> (flank ±$(flank_bp) bp) &bull; PeakName: <strong>$resolved_peakname</strong> &bull; Category: <strong>$resolved_category</strong></p>
          <div class="badge-bar">
            <span class="badge badge-k">Optimal k* = $opt_k</span>
            <span class="badge badge-reads">$total_reads Spanning Reads</span>
            <span class="badge badge-ctx">Category: $resolved_category</span>
            $(!isempty(annot.nearest_gene) ? "<span class='badge badge-gene'>Gene: $(annot.nearest_gene)</span>" : "")
            $(!isempty(annot.context) ? "<span class='badge badge-ctx'>$(annot.context)</span>" : "")
            $(annot.tss_distance != 0 ? "<span class='badge badge-tss'>$tss_dist_str</span>" : "")
            $(has_open_cluster ? "<span class='badge badge-open'>Accessible Sub-state Detected</span>" : "")
            $(has_strand_skew ? "<span class='badge badge-alert'>Strand Skew (p < 0.01)</span>" : "")
            $(has_hap_skew ? "<span class='badge badge-alert'>Allelic / Hap Skew (p < 0.01)</span>" : "")
            $(has_sample_skew ? "<span class='badge badge-alert'>Sample Skew (p < 0.01)</span>" : "")
          </div>
        </header>

        <!-- 1. Model Selection & Summary -->
        <section class="report-section">
          <h3>1. Model Selection & Cluster Overview</h3>
          <p class="section-desc">Iterated local search (ILS) hill climbing k-means convergence, ΔBIC criterion, and cluster occupancy.</p>
          <div class="figure-wrapper">
            <img src="$(assets.summary)" alt="Cluster Sizes and Model BIC Curve" onclick="openLightbox(this.src)">
          </div>
          <h4>Cluster Properties Summary</h4>
          $clustertable_html
        </section>

        <!-- 2. Single-Molecule Tracks -->
        <section class="report-section">
          <h3>2. Single-Molecule Footprint & Methylation Tracks</h3>
          <p class="section-desc">Synchronized single-molecule channels: Stencil Feature Heatmap (M), FIRE/Linker/Nucleosome Footprints, 6mA Accessibility, and 5mC Methylation.</p>
          
          $(!isnothing(assets.genemodels) ? """
          <div style="margin-bottom: 0.75rem;">
            <div style="font-size: 0.85rem; font-weight: 600; color: var(--text-muted); margin-bottom: 0.35rem;">Gene Models (Shaded box indicates clustered peak target: $chrom:$peak_start-$peak_stop)</div>
            <div class="figure-wrapper">
              <img src="$(assets.genemodels)" alt="Overlapping Gene Models" onclick="openLightbox(this.src)">
            </div>
          </div>
          """ : "")

          <div class="figure-wrapper">
            <img src="$(assets.tracks)" alt="Single-Molecule Tracks" onclick="openLightbox(this.src)">
          </div>
        </section>

        <!-- 3. Centroid Profiles -->
        <section class="report-section">
          <h3>3. Cluster Centroid Mean Profiles</h3>
          <p class="section-desc">Average per-cluster accessibility and methylation profiles highlighting chromatin state signatures.</p>
          <div class="figure-wrapper">
            <img src="$(assets.profiles)" alt="Cluster Centroid Profiles" onclick="openLightbox(this.src)">
          </div>
        </section>

        <!-- 4. Nucleosome Architecture & Spectral NRL -->
        <section class="report-section">
          <h3>4. Nucleosome Phasing & Repeat Length (NRL)</h3>
          <p class="section-desc">Single-molecule nucleosome repeat length distributions (FIRE and 6mA) alongside cluster centroid autocorrelation wave fits and phasing summary metrics.</p>
          <div style="font-size: 0.85rem; font-weight: 600; color: var(--text-muted); margin-bottom: 0.35rem;">Single-Molecule Read-Level NRL Distributions (FIRE vs 6mA)</div>
          <div class="figure-wrapper">
            <img src="$(assets.nrl_box)" alt="Read-Level NRL Boxplots (FIRE & 6mA)" onclick="openLightbox(this.src)">
          </div>
          <div style="font-size: 0.85rem; font-weight: 600; color: var(--text-muted); margin-top: 1rem; margin-bottom: 0.35rem;">Cluster Centroid Spectral Autocorrelation & Harmonic Cosine Fits (FIRE & 6mA)</div>
          <div class="figure-wrapper">
            <img src="$(assets.nrl_spec)" alt="Centroid Spectral NRL Fits" onclick="openLightbox(this.src)">
          </div>
          <h4>Spectral NRL & Phasing Calls (FIRE vs 6mA)</h4>
          $nrl_table_html
        </section>

        <!-- 5. Assignment Uncertainty -->
        <section class="report-section">
          <h3>5. Assignment Uncertainty & Transition Dynamics</h3>
          <p class="section-desc">Boltzmann-Gibbs assignment entropy H across k values and single-read assignment uncertainty distributions.</p>
          <div class="figure-wrapper">
            <img src="$(assets.entropy)" alt="Read Assignment Uncertainty" onclick="openLightbox(this.src)">
          </div>
        </section>

        <!-- 6. 5mC Epigenetics -->
        <section class="report-section">
          <h3>6. Endogenous 5mC CpG Epigenetics</h3>
          <p class="section-desc">Binomial differential 5mC methylation volcano plot, cluster-resolved rate map, and single-molecule differential sites.</p>
          <div class="figure-wrapper">
            <img src="$(assets.cpg)" alt="5mC Methylation Differences" onclick="openLightbox(this.src)">
          </div>
          <h4>Differentially Methylated Sites</h4>
          $cpg_table_html
        </section>

        <!-- 7A. Technical Strand Validation -->
        <section class="report-section">
          <h3>7A. Technical Strand Diagnostics (Sequencing Bias Control)</h3>
          <p class="section-desc">Evaluation of cluster distribution across positive (+) vs negative (-) sequencing strands to detect PCR or sequencing artifacts.</p>
          <div class="figure-wrapper">
            <img src="$(assets.strand)" alt="Strand Bias Diagnostics" onclick="openLightbox(this.src)">
          </div>
          <h4>Strand Fisher's Exact Test Results</h4>
          $strand_table_html
        </section>

        <!-- 7B. Allelic / Haplotype Validation -->
        <section class="report-section">
          <h3>7B. Allele-Specific & Haplotype Diagnostics (Genomic Imprinting Control)</h3>
          <p class="section-desc">Assessment of cluster enrichment across phased parental haplotypes (Hap 1 vs Hap 2) to identify allele-specific chromatin states.</p>
          <div class="figure-wrapper">
            <img src="$(assets.hap)" alt="Haplotype Bias Diagnostics" onclick="openLightbox(this.src)">
          </div>
          <h4>Haplotype Fisher's Exact Test Results</h4>
          $hap_table_html
        </section>

        <!-- 7C. Sample Representation Validation -->
        <section class="report-section">
          <h3>7C. Sample Representation Diagnostics (Sample / Batch Control)</h3>
          <p class="section-desc">Assessment of cluster enrichment across individual sequencing libraries or biological samples to detect batch effects or state-specific sample representation.</p>
          <div class="figure-wrapper">
            <img src="$(assets.sample)" alt="Sample Bias Diagnostics" onclick="openLightbox(this.src)">
          </div>
          <h4>Sample Representation Fisher's Exact Test Results</h4>
          $sample_table_html
        </section>

      </main>

      <!-- Lightbox Modal -->
      <div id="lightboxModal" class="lightbox-modal" onclick="closeLightbox()">
        <img id="lightboxImg" src="">
      </div>

      <script>
        function openLightbox(src) {
          document.getElementById('lightboxImg').src = src;
          document.getElementById('lightboxModal').style.display = 'flex';
        }
        function closeLightbox() {
          document.getElementById('lightboxModal').style.display = 'none';
        }
        document.addEventListener('keydown', function(e) {
          if (e.key === 'Escape') closeLightbox();
        });
      </script>
    </body>
    </html>
    """

    report_path = joinpath(outdir, filename)
    open(report_path, "w") do io
        write(io, html_content)
    end
    
    smglog("Cluster report written: ", report_path)
    return report_path
end
