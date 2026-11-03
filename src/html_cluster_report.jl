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
    px_per_unit=2
)
    mkpath(asset_dir)
    
    # 1. Summary: Cluster sizes, ΔBIC, Schwarz weights (SVG @ 1600x350)
    f_sum = plotclustersummary(clustering; figure=(; size=(1600, 350)))
    summary_path = joinpath(asset_dir, "summary.svg")
    save(summary_path, f_sum)
    
    # 2. Gene Model track (separate figure, spanning ±15kb flanking context with dual loc/peakloc shading)
    has_genes = false
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
    
    # 5. NRL: Two separate figures to prevent any collision between the boxplot and the AC grid
    k = clustering.optimal_clustering.optimal_k
    f_nrl_box = Figure(size=(550, 340))
    boxplot_nrl(clustering; f=f_nrl_box, axx=1, axy=1, direction=:vertical)
    save(joinpath(asset_dir, "nrl_boxplot.svg"), f_nrl_box)

    f_nrl_spec = Figure(size=(1100, max(260, 200 * cld(k, 4))))
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
    plot_5mC_volcano(clustering; f=f_5mc, axx=1, axy=1)
    plot_5mC_diff_heat(clustering; f=f_5mc, axx=3, axy=1)
    
    sig_df = filter(r -> r.pvalue < 0.05, clustering.modclustertable)
    sig_pos = unique(sig_df.pos)
    levels = sortperm(sortperm(clustering.optimal_clustering.optimal_model.assignments))
    if !isempty(sig_pos)
        ax_mod = Axis(f_5mc[1, 5], title="Single-Molecule 5mC (Diff Sites in Red)", xgridvisible=false, ygridvisible=false)
        # Background: all CpGs in subtle gray
        modplot(clustering.chrom, clustering.loc, clustering.reads, levels; 
                mod=:mod_5mC, f=f_5mc, ax=ax_mod, markersize=3, plotline=false, color=:grey75)
        # Foreground: differential CpGs in vivid crimson
        modplot(clustering.chrom, clustering.loc, clustering.reads, levels; 
                mod=:mod_5mC, f=f_5mc, ax=ax_mod, validpos=(Set(sig_pos), Set{Int}()), markersize=4.5, plotline=false, color=:crimson)
        clusterlines!(clustering, ax=ax_mod, plotlines=true, setticks=true, color=:black)
    else
        ax_mod = Axis(f_5mc[1, 5], title="Single-Molecule 5mC (No sig sites)")
        modplot(clustering.chrom, clustering.loc, clustering.reads, levels; 
                mod=:mod_5mC, f=f_5mc, ax=ax_mod, markersize=3, plotline=false, color=:grey75)
        clusterlines!(clustering, ax=ax_mod, plotlines=true, setticks=true, color=:black)
    end
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
    label="", 
    filename="report.html", 
    catalog_path=nothing
)
    update_clustering_nrl!(clustering; phase_threshold=phase_threshold)
    mkpath(outdir)
    asset_dir = joinpath(outdir, "assets")
    assets = save_cluster_report_assets(clustering, asset_dir; 
        exondata=exondata, 
        longest_only=longest_only, 
        gene_flank=gene_flank,
        phase_threshold=phase_threshold
    )
    
    chrom = clustering.chrom
    peak_start = first(clustering.peakloc)
    peak_stop  = last(clustering.peakloc)
    loc_start  = first(clustering.loc)
    loc_stop   = last(clustering.loc)
    flank_bp   = peak_start - loc_start
    total_reads = length(clustering.reads)
    opt_k      = clustering.optimal_clustering.optimal_k
    
    # Annotate locus context if exondata provided
    annot = if !isnothing(exondata)
        annotate_locus(exondata, chrom, clustering.loc)
    else
        (; nearest_gene="Unannotated", nearest_gene_id="", gene_type="", tss_distance=0, context="Genomic Region", transcripts=TranscriptModel[])
    end
    
    # Total reads per cluster
    cluster_counts = [count(clustering.readtable.cluster .== k) for k in 1:opt_k]
    
    # Check for biologically notable flags
    has_open_cluster = any(clustering.clustertable.proportion_has_fire .> 0.25)
    has_closed_cluster = any(clustering.clustertable.is_closed)
    
    # Strand Fisher tests per cluster
    has_strand_skew = false
    strand_table_rows = []
    if opt_k > 1 && hasproperty(clustering.readtable, :strand) && length(unique(clustering.readtable.strand)) > 1
        for k in 1:opt_k
            is_k = clustering.readtable.cluster .== k
            is_pos = clustering.readtable.strand .== true
            ft = fisher(is_k, is_pos)
            n_pos = count(is_k .& is_pos)
            n_neg = count(is_k .& .!is_pos)
            tot_k = cluster_counts[k]
            pct_pos = tot_k > 0 ? round(n_pos / tot_k * 100, digits=1) : 0.0
            pv = ft.fisher_bpv
            if !isnan(pv) && pv < 0.01
                has_strand_skew = true
            end
            or_str = isinf(ft.or) ? "Inf" : string(round(ft.or, digits=2))
            ci_low = isinf(ft.logor_se) ? "—" : string(round(exp(log(max(1e-6, ft.or)) - 1.96 * ft.logor_se), digits=2))
            ci_high = isinf(ft.logor_se) ? "—" : string(round(exp(log(max(1e-6, ft.or)) + 1.96 * ft.logor_se), digits=2))
            push!(strand_table_rows, (;
                cluster = k,
                tot = tot_k,
                pos = n_pos,
                neg = n_neg,
                pct_pos = pct_pos,
                or = or_str,
                ci = "$ci_low – $ci_high",
                pv = isnan(pv) ? "—" : string(round(pv, sigdigits=3)),
                is_skew = (!isnan(pv) && pv < 0.01)
            ))
        end
    end
    
    # Haplotype Fisher tests per cluster
    has_hap_skew = false
    hap_table_rows = []
    if opt_k > 1 && hasproperty(clustering.readtable, :hap) && any(clustering.readtable.hap .== 1) && any(clustering.readtable.hap .== 2)
        hap_mask = (clustering.readtable.hap .== 1) .| (clustering.readtable.hap .== 2)
        if count(hap_mask) >= 10
            sub_clust = clustering.readtable.cluster[hap_mask]
            sub_hap = clustering.readtable.hap[hap_mask] .== 1
            for k in 1:opt_k
                is_k = sub_clust .== k
                ft = fisher(is_k, sub_hap)
                tot_k = cluster_counts[k]
                n_h1 = count(clustering.readtable.cluster .== k .&& clustering.readtable.hap .== 1)
                n_h2 = count(clustering.readtable.cluster .== k .&& clustering.readtable.hap .== 2)
                n_h0 = count(clustering.readtable.cluster .== k .&& clustering.readtable.hap .== 0)
                phased_tot = n_h1 + n_h2
                pct_h1 = phased_tot > 0 ? round(n_h1 / phased_tot * 100, digits=1) : 0.0
                pv = ft.fisher_bpv
                if !isnan(pv) && pv < 0.01
                    has_hap_skew = true
                end
                or_str = isinf(ft.or) ? "Inf" : string(round(ft.or, digits=2))
                ci_low = isinf(ft.logor_se) ? "—" : string(round(exp(log(max(1e-6, ft.or)) - 1.96 * ft.logor_se), digits=2))
                ci_high = isinf(ft.logor_se) ? "—" : string(round(exp(log(max(1e-6, ft.or)) + 1.96 * ft.logor_se), digits=2))
                push!(hap_table_rows, (;
                    cluster = k,
                    tot = tot_k,
                    h1 = n_h1,
                    h2 = n_h2,
                    h0 = n_h0,
                    pct_h1 = pct_h1,
                    or = or_str,
                    ci = "$ci_low – $ci_high",
                    pv = isnan(pv) ? "—" : string(round(pv, sigdigits=3)),
                    is_skew = (!isnan(pv) && pv < 0.01)
                ))
            end
        end
    end
    
    # Sample Fisher tests per cluster
    has_sample_skew = false
    sample_table_rows = []
    samples = hasproperty(clustering.readtable, :sample) ? unique(clustering.readtable.sample) : String[]
    if opt_k > 1 && length(samples) >= 2
        for k in 1:opt_k
            best_pv = 1.0
            best_s = first(samples)
            best_ft = nothing
            best_n_s = 0
            
            for s in samples
                is_k = clustering.readtable.cluster .== k
                is_s = clustering.readtable.sample .== s
                ft = fisher(is_k, is_s)
                if !isnan(ft.fisher_bpv) && ft.fisher_bpv <= best_pv
                    best_pv = ft.fisher_bpv
                    best_s = s
                    best_ft = ft
                    best_n_s = count(is_k .& is_s)
                end
            end
            
            if !isnothing(best_ft)
                pv = best_ft.fisher_bpv
                if !isnan(pv) && pv < 0.01
                    has_sample_skew = true
                end
                tot_k = cluster_counts[k]
                pct_s = tot_k > 0 ? round(best_n_s / tot_k * 100, digits=1) : 0.0
                or_str = isinf(best_ft.or) ? "Inf" : string(round(best_ft.or, digits=2))
                ci_low = isinf(best_ft.logor_se) ? "—" : string(round(exp(log(max(1e-6, best_ft.or)) - 1.96 * best_ft.logor_se), digits=2))
                ci_high = isinf(best_ft.logor_se) ? "—" : string(round(exp(log(max(1e-6, best_ft.or)) + 1.96 * best_ft.logor_se), digits=2))
                
                push!(sample_table_rows, (;
                    cluster = k,
                    tot = tot_k,
                    top_sample = string(best_s),
                    count_sample = best_n_s,
                    pct_sample = pct_s,
                    or = or_str,
                    ci = "$ci_low – $ci_high",
                    pv = isnan(pv) ? "—" : string(round(pv, sigdigits=3)),
                    is_skew = (!isnan(pv) && pv < 0.01)
                ))
            end
        end
    end
    
    # 1. Format compact cluster properties table HTML
    ct_buf = IOBuffer()
    println(ct_buf, """
    <div style="overflow-x: auto;">
    <table class="compact-table">
    <thead>
      <tr>
        <th>Cluster</th>
        <th>Reads</th>
        <th>Mean Cost</th>
        <th>Length (bp)</th>
        <th>FIRE (%)</th>
        <th>6mA (%)</th>
        <th>5mC (%)</th>
        <th>FIRE NRL (bp)</th>
        <th>6mA NRL (bp)</th>
        <th>Phasing</th>
        <th>Status</th>
        <th>Strand Bias</th>
        <th>Allelic Bias</th>
        <th>Sample Bias</th>
      </tr>
    </thead>
    <tbody>
    """)
    for (k, row) in enumerate(eachrow(clustering.clustertable))
        fire_nrl_str = ismissing(row.nrl) ? "—" : string(row.nrl)
        nrl_6mA_val = hasproperty(row, :nrl_6mA) ? row.nrl_6mA : missing
        nrl_6mA_str = ismissing(nrl_6mA_val) ? "—" : string(nrl_6mA_val)
        
        status_pill = row.is_closed ? 
            "<span class='status-pill status-closed'>Periodic</span>" : 
            "<span class='status-pill status-open'>Open</span>"
            
        strand_badge = "—"
        if !isempty(strand_table_rows) && k <= length(strand_table_rows)
            st = strand_table_rows[k]
            strand_badge = st.is_skew ? "<span class='status-pill status-alert'>p=$(st.pv)</span>" : "<span class='status-dim'>p=$(st.pv)</span>"
        elseif hasproperty(row, :strand_pvalue) && !isnan(row.strand_pvalue)
            pv = row.strand_pvalue
            strand_badge = (pv < 0.01) ? "<span class='status-pill status-alert'>p=$(round(pv, sigdigits=3))</span>" : "<span class='status-dim'>p=$(round(pv, sigdigits=3))</span>"
        end

        hap_badge = "—"
        if !isempty(hap_table_rows) && k <= length(hap_table_rows)
            ht = hap_table_rows[k]
            hap_badge = ht.is_skew ? "<span class='status-pill status-alert'>p=$(ht.pv)</span>" : "<span class='status-dim'>p=$(ht.pv)</span>"
        elseif hasproperty(row, :hap_pvalue) && !isnan(row.hap_pvalue)
            pv = row.hap_pvalue
            hap_badge = (pv < 0.01) ? "<span class='status-pill status-alert'>p=$(round(pv, sigdigits=3))</span>" : "<span class='status-dim'>p=$(round(pv, sigdigits=3))</span>"
        end

        sample_badge = "—"
        if !isempty(sample_table_rows) && k <= length(sample_table_rows)
            s_t = sample_table_rows[k]
            sample_badge = s_t.is_skew ? "<span class='status-pill status-alert'>p=$(s_t.pv)</span>" : "<span class='status-dim'>p=$(s_t.pv)</span>"
        elseif hasproperty(row, :sample_pvalue) && !isnan(row.sample_pvalue)
            pv = row.sample_pvalue
            sample_badge = (pv < 0.01) ? "<span class='status-pill status-alert'>p=$(round(pv, sigdigits=3))</span>" : "<span class='status-dim'>p=$(round(pv, sigdigits=3))</span>"
        end

        c_reads = hasproperty(row, :totalreads) ? row.totalreads : (hasproperty(row, :count) ? row.count : (k <= length(cluster_counts) ? cluster_counts[k] : row.counts_sample))

        phasing_fire = round(row.phasing_score, digits=3)
        phasing_str = string(phasing_fire)
        if hasproperty(row, :phasing_score_6mA) && !isnan(row.phasing_score_6mA)
            phasing_str *= " / " * string(round(row.phasing_score_6mA, digits=3))
        end

        println(ct_buf, """
        <tr>
          <td><strong>$(row.cluster)</strong></td>
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
        </tr>
        """)
    end
    println(ct_buf, "</tbody></table></div>")
    clustertable_html = String(take!(ct_buf))
    
    # 2. Format Strand Fisher table HTML
    strand_table_html = if !isempty(strand_table_rows)
        sbuf = IOBuffer()
        println(sbuf, """
        <div style="overflow-x: auto; margin-top: 1rem;">
        <table class="compact-table">
        <thead>
          <tr>
            <th>Cluster</th>
            <th>Total Reads</th>
            <th>(+) Strand</th>
            <th>(-) Strand</th>
            <th>% (+) Strand</th>
            <th>Odds Ratio</th>
            <th>95% CI</th>
            <th>Fisher p-value</th>
            <th>Assessment</th>
          </tr>
        </thead>
        <tbody>
        """)
        for r in strand_table_rows
            eval_pill = r.is_skew ? 
                "<span class='status-pill status-alert'>Strand Skew (p < 0.01)</span>" : 
                "<span class='status-pill status-open'>Balanced</span>"
            println(sbuf, """
            <tr>
              <td><strong>$(r.cluster)</strong></td>
              <td>$(r.tot)</td>
              <td>$(r.pos)</td>
              <td>$(r.neg)</td>
              <td>$(r.pct_pos)%</td>
              <td>$(r.or)</td>
              <td>$(r.ci)</td>
              <td>$(r.pv)</td>
              <td>$eval_pill</td>
            </tr>
            """)
        end
        println(sbuf, "</tbody></table></div>")
        String(take!(sbuf))
    else
        "<p class='text-dim'><em>Strand metadata not available or locus has only 1 cluster.</em></p>"
    end

    # 3. Format Haplotype Fisher table HTML
    hap_table_html = if !isempty(hap_table_rows)
        hbuf = IOBuffer()
        println(hbuf, """
        <div style="overflow-x: auto; margin-top: 1rem;">
        <table class="compact-table">
        <thead>
          <tr>
            <th>Cluster</th>
            <th>Total Reads</th>
            <th>Hap 1</th>
            <th>Hap 2</th>
            <th>Unphased (0)</th>
            <th>Hap 1 / (Hap 1 + 2)</th>
            <th>Odds Ratio</th>
            <th>95% CI</th>
            <th>Fisher p-value</th>
            <th>Assessment</th>
          </tr>
        </thead>
        <tbody>
        """)
        for r in hap_table_rows
            eval_pill = r.is_skew ? 
                "<span class='status-pill status-alert'>Allelic Skew (p < 0.01)</span>" : 
                "<span class='status-pill status-open'>Biallelic</span>"
            println(hbuf, """
            <tr>
              <td><strong>$(r.cluster)</strong></td>
              <td>$(r.tot)</td>
              <td>$(r.h1)</td>
              <td>$(r.h2)</td>
              <td>$(r.h0)</td>
              <td>$(r.pct_h1)%</td>
              <td>$(r.or)</td>
              <td>$(r.ci)</td>
              <td>$(r.pv)</td>
              <td>$eval_pill</td>
            </tr>
            """)
        end
        println(hbuf, "</tbody></table></div>")
        String(take!(hbuf))
    else
        "<p class='text-dim'><em>Haplotype phase information not available or insufficient heterozygous SNVs.</em></p>"
    end

    # 4. Format Sample Representation Fisher table HTML
    sample_table_html = if !isempty(sample_table_rows)
        sbuf = IOBuffer()
        println(sbuf, """
        <div style="overflow-x: auto; margin-top: 1rem;">
        <table class="compact-table">
        <thead>
          <tr>
            <th>Cluster</th>
            <th>Total Reads</th>
            <th>Top Enriched Sample</th>
            <th>Sample Reads</th>
            <th>Sample Share (%)</th>
            <th>Odds Ratio</th>
            <th>95% CI</th>
            <th>Fisher p-value</th>
            <th>Assessment</th>
          </tr>
        </thead>
        <tbody>
        """)
        for r in sample_table_rows
            eval_pill = r.is_skew ? 
                "<span class='status-pill status-alert'>Sample Bias (p < 0.01)</span>" : 
                "<span class='status-pill status-open'>Balanced</span>"
            println(sbuf, """
            <tr>
              <td><strong>$(r.cluster)</strong></td>
              <td>$(r.tot)</td>
              <td><code>$(r.top_sample)</code></td>
              <td>$(r.count_sample)</td>
              <td>$(r.pct_sample)%</td>
              <td>$(r.or)</td>
              <td>$(r.ci)</td>
              <td>$(r.pv)</td>
              <td>$eval_pill</td>
            </tr>
            """)
        end
        println(sbuf, "</tbody></table></div>")
        String(take!(sbuf))
    else
        "<p class='text-dim'><em>Sample representation metadata not available or locus has only 1 cluster / 1 sample.</em></p>"
    end

    # Format top differentially methylated CpGs table
    mct_buf = IOBuffer()
    sig_cpgs = filter(r -> r.pvalue < 0.05, clustering.modclustertable)
    if nrow(sig_cpgs) > 0
        println(mct_buf, """
        <div style="overflow-x: auto;">
        <table class="compact-table">
        <thead>
          <tr>
            <th>Pos</th>
            <th>Cluster</th>
            <th>Cluster Rate</th>
            <th>Background Rate</th>
            <th>Coverage</th>
            <th>p-value (Binomial)</th>
          </tr>
        </thead>
        <tbody>
        """)
        for row in eachrow(sig_cpgs[1:min(15, nrow(sig_cpgs)), :])
            cov = hasproperty(row, :total) ? row.total : (hasproperty(row, :depth) ? row.depth : (row.mod + row.unmod))
            println(mct_buf, """
            <tr>
              <td><code>$(row.pos)</code></td>
              <td>$(row.cluster)</td>
              <td>$(round(row.rate * 100, digits=1))%</td>
              <td>$(round(row.total_rate * 100, digits=1))%</td>
              <td>$cov</td>
              <td>$(string(round(row.pvalue, sigdigits=3)))</td>
            </tr>
            """)
        end
        println(mct_buf, "</tbody></table></div>")
        if nrow(sig_cpgs) > 15
            println(mct_buf, "<small style='color: var(--text-dim);'>Showing top 15 of $(nrow(sig_cpgs)) significant CpG sites (p < 0.05)</small>")
        end
    else
        println(mct_buf, "<p class='text-dim'><em>No single-CpG sites reached p < 0.05 differential methylation.</em></p>")
    end
    cpg_table_html = String(take!(mct_buf))
    
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
          $(!isnothing(catalog_path) ? """<div style="margin-bottom: 0.85rem;"><a href="$(catalog_path)" style="display: inline-flex; align-items: center; gap: 0.4rem; color: #0284c7; text-decoration: none; font-size: 0.85rem; font-weight: 600; padding: 0.3rem 0.75rem; background: #f0f9ff; border: 1px solid #bae6fd; border-radius: 6px; transition: all 0.15s ease;">&larr; Back to Catalog Dashboard</a></div>""" : "")
          <h1>Single-Molecule Cluster Report</h1>
          <p>Locus: <strong>$chrom:$loc_start-$loc_stop</strong> (Peak: <code>$chrom:$peak_start-$peak_stop</code>, flank ±$(flank_bp) bp)</p>
          <div class="badge-bar">
            <span class="badge badge-k">Optimal k* = $opt_k</span>
            <span class="badge badge-reads">$total_reads Spanning Reads</span>
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
          <p class="section-desc">Single-molecule nucleosome repeat length distributions alongside cluster centroid autocorrelation wave fits (FIRE and 6mA channels).</p>
          <div class="two-col-grid">
            <div>
              <div style="font-size: 0.85rem; font-weight: 600; color: var(--text-muted); margin-bottom: 0.35rem;">Single-Molecule NRL Distribution</div>
              <div class="figure-wrapper">
                <img src="$(assets.nrl_box)" alt="NRL Boxplot per Cluster" onclick="openLightbox(this.src)">
              </div>
            </div>
            <div>
              <div style="font-size: 0.85rem; font-weight: 600; color: var(--text-muted); margin-bottom: 0.35rem;">Spectral Autocorrelation & Harmonic Cosine Fits (FIRE & 6mA)</div>
              <div class="figure-wrapper">
                <img src="$(assets.nrl_spec)" alt="Centroid Spectral NRL Fits" onclick="openLightbox(this.src)">
              </div>
            </div>
          </div>
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
