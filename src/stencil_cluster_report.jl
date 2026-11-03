#=
    stencil_cluster_report.jl: CairoMakie and AlgebraOfGraphics reporting suite 
    for single-molecule chromatin fiber clustering (smgstats cluster)
=#

using AlgebraOfGraphics
using CairoMakie
using DataFrames
using Statistics
using HypothesisTests
using SwarmMakie

## 

"""
    fireplot(chrom, loc, readdata, leveldata; kwargs...)

Plots single-molecule FIRE accessible elements (red), sub-nucleosomal linker patches (dodgerblue),
and nucleosome footprints (lightgrey) as horizontal track polygons across `loc`.
"""
function fireplot(chrom, loc, readdata, leveldata; hap=nothing, strand=nothing, 
                  f=Figure(), axv=1, axh=1, 
                  ax=Axis(f[axv, axh], xgridvisible=false, ygridvisible=false, xticks=([], [])), 
                  ymax=maximum(leveldata), showlegend=true, expand_fire=0.0, rasterize=false)

    if !isnothing(hap) && !isnothing(strand)
        ind = [(r.haplotype == hap) && (r.strand == strand) for r in readdata]
        reads = @view readdata[ind]
        levels = @view leveldata[ind]
    elseif !isnothing(hap)
        hapind = [r.haplotype == hap for r in readdata]
        reads = @view readdata[hapind]
        levels = @view leveldata[hapind]
    elseif !isnothing(strand)
        strandind = [r.strand == strand for r in readdata]
        reads = @view readdata[strandind]
        levels = @view leveldata[strandind]
    else
        reads = readdata
        levels = leveldata
    end

    if isempty(reads)
        hidespines!(ax, :l, :t, :r, :b)
        xlims!(ax, first(loc), last(loc))
        return f
    end

    lr = mapreduce((l, fr) -> [Rect(start, l + (1 - 0.9) / 2, stop - start, 0.9) for (start, stop, q) in fr.msps if q < 0.9 * 255], vcat, levels, reads)
    mr = mapreduce((l, fr) -> [Rect(start - (stop - start) * expand_fire, l + (1 - 0.9) / 2, (stop - start) * (1 + expand_fire), 0.9) for (start, stop, q) in fr.msps if q >= 0.9 * 255], vcat, levels, reads)
    nr = mapreduce((l, fr) -> [Rect(start, l + (1 - 0.6) / 2, stop - start, 0.6) for (start, stop) in fr.nucs], vcat, levels, reads)

    !isempty(nr) && poly!(ax, nr, color=:lightgrey, label="Nucleosome", rasterize=rasterize)
    
    # @show "one nuc"
    # nr = [Rect(fr.nucs[1][1], l + (1 - 0.6) / 2, fr.nucs[end][2] - fr.nucs[1][1], 0.6) for (fr, l) in zip(reads, levels) if !isempty(fr.nucs)]
    # !isempty(nr) && poly!(ax, nr, color=:lightgrey, label="Nucleosome")

    !isempty(lr) && poly!(ax, lr, color=:dodgerblue, label="Linker", rasterize=rasterize)
    !isempty(mr) && poly!(ax, mr, color=:red, label="FIRE", rasterize=rasterize)

    hidespines!(ax, :l, :t, :r, :b)
    ylims!(ax, ymax, 0)
    showlegend && axislegend(ax, framevisible=false, position=:rc)
    xlims!(ax, first(loc), last(loc))
    f
end


function getmods(read, modfield; unify_gc_strand=true)
    pos = getproperty(read, modfield)
    if unify_gc_strand && !read.strand
        return pos .- 1
    else
        return pos
    end
end

"""
    modplot(chrom, loc, readdata, leveldata; kwargs...)

Plots single-base modification calls (6mA, 5mC, 5hmC) as scatter points along read horizontal levels.
"""
function modplot(chrom, loc, readdata, leveldata; mod=:mod_6mA, unmod=false, plotline=false, 
                 unmodinfire=true, color=nothing, hap=nothing, sample=nothing, showline=true, 
                 validpos=(Set{Int}(), Set{Int}()), strand=nothing, 
                 f=Figure(), axv=1, axh=1, 
                 ax=Axis(f[axv, axh], xgridvisible=false, ygridvisible=false, xticks=([], [])), 
                 ymax=maximum(leveldata), markersize=4, markeralpha=1.0, depthfilter=5, unify_gc_strand=true,
                 labelfun=x -> replace(basename(x), r".bam|.cram|.haplotagged" => ""))

    ind = trues(length(readdata))
    if !isnothing(hap)
        ind .&= [r.haplotype == hap for r in readdata]
    end
    if !isnothing(strand)
        ind .&= [r.strand == strand for r in readdata]
    end
    if !isnothing(sample)
        ind .&= [labelfun(r.label) == sample for r in readdata]
    end
    reads = @view readdata[ind]
    levels = @view leveldata[ind]
    
    if isempty(reads)
        hidespines!(ax, :l, :t, :r, :b)
        xlims!(ax, first(loc), last(loc))
        return f
    end

    sym_mod = Symbol(mod)
    if sym_mod == :mod_6mA
        modfield = ifelse(unmod, :unmods_6mA, :mods_6mA)
        ufield   = :unmods_6mA
        c = :steelblue
        unify_gc_strand = false ### prevent subtraction of 1 on negative strand
    elseif sym_mod == :mod_5mC
        modfield = ifelse(unmod, :unmods_5mC, :mods_5mC)
        ufield   = :unmods_5mC
        c = :green
    elseif sym_mod == :mod_5hmC
        modfield = ifelse(unmod, :unmods_5hmC, :mods_5hmC)
        ufield   = :unmods_5hmC
        c = :orange
    else
        error("Mod: $mod not recognised")
    end
    
    if !isnothing(color)
        c = color
    end

    if unmod && unmodinfire
        data, datalevels = filtermodfire(reads, levels, modfield)
    else 
        data = mapreduce(r -> getmods(r, modfield, unify_gc_strand=unify_gc_strand), vcat, reads)
        datalevels = mapreduce((r, l) -> fill(l, length(getproperty(r, modfield))), vcat, reads, levels)
    end
    
    if !isempty(validpos[1])
        ind = (data .∈ Ref(validpos[1])) .| (data .∈ Ref(validpos[2]))
        data = data[ind]
        datalevels = datalevels[ind]
    end

    if depthfilter > 0 && !isempty(data)
        pile = zeros(Int, length(loc))
        undata = mapreduce(r -> getmods(r, ufield; unify_gc_strand=unify_gc_strand), vcat, reads)
        
        for d in data
            p = d - first(loc) + 1
            if 1 <= p <= length(loc)
                pile[p] += 1
            end
        end
        for d in undata
            p = d - first(loc) + 1
            if 1 <= p <= length(loc)
                pile[p] += 1
            end
        end
        
        ind = falses(length(data))
        for (k, d) in enumerate(data)
            p = d - first(loc) + 1
            if (1 <= p <= length(loc)) && (pile[p] >= depthfilter)
                ind[k] = true
            end
        end
        data = data[ind]
        datalevels = datalevels[ind]
    end

    if !isempty(data)
        if unmod
            scatter!(ax, data, datalevels, marker=:rect, color=:orange, markersize=markersize, alpha=markeralpha)
        else
            scatter!(ax, data, datalevels, marker=:rect, color=c, markersize=markersize, alpha=markeralpha)
        end
    end
    
    if plotline
        for (r, l) in zip(reads, levels)
            lines!(ax, [r.lp, r.rp], [l, l], color=c)
        end
    end

    hidespines!(ax, :l, :t, :r, :b)
    ylims!(ax, ymax, 0)
    # axislegend(ax, framevisible=false, position=:rc)
    xlims!(ax, first(loc), last(loc))
    f
end

function filtermodfire(reads, levels, modfield)
    data = Int[]
    datalevels = Int[]
    
    for (r, l) in zip(reads, levels)    
        for (s, e, q) in r.msps
            (q > 0.8*255) || continue
            for p in getproperty(r, modfield)
                if s <= p <= e
                    push!(data, p)
                    push!(datalevels, l)
                end
            end
        end        
    end
    data, datalevels
    
end

function pilefire!(loc, read, fgpile, bgpile)
    for (s, e, q) in read.msps
        if q > 0.9 * 255
            f_ind = intersect(s:e, loc) .- first(loc) .+ 1
            @views fgpile[f_ind] .+= 1
        end
    end
    bgpile .+= 1
end

function pilemsp!(loc, read, fgpile, bgpile)
    for (s, e, q) in read.msps
        if q <= 0.9 * 255
            f_ind = intersect(s:e, loc) .- first(loc) .+ 1
            @views fgpile[f_ind] .+= 1
        end
    end
    bgpile .+= 1
end

function pilenuc!(loc, read, fgpile, bgpile)
    for (s, e) in read.nucs
        f_ind = intersect(s:e, loc) .- first(loc) .+ 1
        @views fgpile[f_ind] .+= 1
    end
    bgpile .+= 1
end

function pilemod!(loc, read, fgpile, bgpile; modfield=:mods_6mA, unmodfield=:unmods_6mA, unmod=false, w=1)
    @inbounds for p in getproperty(read, modfield)
        r = p - first(loc) + 1
        for i = -w:w
            q = r + i
            if 1 <= q <= length(loc)
                !unmod && (fgpile[q] += 1)
                bgpile[q] += 1
            end
        end
    end
    @inbounds for p in getproperty(read, unmodfield)
        r = p - first(loc) + 1
        for i = -w:w
            q = r + i
            if 1 <= q <= length(loc)
                unmod && (fgpile[q] += 1)
                bgpile[q] += 1
            end
        end
    end
end

@inline pilemod_6mA!(loc, read, fgpile, bgpile)  = pilemod!(loc, read, fgpile, bgpile; modfield=:mods_6mA, unmodfield=:unmods_6mA)
@inline pilemod_5mC!(loc, read, fgpile, bgpile)  = pilemod!(loc, read, fgpile, bgpile; modfield=:mods_5mC, unmodfield=:unmods_5mC)
@inline pilemod_5hmC!(loc, read, fgpile, bgpile) = pilemod!(loc, read, fgpile, bgpile; modfield=:mods_5hmC, unmodfield=:unmods_5hmC)

pilelabel(::typeof(pilefire!))     = "FIRE"
pilelabel(::typeof(pilemsp!))      = "MSP/Linker"
pilelabel(::typeof(pilenuc!))      = "Nucleosome"
pilelabel(::typeof(pilemod!))      = "Mod"
pilelabel(::typeof(pilemod_6mA!))  = "6mA"
pilelabel(::typeof(pilemod_5mC!))  = "5mC"
pilelabel(::typeof(pilemod_5hmC!)) = "5hmC"

function clustermean(clustering; pilefun=pilefire!)
    k = clustering.optimal_clustering.optimal_k  
    n = length(clustering.loc)
    fgcenters = zeros(n, k)
    bgcenters = zeros(n, k)
    
    for (c, r) in zip(clustering.optimal_clustering.optimal_model.assignments, clustering.reads)
        pilefun(clustering.loc, r, view(fgcenters, :, c), view(bgcenters, :, c))
    end
    
    # Avoid divide-by-zero
    @inbounds for i in eachindex(bgcenters)
        if bgcenters[i] == 0.0
            bgcenters[i] = 1.0
        end
    end
    return fgcenters ./ bgcenters
end

function clustermeandf(clustering; pilefun=pilefire!)
    cmean = clustermean(clustering, pilefun=pilefun)
    k = clustering.optimal_clustering.optimal_k
    n_loc = length(clustering.loc)
    df = DataFrame(
        Loc = repeat(collect(clustering.loc), k),
        Cluster = repeat(1:k, inner=n_loc),
        ClusterMean = vec(cmean),
        Label = pilelabel(pilefun)
    )
    return df
end

clustermeans(clustering; pilefuns = [pilefire!, pilemsp!, pilenuc!, pilemod_6mA!, pilemod_5mC!]) = mapreduce(f -> clustermeandf(clustering, pilefun=f), vcat, pilefuns)

function plotclustertotals(clustering; f=Figure(), ax=Axis(f[1, 1]))
    k = clustering.optimal_clustering.optimal_k
    counts = clustering.optimal_clustering.optimal_model.counts
    ax.xticks = (1:k)
    barplot!(ax, 1:k, counts, bar_labels=:y, label_formatter=y -> string(Int(y)), color=:steelblue)
    ylims!(ax, 0, 1.15 * maximum(counts))
    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.xlabel = "Cluster"
    ax.ylabel = "Read Count"
    hidespines!(ax, :t, :r)
    return f
end

function plotbics(clustering; f=Figure(), ax=Axis(f[1, 1]), plotdeltabic=true)
    ks = clustering.optimal_clustering.ks
    bics = copy(clustering.optimal_clustering.bics)
    if plotdeltabic
        bics .-= minimum(bics)
    end
    opt_k = clustering.optimal_clustering.optimal_k
    
    scatterlines!(ax, ks, bics, color=:navy, markercolor=:navy, label=plotdeltabic ? "ΔBIC" : "BIC")
    vlines!(ax, opt_k, color=:crimson, linestyle=:dash, label="Optimal k* = $opt_k")
    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.xlabel = "Candidate k"
  
    hidespines!(ax, :t, :r)
    axislegend(ax, framevisible=false, position=:rb)
    return f
end

softmaxtemp(clustering, temp) = (temp == :numreads) ? Float64(length(clustering.reads)) : 1.0
softmaxtemp(clustering, temp::T) where {T <: Number} = Float64(temp)

function plotbicscharwzweights(clustering; f=Figure(), ax=Axis(f[1, 1]), temp=:numreads)
    ks = clustering.optimal_clustering.ks
    bics = clustering.optimal_clustering.bics .- minimum(clustering.optimal_clustering.bics)
    
    T = softmaxtemp(clustering, temp)
    expbics = exp.(-bics ./ (2.0 * T))
    expbics ./= sum(expbics)
    
    barplot!(ax, ks, expbics, color=:teal)
    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.xlabel = "Candidate k"
    ax.ylabel = "Posterior Weight P(k)"
    hidespines!(ax, :t, :r)
    ylims!(ax, 0, 1.05)
    return f
end

"""
    plotclustersummary(clustering; figure=(; size=(1400, 300)))

Renders a 3-panel figure with Cluster Totals, ΔBIC Curve, and Model Posterior Weights.
"""
function plotclustersummary(clustering; figure=(; size=(1600, 350)))
    f = Figure(; figure...)
    ax1 = Axis(f[1, 1], title="Cluster Sizes")
    plotclustertotals(clustering, f=f, ax=ax1)
    ax2 = Axis(f[1, 2], title="Model Selection (BIC)")
    plotbics(clustering, f=f, ax=ax2)
    ax3 = Axis(f[1, 3], title="Posterior Weights P(k)")
    plotbicscharwzweights(clustering, f=f, ax=ax3)
    return f
end

function clusterlines!(clustering; ax, plotlines=true, setticks=true, color=:white)
    counts = clustering.optimal_clustering.optimal_model.counts
    plotlines && hlines!(ax, cumsum(counts), color=color, linewidth=1.5)
    if setticks
        k = clustering.optimal_clustering.optimal_k
        ax.yticks = (cumsum(counts) .- counts ./ 2, string.(1:k))
    end
end

function featureheatmap(clustering; plotlines=true, setticks=true, f=Figure(), ax=Axis(f[1, 1]))
    perm = sortperm(clustering.optimal_clustering.optimal_model.assignments)
    heatmap!(ax, clustering.CM.ldx, 1:size(clustering.CM.M, 2), clustering.CM.M[:, perm], colormap=:viridis)
    clusterlines!(clustering, ax=ax, plotlines=plotlines, setticks=setticks)

    f
end

function plot_mean_assignment_uncertainty(clustering; f=Figure(), ax=Axis(f[1, 1]))
    ks = clustering.optimal_clustering.ks
    mean_uncertainty = exp.(clustering.optimal_clustering.mean_assignment_entropy)
    max_k = ks[argmax(mean_uncertainty)]
    
    scatterlines!(ax, ks, mean_uncertainty, color=:darkorange, label="exp(Mean Entropy)")
    vlines!(ax, max_k, color=:gray30, linestyle=:dash, label="Max Uncertainty (k = $max_k)")
    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.xlabel = "Candidate k"
    ax.ylabel = "Mean Assignment Uncertainty"
    axislegend(ax, framevisible=false, position=:rb)
    return f
end

function boxplot_read_uncertainty(clustering ; f = Figure(), axx=1, axy=1, direction=:horizontal)
    
     plt = data(clustering.readtable) * 
            mapping(:cluster => string, :max_assignment_entropy => exp => "Assignment Uncertainty", color=:cluster => string) * 
            (visual(BoxPlot, orientation = direction, show_outliers=false) + 
             visual(Beeswarm, direction = direction == :horizontal ? :x : :y, algorithm=:pseudorandom, strokecolor=:black, strokewidth=1))

    draw!(f[axy, axx], plt, axis=(; xgridvisible=false, ygridvisible=false, yreversed = direction == :horizontal))
    f
end

function boxplot_nrl(clustering; f=Figure(), axx=1, axy=1, direction=:vertical)
    valid_nrls = filter(r -> !ismissing(r.nrl), clustering.readtable)
    if isempty(valid_nrls)
        ax = Axis(f[axy, axx], title="No valid NRL calls")
        hidespines!(ax, :l, :t, :r, :b)
        ax.xticks = ([], [])
        ax.yticks = ([], [])
        return f
    end
    
     plt = data(valid_nrls) * mapping(:cluster => string, :nrl, color=:cluster => string) * 
          (visual(BoxPlot, orientation = direction, show_outliers=false) + 
           visual(Beeswarm, direction = direction == :horizontal ? :x : :y, algorithm=:pseudorandom, strokecolor=:black, strokewidth=1))

    draw!(f[axy, axx], plt, axis=(; xgridvisible=false, ygridvisible=false, yreversed = direction == :horizontal))
    f
end


function plot_5mC_volcano(clustering; f=Figure(), axx=1, axy=1, weight=0.05, pv_thresh=0.025)
    if nrow(clustering.modclustertable) == 0
        ax = Axis(f[axy, axx], title="No 5mC sites")
        hidespines!(ax, :l, :t, :r, :b)
        return f
    end

    if !hasproperty(clustering.modclustertable, :log2FC)
        clustering.modclustertable.log2FC = log2.(clustering.modclustertable.rate .+ weight) .- log2.(clustering.modclustertable.total_rate .+ weight)
    end
    
    # Non-significant points in gray, significant points colored by cluster
    sig_mask = clustering.modclustertable.pvalue .< pv_thresh
    df_sig = view(clustering.modclustertable, sig_mask, :)
    df_nonsig = view(clustering.modclustertable, .!sig_mask, :)
    
    plt_nonsig = data(df_nonsig) * mapping(:log2FC, :pvalue => (x -> -log10(x))) * visual(Scatter, color=(:gray75, 0.6), markersize=6)
    plt_sig = data(df_sig) * mapping(:log2FC, :pvalue => (x -> -log10(x)), color=:cluster => string) * visual(Scatter, markersize=8)
    
    plt = isempty(df_sig) ? plt_nonsig : (plt_nonsig + plt_sig)
    pd = draw!(f[axy, axx], plt, axis=(; xgridvisible=false, ygridvisible=false, xlabel="log2 Fold-Change (Cluster vs Background)", ylabel="-log10 p-value"))
    if !isempty(df_sig)
        legend!(f[axy, maximum(axx) + 1], pd)
    end
    
    # Guideline markers
    ax_idx = findlast(x -> x isa Axis, f.content)
    if !isnothing(ax_idx)
        ax_obj = f.content[ax_idx]
        hlines!(ax_obj, -log10(pv_thresh), color=:black, linestyle=:dash)
        vlines!(ax_obj, 0, color=:black, linestyle=:dash)
        valid_fcs = filter(!isnan, clustering.modclustertable.log2FC)
        fcmax = isempty(valid_fcs) ? 2.0 : maximum(abs, valid_fcs) * 1.15
        xlims!(ax_obj, -fcmax, fcmax)
        ylims!(ax_obj, 0, nothing)
    end
    return f
end

function plot_5mC_diff_heat(clustering; f=Figure(), axx=1, axy=1, pv_thresh=0.025)
    df = filter(r -> r.pvalue < pv_thresh, clustering.modclustertable)
    if nrow(df) == 0
        ax = Axis(f[axy, axx], title="No significant 5mC sites (p < $pv_thresh)")
        hidespines!(ax, :l, :t, :r, :b)
        return f
    end

    opt_k = clustering.optimal_clustering.optimal_k
    plt = data(df) * 
          mapping(:pos, :cluster, color=(:rate, :total_rate) => ((r, tr) -> log2((r + 0.05) / (tr + 0.05)))) * 
          visual(Scatter, marker=:vline, markersize=25)

    pd = draw!(f[axy, axx], plt, 
               scales(Color=(; colorrange=(-2, 2), colormap=cgrad([:blue, :white, :red]))), 
               axis=(; yticks=1:opt_k, xgridvisible=false, ygridvisible=false, xlabel="Genomic Position", ylabel="Cluster"))
    xlims!(first(clustering.loc), last(clustering.loc))
    colorbar!(f[axy, maximum(axx) + 1], pd)
    return f
end

function strand_diagnostics(clustering)
    opt_k = clustering.optimal_clustering.optimal_k
    f = Figure(size=(1200, 500))
    levels = sortperm(sortperm(clustering.optimal_clustering.optimal_model.assignments))
    
    if opt_k <= 1
        ax_msg = Axis(f[1, 1], xgridvisible=false, ygridvisible=false)
        hidespines!(ax_msg, :t, :r, :b, :l)
        hidexdecorations!(ax_msg)
        hideydecorations!(ax_msg)
        text!(ax_msg, 0.5, 0.5, text="Homogeneous population (k = 1):\nAcross-cluster strand skew test not applicable.", align=(:center, :center), fontsize=13, color=:grey40)
        
        ax_mod = Axis(f[1, 2], xgridvisible=false, ygridvisible=false, title="6mA (+ vs - Strand)")
        modplot(clustering.chrom, clustering.loc, clustering.reads, levels; mod=:mod_6mA, f=f, ax=ax_mod, unmod=false, strand=false, color=:crimson)
        modplot(clustering.chrom, clustering.loc, clustering.reads, levels; mod=:mod_6mA, f=f, ax=ax_mod, unmod=false, strand=true, color=:steelblue)
        clusterlines!(clustering, ax=ax_mod, plotlines=true, setticks=true, color=:black)
        return f
    end

    cldf = [(cluster=k, fisher(clustering.readtable.cluster .== k, clustering.readtable.strand)...) for k = 1:opt_k] |> DataFrame

    ax = Axis(f[1, 1], xgridvisible=false, ygridvisible=false, yticks=1:opt_k, xlabel="log(Odds Ratio) [+ / -]", ylabel="Cluster")
    plt = data(cldf) * 
          (mapping(:or => log, :cluster => nonnumeric, :logor_se => x -> 1.96 * x, color=:fisher_bpv => (x -> -log10(x))) * visual(Errorbars, direction=:x) + 
           mapping(:or => log, :cluster => nonnumeric, color=:fisher_bpv => (x -> -log10(x)), marker=:fisher_bpv => (x -> x < 0.05) => "Fisher p < 0.05") * visual(Scatter))
    d = draw!(ax, plt, scales(Color=(; colormap=:blues)))
    colorbar!(f[1, 2], d, label="-log10(p)")
    legend!(f[1, 3], d)
    vlines!(ax, 0, color=:black, linestyle=:dash)
    hidespines!(ax, :t, :r)
    ax.yreversed = true
    
    valid_ors = filter(!isinf, abs.(log.(cldf.or)) .+ cldf.logor_se .* 1.96)
    valid_ors = filter(!isnan, valid_ors)
    xmax = isempty(valid_ors) ? 2.0 : max(1.0, 1.15 * maximum(valid_ors))
    xlims!(ax, -xmax, xmax)
    
    ax_mod = Axis(f[1, 4], xgridvisible=false, ygridvisible=false, title="6mA (+ vs - Strand)")
    modplot(clustering.chrom, clustering.loc, clustering.reads, levels; mod=:mod_6mA, f=f, ax=ax_mod, unmod=false, strand=false, color=:crimson)
    modplot(clustering.chrom, clustering.loc, clustering.reads, levels; mod=:mod_6mA, f=f, ax=ax_mod, unmod=false, strand=true, color=:steelblue)
    clusterlines!(clustering, ax=ax_mod, plotlines=true, setticks=true, color=:black)
    
    f
end

function hap_diagnostics(clustering)
    opt_k = clustering.optimal_clustering.optimal_k
    f = Figure(size=(1200, 500))
    levels = sortperm(sortperm(clustering.optimal_clustering.optimal_model.assignments))

    if opt_k <= 1
        ax_msg = Axis(f[1, 1], xgridvisible=false, ygridvisible=false)
        hidespines!(ax_msg, :t, :r, :b, :l)
        hidexdecorations!(ax_msg)
        hideydecorations!(ax_msg)
        text!(ax_msg, 0.5, 0.5, text="Homogeneous population (k = 1):\nAcross-cluster haplotype skew test not applicable.", align=(:center, :center), fontsize=13, color=:grey40)
        
        ax_mod = Axis(f[1, 2], xgridvisible=false, ygridvisible=false, title="6mA by Haplotype (0=Unphased, 1=Red, 2=Blue)")
        modplot(clustering.chrom, clustering.loc, clustering.reads, levels; mod=:mod_6mA, f=f, ax=ax_mod, unmod=false, hap=1, color=:crimson)
        modplot(clustering.chrom, clustering.loc, clustering.reads, levels; mod=:mod_6mA, f=f, ax=ax_mod, unmod=false, hap=2, color=:steelblue)
        modplot(clustering.chrom, clustering.loc, clustering.reads, levels; mod=:mod_6mA, f=f, ax=ax_mod, unmod=false, hap=0, color=:black)
        clusterlines!(clustering, ax=ax_mod, plotlines=true, setticks=true, color=:black)
        return f
    end

    cldf = [(cluster=k, hap=h, fisher(clustering.readtable.cluster .== k, clustering.readtable.hap .== h)...) for k = 1:opt_k, h = [0, 1, 2]] |> vec |> DataFrame
    
    ax = Axis(f[1, 1], xgridvisible=false, ygridvisible=false, yticks=1:opt_k, xlabel="log(Odds Ratio) [Hap vs Other]", ylabel="Cluster")
    
    plt = data(cldf) * 
          (mapping(:or => log, (:cluster, :hap) => (c, h) -> c + (h - 1) * 0.22, :logor_se => x -> 1.96 * x, color=:hap => nonnumeric) * visual(Errorbars, direction=:x) + 
           mapping(:or => log, (:cluster, :hap) => (c, h) -> c + (h - 1) * 0.22, color=:hap => nonnumeric, marker=:fisher_bpv => (x -> x < 0.05) => "Fisher p < 0.05") * visual(Scatter))
    d = draw!(ax, plt, scales(Color=(; palette=[:black, :crimson, :steelblue])))
    legend!(f[1, 2], d)
    vlines!(ax, 0, color=:black, linestyle=:dash)
    hidespines!(ax, :t, :r)
    ax.yreversed = true
    
    valid_ors = filter(!isinf, abs.(log.(cldf.or)) .+ cldf.logor_se .* 1.96)
    valid_ors = filter(!isnan, valid_ors)
    xmax = isempty(valid_ors) ? 2.0 : max(1.0, 1.15 * maximum(valid_ors))
    xlims!(ax, -xmax, xmax)
    
    ax_mod = Axis(f[1, 3], xgridvisible=false, ygridvisible=false, title="6mA by Haplotype (0=Unphased/black, 1=Red, 2=Blue)")
    modplot(clustering.chrom, clustering.loc, clustering.reads, levels; mod=:mod_6mA, f=f, ax=ax_mod, unmod=false, hap=1, color=:crimson)
    modplot(clustering.chrom, clustering.loc, clustering.reads, levels; mod=:mod_6mA, f=f, ax=ax_mod, unmod=false, hap=2, color=:steelblue)
    modplot(clustering.chrom, clustering.loc, clustering.reads, levels; mod=:mod_6mA, f=f, ax=ax_mod, unmod=false, hap=0, color=:black)
    clusterlines!(clustering, ax=ax_mod, plotlines=true, setticks=true, color=:black)
    
    f
end

"""
    sample_diagnostics(clustering)

Generates a diagnostic figure evaluating sample representation across clusters.
Left panel: Forest plot of sample log(Odds Ratio) and 95% CI from Fisher's exact test.
Right panel: Single-molecule 6mA tracks colored by sample.
Handles k = 1 or single-sample datasets gracefully.
"""
function sample_diagnostics(clustering)
    opt_k = clustering.optimal_clustering.optimal_k
    samples = hasproperty(clustering.readtable, :sample) ? unique(clustering.readtable.sample) : String[]
    n_samples = length(samples)
    levels = sortperm(sortperm(clustering.optimal_clustering.optimal_model.assignments))
    
    f = Figure(size=(1200, 500))
    
    if opt_k <= 1 || n_samples < 2
        ax_msg = Axis(f[1, 1], xgridvisible=false, ygridvisible=false)
        hidespines!(ax_msg, :t, :r, :b, :l)
        hidexdecorations!(ax_msg)
        hideydecorations!(ax_msg)
        msg = opt_k <= 1 ? "Homogeneous population (k = 1):\nAcross-cluster sample skew test not applicable." :
                           "Single sample dataset ($(first(samples))):\nAcross-sample skew test not applicable."
        text!(ax_msg, 0.5, 0.5, text=msg, align=(:center, :center), fontsize=13, color=:grey40)
        
        ax_mod = Axis(f[1, 2], xgridvisible=false, ygridvisible=false, title="6mA Reads by Sample")
        modplot(clustering.chrom, clustering.loc, clustering.reads, levels; mod=:mod_6mA, f=f, ax=ax_mod, unmod=false, strand=nothing, color=:steelblue)
        clusterlines!(clustering, ax=ax_mod, plotlines=true, setticks=true, color=:black)
        return f
    end

    # Multi-cluster, multi-sample: Compute Fisher exact tests
    
    cldf = [(cluster=k, sample=s, fisher(clustering.readtable.cluster .== k, clustering.readtable.sample .== s)...) for k in 1:opt_k, s in samples] |> vec |> DataFrame
    
    
    ax = Axis(f[1, 1], xgridvisible=false, ygridvisible=false, yticks=1:opt_k, xlabel="log(Odds Ratio) [Sample vs Other]", ylabel="Cluster")
    
    palette = [:steelblue, :crimson, :darkorange, :forestgreen, :purple, :teal, :magenta]
    sample_palette = palette[1:min(length(palette), n_samples)]
    
    sample_dict = Dict(s => i for (i, s) in enumerate(samples))
    
    plt = data(cldf) * 
          (mapping(:or => log, (:cluster, :sample) => ((c, s) -> c + (sample_dict[s] - (n_samples + 1)/2)*(0.4/n_samples)), :logor_se => (x -> 1.96 * x), color=:sample => nonnumeric) * visual(Errorbars, direction=:x) + 
           mapping(:or => log, (:cluster, :sample) => ((c, s) -> c + (sample_dict[s] - (n_samples + 1)/2)*(0.4/n_samples)), color=:sample => nonnumeric, marker=:fisher_bpv => (x -> x < 0.05) => "Fisher p < 0.05") * visual(Scatter))
    d = draw!(ax, plt, scales(Color=(; palette=sample_palette)))
    legend!(f[1, 2], d)
    vlines!(ax, 0, color=:black, linestyle=:dash)
    hidespines!(ax, :t, :r)
    ax.yreversed = true
    
    valid_ors = filter(!isinf, abs.(log.(cldf.or)) .+ cldf.logor_se .* 1.96)
    valid_ors = filter(!isnan, valid_ors)
    xmax = isempty(valid_ors) ? 2.0 : max(1.0, 1.15 * maximum(valid_ors))
    xlims!(ax, -xmax, xmax)
    
    ax_mod = Axis(f[1, 3], xgridvisible=false, ygridvisible=false, title="6mA Single-Molecule Reads by Sample")
    for (i, s) in enumerate(samples)
        c = sample_palette[mod1(i, length(sample_palette))]
        modplot(clustering.chrom, clustering.loc, clustering.reads, levels; mod=:mod_6mA, f=f, ax=ax_mod, unmod=false, strand=nothing, sample=s, color=c)
    end
    clusterlines!(clustering, ax=ax_mod, plotlines=true, setticks=true, color=:black)
    
    f
end

"""
    plot_genemodel_track(clustering, exondata; 
                         f=nothing, 
                         gene_flank::Int=15_000, 
                         longest_only::Bool=true)

Renders a dedicated gene model track spanning the clustered locus with flanking genomic context (default ±15 kb),
shading both the full fiber footprint window (`clustering.loc`) and target peak (`clustering.peakloc`).
Dynamically scales figure height based on packed transcript rows.
"""
function plot_genemodel_track(clustering, exondata; f=nothing, gene_flank::Int=15_000, longest_only::Bool=true)
    view_loc = max(1, first(clustering.loc) - gene_flank):(last(clustering.loc) + gene_flank)
    
    # 1. Fetch transcript models
    models = TranscriptModel[]
    if !isnothing(exondata)
        annot = annotate_locus(exondata, clustering.chrom, view_loc)
        models = filter_gene_models(annot.transcripts; longest_only=longest_only)
    end
    
    # Calculate row packing and dynamic height if f is not provided
    buffer = max(1000, div(length(view_loc), 50))
    rows = pack_gene_rows(models; buffer=buffer)
    n_rows = isempty(rows) ? 1 : maximum(rows)
    fig_h = max(130, 32 * n_rows + 55)
    
    fig = isnothing(f) ? Figure(size=(1600, fig_h)) : f
    ax = Axis(fig[1, 1], xgridvisible=false, ygridvisible=false, 
              xlabel="Genomic Position ($(clustering.chrom):$(first(view_loc))-$(last(view_loc)))")
    
    # 2. Shaded area for full fiber clustering window (loc)
    loc_s = first(clustering.loc)
    loc_e = last(clustering.loc)
    vspan!(ax, loc_s, loc_e, color=(:dodgerblue, 0.08))
    vlines!(ax, [loc_s, loc_e], color=(:dodgerblue, 0.35), linestyle=:dot, linewidth=1)
    
    # 3. Shaded area for target regulatory peak (peakloc)
    peak_s = first(clustering.peakloc)
    peak_e = last(clustering.peakloc)
    if peak_s != loc_s || peak_e != loc_e
        vspan!(ax, peak_s, peak_e, color=(:royalblue, 0.22))
        vlines!(ax, [peak_s, peak_e], color=(:royalblue, 0.6), linestyle=:dash, linewidth=1.2)
    end
    
    # 4. Plot gene models
    if !isempty(models)
        plot_genemodels!(ax, models, view_loc; longest_only=false, rows=rows)
    else
        ax.title = "No overlapping gene models in $(clustering.chrom):$(first(view_loc))-$(last(view_loc))"
    end
    
    xlims!(ax, first(view_loc), last(view_loc))
    hidespines!(ax, :t, :r)
    return fig
end

"""
    plot_cluster_tracks(clustering; f=Figure(size=(1600, 650)))

Builds the synchronized multi-channel single-molecule fiber track plot:
- Read assignment uncertainty beeswarm.
- Stencil Feature Heatmap (M).
- FIRE / Linker / Nucleosome track blocks.
- 6mA Accessibility scatter track.
- 5mC Methylation scatter track.
"""
function plot_cluster_tracks(clustering; f=Figure(size=(1600, 650)))
    levels = sortperm(sortperm(clustering.optimal_clustering.optimal_model.assignments))
    axs = Axis[]
    
    # 1. Uncertainty beeswarm on left
    boxplot_read_uncertainty(clustering; f=f, axx=1, axy=1, direction=:horizontal)
    
    # 2. Heatmap
    ax_feat = Axis(f[1, 2], xgridvisible=false, ygridvisible=false, title="Feature Heatmap", yreversed=true)
    featureheatmap(clustering, ax=ax_feat)
    push!(axs, ax_feat)

    # 3. FIRE elements
    ax_fire = Axis(f[1, 3], xgridvisible=false, ygridvisible=false, title="FIRE / Nucleosomes")
    fireplot(clustering.chrom, clustering.loc, clustering.reads, levels; f=f, ax=ax_fire, showlegend=false, rasterize=true)
    clusterlines!(clustering, ax=ax_fire, plotlines=true, setticks=true, color=:black)
    push!(axs, ax_fire)

    # 4. 6mA
    ax_6mA = Axis(f[1, 4], xgridvisible=false, ygridvisible=false, title="6mA Accessibility")
    modplot(clustering.chrom, clustering.loc, clustering.reads, levels; mod=:mod_6mA, f=f, ax=ax_6mA, unmod=false)
    clusterlines!(clustering, ax=ax_6mA, plotlines=true, setticks=true, color=:black)
    push!(axs, ax_6mA)

    # 5. 5mC
    ax_5mC = Axis(f[1, 5], xgridvisible=false, ygridvisible=false, title="5mC Methylation")
    modplot(clustering.chrom, clustering.loc, clustering.reads, levels; mod=:mod_5mC, f=f, ax=ax_5mC)
    clusterlines!(clustering, ax=ax_5mC, plotlines=true, setticks=true, color=:black)
    push!(axs, ax_5mC)

    for (k, a) in enumerate(axs)
        if k > 1
            a.yticks = ([], [])
        end
    end
    linkyaxes!(axs...)
    return f
end

"""
    plot_cluster_profiles(clustering; f=Figure(size=(1600, 750)))

Renders the 3 overlay profile charts comparing cluster centroids:
- FIRE vs 6mA
- Nucleosome vs 6mA
- 5mC vs FIRE
Uses compact row gaps and a single header row for column titles.
"""
function plot_cluster_profiles(clustering; f=Figure(size=(1600, 750)))
    cldf = clustermeans(clustering)
    
    plt_fire_6mA = data(filter(r -> r.Label in ["FIRE", "6mA"], cldf)) * 
                   mapping(:Loc, :ClusterMean, row=:Cluster => string, color=:Label => presorted) * 
                   visual(BarPlot, alpha=0.9)
    plt_nuc_6mA  = data(filter(r -> r.Label in ["Nucleosome", "6mA"], cldf)) * 
                   mapping(:Loc, :ClusterMean, row=:Cluster => string, color=:Label => presorted) * 
                   visual(BarPlot, alpha=0.9)
    plt_5mC_6mA  = data(filter(r -> r.Label in ["5mC", "FIRE"], cldf)) * 
                   mapping(:Loc, :ClusterMean, row=:Cluster => string, color=:Label => presorted) * 
                   visual(BarPlot, alpha=0.9)

    gl_fire = GridLayout(f[1, 1])
    gl_nuc  = GridLayout(f[1, 2])
    gl_5mc  = GridLayout(f[1, 3])

    lf = draw!(gl_fire, plt_fire_6mA; axis=(; xgridvisible=false, ygridvisible=false))
    cf = draw!(gl_nuc,  plt_nuc_6mA;  axis=(; xgridvisible=false, ygridvisible=false))
    rf = draw!(gl_5mc,  plt_5mC_6mA;  axis=(; xgridvisible=false, ygridvisible=false))

    rowgap!(gl_fire, 4)
    rowgap!(gl_nuc, 4)
    rowgap!(gl_5mc, 4)

    Label(f[1, 1, Top()], "FIRE vs 6mA", font=:bold, fontsize=13, padding=(0, 0, 8, 0))
    Label(f[1, 2, Top()], "Nucleosome vs 6mA", font=:bold, fontsize=13, padding=(0, 0, 8, 0))
    Label(f[1, 3, Top()], "5mC vs FIRE", font=:bold, fontsize=13, padding=(0, 0, 8, 0))

    legend!(f[2, 1], lf, orientation=:horizontal, tellheight=true)
    legend!(f[2, 2], cf, orientation=:horizontal, tellheight=true)
    legend!(f[2, 3], rf, orientation=:horizontal, tellheight=true)
    rowgap!(f.layout, 8)
    return f
end

"""
    plot_spectral_nrl(clustering; f=nothing, nrl_feature=:fire, phase_threshold=0.035)

Plots the detrended autocorrelation signal and best-fit cosine wave for each cluster centroid,
comparing both FIRE and 6mA derived spectral NRLs.
"""
function plot_spectral_nrl(clustering; f=nothing, nrl_feature=:fire, phase_threshold=0.035)
    k = clustering.optimal_clustering.optimal_k
    ncols = min(k, 4)
    nrows = cld(k, ncols)
    
    fig = isnothing(f) ? Figure(size=(1100, max(260, 220 * nrows))) : f
    
    has_fire = :fire in clustering.CM.features
    has_6mA  = :mod_6mA in clustering.CM.features
    
    FM_fire = has_fire ? featureblock(clustering.CM, :fire, M=clustering.optimal_clustering.optimal_model.centers) : nothing
    FM_6mA  = has_6mA  ? featureblock(clustering.CM, :mod_6mA, M=clustering.optimal_clustering.optimal_model.centers) : nothing
    
    ref_FM = !isnothing(FM_fire) ? FM_fire : FM_6mA
    n_lags = size(ref_FM, 1) - 1
    lags = clustering.CM.ldx[1:n_lags] .- first(clustering.CM.ldx)
    
    A_fire = !isnothing(FM_fire) ? autocor(FM_fire, 1:n_lags) : nothing
    A_6mA  = !isnothing(FM_6mA)  ? autocor(FM_6mA,  1:n_lags) : nothing
    
    detrended_fire = Vector{Float64}(undef, n_lags)
    detrended_6mA  = Vector{Float64}(undef, n_lags)

    for c in 1:k
        row = cld(c, ncols)
        col = mod1(c, ncols)
        
        res_fire = !isnothing(A_fire) ? characterise_nucleosomes_spectral!(A_fire[:, c], lags, detrended=detrended_fire, phase_threshold=phase_threshold) : nothing
        res_6mA  = !isnothing(A_6mA)  ? characterise_nucleosomes_spectral!(A_6mA[:, c],  lags, detrended=detrended_6mA,  phase_threshold=phase_threshold) : nothing
        
        title_parts = String[]
        if !isnothing(res_fire)
            push!(title_parts, res_fire.is_closed ? "FIRE: $(res_fire.nrl) bp" : "FIRE: open")
        end
        if !isnothing(res_6mA)
            push!(title_parts, res_6mA.is_closed ? "6mA: $(res_6mA.nrl) bp" : "6mA: open")
        end
        panel_title = "Cluster $c (" * join(title_parts, ", ") * ")"
        
        ax = Axis(fig[row, col], xgridvisible=false, ygridvisible=false,
                  title=panel_title, xlabel="Lag (bp)", ylabel="Detrended AC")
        
        if !isnothing(res_fire)
            lines!(ax, lags, detrended_fire, color=(:dodgerblue, 0.85), linewidth=1.5, label="FIRE AC")
            if res_fire.is_closed
                wave_fire = (2.0 * res_fire.phasing_score) .* cos.(2.0 * π .* lags ./ res_fire.nrl)
                lines!(ax, lags, wave_fire, color=:blue, linestyle=:dash, linewidth=1.5, label="FIRE $(res_fire.nrl) bp")
                vlines!(ax, res_fire.nrl, color=:blue, linestyle=:dot)
            end
        end
        
        if !isnothing(res_6mA)
            lines!(ax, lags, detrended_6mA, color=(:darkorange, 0.85), linewidth=1.5, label="6mA AC")
            if res_6mA.is_closed
                wave_6mA = (2.0 * res_6mA.phasing_score) .* cos.(2.0 * π .* lags ./ res_6mA.nrl)
                lines!(ax, lags, wave_6mA, color=:crimson, linestyle=:dash, linewidth=1.5, label="6mA $(res_6mA.nrl) bp")
                vlines!(ax, res_6mA.nrl, color=:crimson, linestyle=:dot)
            end
        end
        
        hidespines!(ax, :t, :r)
        axislegend(ax, framevisible=false, position=:rt, labelsize=8)
    end
    return fig
end