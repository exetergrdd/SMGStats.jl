module SMGStats

import YAML

using AlgebraOfGraphics
using CSV
using CairoMakie
using ColorSchemes
using CodecZlib
using DataFrames
using Dates
using GenomicFeatures
using Markdown
using OnlineStats
using ProgressMeter
using SMGReader
using Statistics
using StatsBase
using Random
using Clustering
using Distances
using JLD2
using ImageFiltering

using HypothesisTests
using SwarmMakie


export firestats, unicodeplot, displayplots, writeallstats, readstats, readmultistats, plotstat, plotmultistat, htmlreport, htmlmultireport, statname, samplesummary,
    ChromStat, ModCrossCor, ModHist, ModMetaHist, ModRate, NucMSPLenHist, NucMSPRate, ReadLengthCounter, MetaPlot, NucMonoDiRatio,
    @smgstats,
    metaplot, loadregions,
    loadmultichromintervals, loadstatconfig,
    runstats,
    combinereports,
    modtable,
    coaccess, coaccess_wrapper,
    countfire, countfire_wrapper, countfire_sweep,
    loadbed, bedintervals, buildintervalset, bedchromintervals, buildchromintervalset,
    differentialaccessibility,
    smglog,
    ModAllFireRead, ModGroupFileReader, collectstencilreads, read_depth_table, featurematrix, select_k_kmeans_bic,
    default_cluster_order, cluster_order, reorder_clusters!,
    load_exondata, annotate_locus, filter_gene_models, pack_gene_rows, plot_genemodels!,
    plotclustersummary, plot_genemodel_track, plot_cluster_tracks, plot_cluster_profiles, plot_spectral_nrl,
    plot_5mC_volcano, plot_5mC_diff_heat, strand_diagnostics, hap_diagnostics, sample_diagnostics,
    render_cluster_report, render_cluster_guide_html, LocusSummary, extract_locus_summary, build_cluster_catalog, load_clustering




"""
    smglog(msg)

Prints a timestamped, grep-friendly log line:

    [2025-02-12T14:22:03] [SMGStats] <msg>
"""
function smglog(msg...)
    msg = string(msg...)
    ts = Dates.format(Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS")
    println("[$ts] [SMGStats] $msg")
end


# include("unicode_tools.jl")
include("utils.jl")
include("bitpacked.jl")
include("regions.jl")
include("stats/stats.jl")
include("report.jl")
include("multireport.jl")
include("summarystats.jl")

include("coaccess.jl")
include("diff.jl")

include("countfires.jl")
include("collectreads.jl")
include("stencilcluster.jl")
include("exondata.jl")
include("stencil_cluster_report.jl")
include("html_cluster_report.jl")
include("report_catalog.jl")
# include("countfires_sweep.jl")

include("main.jl")

end
