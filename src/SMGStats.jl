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
using Clustering
using Distances


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
    ModAllFireRead, ModGroupFileReader, collectstencilreads, read_depth_table, featurematrix, select_kmeans_k_bic




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
# include("countfires_sweep.jl")

include("main.jl")

end
