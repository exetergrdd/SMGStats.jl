module SMGStats

import YAML

using AlgebraOfGraphics
using CSV
using CairoMakie
using ColorSchemes
using DataFrames
using Dates
using GenomicFeatures
using OnlineStats
using ProgressMeter
using SMGReader
using Statistics
using StatsBase

using ImageFiltering

export firestats, unicodeplot, displayplots, writeallstats, readstats, plotstat, htmlreport, statname, samplesummary,
        ChromStat, ModCrossCor, ModHist, ModMetaHist, ModRate, NucMSPLenHist, NucMSPRate, ReadLengthCounter, MetaPlot, NucMonoDiRatio,
        @smgstats, 
        metaplot, loadregions,
        loadmultichromintervals, loadstatconfig,
        runstats,
        combinereports
        



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
include("summarystats.jl")
include("multireport.jl")
include("main.jl")

end
