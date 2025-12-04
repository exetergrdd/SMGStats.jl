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

export firestats, unicodeplot, displayplots, writeallstats, readstats, plotstat, htmlreport, statname, samplesummary,
        ChromStat, ModCrossCor, ModHist, ModMetaHist, ModRate, NucMSPLenHist, NucMSPRate, ReadLengthCounter, MetaPlot,
        @smgstats, 
        metaplot, loadregions,
        loadmultichromintervals, loadstatconfig,
        runstats
        


# include("unicode_tools.jl")
include("utils.jl")
include("bitpacked.jl")
include("regions.jl")
include("stats/stats.jl")
include("report.jl")
include("summarystats.jl")
include("main.jl")
end
