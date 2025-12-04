module SMGStats

using SMGReader
using Statistics
using StatsBase
using OnlineStats
using ProgressMeter
# using UnicodePlots
using DataFrames, CSV
using GenomicFeatures

import YAML

using CairoMakie, AlgebraOfGraphics, ColorSchemes

export firestats, unicodeplot, displayplots, writeallstats, readstats, plotstat, htmlreport, statname, samplesummary,
        ChromStat, ModCrossCor, ModHist, ModMetaHist, ModRate, NucMSPLenHist, NucMSPRate, ReadLengthCounter, MetaPlot,
        @smgstats, 
        metaplot, loadregions,
        loadmultichromintervals, loadstatconfig
        


# include("unicode_tools.jl")
include("utils.jl")
include("bitpacked.jl")
include("regions.jl")
include("stats/stats.jl")
include("report.jl")
include("summarystats.jl")
end
