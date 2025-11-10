module SMGStats

using SMGReader
using Statistics
using StatsBase
using OnlineStats
using ProgressMeter
# using UnicodePlots
using DataFrames, CSV

using CairoMakie, AlgebraOfGraphics, ColorSchemes

export firestats, unicodeplot, displayplots, writeallstats, readstats, plotstat, htmlreport, statname, samplesummary,
        ChromStat, ModCrossCor, ModHist, ModMetaHist, ModRate, NucMSPLenHist, NucMSPRate, ReadLengthCounter,
        @smgstats
        


# include("unicode_tools.jl")
include("utils.jl")
include("bitpacked.jl")
include("stats/stats.jl")
include("report.jl")
end
