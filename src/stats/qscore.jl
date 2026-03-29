
############################# mod rate #############################
struct QscoreHist <: RecordStat
    n::Int
    stat::OnlineStats.Hist{Float64}
end
QscoreHist( ; n=120) = QscoreHist(n, OnlineStats.Hist(range(0, 60, length=n+1)))

recordupdates(::Type{QscoreHist}) = RecordUpdates
modupdates(::Type{QscoreHist}) = nothing
postmodupdates(::Type{QscoreHist}) = nothing

instantiate(::Type{QscoreHist}, config) = QscoreHist()
instantiate(::Type{QscoreHist}, reader, mods) = QscoreHist()
statname(::Type{QscoreHist}) = "Mean Qscore Histogram"

compatible(::Type{QscoreHist}, ::Type{AuxMapMod}) = false
compatible(::Type{QscoreHist}, ::Type{AuxMapModFire}) = false
compatible(::Type{QscoreHist}, ::Type{AuxMapModFiberTools}) = false

@inline function updaterecord!(stat::QscoreHist, record::BamRecord, recorddata)
    qs = basecallqscore(record, recorddata)
    fit!(stat.stat, Float64(qs))
end

statdf(stat::QscoreHist) =  DataFrame(Qscore=edges_to_bins(stat.stat.edges), count=stat.stat.counts)


function writestats(stat::QscoreHist, path::String, file="qscorehist.tsv.gz")
    filepath = joinpath(path, file)
    CSV.write(filepath, statdf(stat), delim='\t', compress=endswith(file, ".gz"))
end

function plotstat(::Type{QscoreHist}, df::DataFrame)
    
    df.proportion = df.count/sum(df.count)

    f  = Figure(size=(1000, 500))
    ax = Axis(f[1, 1], xgridvisible=false, ygridvisible=false, ylabel="Proportion")


    barplot!(df.Qscore, df.proportion, color=:darkgrey, width=0.5)

    w = Weights(df.count)
    mean_qscore = mean(df.Qscore, w)
    std_qscore = std(df.Qscore, w)
    
    vlines!([mean_qscore], color=:red, label=string("Mean = ", tt(mean_qscore), " ± ", tt(std_qscore)))

    axislegend(ax)    
    xlims!(0, 60)

    f
end

