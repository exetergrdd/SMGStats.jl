############################# nuc mono di ratio size #############################
struct NucMonoDiRatio{F} <: RecordStat
    stat::OnlineStats.Hist{Float64}
    trans::F
end
NucMonoDiRatio(bins=range(-5, 20, length=100), trans=log2) = NucMonoDiRatio(OnlineStats.Hist(bins), trans)


instantiate(::Type{NucMonoDiRatio}, config) = NucMonoDiRatio()



recordupdates(::Type{<:NucMonoDiRatio})  = RecordUpdates
modupdates(::Type{<:NucMonoDiRatio})     = nothing
postmodupdates(::Type{<:NucMonoDiRatio}) = nothing

statname(::Type{NucMonoDiRatio}) = "Nucleosome and Mono to Di (or above) ratio"
### stat not compatible with stats missing the firefields
compatible(::Type{NucMonoDiRatio}, ::Type{AuxMapMod}) = false


@inline function updaterecord!(stat::NucMonoDiRatio, record::BamRecord, recorddata, mono_thresh=225)

    total_mono = 0
    total_di = 0
    @inbounds for nl in SMGReader.firenuclen(record, recorddata.auxmap)
        if nl < mono_thresh
            total_mono += 1
        else
            total_di += 1
        end
    end
    fit!(stat.stat, stat.trans(total_mono./total_di))
end

statdf(stat::NucMonoDiRatio) =  DataFrame(Ratio=edges_to_bins(stat.stat.edges), count=stat.stat.counts, fun=stat.trans)

function writestats(stat::NucMonoDiRatio, path::String, file="nuc_mono_di_ratio.tsv.gz")
    filepath = joinpath(path, file)
    CSV.write(filepath, statdf(stat), delim='\t', compress=endswith(file, ".gz"))
end

function plotstat(::Type{NucMonoDiRatio}, df::DataFrame, scale=1)
    df.proportion = df.count/sum(df.count)
    
    f = Figure(size=(1000, 500))
    ax = Axis(f[1, 1], xlabel=string(df.fun[1], " ratio"), xgridvisible=false, ygridvisible=false)
    barplot!(df.Ratio, df.proportion, color=:darkgrey)
    
    # plt = data(df) * mapping(:Ratio, :proportion) * visual(BarPlot)
    # f = draw(plt, figure=(; size=(1000, 500)), axis= (; xgridvisible=false, ygridvisible=false))

    w = Weights(df.count)
    mean_ratio = mean(df.Ratio, w)
    std_ratio = std(df.Ratio, w)
    

    vlines!([mean_ratio], color=:red, label=string("Mean = ", tt(mean_ratio), " ± ", tt(std_ratio), ""))
    axislegend()
    f
end

