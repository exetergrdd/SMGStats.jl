
############################# mod rate #############################
struct ModRate <: RecordStat
    n::Int
    thresh::Float64
    stat::Dict{Modification, KHist{Float64}}
    totalmods::Dict{Modification, Int}
end
ModRate(mods ; n=500, thresh=0.9*255) = ModRate(n, thresh, Dict(m => KHist(n) for m in mods), Dict(m => 0 for m in mods))

recordupdates(::Type{ModRate}) = nothing
modupdates(::Type{ModRate}) = ModUpdates
postmodupdates(::Type{ModRate}) = PostModUpdates

instantiate(::Type{ModRate}, reader, mods) = ModRate(mods)

statname(::Type{ModRate}) = "Modification Rate Histogram"

@inline function updatemod!(stat::ModRate, mod::ModificationInfo, record::BamRecord, recorddata)
    if mod.prob > stat.thresh
        stat.totalmods[mod.mod] = stat.totalmods[mod.mod] + 1
    end
end

@inline function updatepostmod!(stat::ModRate, record::BamRecord, recorddata)
    for (mod, total) in stat.totalmods
        fit!(stat.stat[mod], total/querylength(record))
    end
    ## reset for next read
    for mod in keys(stat.totalmods)
        stat.totalmods[mod] = 0
    end
end



# function update!(stat::ModRate, reader::HTSFileReader, record::BamRecord, recorddata)
#     totalmods = Dict{Modification, Int}()
#     for mi in ModIterator(record, recorddata)
#         (mi.prob > stat.thresh)   || continue
#         haskey(totalmods, mi.mod) || (totalmods[mi.mod] = 0)
#         totalmods[mi.mod] = totalmods[mi.mod] + 1
#     end
#     for (mod, total) in totalmods
#         haskey(stat.stat, mod) || (stat.stat[mod] = KHist(stat.n))
#         fit!(stat.stat[mod], total/querylength(record))
#     end
# end
# unicodeplot(stat::ModRate) = [lineplot(first.(kh.bins), last.(kh.bins), title=string(mod)) for (mod, kh) in stat.stat]

function writestats(stat::ModRate, path::String, file="modrate.tsv.gz")
    filepath = joinpath(path, file)
    df = mapreduce(((k, v), ) -> DataFrame(mod=k, rate=first.(v.bins), count=last.(v.bins)), vcat, collect(stat.stat))
    CSV.write(filepath, df, delim='\t', compress=endswith(file, ".gz"))
end

function plotstat(::Type{ModRate}, df::DataFrame, scale=1)
    transform!(groupby(df, :mod), :count => (c -> c/sum(c)) => :proportion)

    ### quantiles for boxplot
    qdf = combine(groupby(df, :mod), df -> quantile_from_freq(df.rate, df.proportion))
    sort!(qdf, :mod)
    # display(qdf)
    n = length(unique(df.mod))
    f  = Figure(size=(1000, 400))
    ### draw boxplots from qdf
    ax = Axis(f[1:2, 1:n], xgridvisible=false, ygridvisible=false, ylabel="mod/bp", xticks=(1:n, qdf.mod))
    quantileboxplot!(qdf)

    plt = data(df) * mapping(:rate, :proportion, color=:mod) * mapping(col=:mod) * visual(Lines)
    draw!(f[1, (1:n) .+ n], plt, axis= (; xgridvisible=false, ygridvisible=false))
    plt_overlay = data(df) * mapping(:rate, :proportion, color=:mod) * visual(Lines)
    ag = draw!(f[2, (1:n) .+ n], plt_overlay, axis=(; xlabel="mod/bp", xgridvisible=false, ygridvisible=false))
    legend!(f[1:2, 2*n + 1], ag)

    f
end
