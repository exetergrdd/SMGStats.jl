
############################# mod rate #############################
struct ModRate <: RecordStat
    n::Int
    thresh::Float64
    stat::Dict{Modification, OnlineStats.Hist{Float64}}
    totalmods::Dict{Modification, Int}
end
ModRate(mods ; n=10_000, thresh=0.9*255) = ModRate(n, thresh, Dict(m => OnlineStats.Hist(range(0, 1, length=n+1)) for m in mods), Dict(m => 0 for m in mods))

recordupdates(::Type{ModRate}) = nothing
modupdates(::Type{ModRate}) = ModUpdates
postmodupdates(::Type{ModRate}) = PostModUpdates

instantiate(::Type{ModRate}, config) = ModRate(config.mods)

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


statdf(stat::ModRate) =  mapreduce(((k, v), ) -> DataFrame(mod=k, rate=edges_to_bins(v.edges), count=v.counts), vcat, collect(stat.stat))


function writestats(stat::ModRate, path::String, file="modrate.tsv.gz")
    filepath = joinpath(path, file)
    CSV.write(filepath, statdf(stat), delim='\t', compress=endswith(file, ".gz"))
end

function plotstat(::Type{ModRate}, df::DataFrame, scale=1, smooth_sig=10, axis_thresh=0.99999)
    transform!(groupby(df, :mod), :count => (c -> c/sum(c)) => :proportion)
    transform!(groupby(df, :mod), :proportion => cumsum => :cumsum, :proportion => (x -> kernelsmooth(x, smooth_sig)) => :smoothprop)

    n = length(unique(df.mod))
    xlt = combine(groupby(df[df.cumsum .>= axis_thresh, :], :mod), df -> df[argmin(df.cumsum), :]).rate |> maximum

    f  = Figure(size=(1000, 750))
    #### viobox
    ax = Axis(f[1, 1:n], xgridvisible=false, ygridvisible=false, ylabel="mod/bp")
    vioboxdist!(df, :mod, :rate, :proportion)
    ylims!(0, xlt)

    ### per mod dist
    
    plt = data(df) * mapping(:rate, direct(0), :smoothprop, color=:mod) * mapping(col=:mod) * visual(Band)
    draw!(f[2, (1:2*n)], plt, axis= (; xgridvisible=false, ygridvisible=false), facet=(; linkyaxes=false))
    xlims!(0, xlt)
    plt_overlay = data(df) * ((mapping(:rate, direct(0), :smoothprop, color=:mod) * visual(Band, alpha=0.5)) +  
                              (mapping(:rate, :smoothprop, color=:mod) * visual(Lines))) 
    ag = draw!(f[1, (1:n) .+ n], plt_overlay, axis=(; xlabel="mod/bp", xgridvisible=false, ygridvisible=false))
    legend!(f[1:2, 2*n + 1], ag)
    xlims!(0, xlt)
    f
end

