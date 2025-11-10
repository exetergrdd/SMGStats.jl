

############################# mod meta read size #############################

struct ModMetaHist <: RecordStat
    n::Int
    thresh::Float64
    data::Matrix{Int}
end
ModMetaHist(mods; n=1000, thresh=0.9*255) = ModMetaHist(n, thresh, zeros(n, length(instances(Modification))))
instantiate(::Type{ModMetaHist}, reader, mods) = ModMetaHist(mods)


recordupdates(::Type{ModMetaHist})  = nothing
modupdates(::Type{ModMetaHist})     = ModUpdates
postmodupdates(::Type{ModMetaHist}) = nothing


statname(::Type{ModMetaHist}) = "Modification Distribution along read"
@inline function updatemod!(stat::ModMetaHist, mod::ModificationInfo, record::BamRecord, recorddata)
    if mod.prob > stat.thresh
        i = Int(floor(stat.n*mod.pos/querylength(record)))
        i = max(i, 1)
        @inbounds stat.data[i, Int(mod.mod) + 1] += 1
    end
end



function writestats(stat::ModMetaHist, path::String, file="mod_meta_hist.tsv.gz")
    filepath = joinpath(path, file)
    ind = vec(sum(stats.data, dims=1) .> 0)
    counts = stat.data[:, ind]
    mods = instances(Modification)[ind]

    df = mapreduce((mod, counts) -> DataFrame(mod=mod, bin=range(0, 1, length=stat.n), count=counts), vcat, mods, eachcol(counts))

    CSV.write(filepath, df, delim='\t', compress=endswith(file, ".gz"))
end

function plotstat(::Type{ModMetaHist}, df::DataFrame, scale=1)

    plt = data(df) * mapping(:bin, :count, color=:mod) * mapping(col=:mod) * visual(Lines)
    f = draw(plt, axis= (; xgridvisible=false, ygridvisible=false))
    plt_overlay = data(df) * mapping(:bin, :count, color=:mod) * visual(Lines)
    n = length(unique(df.mod))
    
    draw!(f.figure[2, 1:n], plt_overlay, axis=(; xlabel="mod/bp", xgridvisible=false, ygridvisible=false))
    f
end
