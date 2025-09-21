

############################# mod meta read size #############################
struct ModMetaHist <: RecordStat
    n::Int
    thresh::Float64
    data::Dict{Modification, Vector{Int}}
end
ModMetaHist(n=1000) = ModMetaHist(n, 0.9*255, Dict{Modification, Vector{Int}}())

statname(::Type{ModMetaHist}) = "Modification Distribution along read"

function update!(stat::ModMetaHist, reader::HTSFileReader, record::BamRecord, recorddata)
    for mi in ModIterator(record, recorddata)
        (mi.prob > stat.thresh)   || continue
        haskey(stat.data, mi.mod) || (stat.data[mi.mod] = zeros(Int, stat.n))

        i = Int(floor(stat.n*mi.pos/querylength(record)))
        i = max(i, 1)
        stat.data[mi.mod][i] += 1
    end
end
# unicodeplot(stat::ModMetaHist) = ...


function writestats(stat::ModMetaHist, path::String, file="mod_meta_hist.tsv.gz")
    filepath = joinpath(path, file)
    df = mapreduce(((mod, counts), ) -> DataFrame(mod=mod, bin=range(0, 1, length=stat.n), count=counts), vcat, collect(stat.data))
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
