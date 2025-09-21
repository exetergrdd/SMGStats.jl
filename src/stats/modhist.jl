
############################# modifcation hist #############################
struct ModHist <: RecordStat
    data::Dict{Modification, Vector{Int}}
end
ModHist() = ModHist(Dict{Modification, Vector{Int}}())

statname(::Type{ModHist}) = "Modification histogram" 

function update!(stat::ModHist, reader::HTSFileReader, record::BamRecord, recorddata)
    for mi in ModIterator(record, recorddata)
        haskey(stat.data, mi.mod) || (stat.data[mi.mod] = zeros(Int, 256))
        stat.data[mi.mod][mi.prob + 1] += 1
    end
end

function writestats(stat::ModHist, path::String, file="modhist.tsv.gz")
    filepath = joinpath(path, file)
    df = DataFrame(Dict(string(m) => v for (m, v) in stat.data))
    modlabels = names(df)
    df.ml = 0:255
    df.prob = df.ml./255
    df = df[!, ["ml" ; "prob" ; modlabels]]
    # df = mapreduce((mod, vec) -> DataFrame(mod=mod, ml=ml, prob=prob, count=count), stat.data)
    CSV.write(filepath, df, delim='\t', compress=endswith(file, ".gz"))
end
function plotstat(::Type{ModHist}, df::DataFrame, scale=1e+6, thresh=0.9)
    sdf = stack(df, names(df, r"mod"))
    sdf.ismod = sdf.prob .> thresh
    transform!(groupby(sdf, :variable), :value => (v -> v/sum(v)) => :proportion)


    mdf = combine(groupby(sdf, [:variable, :ismod]), :proportion => sum => :proportion)
    mdf = mdf[mdf.ismod, :]

   

    f = Figure(; size=(1000, 400))
    prate = data(mdf) * mapping(:variable, :proportion, color=:variable) * visual(BarPlot, bar_labels=:y)
    draw!(f[1, 1], prate, axis= (; xgridvisible=false, ygridvisible=false, ylabel="Proportion modified", xlabel="Modification"))
    ylims!(-0.005, 1.2*maximum(mdf.proportion))

    plt = data(sdf) * mapping(:prob, :proportion , color=:variable) * mapping(row=:variable) * visual(BarPlot)
    ag = draw!(f[1, 2], plt, axis= (; xgridvisible=false, ygridvisible=false, ylabel="Mods / million", xlabel="Mod Probability"))
    legend!(f[1, 3], ag)
    f
end
