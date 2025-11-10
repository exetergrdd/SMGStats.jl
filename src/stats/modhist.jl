
############################# modifcation hist #############################


struct ModHist <: RecordStat
    data::Matrix{Int}
end
ModHist() = ModHist(zeros(256, length(instances(Modification))))


instantiate(::Type{ModHist}, reader, mods) = ModHist(zeros(256, length(instances(Modification))))
recordupdates(::Type{ModHist}) = nothing
modupdates(::Type{ModHist}) = ModUpdates
postmodupdates(::Type{ModHist}) = nothing


statname(::Type{ModHist}) = "Modification histogram" 

@inline updatemod!(stat::ModHist, mod::ModificationInfo, record::BamRecord, recorddata) = stat.data[mod.prob + 1, Int(mod.mod) + 1] += 1

function writestats(stat::ModHist, path::String, file="modhist.tsv.gz")
    filepath = joinpath(path, file)
    ind = sum(stat.data, dims=1) .> 0
    df = [DataFrame(ml=0:255, prob=(0:255)/255) DataFrame(stat.data[:, ind], string.(instances(Modification)[ind]))]
    
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
