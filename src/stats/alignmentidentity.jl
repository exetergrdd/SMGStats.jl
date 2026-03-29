
############################# mod rate #############################
struct AlignmentIdentityHist <: RecordStat
    n::Int
    stat::OnlineStats.Hist{Float64}
end
AlignmentIdentityHist( ; n=2000) = AlignmentIdentityHist(n, OnlineStats.Hist(range(0, 1, length=n+1)))

recordupdates(::Type{AlignmentIdentityHist}) = RecordUpdates
modupdates(::Type{AlignmentIdentityHist}) = nothing
postmodupdates(::Type{AlignmentIdentityHist}) = nothing

instantiate(::Type{AlignmentIdentityHist}, config) = AlignmentIdentityHist()
instantiate(::Type{AlignmentIdentityHist}, reader, mods) = AlignmentIdentityHist()
statname(::Type{AlignmentIdentityHist}) = "Alignment Identity Histogram"


compatible(::Type{AlignmentIdentityHist}, ::Type{AuxMapMod}) = false
compatible(::Type{AlignmentIdentityHist}, ::Type{AuxMapModFire}) = false
compatible(::Type{AlignmentIdentityHist}, ::Type{AuxMapModFiberTools}) = false

@inline function updaterecord!(stat::AlignmentIdentityHist, record::BamRecord, recorddata)
    ai = alignmentidentity(record, recorddata)
    fit!(stat.stat, Float64(ai))
end

statdf(stat::AlignmentIdentityHist) =  DataFrame(Identity=edges_to_bins(stat.stat.edges), count=stat.stat.counts)


function writestats(stat::AlignmentIdentityHist, path::String, file="alignmentidentityhist.tsv.gz")
    filepath = joinpath(path, file)
    CSV.write(filepath, statdf(stat), delim='\t', compress=endswith(file, ".gz"))
end

function plotstat(::Type{AlignmentIdentityHist}, df::DataFrame)
    
    df.proportion = df.count/sum(df.count)

    f  = Figure(size=(1000, 500))
    ax = Axis(f[1, 1], xgridvisible=false, ygridvisible=false, ylabel="Proportion")


    # @show maximum(df.proportion)


    barplot!(df.Identity, df.proportion, color=:darkgrey, width=1/size(df, 1))

    w = Weights(df.count)
    mean_identity = mean(df.Identity, w)
    std_identity = std(df.Identity, w)
    
    vlines!([mean_identity], color=:red, label=string("Mean = ", tt(mean_identity), " ± ", tt(std_identity)))

    axislegend(ax, position=:lt)    
    xlims!(mean_identity - 3*std_identity, 1)

    f
end

