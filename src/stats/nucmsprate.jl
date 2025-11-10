
############################# nuc msp rate #############################
struct NucMSPRate <: RecordStat
    nuc::KHist{Float64}
    msp::KHist{Float64}
end
NucMSPRate(n=500) = NucMSPRate(KHist(n), KHist(n))
instantiate(::Type{NucMSPRate}, reader, mods) = NucMSPRate()

recordupdates(::Type{<:NucMSPRate})  = RecordUpdates
modupdates(::Type{<:NucMSPRate})     = nothing
postmodupdates(::Type{<:NucMSPRate}) = nothing


statname(::Type{NucMSPRate}) = "Nucleosome and MSP rate historgram"

@inline function updaterecord!(stat::NucMSPRate, record::BamRecord, recorddata)
   
    nn = length(SMGReader.firenucpos(record, recorddata))
    nm = length(SMGReader.firemsppos(record, recorddata))

    fit!(stat.nuc, nn/querylength(record))
    fit!(stat.msp, nm/querylength(record))
end


function writestats(stat::NucMSPRate, path::String, file="nuc_msp_rate.tsv.gz")
    filepath = joinpath(path, file)

    features = [fill("nuc", length(stat.nuc.bins)) ; fill("msp", length(stat.msp.bins))]
    rates  = [first.(stat.nuc.bins) ; first.(stat.msp.bins)]
    counts = [last.(stat.nuc.bins) ; last.(stat.msp.bins)]

    df = DataFrame(Feature=features, Rate=rates, Count=counts)
    CSV.write(filepath, df, delim='\t', compress=endswith(file, ".gz"))
end


function plotstat(::Type{NucMSPRate}, df::DataFrame, scale=1)
    # plt = data(df) * mapping(:Rate, :Count, color=:Feature) * mapping(col=:Feature) * visual(Lines)
    # f = draw(plt, axis= (; xgridvisible=false, ygridvisible=false))
    # plt_overlay = data(df) * mapping(:Rate, :Count, color=:Feature) * visual(Lines)
    # n = length(unique(df.Feature))
    
    # draw!(f.figure[2, 1:n], plt_overlay, axis=(; xlabel="feature/bp", xgridvisible=false, ygridvisible=false))
    # f


    #####
    transform!(groupby(df, :Feature), :Count => (c -> c/sum(c)) => :proportion)

    ### quantiles for boxplot
    qdf = combine(groupby(df, :Feature, sort=true), df -> quantile_from_freq(df.Rate, df.proportion))
    
    
    n = length(unique(df.Feature))
    f  = Figure(size=(1000, 400))
    ### draw boxplots from qdf
    ax = Axis(f[1:2, 1:n], xgridvisible=false, ygridvisible=false, ylabel="Feature/bp", xticks=(1:n, qdf.Feature))
    quantileboxplot!(qdf)

    plt = data(df) * mapping(:Rate, :proportion, color=:Feature) * mapping(col=:Feature) * visual(Lines)
    draw!(f[1, (1:n) .+ n], plt, axis= (; xgridvisible=false, ygridvisible=false))
    plt_overlay = data(df) * mapping(:Rate, :proportion, color=:Feature) * visual(Lines)
    ag = draw!(f[2, (1:n) .+ n], plt_overlay, axis=(; xlabel="Feature/bp", xgridvisible=false, ygridvisible=false))
    legend!(f[1:2, 2*n + 1], ag)

    f
end
