
############################# nuc msp rate #############################
struct NucMSPRate <: RecordStat
    nuc::KHist{Float64}
    msp::KHist{Float64}
end
NucMSPRate(n=150) = NucMSPRate(KHist(n), KHist(n))
instantiate(::Type{NucMSPRate}, config) = NucMSPRate()

instantiate(::Type{NucMSPRate}, reader, mods) = NucMSPRate()

recordupdates(::Type{<:NucMSPRate})  = RecordUpdates
modupdates(::Type{<:NucMSPRate})     = nothing
postmodupdates(::Type{<:NucMSPRate}) = nothing


statname(::Type{NucMSPRate}) = "Nucleosome and MSP rate histogram"
### stat not compatible with stats missing the firefields
compatible(::Type{NucMSPRate}, ::Type{AuxMapMod}) = false

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


function plotstat(::Type{NucMSPRate}, df::DataFrame, scale=1, smooth_sig=10, axis_thresh=0.999)
  
    transform!(groupby(df, :Feature), :Count => (c -> c/sum(c)) => :proportion)
    transform!(groupby(df, :Feature), :proportion => cumsum => :cumsum, :proportion => (x -> kernelsmooth(x, smooth_sig)) => :smoothprop)

    n = length(unique(df.Feature))
    xlt_thresh = combine(groupby(df[df.cumsum .> axis_thresh, :], :Feature), df -> df[argmin(df.cumsum), :]).Rate |> maximum
    xlt_med = combine(groupby(df[df.cumsum .> 0.5, :], :Feature), df -> df[argmin(df.cumsum), :]).Rate |> maximum

    xlt = max(xlt_thresh, 2*xlt_med)

    f  = Figure(size=(1000, 500))
    #### viobox
    ax = Axis(f[1:2, 1:n], xgridvisible=false, ygridvisible=false, ylabel="Feature/bp")
    vioboxdist!(df, :Feature, :Rate, :proportion)
    ylims!(0, xlt)

    ### per Feature dist
    
    plt = data(df) * mapping(:Rate, direct(0), :smoothprop, color=:Feature) * mapping(col=:Feature) * visual(Band)
    draw!(f[1, (1:n) .+ n], plt, axis= (; xgridvisible=false, ygridvisible=false))
    xlims!(0, xlt)
    plt_overlay = data(df) * ((mapping(:Rate, direct(0), :smoothprop, color=:Feature) * visual(Band, alpha=0.5)) +  
                              (mapping(:Rate, :smoothprop, color=:Feature) * visual(Lines))) 
    ag = draw!(f[2, (1:n) .+ n], plt_overlay, axis=(; xlabel="Feature/bp", xgridvisible=false, ygridvisible=false))
    legend!(f[1:2, 2*n + 1], ag)
    xlims!(0, xlt)
    f
end
