

############################# nuc msp size #############################
struct NucMSPLenHist <: RecordStat
    nuc::Vector{Int}
    msp::Vector{Int}
end
NucMSPLenHist(n=1000) = NucMSPLenHist(zeros(n), zeros(n))

instantiate(::Type{NucMSPLenHist}, reader, mods) = NucMSPLenHist()

recordupdates(::Type{<:NucMSPLenHist})  = RecordUpdates
modupdates(::Type{<:NucMSPLenHist})     = nothing
postmodupdates(::Type{<:NucMSPLenHist}) = nothing

statname(::Type{NucMSPLenHist}) = "Nucleosome and MSP length histogram"

@inline function updaterecord!(stat::NucMSPLenHist, record::BamRecord, recorddata)

    @inbounds for nl in SMGReader.firenuclen(record, recorddata.auxmap)
        if 1 ≤ nl ≤ length(stat.nuc)
            stat.nuc[nl] += 1
        end
    end

    @inbounds for al in SMGReader.firemsplen(record, recorddata.auxmap)
        if 1 ≤ al ≤ length(stat.nuc)
            stat.msp[al] += 1
        end
    end
end
# unicodeplot(stat::NucMSPLenHist) = [barplot(1:length(stat.nuc), stat.nuc, title="Nucleosome Length Hist"),
                                    # barplot(1:length(stat.msp), stat.msp, title="MSP Length HIst")]

function writestats(stat::NucMSPLenHist, path::String, file="nuc_msp_length_hist.tsv.gz")
    filepath = joinpath(path, file)

    features = [fill("nuc", length(stat.nuc)) ; fill("msp", length(stat.msp))]
    lengths  = [1:length(stat.nuc) ; 1:length(stat.msp)]
    counts   = [stat.nuc ; stat.msp]

    df = DataFrame(Feature=features, Length=lengths, Count=counts)
    CSV.write(filepath, df, delim='\t', compress=endswith(file, ".gz"))
end

function plotstat(::Type{NucMSPLenHist}, df::DataFrame, scale=1)
    df = df[df.Length .> 1, :]
    plt = data(df) * mapping(:Length, :Count, color=:Feature) * mapping(col=:Feature) * visual(Lines)
    f = draw(plt, axis= (; xgridvisible=false, ygridvisible=false))
    plt_overlay = data(df) * mapping(:Length, :Count, color=:Feature) * visual(Lines)
    n = length(unique(df.Feature))
    
    draw!(f.figure[2, 1:n], plt_overlay, axis=(; xgridvisible=false, ygridvisible=false))
    f
end
