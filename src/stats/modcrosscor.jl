
############################# mod cross corr #############################
# mutable struct ModCrossCor{L} <: RecordStat
#     lags::L
#     t::Int
#     data::Dict{Tuple{Modification, Modification}, Vector{Float64}}
#     buffer::Dict{Modification, VectorBuffer{Bool}}
#     result::Vector{Float64}
# end
# ModCrossCor(lags=1:10:360, mods=(mod_6mA, mod_5mC, mod_5hmC)) = ModCrossCor(lags, 0,
#                                  Dict{Tuple{Modification, Modification}, Vector{Float64}}(),
#                                  Dict(m => VectorBuffer{Bool}(SMGReader.READ_BUFFER_LENGTH) for m in mods), Vector{Float64}(undef, length(lags)))

mutable struct ModCrossCor{L} <: RecordStat
    lags::L       # lags for cross correlation
    length::Int   # length of read where correlation is performed multiple of 64 currently 1024
    offset::Int   # start position within read for correlation
    t::Int        # total reads
    data::Dict{Tuple{Modification, Modification}, Vector{Float64}}  ### the sum of cross cors
    buffer::Dict{Modification, BitPacked{UInt64}}   ### prealloc buffer to store mod bitvectors
    result::Vector{Float64}                 ### preallocate vector for crosscor
end


function ModCrossCor(; lags=[1:10; 10:10:280], seqlength=1024, offset=180, mods=(mod_6mA, mod_5mC, mod_5hmC)) 
    
    totalreads = 0
    data = Dict((modA, modB) => zeros(length(lags)) for modA in mods, modB in mods if modA <= modB)
    buffer = Dict(m => BitPacked{UInt64}(1024) for m in mods)
    result = Vector{Float64}(undef, length(lags))

    ModCrossCor(lags, seqlength, offset, totalreads, data, buffer, result)
end


recordupdates(::Type{<:ModCrossCor}) = RecordUpdates
modupdates(::Type{<:ModCrossCor}) = ModUpdates
postmodupdates(::Type{<:ModCrossCor}) = PostModUpdates

instantiate(::Type{ModCrossCor}, reader, mods) = ModCrossCor(mods=mods)

statname(::Type{ModCrossCor}) = "Modification Auto/Cross Correlation"

@inline function updaterecord!(stat::ModCrossCor, record, recorddata)
    for v in values(stat.buffer)
        v .= false
    end
end
@inline function updatemod!(stat::ModCrossCor, mod::ModificationInfo, record::BamRecord, recorddata)
    thresh = 0.9*255
    @inbounds if  (querylength(record) >= (stat.length + stat.offset)) && (mod.prob > thresh) && (stat.offset <= mod.pos <= (stat.offset + stat.length - 1))
        i = mod.pos - stat.offset + 1
        stat.buffer[mod.mod][i] = true
    end
    nothing
end


@inline function updatepostmod!(stat::ModCrossCor, record::BamRecord, recorddata)
    if querylength(record) < stat.length + stat.offset
        return nothing
    end


    @inbounds for modA in keys(stat.buffer)
        for modB in keys(stat.buffer)
            if modA <= modB
                k = (modA, modB)
                crosscor!(stat.result, stat.buffer[modA], stat.buffer[modB], stat.lags)
                if !any(isnan, stat.result)
                    stat.data[k] .+= stat.result
                end
            end
        end
    end
    stat.t += 1
    nothing
end


function writestats(stat::ModCrossCor, path::String, file="modcrosscor.tsv.gz")
    filepath = joinpath(path, file)

    df = DataFrame(ModA=String[], ModB=String[], Lag=Int[], CrossCor=Float64[])
    for ((modA, modB), cc) in stat.data
        append!(df, DataFrame(ModA=string(modA), ModB=string(modB), Lag=stat.lags, CrossCor=cc/stat.t))
    end
    
    CSV.write(filepath, df, delim='\t', compress=endswith(file, ".gz"))
end
function plotstat(::Type{ModCrossCor}, df::DataFrame)
    # f = Figure()
    # ax = Axis(f[1, 1], xgridvisible=false, ygridvisible=false, ylabel="Correlation")
    df.Label = string.(replace.(df.ModA, "mod_" => ""), "\n", replace.(df.ModB, "mod_" => ""))

    plt = data(df) * mapping(:Lag, :CrossCor, color=:Label) * mapping(col=:Label) * (visual(Lines) + visual(Scatter))
    f = draw(plt; figure=(; size=(1000, 500)), legend=(; position=:bottom))
    f
end
