
############################# read length counter #############################
struct ReadLengthCounter <: RecordStat
    data::Dict{Int, Int}
end
ReadLengthCounter() = ReadLengthCounter(Dict{Int, Int}())

instantiate(::Type{ReadLengthCounter}, config) = ReadLengthCounter()

instantiate(::Type{ReadLengthCounter}, reader, mods) = ReadLengthCounter()

recordupdates(::Type{ReadLengthCounter}) = RecordUpdates
modupdates(::Type{ReadLengthCounter}) = nothing
postmodupdates(::Type{ReadLengthCounter}) = nothing

statname(::Type{ReadLengthCounter}) = "Read length histogram"

@inline function updaterecord!(stat::ReadLengthCounter, record::BamRecord, recorddata)
    l = querylength(record)
    stat.data[l] = get(stat.data, l, 0) + 1
    nothing
end


function writestats(stat::ReadLengthCounter, path::String, file="readlengthcounts.tsv.gz")
    filepath = joinpath(path, file)
    readlengths = sort(collect(keys(stat.data)))
    counts = getindex.(Ref(stat.data), readlengths)
    df = DataFrame(readlength=readlengths, count=counts)
    CSV.write(filepath, df, delim='\t', compress=endswith(file, ".gz"))
end
function plotstat(::Type{ReadLengthCounter}, data::DataFrame, xscale=1000, yscale=1000, binsize=500)
    f = Figure()
    ax = Axis(f[1, 1], xgridvisible=false, ygridvisible=false, ylabel="Counts (k)")

    data.bins = floor.(data.readlength/binsize)
    gd = combine(groupby(data, :bins), :count => sum => :count, :readlength => mean => :readlength)

    barplot!(gd.readlength/xscale, gd.count/yscale, color=:lightgrey, width=0.5)

    w = Weights(data.count)
    mean_read_length = mean(data.readlength, w)
    std_read_length = std(data.readlength, w)
    
    mean_read_length /= xscale
    std_read_length  /= xscale
    
    vlines!([mean_read_length], color=:red, label=string("Mean = ", tt(mean_read_length), " ± ", tt(std_read_length), "kb"))

    axislegend(ax)

    
    xlims!(0, mean_read_length + 6*std_read_length)
    f
end
