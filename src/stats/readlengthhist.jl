
############################# read length hist #############################
struct ReadLengthHist{S} <: RecordStat
    stat::S
end
ReadLengthHist() = ReadLengthHist(KHist(250))

statname(::Type{ReadLengthHist}) = "Read length histogram"

update!(stat::ReadLengthHist, reader::HTSFileReader, record::BamRecord, recorddata) = fit!(stat.stat, log10(Float64(querylength(record))))

function writestats(stat::ReadLengthHist, path::String, file="readlengthhist.tsv.gz")
    filepath = joinpath(path, file)
    readlengths = first.(stat.stat.bins)
    counts = last.(stat.stat.bins)
    df = DataFrame(readlength=readlengths, count=counts)
    CSV.write(filepath, df, delim='\t', compress=endswith(file, ".gz"))
end
function plotstat(::Type{ReadLengthHist}, data::DataFrame, scale=1, xtrans=x -> x)
    f = Figure()
    ax = Axis(f[1, 1], xgridvisible=false, ygridvisible=false, ylabel="Counts (k)")
    lines!(xtrans.(data.readlength), data.count, color=:black)
    mean_read_length = (xtrans.(data.readlength)'*(data.count/sum(data.count)))
    vlines!([mean_read_length/scale], color=:red, label=string("Mean = ", round(mean_read_length/scale, digits=2), " kb"))
    # Legend(f[1, 2], ax)
    axislegend(ax)
    if xtrans == identity
        xmax = 1.1*maximum(data.readlength)
    else
        xmax = min(1.1*xtrans(maximum(data.readlength)), 50_000)
    end
    # xlims!(0, xmax)
    
    # xlims!(0, 100)
    f
end
