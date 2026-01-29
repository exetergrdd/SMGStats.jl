### functions for coaccessibility of fire peaks

### loading fire intervals


struct CoaccessContingency
    a::Int32
    b::Int32
    c::Int32
    d::Int32
end
CoaccessContingency() = CoaccessContingency(zero(Int32), zero(Int32), zero(Int32), zero(Int32))

struct IntervalPair{T}
    ivA::T
    ivB::T
end


function updatecontigency(c, fA, fB)
    fA && fB && return CoaccessContingency(c.a + 1, c.b, c.c, c.d)
    fA && .!fB && return CoaccessContingency(c.a, c.b + 1, c.c, c.d)
    .!fA && fB && return CoaccessContingency(c.a, c.b, c.c + 1, c.d)
    .!fA && .!fB && return CoaccessContingency(c.a, c.b, c.c, c.d + 1)
end


function coaccess(file, intervals::GenomicIntervalCollection{GenomicInterval{T}}; total_reads=1000, fire_backing=5_000) where T

    reader = open(HTSFileReader, file)
    recorddata = StencillingData(AuxMapModFire())

    ### fireinteresctions is a buffer of those intervals each read intersects with
    fireintersections = VectorBuffer{GenomicInterval{T}}(fire_backing)
    ### a buffer of bool to determine whether a read has a fire element that intersects with the interval sets
    intersects = VectorBuffer{Bool}(fire_backing)

    codict = Dict{IntervalPair,CoaccessContingency}()

    t = 0
    for record in eachrecord(reader)
        validflag(record) || continue
        processread!(record, recorddata) || continue
        ref = refname(reader, record)
        readinterval = GenomicInterval(ref, SMGReader.leftposition(record) + 1, SMGReader.rightposition(recorddata)) ## convert to 1-based inclusive
        fi = 0
        setlength!(fireintersections, length(fireintersections.data))
        for fiv in eachoverlap(intervals, readinterval)
            fi += 1
            fireintersections[fi] = fiv
        end
        setlength!(fireintersections, fi)
        setlength!(intersects, fi)
        intersects .= false

        for (s, l, q) in firemsps(record, recorddata)
            if q >= 0.9 * 255
                genstart, genstop = genomecoords(s, l, record, recorddata) ### this could be improved if s+l lands in a coord that doesn't map this doesn't search back for a match
                if !iszero(genstart) && !iszero(genstop)
                    # giv = genstart:genstop
                    @inbounds @simd for i = eachindex(fireintersections)
                        intersects[i] |= (genstart <= fireintersections[i].last) && (genstop >= fireintersections[i].first)
                    end
                end
            end
        end

        for i = 1:(length(intersects)-1)
            for j = (i+1):length(intersects)
                ip = IntervalPair(fireintersections[i], fireintersections[j])
                if !haskey(codict, ip)
                    codict[ip] = CoaccessContingency()
                end
                codict[ip] = updatecontigency(codict[ip], intersects[i], intersects[j])
            end
        end

        t += 1
        (t == total_reads) && break
    end

    close(reader)

    codict
end


function coaccessdf(cod, calcfisher=false, n_thresh=10)
    nts = [(chromA=seqname(ip.ivA), startA=ip.ivA.first, stopA=ip.ivA.last, chromB=seqname(ip.ivB), startB=ip.ivB.first, stopB=ip.ivB.last, a=cc.a, b=cc.b, c=cc.c, d=cc.d) for (ip, cc) in cod]
    df = DataFrame(nts)
    df.a = Int.(df.a)
    df.b = Int.(df.b)
    df.c = Int.(df.c)
    df.d = Int.(df.d)
    df.Total = df.a .+ df.b .+ df.c .+ df.d

    if calcfisher
        df = coaccessfisher(df, n_thresh)
    end

    df
end

function coaccessfisher(codf, n_thresh=10)
    df = codf[codf.Total .> n_thresh, :]

    # df = @subset(codf, :Total .> n_thresh)
    ft = FisherExactTest.(df.a, df.b, df.c, df.d)
    df.or = [f.ω for f in ft]
    df.logor_se = sqrt.((1.0 ./ df.a) .+ (1.0 ./ df.b) .+ (1.0 ./ df.c) .+ (1.0 ./ df.d))
    df.lpv = pvalue.(ft, tail=:left)
    df.bpv = pvalue.(ft, tail=:both)
    df.rpv = pvalue.(ft, tail=:right)

    df
end
###
catcoaccess(dfs) = combine(groupby(reduce(vcat, dfs), [:chromA, :startA, :stopA, :chromB, :startB, :stopB]), :a => sum => :a, :b => sum => :b, c => :c, :d => :d)

# cat /Users/ndlo201/julia/dev/SMGStats/.vscode/settings.json

"""
    coaccess_wrapper(files, peakfile, output; chromlabel=:chrom, startlabel=:start, stoplabel=:stop)

Evaluates coaccessibility for a set of BAM files against a set of peaks.
Combines results and writes to output.
"""
function coaccess_wrapper(files::Vector{String}, peakfile::String, output::String; chromlabel=:chrom, startlabel=:start, stoplabel=:stop, n_thresh=10, calcfisher=true)

    # Load intervals using robust loadbed -> bedintervals pipeline
    bdf = loadbed(peakfile)
    intervals = bedintervals(bdf; chromlabel=Symbol(chromlabel), startlabel=Symbol(startlabel), stoplabel=Symbol(stoplabel))

    smglog("Loaded $(length(intervals)) intervals from $peakfile")

    dfs = DataFrame[]
    for file in files
        smglog("Processing coaccessibility for: $file")
        # Assuming coaccess returns a Dict
        cod = coaccess(file, intervals)
        # Convert to DataFrame
        df = coaccessdf(cod, false)
        push!(dfs, df)
    end


    combineddf = catcoaccess(dfs)

    if calcfisher
        combineddf = coaccessfisher(combineddf, n_thresh)
    end

    smglog("Writing results to: $output")
    if endswith(output, ".gz")
        CSV.write(output, combineddf, delim='\t', compress=true)
    else
        CSV.write(output, combineddf, delim='\t')
    end
end

