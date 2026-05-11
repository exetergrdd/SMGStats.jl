



"""
    countfire(files, peakfile, output; chromlabel=:chrom, startlabel=:start, stoplabel=:stop)

Counts fires in specified peaks for one or more BAM files.
Writes results to output file.
"""
function countfire(files::Vector{String}, peakfile::String, output::String; chromlabel=:chrom, startlabel=:start, stoplabel=:stop)

    beddf = loadbed(peakfile)
    intervals = bedintervals(bdf; chromlabel=Symbol(chromlabel), startlabel=Symbol(startlabel), stoplabel=Symbol(stoplabel))

    smglog("Loaded $(length(intervals)) intervals from $peakfile")


    # Ensure result_df has the standard columns if we want to write valid bed-like output
    # But for analysis, just appending counts is cleaner.

    # Just to be safe, let's process each file

    for file in files
        smglog("Counting fires for: $file")
        fc, rc = countfire(file, intervals)

        samplename = basename(file)
        beddf[!, Symbol(samplename, "_fire")] = fc
        beddf[!, Symbol(samplename, "_read")] = rc
    end

    smglog("Writing results to: $output")
    if endswith(output, ".gz")
        CSV.write(output, beddf, delim='\t', compress=true)
    else
        CSV.write(output, beddf, delim='\t')
    end

end


function countfire(file, intervals::GenomicIntervalCollection{GenomicInterval{T}}; total_reads=-1, fire_backing=5_000) where T

    reader = open(HTSFileReader, file)
    recorddata = StencillingData(AuxMapModFire())

    ### fireinteresctions is a buffer of those intervals each read intersects with
    fireintersections = VectorBuffer{GenomicInterval{T}}(fire_backing)
    intersects = VectorBuffer{Bool}(fire_backing)


    fire_counts = zeros(Int(length(intervals)))
    read_counts = zeros(Int(length(intervals)))

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
            k = GenomicFeatures.metadata(fiv)
            read_counts[k] += 1
        end
        setlength!(fireintersections, fi)
        setlength!(intersects, fi)
        intersects .= false
        ### ensure that each fire is only counted once per read
        for (s, l, q) in firemsps(record, recorddata)
            if q >= 0.9 * 255
                genstart, genstop = genomecoords(s, l, record, recorddata)
                if !iszero(genstart) && !iszero(genstop)

                    @inbounds @simd for i = eachindex(fireintersections)
                        if (genstart <= fireintersections[i].last) && (genstop >= fireintersections[i].first)
                            # k = GenomicFeatures.metadata(fireintersections[i])
                            # fire_counts[k] += 1
                            intersects[i] = true
                        end
                    end
                end
            end
        end

        for i = eachindex(intersects)
            if intersects[i]
                k = GenomicFeatures.metadata(fireintersections[i])
                fire_counts[k] += 1
            end
        end

        t += 1
        (t == total_reads) && break
    end
    close(reader)

    fire_counts, read_counts

end






function countmods(file, mod, intervals::GenomicIntervalCollection{GenomicInterval{T}}; total_reads=-1, fire_backing=5_000) where T

    reader = open(HTSFileReader, file)
    recorddata = StencillingData(AuxMapModFire())

    ### fireinteresctions is a buffer of those intervals each read intersects with
    fireintersections = VectorBuffer{GenomicInterval{T}}(fire_backing)


    mod_counts   = zeros(Int(length(intervals)))
    unmod_counts = zeros(Int(length(intervals)))

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
            k = GenomicFeatures.metadata(fiv)
            # read_counts[k] += 1
        end
        setlength!(fireintersections, fi)

        for mi in ModIterator(record, recorddata)
            (mi.mod == mod) || continue
            genomepos = genomecoords(mi.pos, record, recorddata)
            @inbounds @simd for i = eachindex(fireintersections)
                if (genomepos <= fireintersections[i].last) && (genomepos >= fireintersections[i].first)
                    k = GenomicFeatures.metadata(fireintersections[i])
                    if mi.prob >= 0.9 * 255
                        mod_counts[k] += 1
                    else
                        unmod_counts[k] += 1
                    end
                end
            end


        end

        t += 1
        (t == total_reads) && break
    end
    close(reader)

    mod_counts, unmod_counts

end



function countfiremods(file, intervals::GenomicIntervalCollection{GenomicInterval{T}}; total_reads=-1, fire_backing=5_000) where T

    reader = open(HTSFileReader, file)
    recorddata = StencillingData(AuxMapModFire())

    ### fireinteresctions is a buffer of those intervals each read intersects with
    fireintersections = VectorBuffer{GenomicInterval{T}}(fire_backing)
    intersects = VectorBuffer{Bool}(fire_backing)

    fields = [:fire, :mod_6mA, :mod_5mC, :mod_5hmC]
    fg_counts = Dict(f => zeros(Int, length(intervals)) for f in fields)
    total_counts = Dict(f => zeros(Int, length(intervals)) for f in fields)

    

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
            k = GenomicFeatures.metadata(fiv)
            total_counts[:fire][k] += 1
        end
        setlength!(fireintersections, fi)
        setlength!(intersects, fi)
        intersects .= false
        ### ensure that each fire is only counted once per read
        for (s, l, q) in firemsps(record, recorddata)
            if q >= 0.9 * 255
                genstart, genstop = genomecoords(s, l, record, recorddata)
                if !iszero(genstart) && !iszero(genstop)

                    @inbounds @simd for i = eachindex(fireintersections)
                        if (genstart <= fireintersections[i].last) && (genstop >= fireintersections[i].first)
                            # k = GenomicFeatures.metadata(fireintersections[i])
                            # fire_counts[k] += 1
                            intersects[i] = true
                        end
                    end
                end
            end
        end

        for mi in ModIterator(record, recorddata)
            
            genomepos = genomecoords(mi.pos, record, recorddata)
            
            @inbounds @simd for i = eachindex(fireintersections)
                if (genomepos <= fireintersections[i].last) && (genomepos >= fireintersections[i].first)
                    k = GenomicFeatures.metadata(fireintersections[i])
                    if mi.prob >= 0.9 * 255
                        fg_counts[Symbol(mi.mod)][k] += 1
                    end
                    total_counts[Symbol(mi.mod)][k] += 1
                end
            end
        end

        for i = eachindex(intersects)
            if intersects[i]
                k = GenomicFeatures.metadata(fireintersections[i])
                fg_counts[:fire][k] += 1
            end
        end

        t += 1
        (t == total_reads) && break
    end
    close(reader)

    fg_counts, total_counts

end
