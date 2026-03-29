
mutable struct ActiveSweep{S}
    intervals::ChromIntervalCollection{S}
    tid::Int
    regionpointer::Int
    firstintersection::Int
    lastintersection::Int
end

function instantiate_sweep(intervals::ChromIntervalCollection{String}, reader)
    # Convert String keys to TID keys using existing regions.jl function
    tid_intervals = chromintervals_to_tid(intervals, reader)
    ActiveSweep(tid_intervals, -1, 1, typemax(Int), typemin(Int))
end

# If already TID keys (Int), just wrap
function instantiate_sweep(intervals::ChromIntervalCollection{Int}, reader)
    ActiveSweep(intervals, -1, 1, typemax(Int), typemin(Int))
end


function countfire_sweep(file, intervals::ChromIntervalCollection{S}; total_reads=-1, fire_backing=5_000) where S

    reader = open(HTSFileReader, file)
    recorddata = StencillingData(AuxMapModFire())

    sweep = instantiate_sweep(intervals, reader)
    
    # Calculate total regions. 
    # ChromIntervalCollection doesn't store total count directly except implicitly.
    # We can iterate to find max regionindex or assume user passed max index.
    # Regions are indexed 1..N in bedchromintervals.
    # Finding max index:
    max_index = 0
    for v in values(intervals.intervals)
        for ci in v
            if ci.regionindex > max_index
                max_index = ci.regionindex
            end
        end
    end

    fire_counts = zeros(Int, max_index)
    read_counts = zeros(Int, max_index)

    t = 0
    for record in eachrecord(reader)
        validflag(record) || continue
        
        tid = record.core.tid
        # skip if tid not in intervals
        if !haskey(sweep.intervals.intervals, tid)
            continue
        end

        processread!(record, recorddata) || continue
        
        lp = SMGReader.leftposition(record) + 1
        rp = SMGReader.rightposition(recorddata)
        
        # Sorted sweep logic
        if sweep.tid != tid
            sweep.tid = tid
            sweep.regionpointer = 1
        end

        reg = sweep.intervals.intervals[tid]
        numregions = length(reg)
        sweep.firstintersection = typemax(Int)
        sweep.lastintersection = typemin(Int)

        for i = sweep.regionpointer:numregions
            ci = reg[i]

            if ci.stop < lp # interval finishes before read starts
                sweep.regionpointer += 1
            elseif rp < ci.start # read finishes before interval starts
                break
            else
                # overlaps
                if sweep.firstintersection == typemax(Int)
                    sweep.firstintersection = i
                end
                sweep.lastintersection = i
            end
        end


        # If we found intersections
        if sweep.firstintersection <= sweep.lastintersection
            
            for ri in sweep.firstintersection:sweep.lastintersection
                iv = reg[ri]
                k = iv.regionindex
                read_counts[k] += 1
            end

            # Now check for fire MSPs
            for (s, l, q) in firemsps(record, recorddata)
                if q >= 0.9 * 255
                    genstart, genstop = genomecoords(s, l, record, recorddata)
                    if !iszero(genstart) && !iszero(genstop)
                        
                        for ri in sweep.firstintersection:sweep.lastintersection
                            iv = reg[ri]
                            if (genstart <= iv.stop) && (genstop >= iv.start)
                                k = iv.regionindex
                                fire_counts[k] += 1
                            end
                        end
                    end
                end
            end
        end

        t += 1
        (t == total_reads) && break
    end

    close(reader)
    fire_counts, read_counts
end
