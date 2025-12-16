

function modtable(file, chrom, loc)

    reader = open(HTSFileReader, file)
    recorddata = StencillingData(AuxMapMod())

    reads = String[]
    tid = Int[]
    pos = Int[]
    strand = String[]
    prob = Int[]
    genomepos = Int[]
    mod = Modification[]

    for record in eachintersection(reader, chrom, loc)
        validflag(record) || continue
        processread!(record, recorddata) || continue

        for mi in ModIterator(record, recorddata)
            gp = recorddata.alignmap[mi.pos]

            push!(reads, qname(record))
            push!(tid, record.core.tid)
            push!(pos, mi.pos)
            push!(strand, ifelse(mi.strand, "+", "-"))
            push!(prob, mi.prob)
            push!(genomepos, gp)
            push!(mod, mi.mod)
            
        end
    end
    tid_to_chrom = referencedict(reader)
    chroms = getindex.(Ref(tid_to_chrom), tid)
    DataFrame((; reads, mod, pos, strand, prob, chroms, genomepos))
end



### metaplots heatmaps etc.

# iv(s, e) = s:e
# function loadregions(file, dw=1500)
#     df = CSV.read(file, DataFrame)
#     strand = "strand" ∈ names(df) ? df.strand : fill("+", size(df, 1))
#     String.(df.chrom), [div(s + e, 2) .+ (-dw:dw) for (s, e) in zip(df.start, df.stop)], strand
# end


# function metaplot(bedfile, file, dw=1500, thresh = 0.9*255)
   
#     chroms, locs, strands = loadregions(bedfile)

#     # n = length(first(locs))
#     xp = -dw:dw
#     m = length(xp)
#     fg = zeros(m)
#     bg = zeros(m)
#     reader = open(HTSFileReader, file)
#     recorddata = StencillingData(AuxMapMod())


#     @showprogress for (chrom, loc, strand) in zip(chroms, locs, strands)
#         for record in eachintersection(reader, chrom, loc)
#             processread!(record, recorddata)
#             for mod in ModIterator(record, recorddata)
#                 genomepos = recorddata.alignmap[mod.pos]
#                 if !iszero(genomepos)
#                     if strand == "+"
#                         i = genomepos - first(loc) + 1
#                     else
#                         i = last(loc) - genomepos + 1
#                     end
#                     if 1 <= i <= m
#                         bg[i] += 1
#                         if mod.prob > thresh
#                             fg[i] += 1
#                         end
#                     end
#                 end
#             end
#         end
#     end
#     close(reader)
#     fg, bg
# end