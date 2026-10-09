struct ModAllFireRead
    chrom::String
    lp::Int
    rp::Int
    strand::Bool
    mods_6mA::Vector{Int}
    mods_5mC::Vector{Int}
    mods_5hmC::Vector{Int}
    unmods_6mA::Vector{Int}
    unmods_5mC::Vector{Int}
    unmods_5hmC::Vector{Int}
    nucs::Vector{NTuple{2,Int}}
    msps::Vector{NTuple{3,Int}}
    haplotype::UInt8
    label::String
end


mutable struct ModGroupFileReader{T}
    readers::Vector{HTSFileReader}
    labels::Vector{String}
    recorddata::T
    buff_mods_6mA::VectorBuffer{Int}
    buff_mods_5mC::VectorBuffer{Int}
    buff_mods_5hmC::VectorBuffer{Int}
    buff_umat_6mA::VectorBuffer{Int}
    buff_umat_5mC::VectorBuffer{Int}
    buff_umat_5hmC::VectorBuffer{Int}
    
end

@inline Base.open(::Type{ModGroupFileReader}, path::AbstractString; label=path, am=AuxMapModFire, init_buffer_length = 10_000) = opengroupreader([path], labels=[label], am=am, init_buffer_length=init_buffer_length)
@inline Base.open(::Type{ModGroupFileReader}, paths::Vector{T}; labels=paths, am=AuxMapModFire, init_buffer_length = 10_000) where {T} = opengroupreader(paths, labels=labels, am=am, init_buffer_length=init_buffer_length)

function opengroupreader(files::Vector{T}; labels=["" for f in files], am=AuxMapModFire, init_buffer_length = 10_000) where {T}

    readers = [open(HTSFileReader, file) for file in files]

    if all(isempty, labels)
        labels = files
    end

    recorddata = StencillingData(am())
    buff_mods_6mA  = VectorBuffer{Int}(init_buffer_length)
    buff_mods_5mC  = VectorBuffer{Int}(init_buffer_length)
    buff_mods_5hmC = VectorBuffer{Int}(init_buffer_length)
    buff_umat_6mA  = VectorBuffer{Int}(init_buffer_length)
    buff_umat_5mC  = VectorBuffer{Int}(init_buffer_length)
    buff_umat_5hmC = VectorBuffer{Int}(init_buffer_length)

    ModGroupFileReader(readers, labels, recorddata, buff_mods_6mA, buff_mods_5mC, buff_mods_5hmC, buff_umat_6mA, buff_umat_5mC, buff_umat_5hmC)
end

Base.close(reader::ModGroupFileReader) = close.(reader.readers)


function collectstencilreads(file, chrom::AbstractString, loc::UnitRange{Int}; byhap=false, bystrand=false, strictintersect=0, strictfirefeatures=0, mod_6mA_thresh=0.1, assign_no_hap=false, label::String="")
    reader = open(ModGroupFileReader, file)
    reads = collectstencilreads(chrom, loc, reader, byhap=byhap, bystrand=bystrand, strictintersect=strictintersect, strictfirefeatures=strictfirefeatures, mod_6mA_thresh=mod_6mA_thresh, assign_no_hap=assign_no_hap, label=label)
    close(reader)
    reads
end
    


function collectstencilreads(chrom::AbstractString, loc::UnitRange{Int}, reader::ModGroupFileReader ; byhap=false, bystrand=false, strictintersect=0, strictfirefeatures=0, mod_6mA_thresh=0.1, assign_no_hap=false, label::String="")

    reads = Vector{ModAllFireRead}()

    prob_threshold = round(Int, 0.9 * 255)
    recorddata = reader.recorddata

    for (filereader, readerlabel) in zip(reader.readers, reader.labels)
        for (k, record) in enumerate(eachintersection(filereader, chrom, loc))
            validflag(record) || continue
            processread!(record, recorddata) || continue

            if strictintersect > 0
                lp = SMGReader.leftposition(record) + 1
                rp = SMGReader.rightposition(recorddata)
                ((lp + strictintersect < first(loc)) && (last(loc) < rp + strictintersect)) || continue
            end

            firegenomecoords = firegenome(record, recorddata)
            nucs = collect(firegenomecoords(firenucs(record, recorddata)))
            msps = collect(firegenomecoords(firemsps(record, recorddata)))

            if strictfirefeatures > 0
                featureleft = false
                featureright = false
                
                for (s, e) in nucs
                    featureleft  |= (s - strictfirefeatures) < first(loc)
                    featureright |= (e + strictfirefeatures) > last(loc)
                end
                if !(featureleft && featureright)
                    for (s, e, q) in msps
                        featureleft  |= (s - strictfirefeatures) < first(loc)
                        featureright |= (e + strictfirefeatures) > last(loc)     
                    end
                end
                !(featureleft && featureright) && continue
            end



            empty!(reader.buff_mods_6mA)
            empty!(reader.buff_mods_5mC)
            empty!(reader.buff_mods_5hmC)
            empty!(reader.buff_umat_6mA)
            empty!(reader.buff_umat_5mC)
            empty!(reader.buff_umat_5hmC)
            
            total_6mA_mod = 0
            total_6mA_unmod = 0

            for mi in ModIterator(record, recorddata)
                genpos = genomecoords(mi.pos, record, recorddata, onebased=true)
                if genpos > 1
                    if mi.prob > prob_threshold
                        mi.mod == mod_6mA  && push!(reader.buff_mods_6mA, genpos)
                        mi.mod == mod_5mC  && push!(reader.buff_mods_5mC, genpos)
                        mi.mod == mod_5hmC && push!(reader.buff_mods_5hmC, genpos)
                        (mi.mod == mod_6mA) && (total_6mA_mod += 1)
                    else
                        mi.mod == mod_6mA  && push!(reader.buff_umat_6mA, genpos)
                        mi.mod == mod_5mC  && push!(reader.buff_umat_5mC, genpos)
                        mi.mod == mod_5hmC && push!(reader.buff_umat_5hmC, genpos)
                        (mi.mod == mod_6mA) && (total_6mA_unmod += 1)
                    end
                end
            end
            
            total_6mA_sites = total_6mA_mod + total_6mA_unmod
            frac_6mA_mod = total_6mA_sites > 0 ? (total_6mA_mod / total_6mA_sites) : 0.0
            frac_6mA_mod > mod_6mA_thresh && continue ### ignore if overmethylated

            hap = haplotype(record, recorddata, ifelse(assign_no_hap, ifelse(iseven(k), 0x01, 0x02), 0x00))

            fr = ModAllFireRead(
                chrom,
                SMGReader.leftposition(record), 
                SMGReader.rightposition(recorddata), 
                ispositive(record), 
                copy(reader.buff_mods_6mA), 
                copy(reader.buff_mods_5mC), 
                copy(reader.buff_mods_5hmC), 
                copy(reader.buff_umat_6mA), 
                copy(reader.buff_umat_5mC), 
                copy(reader.buff_umat_5hmC), 
                nucs, 
                msps, 
                hap,
                string(label, readerlabel)
            )
            push!(reads, fr)
        end
        
        
        sort!(reads, by = fr -> fr.lp - fr.rp)
    end
    
    return reads
end