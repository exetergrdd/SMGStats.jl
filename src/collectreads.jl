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

const _GENERIC_BAM_STEMS = Set(["", "basecall", "reads", "calls", "aligned", "alignments", "output", "bam", "cram", "pass", "fail"])
const _GENERIC_DIR_NAMES = r"^(dorado.*|epi2me.*|fire.*|results|bam|cram|alignment|alignments|output|data|pass|fail)$"i
const _PIPELINE_SUFFIX_REGEX = r"(\.(bam|cram|haplotagged|fire|ftnucs|calledmods|mods|sort|sorted|hg38|chm13|grch38|t2t|mm10|mm39))+$"i

"""
    infer_sample_label(path::AbstractString) -> String

Extracts an informative sample identifier from a BAM/CRAM file path or label string.
Recognises standard ONT / Dorado directory layouts:
- `*/sample_id/sample_id_rep/run_ind/dorado*/basecall*.bam` -> `sample_id_rep`
- `*/sample_id/sample_id_rep/run_ind/dorado*/epi2me_*/results/sample_id/sample_id*.bam` -> `sample_id_rep`
- `*/sample_id/sample_id_rep/run_ind/dorado*/fire_*/results/sample_id/sample_id*.cram` -> `sample_id_rep`
- `*/results/sample_id/*.bam` -> `sample_id`
Falls back to stripping standard pipeline suffixes (`.bam`, `.cram`, `.haplotagged`, `.fire`, `.ftnucs`,
`.calledmods`, `.sort`, `.hg38`, `basecall.`), and if the filename is generic (`basecall*.bam`),
uses the nearest informative parent directory name.
"""
function infer_sample_label(path::AbstractString)::String
    s = strip(String(path))
    isempty(s) && return "sample"

    norm_s = replace(s, '\\' => '/')
    parts = split(norm_s, '/'; keepempty=false)
    isempty(parts) && return "sample"

    # Already a clean label (no directory separators and no BAM/CRAM extension)
    if length(parts) == 1 && !occursin(r"\.(bam|cram)$"i, parts[1])
        return String(parts[1])
    end

    # 1. Check for `dorado*` directory in path:
    #    */sample_id/sample_id_rep/run_ind/dorado*/...
    if length(parts) >= 2
        d_idx = findlast(p -> occursin(r"^dorado"i, p), @view parts[1:end-1])
        if !isnothing(d_idx)
            if d_idx >= 3
                return String(parts[d_idx - 2])
            elseif d_idx == 2
                # Check if there is a downstream results/<sample_id>/<file> before falling back to parts[1]
                r_idx = findlast(p -> lowercase(p) == "results", @view parts[d_idx+1:end-1])
                if !isnothing(r_idx) && (d_idx + r_idx + 1) <= length(parts) - 1
                    return String(parts[d_idx + r_idx + 1])
                end
                return String(parts[1])
            end
        end

        # 2. Check for `.../results/<sample_id>/<file>` layout without dorado* prefix
        r_idx = findlast(p -> lowercase(p) == "results", @view parts[1:end-2])
        if !isnothing(r_idx)
            return String(parts[r_idx + 1])
        end
    end

    # 3. Clean filename stem
    fname = String(parts[end])
    stem = replace(fname, _PIPELINE_SUFFIX_REGEX => "")
    stem = replace(stem, r"^basecall\.?"i => "")
    if !isempty(stem) && !(lowercase(stem) in _GENERIC_BAM_STEMS)
        return stem
    end

    # 4. Filename is generic (e.g. basecall.hg38.sort.calledmods.ftnucs.fire.haplotagged.bam):
    #    walk up parent directories to find the first informative directory name
    for idx in (length(parts) - 1):-1:1
        dir_name = String(parts[idx])
        if !occursin(_GENERIC_DIR_NAMES, dir_name)
            return dir_name
        end
    end

    return length(parts) >= 2 ? String(parts[end - 1]) : (isempty(stem) ? fname : stem)
end

"""
    parse_sample_arg(samples_str::Union{Nothing, AbstractString})

Parses a comma-separated `--samples` CLI string into either:
- `nothing` (if empty or `nothing`),
- a `Vector{Pair{String,String}}` if `pattern=label` syntax is used (e.g. `"050526=endo_as1,200126=endo_as2"`), or
- a `Vector{String}` of positional sample labels (e.g. `"endo_as1,endo_as2,endo_wgs"`).
"""
function parse_sample_arg(samples_str::Union{Nothing, AbstractString})
    isnothing(samples_str) && return nothing
    s = strip(String(samples_str))
    isempty(s) && return nothing
    tokens = [ String(strip(t)) for t in split(s, ',') if !isempty(strip(t)) ]
    isempty(tokens) && return nothing
    if all(t -> occursin('=', t), tokens)
        pairs_vec = Pair{String, String}[]
        for t in tokens
            k, v = split(t, '='; limit=2)
            push!(pairs_vec, String(strip(k)) => String(strip(v)))
        end
        return pairs_vec
    end
    return tokens
end

function _disambiguate_colliding_labels!(labels::Vector{String}, raw_paths::AbstractVector{<:AbstractString})
    n = length(labels)
    n <= 1 && return labels

    # Group indices by label for distinct raw_paths
    groups = Dict{String, Vector{Int}}()
    for i in 1:n
        push!(get!(groups, labels[i], Int[]), i)
    end

    for (base_lbl, idxs) in groups
        distinct_paths = unique([String(raw_paths[i]) for i in idxs])
        length(distinct_paths) <= 1 && continue

        # Try disambiguating using the nearest differing parent directory component (e.g. run_ind)
        split_map = Dict{String, Vector{SubString{String}}}(
            p => split(replace(p, '\\' => '/'), '/'; keepempty=false) for p in distinct_paths
        )
        candidate_map = Dict{String, String}()
        min_len = minimum(length(v) for v in values(split_map))
        disambiguated = false
        if min_len >= 2
            for offset in 1:(min_len - 1)
                parts_at_offset = [String(split_map[p][end - offset]) for p in distinct_paths]
                if length(unique(parts_at_offset)) == length(distinct_paths) && all(p -> p != base_lbl, parts_at_offset)
                    for (p, diff_part) in zip(distinct_paths, parts_at_offset)
                        short_diff = length(diff_part) > 24 ? first(diff_part, 24) : diff_part
                        candidate_map[p] = "$(base_lbl)_$(short_diff)"
                    end
                    disambiguated = true
                    break
                end
            end
        end

        if !disambiguated || length(unique(values(candidate_map))) < length(distinct_paths)
            for (rank, p) in enumerate(distinct_paths)
                candidate_map[p] = "$(base_lbl)_$(rank)"
            end
        end

        for i in idxs
            labels[i] = candidate_map[String(raw_paths[i])]
        end
    end
    return labels
end

"""
    infer_sample_labels(paths::AbstractVector{<:AbstractString}; samples=nothing) -> Vector{String}

Maps a list of BAM/CRAM file paths to unique sample labels.
- If `samples` is a `Vector{<:AbstractString}` of the same length as `paths`, uses those labels directly.
- If `samples` is a `Vector{Pair{String,String}}` (`pattern => label`), matches each path against `pattern`
  and falls back to `infer_sample_label(path)`.
- Guarantees that distinct input paths in `paths` always map to distinct sample labels.
"""
function infer_sample_labels(paths::AbstractVector{<:AbstractString}; samples=nothing)::Vector{String}
    parsed = samples isa AbstractString ? parse_sample_arg(samples) : samples
    n = length(paths)
    labels = Vector{String}(undef, n)

    if parsed isa AbstractVector{<:AbstractString} && !isempty(parsed)
        if length(parsed) != n
            throw(ArgumentError("Number of provided sample labels ($(length(parsed)): $(join(parsed, ", "))) does not match number of BAM/CRAM files ($n)."))
        end
        for i in 1:n
            labels[i] = String(parsed[i])
        end
    elseif parsed isa AbstractVector{<:Pair}
        for i in 1:n
            p = String(paths[i])
            inferred = infer_sample_label(p)
            matched = nothing
            for (pat, target) in parsed
                if occursin(pat, p) || occursin(pat, inferred)
                    matched = String(target)
                    break
                end
            end
            labels[i] = isnothing(matched) ? inferred : matched
        end
    else
        for i in 1:n
            labels[i] = infer_sample_label(paths[i])
        end
    end

    _disambiguate_colliding_labels!(labels, paths)
    return labels
end

"""
    resolve_read_samples(reads; samples=nothing, labelfun=nothing) -> Vector{String}

Resolves the sample label for each read in `reads` (where each element has a `.label` property),
guaranteeing that distinct raw BAM paths in `r.label` map to unique sample names.
"""
function resolve_read_samples(reads; samples=nothing, labelfun=nothing)::Vector{String}
    n = length(reads)
    n == 0 && return String[]

    raw_labels = String[String(r.label) for r in reads]
    if !isnothing(labelfun) && isnothing(samples)
        return String[String(labelfun(raw)) for raw in raw_labels]
    end

    # Sort unique raw labels so positional/disambiguated mapping is deterministic regardless of read order
    unique_raw = sort!(unique(raw_labels))
    resolved_unique = if samples isa AbstractVector{<:AbstractString} && length(samples) != length(unique_raw)
        # If a locus only has reads from a subset of BAMs and positional samples were given without patterns,
        # infer from paths directly so we never misassign labels
        infer_sample_labels(unique_raw; samples=nothing)
    else
        infer_sample_labels(unique_raw; samples=samples)
    end

    label_map = Dict{String, String}(raw => res for (raw, res) in zip(unique_raw, resolved_unique))
    return String[label_map[raw] for raw in raw_labels]
end


function collectstencilreads(file::AbstractString, chrom::AbstractString, loc::UnitRange{Int}; byhap=false, bystrand=false, strictintersect=0, strictfirefeatures=0, mod_6mA_thresh=0.1, assign_no_hap=false, label::String="", labels=nothing)
    file_label = if !isnothing(labels) && !isempty(labels)
        String(first(labels))
    else
        String(file)
    end
    reader = open(ModGroupFileReader, file; label=file_label)
    reads = collectstencilreads(chrom, loc, reader, byhap=byhap, bystrand=bystrand, strictintersect=strictintersect, strictfirefeatures=strictfirefeatures, mod_6mA_thresh=mod_6mA_thresh, assign_no_hap=assign_no_hap, label=label)
    close(reader)
    reads
end

function collectstencilreads(files::AbstractVector{<:AbstractString}, chrom::AbstractString, loc::UnitRange{Int}; byhap=false, bystrand=false, strictintersect=0, strictfirefeatures=0, mod_6mA_thresh=0.1, assign_no_hap=false, label::String="", labels=nothing)
    file_labels = if !isnothing(labels) && !isempty(labels)
        infer_sample_labels(files; samples=labels)
    else
        String.(files)
    end
    reader = open(ModGroupFileReader, files; labels=file_labels)
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