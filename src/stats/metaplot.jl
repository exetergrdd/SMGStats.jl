############################# metaplot  #############################
mutable struct MetaPlot{S} <: RecordStat
    fg::Vector{Matrix{Int}}
    bg::Vector{Matrix{Int}}
    cic::ChromIntervalCollection{S}
    
    tid::Int
    regionpointer::Int
    firstintersection::Int
    lastintersection::Int
end


function instantiate(::Type{MetaPlot}, config)

    regions = chromintervals_to_tid(loadmultichromintervals(config.MetaPlot), config.reader)
    
    fg = [zeros(Int, 2*w + 1, length(instances(Modification))) for w in regions.dw]
    bg = [zeros(Int, 2*w + 1, length(instances(Modification))) for w in regions.dw]
    
    
    MetaPlot(fg, bg, regions, -1, 0, typemax(Int), typemin(Int))
end


recordupdates(::Type{<:MetaPlot})  = RecordUpdates
modupdates(::Type{<:MetaPlot})     = ModUpdates
postmodupdates(::Type{<:MetaPlot}) = nothing
statname(::Type{MetaPlot}) = "Metaplot"

@inline function updaterecord!(stat::MetaPlot, record::BamRecord, recorddata)
    ### intersect read with intervals 
    lp = SMGReader.leftposition(record) + 1 ### convert to 1-based inclusive
    rp = SMGReader.rightposition(recorddata)
    tid = record.core.tid

    
    if tid ∈ keys(stat.cic.intervals)
        if stat.tid != tid
            stat.tid = tid
            stat.regionpointer = 1
        end
        

        reg = stat.cic.intervals[tid]
        
        numregions = length(reg)
        stat.firstintersection = typemax(Int)
        stat.lastintersection = typemin(Int)

        for i = stat.regionpointer:numregions  
            ci = reg[i]

            if ci.stop < lp
                ## interval finishes before read advance
                stat.regionpointer += 1    
            elseif rp < ci.start
                ## read finishes before region
                break
            else
                # read intersects
                if stat.firstintersection == typemax(Int)
                    stat.firstintersection = i
                end
                stat.lastintersection = i
            end
        end
    else
        stat.tid = -1
        stat.firstintersection = 1
        stat.lastintersection = 0
    end
    nothing
end


@inline function updatemod!(stat::MetaPlot, mod::ModificationInfo, record::BamRecord, recorddata)

    (stat.firstintersection > stat.lastintersection ) && return nothing
    thresh = 0.9*255
    
    genomepos = recorddata.alignmap[mod.pos]
    if !iszero(genomepos)
        for ri in stat.firstintersection:stat.lastintersection
            ci = stat.cic.intervals[stat.tid][ri]
            # @show ci
            bg = stat.bg[ci.group]
            fg = stat.fg[ci.group]
            pile_i = genomepos - ci.start + 1
            
            if ci.strand
                pile_i = genomepos - ci.start + 1
            else
                pile_i = ci.stop - genomepos + 1
            end

            ### strand
            if 1 <= pile_i <= size(bg, 1)
                bg[pile_i, Int(mod.mod) + 1] += 1
                if mod.prob > thresh
                    fg[pile_i, Int(mod.mod) + 1] += 1
                end
            end
        end
    end
    nothing
end

###### as we will be streaming in order change to do sorted scan

function writestats(stat::MetaPlot, path::String, file="metaplot.tsv.gz")
    filepath = joinpath(path, file)

    mods = instances(Modification)

    df = DataFrame(Region=String[], Mod=[], xp=Int[], FG=Int[], BG=Int[])

    for i = 1:length(stat.cic.grouplabels)
        fg = stat.fg[i]
        bg = stat.bg[i]
        modind = vec(sum(bg, dims=1) .> 0)
        dw = stat.cic.dw[i]
        xp = -dw:dw

        mdf = mapreduce((f, b, m) -> DataFrame(Region=stat.cic.grouplabels[i], Mod=m, xp=xp, FG=f, BG=b), vcat, eachcol(fg[:, modind]), eachcol(bg[:, modind]), mods[modind])
        append!(df, mdf)
    end

    CSV.write(filepath, df, delim='\t', compress=endswith(file, ".gz"))
end
function plotstat(::Type{MetaPlot}, df::DataFrame)
    # f = Figure()
    # ax = Axis(f[1, 1], xgridvisible=false, ygridvisible=false, ylabel="Correlation")

    df.R = df.FG./df.BG

    transform!(groupby(df, [:Region, :Mod]), :R => (x -> (x .- mean(x))/std(x)) => :R)
    # df.Label = string.(replace.(df.ModA, "mod_" => ""), "\n", replace.(df.ModB, "mod_" => ""))

    plt = data(df) * mapping(:xp, :R, color=:Mod) * mapping(col=:Mod, row=:Region) * (visual(Scatter, alpha=0.1, marker=:vline) + smooth(span=0.05, degree=4, npoints=3000))
    f = draw(plt; figure=(; size=(1000, 500)), legend=(; position=:bottom))
    f
end
