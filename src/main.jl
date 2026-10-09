using ArgParse

function default_config(genome)

    if gethostname() == "penrose"
        return "/penrose/projects/ont/scripts/smgstats/smgstats.config.$(genome).yaml"
    else
        error("Default not defined for $genomee")
    end

end

statreport(bamfile, genome="hg19"; yaml=default_config(genome), nr=-1) = runstats(bamfile, yaml, nr=nr)

function runstats(file, yaml; nr=-1, statdirfilename=false, outputdir="")

    start = time()

    if !isempty(outputdir)
        statdir = outputdir
    elseif statdirfilename
        statdir = string("stats_", join(split(basename(file), ".")[1:end-1], "."))
    else
        statdir = "stats"
    end


    #### print banner
    banner()
    smglog("File                   : ", file)
    smglog("YAML                   : ", yaml)

    ### this is all a bit unclean htsdata could be inferred from the auxmap
    auxmap = autodetectaux(file)

    if auxmap isa AuxMapModFire
        auxmap = AuxMapModFireQC()
    end
    analysis, statset, config = loadstatconfig(yaml, auxmap)

    ##### auto detect file details
    # htsdata = autodetecthtsdata(file) 
    if analysis == "stencilling"
        htsdata = StencillingData
    elseif analysis == "directrna"
        htsdata = DirectRNA
    else
        error("unrecognised analysis: $analysis")
    end

 

    mods = autodetectmods(file)
    recorddata = htsdata(auxmap)

    smglog("Analysis               : ", htsdata)
    smglog("Detected Aux map       : ", typeof(auxmap))
    smglog("Detected Modifications : ", mods)

    smglog("Config                 : ", yaml)
    smglog("Stats                  : ", replace(string(statset), "Tuple{" => "", "}" => ""))
    smglog("Outputdir              : ", statdir)

    if nr != -1
        smglog("Sampling               : ", nr, " reads")
    end



    if hasproperty(config, :MetaPlot)
        for mp in config.MetaPlot
            smglog("Metaplot               : ", mp.label, "\t: ", mp.width, " : ", mp.file)
        end
    end
    smglog("Calculating stats:")

    stats = calculatestats(file, statset, recorddata, mods; nr=nr, config=config)
    smglog("Writing Stats:")
    writeallstats(stats, bamfile=file, statdir=statdir)
    smglog("Generating Report:")
    stats = readstats(bamfile=file, statdir=statdir)
    reportfile = htmlreport(stats)
    smglog("Written Report    : ", reportfile)
    smglog("Complete in       : ", time() - start, " seconds")
end

function banner()

    ### generated with https://patorjk.com/software/taag Future font
    s = """
    ======================================================================
    ┏━┓┏┳┓┏━╸┏━┓╺┳╸┏━┓╺┳╸┏━┓  ┏┓╻  
    ┗━┓┃┃┃┃╺┓┗━┓ ┃ ┣━┫ ┃ ┗━┓   ┃┃  
    ┗━┛╹ ╹┗━┛┗━┛ ╹ ╹ ╹ ╹ ┗━┛╹┗━┛┗━╸

    """
    println(s)
end


### interface
### smgstats stats -config yaml [output] bam/cram
### smgstats coaccess -peaks bed -o table.tsv.gz bam/cram(s)
### smgstats firecount -peaks bed bam/cram(s)

function parse_commandline()
    s = ArgParseSettings(prog="smgstats", description="SMGStats: Single Molecule Genomics Statistics", preformatted_description=true, help_width=100)

    @add_arg_table! s begin
        "stats"
        action = :command
        help = "Run standard statistics"
        "coaccess"
        action = :command
        help = "Run coaccessibility analysis"
        "firecount"
            action = :command
            help = "Run fire counting"
        "diff"
            action = :command
            help = "Run differential accessibility analysis"
        "cluster"
            action = :command
            help = "Run stencil clustering on genomic regions"
    end

    @add_arg_table! s["stats"] begin
        "--config", "-c"
            help = "Configuration YAML file"
            required = true
        "--output", "-o"
            help = "Output directory (optional)"
            default = ""
        "bamfile"
            help = "BAM file"
            required = true
    end

    @add_arg_table! s["coaccess"] begin
        "--peaks", "-p"
            help = "Peaks BED file"
            required = true
        "--output", "-o"
            help = "Output file (TSV/CSV)"
            required = true
        "--chrom"
            help = "Chromosome column name"
            default = "chrom"
        "--start"
            help = "Start column name"
            default = "start"
        "--stop"
            help = "Stop column name"
            default = "stop"
        "files"
            nargs = '+'
            help = "BAM/CRAM files"
            required = true
    end

    @add_arg_table! s["firecount"] begin
        "--peaks", "-p"
            help = "Peaks BED file"
            required = true
        "--output", "-o"
            help = "Output file (TSV/CSV)"
            required = true
        "--chrom"
            help = "Chromosome column name"
            default = "chrom"
        "--start"
            help = "Start column name"
            default = "start"
        "--stop"
            help = "Stop column name"
            default = "stop"
        "files"
            nargs = '+'
            help = "BAM/CRAM files"
            required = true
    end

    @add_arg_table! s["diff"] begin
        "--peaks", "-p"
            help = "Peaks BED file"
            required = true
        "--output", "-o"
            help = "Output file (TSV/CSV)"
            required = true
        "--groupA", "-a"
            help = "Group A BAM files"
            nargs = '+'
            required = true
        "--groupB", "-b"
            help = "Group B BAM files"
            nargs = '+'
            required = false ### necessary apart from differential haplotype analysis
        "--method", "-m"
            help = "Statistical method: fisher, binomial_pooled, binomial_ref, all"
            default = "all"
        "--chrom"
            help = "Chromosome column name"
            default = "chrom"
        "--start"
            help = "Start column name"
            default = "start"
        "--stop"
            help = "Stop column name"
            default = "stop"
        "--mods"
            help = "Calculate total differential 6mA, 5mC, 5hmC modification over peak regions"
            action = :store_true
        "--hap"
            help = "Perform differential haplotype analysis on Group A files"
            action = :store_true
    end

    @add_arg_table! s["cluster"] begin
        "--peaks", "-p"
            help = "Peaks BED file"
            required = true
        "--output", "-o"
            help = "Output directory"
            required = true
        "--dw"
            help = "Window around peak to include (bp)"
            arg_type = Int
            default = 600
        "--features", "-f"
            help = "Comma separated list of features to cluster (e.g. fire,mod_6mA,mod_5mC,mod_5hmC)"
            default = "fire,mod_6mA"
        "--dx"
            help = "Spatial binning resolution (bp)"
            arg_type = Int
            default = 15
        "--smooth_sigma"
            help = "Gaussian smoothing sigma for the feature matrix"
            arg_type = Float64
            default = 0.0
        "--strictfirefeatures"
            help = "Number of bp that a read's FIRE features must exceed the cluster region, set to >= 1 ensures reads are annotated over cluster region"
            arg_type = Int
            default = 1
        "--unify_gc_strand"
            help = "Unify GC strand (true/false), true subtracts 1 from position of negative strand CpG so 5-CG-3/3-GC-5 align in coordinates"
            arg_type = Bool
            default = true
        "--ks"
            help = "K means cluster numbers to evaluate (comma separated ranges, e.g. '2:20,32:2:60')"
            default = "2:20,32:2:60"
        "--n_ils_restarts"
            help = "Number of iterated local search restarts"
            arg_type = Int
            default = 2
        "--n_restarts"
            help = "Number of random restarts"
            arg_type = Int
            default = 20
        "--random_seed_offset"
            help = "Random seed offset"
            arg_type = Int
            default = 0
        "--epochs"
            help = "Maximum epochs for ILS"
            arg_type = Int
            default = 1000
        "--rejection_neighbourhood_fraction"
            help = "Max consectuive rejections in hill climbing = rejection_neighbourhood_fraction*k*readdepth"
            arg_type = Float64
            default = 0.08
        "--max_rejections"
            help = "Maximum consecutive rejections in hillclimbing, only used in rejection_neighbourhood_fraction = 0"
            arg_type = Int
            default = 150
        "--maxiter"
            help = "Maximum iterations for kmeans"
            arg_type = Int
            default = 300
        "--linkage"
            help = "Hierarchical clustering linkage method (single, complete, average, etc.)"
            default = "single"
        "--logiter"
            help = "Log clustering iterations"
            action = :store_true
        "--label"
            help = "Label for output files"
            default = ""
        "files"
            nargs = '+'
            help = "BAM/CRAM files"
            required = true
    end

    return parse_args(s)
end

function parse_ranges(s::String)
    res = Int[]
    for p in split(s, ',')
        p = strip(p)
        isempty(p) && continue
        if occursin(':', p)
            parts = split(p, ':')
            if length(parts) == 2
                append!(res, parse(Int, parts[1]):parse(Int, parts[2]))
            elseif length(parts) == 3
                append!(res, parse(Int, parts[1]):parse(Int, parts[2]):parse(Int, parts[3]))
            end
        else
            push!(res, parse(Int, p))
        end
    end
    return res
end

# function (@main)(args)
function julia_main()::Cint
    banner()
    parsed_args = parse_commandline()
    cmd = parsed_args["%COMMAND%"]

    if cmd == "stats"
        opts = parsed_args["stats"]
        config = opts["config"]
        outputdir = opts["output"]
        bamfile = opts["bamfile"]
        
        runstats(bamfile, config, outputdir=outputdir)

    elseif cmd == "coaccess"
        opts = parsed_args["coaccess"]
        coaccess_wrapper(opts["files"], opts["peaks"], opts["output"];
                         chromlabel=Symbol(opts["chrom"]),
                         startlabel=Symbol(opts["start"]),
                         stoplabel=Symbol(opts["stop"]))

    elseif cmd == "firecount"
        opts = parsed_args["firecount"]
        countfire(opts["files"], opts["peaks"], opts["output"];
            chromlabel=Symbol(opts["chrom"]),
            startlabel=Symbol(opts["start"]),
            stoplabel=Symbol(opts["stop"]))
    
    elseif cmd == "diff"
        opts = parsed_args["diff"]
        
        # Parse method symbol
        method_str = opts["method"]
        method_sym = if method_str == "fisher"
            :fisher
        elseif method_str == "binomial_pooled"
            :binomial_pooled
        elseif method_str == "binomial_ref"
            :binomial_ref
        elseif method_str == "all"
            :all
        else
            error("Unknown method: $method_str. Options: fisher, binomial_pooled, binomial_ref, all")
        end
        
        bamsB_arg = isnothing(opts["groupB"]) ? String[] : string.(opts["groupB"])
        
        if opts["mods"]
            differentialaccessibility_mod_fire(opts["peaks"], string.(opts["groupA"]), bamsB_arg, opts["output"];
                     diff_by_hap=opts["hap"],
                     method=method_sym,
                     chromlabel=Symbol(opts["chrom"]),
                     startlabel=Symbol(opts["start"]),
                     stoplabel=Symbol(opts["stop"]))
        else
            differentialaccessibility(opts["peaks"], string.(opts["groupA"]), bamsB_arg, opts["output"];
                     method=method_sym,
                     chromlabel=Symbol(opts["chrom"]),
                     startlabel=Symbol(opts["start"]),
                     stoplabel=Symbol(opts["stop"]))
        end
    elseif cmd == "cluster"
        opts = parsed_args["cluster"]
        features = Tuple(Symbol.(strip.(split(opts["features"], ','))))
        ks = parse_ranges(opts["ks"])
        linkage_sym = Symbol(opts["linkage"])
        
        cluster_peaks(opts["peaks"], string.(opts["files"]), opts["output"];
            dw=opts["dw"],
            features=features,
            dx=opts["dx"],
            smooth_sigma=opts["smooth_sigma"],
            strictfirefeatures=opts["strictfirefeatures"],
            unify_gc_strand=opts["unify_gc_strand"],
            ks=ks,
            n_ils_restarts=opts["n_ils_restarts"],
            n_restarts=opts["n_restarts"],
            random_seed_offset=opts["random_seed_offset"],
            epochs=opts["epochs"],
            rejection_neighbourhood_fraction=opts["rejection_neighbourhood_fraction"],
            max_rejections=opts["max_rejections"],
            maxiter=opts["maxiter"],
            linkage=linkage_sym,
            logiter=opts["logiter"],
            label=opts["label"]
        )
    end

    return 0
end
