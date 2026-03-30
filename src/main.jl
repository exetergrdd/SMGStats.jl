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
            required = true
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
    end

    return parse_args(s)
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
        
        differentialaccessibility(opts["peaks"], string.(opts["groupA"]), string.(opts["groupB"]), opts["output"];
                     method=method_sym,
                     chromlabel=Symbol(opts["chrom"]),
                     startlabel=Symbol(opts["start"]),
                     stoplabel=Symbol(opts["stop"]))
    end

    return 0
end
