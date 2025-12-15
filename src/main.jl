
function default_config(genome)

    if gethostname() == "penrose"
        return "/penrose/projects/ont/scripts/smgstats/smgstats.config.$(genome).yaml"
    else
        error("Default not defined for $genomee")
    end

end

statreport(bamfile, genome="hg19"; yaml=default_config(genome), nr = -1) = runstats(bamfile, yaml, nr=nr)

function runstats(file, yaml; nr = -1, statdirfilename=false)

    start = time()

    if statdirfilename
        statdir = string("stats_", join(split(basename(file), ".")[1:end-1], "."))
    else
        statdir = "stats"
    end


    #### print banner
    banner()
    smglog("File                   : ", file)

    ##### auto detect file details
    htsdata = autodetecthtsdata(file)
    auxmap  = autodetectaux(file) 
    mods    = autodetectmods(file)
    recorddata = htsdata(auxmap)

    smglog("Detected Data Type     : ", htsdata)
    smglog("Detected Aux map       : ", typeof(auxmap))
    smglog("Detected Modifications : ", mods)

    statset, config = loadstatconfig(yaml, auxmap)
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

    stats = SMGStats.calculatestats(file, statset, recorddata, mods; nr=nr, config=config);
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
