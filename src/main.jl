"""
    smglog(msg)

Prints a timestamped, grep-friendly log line:

    [2025-02-12T14:22:03] [SMGStats] <msg>
"""
function smglog(msg...)
    msg = string(msg...)
    ts = Dates.format(Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS")
    println("[$ts] [SMGStats] $msg")
end


function runstats(file, yaml; nr = -1)

    start = time()

    statset, config = loadstatconfig(yaml)


    banner()
    smglog("File       : ", file)
    smglog("Config     : ", yaml)
    smglog("Stats      : ", replace(string(statset), "Tuple{" => "", "}" => ""))
    if nr != -1
        smglog("Sampling   : ", nr, " reads")
    end
    if hasproperty(config, :MetaPlot)
        for mp in config.MetaPlot
            smglog("Metaplot   : ", mp.label, "\t: ", mp.width, " : ", mp.file)
        end
    end
    smglog("Calculating FIRE stats:")
    stats = SMGStats.firestats(file, statset; nr=nr, config=config);
    smglog("Writing Stats:")
    writeallstats(stats, bamfile=file)
    smglog("Generating Report:")
    stats = readstats(bamfile=file)
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
