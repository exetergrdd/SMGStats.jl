

function htmlreport(sample, data, path::String, file="report.html")

    
    statnames = sort(collect(keys(data)))

    body = IOBuffer()

    println(body, """
    <header class="container sticky-header">
    <h4>Statistics for sample $sample</h4>
    <hr>
    </header>
    <main class="container">
    <aside>
    <nav>
    <h6>Contents</h6>
    <small>
    <ul>
    """)
    ### Summary
    println(body, "<li><a href=\"#summary\">Summary</a></li>")
    for stat in statnames
        println(body, "<li><a href=\"#", stat, "\">", statname(stat) ,"</a></li>")
    end
    println(body, """
    </ul>
    </small>
    </nav>
    </aside>
    <section>
    """)

    ### summary
    println(body, "<article id=\"summary\">")
    println(body, "<h5>Summary:</h5>")
    samplestats, modstats = samplesummary(data)
    dataframe_to_htmltable(body, samplestats)
    if any(==("N50*"), samplestats.Stat)
        println(body, "<small>*N50 calculation currently not robust, mean read length more representative</small>")
    end
    if !isempty(modstats)
        println(body, "<p></p>")
        println(body, "<p></p>")
        dataframe_to_htmltable(body, modstats)
    end
    println(body, "</article>")

    for stat in statnames
        f = plotstat(String(stat), data[stat])
        println(body, "<article id=\"", stat, "\">")
        println(body, "<h5>", statname(stat), ":</h5>")
        show(body, MIME"image/svg+xml"(), f)
        println(body, "</article>")
    end
    println(body, """
    </section>
    </div>
    </main>
    """)
    # 

    filepath  = joinpath(path, file)
    io = open(filepath, "w")
    # <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/water.css@2/out/water.css">
    print(io, """
    <!DOCTYPE html>
    <html>
    <head>
    <meta charset="utf-8">
    <title>$sample</title>
    
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/@picocss/pico@2/css/pico.min.css">
      <style>
    /* Sidebar layout */
    .sticky-header {
        position: sticky;
        top: 0;
        background: var(--pico-background-color); /* match Pico theme */
        z-index: 1000;
        padding: 0.5rem 0;
    }
    main {
      display: flex;
      flex-wrap: wrap;
      gap: 2rem;
    }
    aside {
      flex: 0 0 220px;       /* fixed width */
      position: sticky;      /* sticks in place */
      top: 4.5rem;             /* spacing from top */
      align-self: flex-start;
      height: fit-content;
    }
    section {
      flex: 1 1 auto;        /* fill the rest */
      min-width: 300px;
    }
    /* Optional: tidy up the nav list */
    aside nav ul {
      margin: 0;
      padding: 0;
      list-style: none;
    }
    article[id] {
        scroll-margin-top: 4.5rem;
    }

 
   </style>
    </head>
    <body>
    $(String(take!(body)))
    </body>
    </html>
    """)
    close(io)
    filepath

end


function dataframe_to_htmltable(io, df)
    println(io, "<table class=\"striped\">") 
    println(io, "<thead><tr>")
    for field in names(df)
        println(io, "<th>", field, "</th>")
    end
    println(io, "</tr></thead>")
    println(io, "<tbody>")
    for row in eachrow(df)
        println(io, "<tr>")
        for field in names(df)
            println(io, "<td>", tt(row[field]), "</td>")
        end
        println(io, "</tr>")
    end
    println(io, "</tbody></table>")
end
tt(x::String) = x
tt(x) = string(x)
tt(x::Float64) = string(round(x, digits=2))