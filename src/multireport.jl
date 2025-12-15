

function combine_iframe_viewer(files::Vector{String}, labels::Vector{String};
                               outfile="viewer.html",
                               title="Nanopore Sample Viewer")

    @assert length(files) == length(labels)

    # Load the template (you can also embed the string directly instead)
    template = """
    <!DOCTYPE html>
    <html>
    <head>
      <meta charset="utf-8">
      <title>$title</title>

      <link rel="stylesheet"
            href="https://cdn.jsdelivr.net/npm/@picocss/pico@2/css/pico.min.css">

      <style>
        body {
          margin: 0;
          display: flex;
          height: 100vh;
          overflow: hidden;
        }
        #sidebar {
          width: 260px;
          border-right: 1px solid #ddd;
          padding: 1rem;
          overflow-y: auto;
          position: sticky;
          top: 0;
        }
        #content {
          flex: 1;
          padding: 0;
          overflow: hidden;
        }
        iframe {
          width: 100%;
          height: 100%;
          border: none;
        }
        #sample-list a {
          display: block;
          padding: 0.25rem 0.5rem;
          text-decoration: none;
        }
        #sample-list a:hover {
          background: var(--pico-muted-border-color);
        }
      </style>

      <script>
        function loadSample(path) {
          document.getElementById("viewer").src = path;
        }
      </script>

    </head>

    <body>
      <div id="sidebar">
        <h4>Samples</h4>
        <aside>
        <nav>
          <ul id="sample-list">
            <!-- LINKS -->
          </ul>
          
        </nav>
        </aside>
      </div>

      <div id="content">
        <iframe id="viewer" src=""></iframe>
      </div>

    </body>
    </html>
    """

    buf = IOBuffer()

    println(buf, template)

    html = String(take!(buf))

    # Build the clickable sidebar list
    links = IOBuffer()
    for (file, label) in zip(files, labels)
        println(links,
            """<li><a href="#" onclick="loadSample('$file')">$label</a></li>"""
        )
    end

    # Insert into placeholder <!-- LINKS -->
    html = replace(html, "<!-- LINKS -->" => String(take!(links)))

    # Write final output
    open(outfile, "w") do io
        write(io, html)
    end

    println("Wrote combined viewer → $outfile")
end
