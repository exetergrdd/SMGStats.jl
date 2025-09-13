
function plotasstring(p)
    iob = IOBuffer()
    show(IOContext(iob, :color => true), MIME("text/plain"), p)
    String(take!(iob))
end


strip_ansi(s) = replace(s, r"\x1B\[[0-9;]*m" => "")
function sidebyside(plots...; col="   │  ")
    # Convert all plots to strings
    strs = [plotasstring(p) for p in plots]
    # Split into lines
    lines = [split(s, '\n') for s in strs]
    # Compute max height
    maxh = maximum(length.(lines))
    # Pad each plot with empty lines
    for l in lines
        while length(l) < maxh
            push!(l, "")
        end
    end
    # Compute display width of each column (strip ANSI first)
    
    widths = [maximum(length.(strip_ansi.(l))) for l in lines]
    # Build combined lines
    buf = IOBuffer()
    for i in 1:maxh
        for (j, l) in enumerate(lines)
            line = l[i]
            vislen = length(strip_ansi(line))
            print(buf, line)
            print(buf, " "^(widths[j] - vislen))
            if j < length(lines)
                print(buf, col)
            end
        end
        println(buf)
    end
    println(String(take!(buf)))
end


function displayplots(ps)

    for p in ps
        if length(p) == 1
            display(p[1])
        else
            sidebyside(p...)
        end
    end
end