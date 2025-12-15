
function quantile_from_freq(x, counts)

    w = Weights(counts)
    q25 = quantile(x, w, 0.25)
    q50 = quantile(x, w, 0.50)
    q75 = quantile(x, w, 0.75)
    iqr = q75 - q25
    lower = max(minimum(x), q25 - 1.5*iqr)
    upper = min(maximum(x), q75 + 1.5*iqr)
    (q25=q25, median=q50, q75=q75, lower=lower, upper=upper)
end

function quantileboxplot!(qdf; offset=0, width=0.6)

    n = size(qdf, 1)
    palette = Makie.current_default_theme().palette[:color].val
    colors = [palette[((i - 1) % length(palette)) + 1] for i in 1:n]

    for (i, row) in enumerate(eachrow(qdf))
        # Box
        pos = i + offset
        poly!(Rect(pos - width/2 , row.q25, width, row.q75 - row.q25), color=colors[i], strokecolor=:black)
        # Median
        lines!([pos - width/2, pos +  width/2], [row.median, row.median], color=:black, linewidth=2)
        # Whiskers
        lines!([pos, pos], [row.q75, row.upper], color=:black)
        lines!([pos, pos], [row.q25, row.lower], color=:black)
        # lines!([i-0.2, i+0.2], [row.upper, row.upper], color=:black)
        # lines!([i-0.2, i+0.2], [row.lower, row.lower], color=:black)
    end
end


kernelsmooth(x::T, σ, pad=NA()) where {T <: AbstractVector} = imfilter(x, KernelFactors.gaussian(σ, Int(ceil(3σ)) + 1), pad)

function vioboxdist!(df, groupfield, xfield, yfield; violinwidth=0.3, boxwidth=0.2, offset=0.15, probcutt=1e-5, scale=:width)

    n = length(unique(df[!, groupfield]))
    palette = Makie.current_default_theme().palette[:color].val
    colors = [palette[((i - 1) % length(palette)) + 1] for i in 1:n]
    m = 1

    if scale == :area
        m = violinwidth/maximum(df[!, yfield])
    end
    labels = String[]
    maxy = 0
    for (k, gdf) in enumerate(groupby(df, groupfield))

        x = gdf[!, xfield]
        y = gdf[!, yfield]
        w = Weights(y)
        q25 = quantile(x, w, 0.25)
        q50 = quantile(x, w, 0.50)
        q75 = quantile(x, w, 0.75)
        iqr = q75 - q25
        lower = max(minimum(x), q25 - 1.5*iqr)
        upper = min(maximum(x), q75 + 1.5*iqr)

        if scale == :width
            m = violinwidth/maximum(y)
        end

        yv = SMGStats.kernelsmooth(y, 10) 
        # yv = y
        ind = yv .>= probcutt
        # lines!(k .- yv[ind]*m, gdf[ind, xfield], color=:black)
        band!(gdf[ind, xfield], k .- yv[ind]*m, k, direction=:y, color=colors[k])
        push!(labels, string(gdf[1, groupfield]))

        ##### box plot
        pos = k + offset
        poly!(Rect(pos - boxwidth/2 , q25, boxwidth, q75 - q25), color=colors[k], strokecolor=:black)
        # Median
        lines!([pos - boxwidth/2, pos +  boxwidth/2], [q50, q50], color=:black, linewidth=2)
        # Whiskers
        lines!([pos, pos], [q75, upper], color=:black)
        lines!([pos, pos], [q25, lower], color=:black)
    end

    ax = current_axis()
    ax.xticks = (1:length(labels), labels)
    xlims!(ax, 0.5, n + 0.5)
end

edges_to_bins(edges) = edges[1:end-1] .+ (edges.step.lo + edges.step.hi)/2
