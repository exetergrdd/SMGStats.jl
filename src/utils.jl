
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

function quantileboxplot!(qdf)

    n = size(qdf, 1)
    palette = Makie.current_default_theme().palette[:color].val
    colors = [palette[((i - 1) % length(palette)) + 1] for i in 1:n]

    for (i, row) in enumerate(eachrow(qdf))
        # Box
        poly!(Rect(i-0.3, row.q25, 0.6, row.q75 - row.q25), color=colors[i], strokecolor=:black)
        # Median
        lines!([i-0.3, i+0.3], [row.median, row.median], color=:black, linewidth=2)
        # Whiskers
        lines!([i, i], [row.q75, row.upper], color=:black)
        lines!([i, i], [row.q25, row.lower], color=:black)
        lines!([i-0.2, i+0.2], [row.upper, row.upper], color=:black)
        lines!([i-0.2, i+0.2], [row.lower, row.lower], color=:black)
    end
end