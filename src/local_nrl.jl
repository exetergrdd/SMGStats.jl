#=
    local_nrl.jl: High-resolution sliding-window local autocorrelation & 2D spectral NRL
    analysis along individual single-molecule chromatin fibers (ModAllFireRead).
=#

struct LocalNRLWorkspace
    raw_1bp::VectorBuffer{Float64}
    binned::VectorBuffer{Float64}
    win_buf::VectorBuffer{Float64}
    weights::VectorBuffer{Float64}
    ac_signal::Vector{Float64}
    detrended::Vector{Float64}
    lags_bp::Vector{Float64}
    periods::Vector{Float64}
    cos_table::Matrix{Float64} # (n_lags, n_periods)
    dx::Int
    max_lag_bins::Int
    detrend_half_w::Int
end

function LocalNRLWorkspace(;
    dx::Int = 5,
    max_lag_bp::Int = 600,
    min_nrl::Real = 140.0,
    max_nrl::Real = 240.0,
    nrl_step::Real = 1.0,
    detrend_bp::Real = 160.0,
    init_bp_capacity::Int = 50_000
)
    max_lag_bins = max(4, fld(max_lag_bp, dx))
    lags_bp = Float64[k * dx for k in 1:max_lag_bins]
    periods = collect(Float64(min_nrl):Float64(nrl_step):Float64(max_nrl))
    n_lags = length(lags_bp)
    n_periods = length(periods)

    cos_table = Matrix{Float64}(undef, n_lags, n_periods)
    @inbounds for j in 1:n_periods
        inv_p = 2.0 * π / periods[j]
        for i in 1:n_lags
            cos_table[i, j] = cos(lags_bp[i] * inv_p)
        end
    end

    detrend_half_w = max(1, round(Int, detrend_bp / dx) ÷ 2)

    return LocalNRLWorkspace(
        VectorBuffer{Float64}(init_bp_capacity),
        VectorBuffer{Float64}(cld(init_bp_capacity, dx)),
        VectorBuffer{Float64}(1024),
        VectorBuffer{Float64}(1024),
        zeros(Float64, n_lags),
        zeros(Float64, n_lags),
        lags_bp,
        periods,
        cos_table,
        dx,
        max_lag_bins,
        detrend_half_w
    )
end

struct LocalNRLReadResult
    read_index::Int
    chrom::String
    lp::Int
    rp::Int
    strand::Bool
    haplotype::UInt8
    label::String
    feature::Symbol
    dx::Int
    window_bp::Int
    step_bp::Int
    centers::Vector{Int}
    periods::Vector{Float64}
    lags_bp::Vector{Float64}
    nrl::Vector{Union{Missing, Int}}
    raw_best_nrl::Vector{Float64}
    ridge_nrl::Vector{Float64}
    phasing_score::Vector{Float64}
    is_closed::BitVector
    peak_type::Vector{Symbol}
    spectrogram::Matrix{Float64}    # (n_periods, n_windows)
    autocor_matrix::Matrix{Float64} # (n_lags, n_windows) detrended local AC
end

struct LocalNRLCohortResult
    loc::UnitRange{Int}
    feature::Symbol
    dx::Int
    window_bp::Int
    step_bp::Int
    centers::Vector{Int}
    periods::Vector{Float64}
    lags_bp::Vector{Float64}
    nrl::Matrix{Union{Missing, Int}}  # (n_windows, n_reads)
    raw_best_nrl::Matrix{Float64}     # (n_windows, n_reads)
    ridge_nrl::Matrix{Float64}        # (n_windows, n_reads)
    phasing_score::Matrix{Float64}    # (n_windows, n_reads)
    is_closed::BitMatrix              # (n_windows, n_reads)
    peak_type::Matrix{Symbol}         # (n_windows, n_reads)
    mean_spectrogram::Matrix{Float64} # (n_periods, n_windows)
    read_results::Vector{LocalNRLReadResult}
end

"""
    extract_read_signal!(raw_buf::VectorBuffer{Float64}, bin_buf::VectorBuffer{Float64}, read::ModAllFireRead, span::UnitRange{Int}; kwargs...)

Extracts a 1D signal for a single `ModAllFireRead` over `span` into `raw_buf` and spatially averages
into `bin_buf` at resolution `dx` without allocating intermediate arrays.
Supports `feature = :fire`, `:mod_6mA`, `:fire_elem`, `:linker`, `:nuc`, or `:mod_5mC`.
"""
function extract_read_signal!(
    raw_buf::VectorBuffer{Float64},
    bin_buf::VectorBuffer{Float64},
    read::ModAllFireRead,
    span::UnitRange{Int};
    feature::Symbol = :mod_6mA,
    dx::Int = 5,
    weight_fire::Float64 = 1.0,
    weight_linker::Float64 = 0.75,
    weight_nuc::Float64 = 0.25,
    weight_mod::Float64 = 1.0
)
    n_bp = length(span)
    start_coord = first(span)
    resize!(raw_buf, n_bp)
    fill!(raw_buf, 0.0)

    if feature == :fire || feature == :fire_elem || feature == :linker
        @inbounds for (s, e, q) in read.msps
            is_fire_elem = q > 0.9 * 255
            if feature == :fire || (feature == :fire_elem && is_fire_elem) || (feature == :linker && !is_fire_elem)
                w_val = is_fire_elem ? weight_fire : weight_linker
                i_start = max(1, s - start_coord + 1)
                i_stop  = min(n_bp, e - start_coord + 1)
                for idx in i_start:i_stop
                    raw_buf[idx] += w_val
                end
            end
        end
    end

    if feature == :fire || feature == :nuc
        @inbounds for (s, e) in read.nucs
            i_start = max(1, s - start_coord + 1)
            i_stop  = min(n_bp, e - start_coord + 1)
            for idx in i_start:i_stop
                raw_buf[idx] += weight_nuc
            end
        end
    end

    if feature == :mod_6mA
        @inbounds for p in read.mods_6mA
            idx = p - start_coord + 1
            if 1 <= idx <= n_bp
                raw_buf[idx] += weight_mod
            end
        end
    elseif feature == :mod_5mC
        @inbounds for p in read.mods_5mC
            p_adj = read.strand ? p : p - 1
            idx = p_adj - start_coord + 1
            if 1 <= idx <= n_bp
                raw_buf[idx] += weight_mod
            end
        end
    elseif feature ∉ (:fire, :fire_elem, :linker, :nuc)
        throw(ArgumentError("Unsupported feature :$feature for local NRL analysis. Use :mod_6mA, :fire, :fire_elem, :linker, :nuc, or :mod_5mC."))
    end

    n_bins = cld(n_bp, dx)
    resize!(bin_buf, n_bins)
    fill!(bin_buf, 0.0)

    if dx == 1
        @inbounds for i in 1:n_bp
            bin_buf[i] = raw_buf[i]
        end
        return n_bins
    end

    @inbounds for b in 1:n_bins
        i1 = (b - 1) * dx + 1
        i2 = min(n_bp, b * dx)
        s = 0.0
        for i in i1:i2
            s += raw_buf[i]
        end
        bin_buf[b] = s / (i2 - i1 + 1)
    end

    return n_bins
end

function _fill_window_weights!(weights::VectorBuffer{Float64}, w_bins::Int, kernel::Symbol)
    resize!(weights, w_bins)
    if kernel == :boxcar || kernel == :rect || w_bins <= 2
        fill!(weights, 1.0)
    elseif kernel == :hann || kernel == :hanning
        # Raised cosine (floor at 0.10 so outer bins still contribute to long-lag pairs)
        denom = Float64(w_bins - 1)
        @inbounds for i in 1:w_bins
            weights[i] = 0.10 + 0.90 * (0.5 * (1.0 - cos(2.0 * π * (i - 1) / denom)))
        end
    elseif kernel == :gaussian
        mid = 0.5 * (w_bins + 1)
        sigma = 0.30 * w_bins
        inv_2s2 = 0.5 / (sigma * sigma)
        @inbounds for i in 1:w_bins
            d = i - mid
            weights[i] = exp(-(d * d) * inv_2s2)
        end
    elseif kernel == :tricube
        mid = 0.5 * (w_bins + 1)
        half = 0.5 * w_bins
        @inbounds for i in 1:w_bins
            u = min(1.0, abs(i - mid) / half)
            t = 1.0 - u * u * u
            weights[i] = 0.10 + 0.90 * (t * t * t)
        end
    else
        throw(ArgumentError("Unknown window kernel :$kernel. Choose from :hann, :boxcar, :gaussian, or :tricube."))
    end
    return weights
end

"""
    _local_weighted_autocor!(ac_signal, win_buf, signal_view, weights, max_lag_bins)

Computes the local weighted autocorrelation of `signal_view` across `1:max_lag_bins` in-place.
Returns `true` if the window has non-zero variance, or `false` if flat/empty.
"""
function _local_weighted_autocor!(
    ac_signal::Vector{Float64},
    win_buf::VectorBuffer{Float64},
    signal_view::SubArray{Float64, 1},
    weights::VectorBuffer{Float64},
    max_lag_bins::Int
)
    w_bins = length(signal_view)
    if w_bins <= max_lag_bins + 2
        fill!(ac_signal, NaN)
        return false
    end

    resize!(win_buf, w_bins)

    # 1. Weighted local mean
    w_sum = 0.0
    wx_sum = 0.0
    @inbounds for i in 1:w_bins
        wi = weights[i]
        xi = signal_view[i]
        w_sum += wi
        wx_sum += wi * xi
    end
    mu = wx_sum / w_sum

    # 2. Mean-center and apply spatial taper
    var_0 = 0.0
    @inbounds for i in 1:w_bins
        z = (signal_view[i] - mu) * sqrt(weights[i])
        win_buf[i] = z
        var_0 += z * z
    end

    if var_0 <= 1e-12
        fill!(ac_signal, NaN)
        return false
    end

    # 3. Unbiased taper-normalized autocorrelation at each lag
    @inbounds for k in 1:max_lag_bins
        n_pairs = w_bins - k
        num = 0.0
        denom_l = 0.0
        denom_r = 0.0
        for i in 1:n_pairs
            zl = win_buf[i]
            zr = win_buf[i + k]
            num += zl * zr
            denom_l += zl * zl
            denom_r += zr * zr
        end
        denom = sqrt(denom_l * denom_r)
        ac_signal[k] = denom > 1e-12 ? (num / denom) : 0.0
    end

    return true
end

function _detrend_ac!(detrended::Vector{Float64}, ac_signal::Vector{Float64}, half_w::Int)
    n = length(ac_signal)
    @inbounds for i in 1:n
        min_idx = max(1, i - half_w)
        max_idx = min(n, i + half_w)
        m = 0.0
        for j in min_idx:max_idx
            m += ac_signal[j]
        end
        detrended[i] = ac_signal[i] - m / (max_idx - min_idx + 1)
    end
    return detrended
end

"""
    _track_nrl_ridge!(ridge_nrl, spectrogram, periods; jump_penalty=0.0015)

Dynamic-programming Viterbi ridge tracker across the 2D local NRL spectrogram `(n_periods, n_windows)`.
Finds a smooth trajectory of NRL regimes along the read that maximizes spectral power while
penalizing abrupt single-step period jumps (`jump_penalty * (L_curr - L_prev)^2`).
"""
function _track_nrl_ridge!(
    ridge_nrl::Vector{Float64},
    spectrogram::Matrix{Float64},
    periods::Vector{Float64};
    jump_penalty::Float64 = 0.0015
)
    n_periods, n_windows = size(spectrogram)
    (n_windows == 0 || n_periods == 0) && return ridge_nrl

    dp_prev = zeros(Float64, n_periods)
    dp_curr = zeros(Float64, n_periods)
    backptr = zeros(Int, n_periods, n_windows)

    @inbounds for p in 1:n_periods
        s = spectrogram[p, 1]
        dp_prev[p] = isfinite(s) ? s : -1e9
        backptr[p, 1] = p
    end

    @inbounds for w in 2:n_windows
        col_has_signal = false
        for p in 1:n_periods
            if isfinite(spectrogram[p, w])
                col_has_signal = true
                break
            end
        end

        for p_curr in 1:n_periods
            s_curr = col_has_signal ? (isfinite(spectrogram[p_curr, w]) ? spectrogram[p_curr, w] : 0.0) : 0.0
            best_val = -Inf
            best_prev = p_curr
            l_curr = periods[p_curr]
            for p_prev in 1:n_periods
                dl = l_curr - periods[p_prev]
                val = dp_prev[p_prev] - jump_penalty * (dl * dl)
                if val > best_val
                    best_val = val
                    best_prev = p_prev
                end
            end
            dp_curr[p_curr] = best_val + s_curr
            backptr[p_curr, w] = best_prev
        end
        copyto!(dp_prev, dp_curr)
    end

    best_idx = argmax(dp_prev)
    @inbounds for w in n_windows:-1:1
        ridge_nrl[w] = isfinite(spectrogram[best_idx, w]) ? periods[best_idx] : NaN
        best_idx = backptr[best_idx, w]
    end

    return ridge_nrl
end

"""
    local_nrl_scan(read::ModAllFireRead; kwargs...) -> LocalNRLReadResult

Computes a sliding-window local autocorrelation and 2D cosine matched-filter NRL spectrogram
along a single `ModAllFireRead`.

# Keyword Arguments
- `loc::Union{Nothing, UnitRange{Int}} = nothing`: Genomic interval to scan (defaults to the full read span `read.lp:read.rp`).
- `feature::Symbol = :mod_6mA`: Feature channel to analyze (`:mod_6mA`, `:fire`, `:fire_elem`, `:linker`, `:nuc`, `:mod_5mC`).
- `dx::Int = 5`: Spatial binning resolution in bp (default `5` bp for high resolution).
- `window_bp::Int = 1000`: Local sliding window width in bp (centered on each evaluated coordinate).
- `step_bp::Int = 25`: Step size in bp between adjacent sliding window centers.
- `max_lag_bp::Int = 600`: Maximum autocorrelation lag in bp (~3 nucleosome repeats).
- `min_nrl::Real = 140.0`, `max_nrl::Real = 240.0`, `nrl_step::Real = 1.0`: Candidate NRL period grid (bp).
- `phase_threshold::Float64 = 0.035`: Minimum phasing power score required to call a window phased (`is_closed = true`).
- `kernel::Symbol = :hann`: Spatial window taper (`:hann`, `:boxcar`, `:gaussian`, `:tricube`).
- `min_coverage_frac::Float64 = 0.85`: Minimum fraction of the window that must overlap the read (`read.lp:read.rp`).
"""
function local_nrl_scan(
    read::ModAllFireRead;
    loc::Union{Nothing, UnitRange{Int}} = nothing,
    read_index::Int = 1,
    feature::Symbol = :mod_6mA,
    dx::Int = 5,
    window_bp::Int = 1000,
    step_bp::Int = 25,
    max_lag_bp::Int = 600,
    min_nrl::Real = 140.0,
    max_nrl::Real = 240.0,
    nrl_step::Real = 1.0,
    phase_threshold::Float64 = 0.035,
    kernel::Symbol = :hann,
    min_coverage_frac::Float64 = 0.85,
    jump_penalty::Float64 = 0.0015,
    workspace::Union{Nothing, LocalNRLWorkspace} = nothing
)
    ws = isnothing(workspace) ? LocalNRLWorkspace(;
        dx=dx, max_lag_bp=max_lag_bp, min_nrl=min_nrl, max_nrl=max_nrl, nrl_step=nrl_step
    ) : workspace

    half_win_bp = window_bp ÷ 2
    target_loc = isnothing(loc) ? ((read.lp + half_win_bp):(read.rp - half_win_bp)) : loc
    if last(target_loc) < first(target_loc)
        target_loc = read.lp:read.rp
    end

    # Extract signal over target_loc padded by half_win_bp on both sides so edge windows have full data
    ext_start = max(read.lp, first(target_loc) - half_win_bp)
    ext_stop  = min(read.rp, last(target_loc) + half_win_bp)
    ext_span  = ext_start:ext_stop

    centers = collect(first(target_loc):step_bp:last(target_loc))
    n_windows = length(centers)
    n_periods = length(ws.periods)
    n_lags    = length(ws.lags_bp)

    nrl_vec          = Vector{Union{Missing, Int}}(missing, n_windows)
    raw_best_nrl_vec = fill(NaN, n_windows)
    ridge_nrl_vec    = fill(NaN, n_windows)
    phasing_vec      = fill(NaN, n_windows)
    closed_vec       = falses(n_windows)
    peak_type_vec    = fill(:none, n_windows)
    spectrogram      = fill(NaN, n_periods, n_windows)
    autocor_matrix   = fill(NaN, n_lags, n_windows)

    if length(ext_span) < max_lag_bp + 4 * dx || n_windows == 0
        return LocalNRLReadResult(
            read_index, read.chrom, read.lp, read.rp, read.strand, read.haplotype, read.label,
            feature, dx, window_bp, step_bp, centers, copy(ws.periods), copy(ws.lags_bp),
            nrl_vec, raw_best_nrl_vec, ridge_nrl_vec, phasing_vec, closed_vec, peak_type_vec,
            spectrogram, autocor_matrix
        )
    end

    n_bins = extract_read_signal!(ws.raw_1bp, ws.binned, read, ext_span; feature=feature, dx=dx)
    binned_data = @view ws.binned.data[1:n_bins]

    w_bins_target = max(ws.max_lag_bins + 4, cld(window_bp, dx))
    _fill_window_weights!(ws.weights, w_bins_target, kernel)

    @inbounds for (w_idx, c_pos) in enumerate(centers)
        win_start_bp = c_pos - half_win_bp
        win_stop_bp  = win_start_bp + window_bp - 1

        # Require sufficient read overlap across the window
        ov_start = max(read.lp, win_start_bp)
        ov_stop  = min(read.rp, win_stop_bp)
        (ov_stop - ov_start + 1 < min_coverage_frac * window_bp) && continue

        b1 = cld(max(1, ov_start - ext_start + 1), dx)
        b2 = min(n_bins, b1 + w_bins_target - 1)
        actual_w_bins = b2 - b1 + 1
        (actual_w_bins <= ws.max_lag_bins + 2) && continue

        if actual_w_bins != length(ws.weights)
            _fill_window_weights!(ws.weights, actual_w_bins, kernel)
        end

        sig_view = @view binned_data[b1:b2]
        valid_ac = _local_weighted_autocor!(ws.ac_signal, ws.win_buf, sig_view, ws.weights, ws.max_lag_bins)
        if actual_w_bins != w_bins_target
            _fill_window_weights!(ws.weights, w_bins_target, kernel)
        end
        valid_ac || continue

        _detrend_ac!(ws.detrended, ws.ac_signal, ws.detrend_half_w)
        for l in 1:n_lags
            autocor_matrix[l, w_idx] = ws.detrended[l]
        end

        # Cosine matched filter across all candidate periods -> spectrogram column
        best_score = -Inf
        best_p_idx = 1
        inv_nlags = 1.0 / n_lags
        for p_idx in 1:n_periods
            s = 0.0
            for l in 1:n_lags
                s += ws.detrended[l] * ws.cos_table[l, p_idx]
            end
            s *= inv_nlags
            spectrogram[p_idx, w_idx] = s
            if s > best_score
                best_score = s
                best_p_idx = p_idx
            end
        end

        best_nrl = ws.periods[best_p_idx]
        raw_best_nrl_vec[w_idx] = best_nrl
        phasing_vec[w_idx] = best_score
        is_closed = best_score > phase_threshold
        closed_vec[w_idx] = is_closed

        if is_closed
            nrl_vec[w_idx] = round(Int, best_nrl)
            idx_1N = searchsortedfirst(ws.lags_bp, best_nrl)
            idx_2N = searchsortedfirst(ws.lags_bp, 2.0 * best_nrl)
            val_1N = (1 <= idx_1N <= n_lags) ? ws.detrended[idx_1N] : 0.0
            val_2N = (1 <= idx_2N <= n_lags) ? ws.detrended[idx_2N] : 0.0
            peak_type_vec[w_idx] = val_2N > val_1N ? :dinucleosome : :mononucleosome
        end
    end

    _track_nrl_ridge!(ridge_nrl_vec, spectrogram, ws.periods; jump_penalty=jump_penalty)

    return LocalNRLReadResult(
        read_index, read.chrom, read.lp, read.rp, read.strand, read.haplotype, read.label,
        feature, dx, window_bp, step_bp, centers, copy(ws.periods), copy(ws.lags_bp),
        nrl_vec, raw_best_nrl_vec, ridge_nrl_vec, phasing_vec, closed_vec, peak_type_vec,
        spectrogram, autocor_matrix
    )
end

"""
    local_nrl_scan(reads::AbstractVector{ModAllFireRead}, loc::Union{Nothing, UnitRange{Int}}=nothing; kwargs...) -> LocalNRLCohortResult

Runs sliding-window local autocorrelation and 2D spectral NRL analysis across a vector of
`ModAllFireRead` molecules from `collectstencilreads`. Reuses a single `VectorBuffer`-backed
`LocalNRLWorkspace` across reads to avoid large intermediate allocations.
"""
function local_nrl_scan(
    reads::AbstractVector{ModAllFireRead},
    loc::Union{Nothing, UnitRange{Int}} = nothing;
    feature::Symbol = :mod_6mA,
    dx::Int = 5,
    window_bp::Int = 1000,
    step_bp::Int = 25,
    max_lag_bp::Int = 600,
    min_nrl::Real = 140.0,
    max_nrl::Real = 240.0,
    nrl_step::Real = 1.0,
    phase_threshold::Float64 = 0.035,
    kernel::Symbol = :hann,
    min_coverage_frac::Float64 = 0.85,
    jump_penalty::Float64 = 0.0015
)
    isempty(reads) && throw(ArgumentError("Cannot run local_nrl_scan on an empty vector of reads."))

    resolved_loc = if !isnothing(loc)
        loc
    else
        min_lp = minimum(r -> r.lp, reads)
        max_rp = maximum(r -> r.rp, reads)
        half_w = window_bp ÷ 2
        (min_lp + half_w < max_rp - half_w) ? ((min_lp + half_w):(max_rp - half_w)) : (min_lp:max_rp)
    end

    ws = LocalNRLWorkspace(;
        dx=dx, max_lag_bp=max_lag_bp, min_nrl=min_nrl, max_nrl=max_nrl, nrl_step=nrl_step,
        init_bp_capacity = length(resolved_loc) + 2 * window_bp + 1024
    )

    centers = collect(first(resolved_loc):step_bp:last(resolved_loc))
    n_windows = length(centers)
    n_reads   = length(reads)
    n_periods = length(ws.periods)

    nrl_mat       = Matrix{Union{Missing, Int}}(missing, n_windows, n_reads)
    raw_best_mat  = fill(NaN, n_windows, n_reads)
    ridge_mat     = fill(NaN, n_windows, n_reads)
    phasing_mat   = fill(NaN, n_windows, n_reads)
    closed_mat    = falses(n_windows, n_reads)
    peak_type_mat = fill(:none, n_windows, n_reads)
    mean_spec     = zeros(Float64, n_periods, n_windows)
    spec_counts   = zeros(Int, n_windows)
    read_results  = Vector{LocalNRLReadResult}(undef, n_reads)

    for (r_idx, r) in enumerate(reads)
        res = local_nrl_scan(
            r;
            loc=resolved_loc,
            read_index=r_idx,
            feature=feature,
            dx=dx,
            window_bp=window_bp,
            step_bp=step_bp,
            max_lag_bp=max_lag_bp,
            min_nrl=min_nrl,
            max_nrl=max_nrl,
            nrl_step=nrl_step,
            phase_threshold=phase_threshold,
            kernel=kernel,
            min_coverage_frac=min_coverage_frac,
            jump_penalty=jump_penalty,
            workspace=ws
        )
        read_results[r_idx] = res

        @inbounds for w in 1:n_windows
            nrl_mat[w, r_idx]       = res.nrl[w]
            raw_best_mat[w, r_idx]  = res.raw_best_nrl[w]
            ridge_mat[w, r_idx]     = res.ridge_nrl[w]
            phasing_mat[w, r_idx]   = res.phasing_score[w]
            closed_mat[w, r_idx]    = res.is_closed[w]
            peak_type_mat[w, r_idx] = res.peak_type[w]

            if isfinite(res.phasing_score[w])
                spec_counts[w] += 1
                for p in 1:n_periods
                    mean_spec[p, w] += res.spectrogram[p, w]
                end
            end
        end
    end

    @inbounds for w in 1:n_windows
        if spec_counts[w] > 0
            inv_c = 1.0 / spec_counts[w]
            for p in 1:n_periods
                mean_spec[p, w] *= inv_c
            end
        else
            for p in 1:n_periods
                mean_spec[p, w] = NaN
            end
        end
    end

    return LocalNRLCohortResult(
        resolved_loc, feature, dx, window_bp, step_bp, centers,
        copy(ws.periods), copy(ws.lags_bp),
        nrl_mat, raw_best_mat, ridge_mat, phasing_mat, closed_mat, peak_type_mat,
        mean_spec, read_results
    )
end


####### VISUALIZATION SUITE FOR LOCAL NRL & 2D SPECTROGRAMS #######

"""
    plot_read_local_nrl(res::LocalNRLReadResult, read::Union{Nothing, ModAllFireRead}=nothing; phase_threshold=0.035, show_autocor=true, f=Figure(size=(1200, show_autocor ? 820 : 620)))
    plot_read_local_nrl(cohort::LocalNRLCohortResult, read_idx::Int, reads::Union{Nothing, AbstractVector{ModAllFireRead}}=nothing; kwargs...)

Visualizes the local sliding-window NRL architecture of a single DNA molecule:
1. **Top Track**: Single-molecule footprint & 6mA calls (nucleosomes in grey, linkers in blue, FIRE elements in red, 6mA calls as ticks).
2. **Middle Panel (2D Local NRL Spectrogram)**: Continuous cosine matched-filter power across candidate NRL periods (`140:240 bp`) along the read, overlaid with the Viterbi ridge trajectory and phased NRL calls.
3. **Phasing Power Track**: Local spectral phasing score along the read vs `phase_threshold`.
4. **Bottom Panel (Optional 2D Local Autocorrelation Map)**: Detrended local autocorrelation `Ã(c, τ)` as a function of genomic position `c` and lag `τ (bp)`.
"""
function plot_read_local_nrl(
    res::LocalNRLReadResult,
    read::Union{Nothing, ModAllFireRead} = nothing;
    phase_threshold::Float64 = 0.035,
    show_autocor::Bool = true,
    show_ridge::Bool = true,
    colormap = :viridis,
    f = Figure(size = (1200, show_autocor ? 850 : 640))
)
    isempty(res.centers) && return f
    x_min = first(res.centers)
    x_max = last(res.centers)
    row_idx = 1

    # 1. Single-Molecule Footprint + 6mA Track (if read struct is provided)
    ax_track = nothing
    if !isnothing(read)
        ax_track = Axis(
            f[row_idx, 1],
            title = "Read #$(res.read_index) ($(res.chrom):$(read.lp)-$(read.rp), strand=$(read.strand ? '+' : '-'), hap=$(read.haplotype)) — Feature: :$(res.feature) (window=$(res.window_bp) bp, step=$(res.step_bp) bp, dx=$(res.dx) bp)",
            ylabel = "Molecule",
            yticks = ([], []),
            xticklabelsvisible = false,
            xgridvisible = false,
            ygridvisible = false
        )
        fireplot(res.chrom, x_min:x_max, [read], [1]; f=f, ax=ax_track, ymax=2, showlegend=false)
        if !isempty(read.mods_6mA)
            m6 = [p for p in read.mods_6mA if x_min <= p <= x_max]
            !isempty(m6) && scatter!(ax_track, m6, fill(1.75, length(m6)), marker=:vline, markersize=8, color=(:navy, 0.75))
        end
        ylims!(ax_track, 2.1, 0.3)
        xlims!(ax_track, x_min, x_max)
        rowsize!(f.layout, row_idx, Fixed(75))
        row_idx += 1
    end

    # 2. 2D Local NRL Spectrogram (Position x Period -> Power)
    spec_title = isnothing(read) ?
        "Read #$(res.read_index) Local NRL Spectrogram (Feature :$(res.feature), window=$(res.window_bp) bp, step=$(res.step_bp) bp, dx=$(res.dx) bp)" :
        "2D Local NRL Spectrogram (Cosine Matched-Filter Power)"
    ax_spec = Axis(
        f[row_idx, 1],
        title = spec_title,
        ylabel = "Candidate NRL (bp)",
        xticklabelsvisible = false,
        xgridvisible = false,
        ygridvisible = false
    )
    hm_spec = heatmap!(ax_spec, res.centers, res.periods, res.spectrogram', colormap=colormap)
    Colorbar(f[row_idx, 2], hm_spec, label = "Spectral Power")

    if show_ridge
        lines!(ax_spec, res.centers, res.ridge_nrl, color=(:white, 0.75), linestyle=:dash, linewidth=1.8, label="Viterbi Ridge")
    end

    phased_x = Float64[]
    phased_y = Float64[]
    for (c, is_cl, best_l) in zip(res.centers, res.is_closed, res.raw_best_nrl)
        if is_cl && isfinite(best_l)
            push!(phased_x, Float64(c))
            push!(phased_y, best_l)
        else
            push!(phased_x, NaN)
            push!(phased_y, NaN)
        end
    end
    lines!(ax_spec, phased_x, phased_y, color=:cyan, linewidth=3.0, label="Phased NRL (> $(phase_threshold))")
    axislegend(ax_spec, position=:rt, framevisible=true, backgroundcolor=(:black, 0.45), labelcolor=:white, labelsize=10)
    xlims!(ax_spec, x_min, x_max)
    ylims!(ax_spec, first(res.periods), last(res.periods))
    row_idx += 1

    # 3. Local Phasing Score Track
    ax_score = Axis(
        f[row_idx, 1],
        ylabel = "Phasing Score",
        xlabel = show_autocor ? "" : "Genomic Position ($(res.chrom))",
        xticklabelsvisible = show_autocor ? false : true,
        xgridvisible = false,
        ygridvisible = false
    )
    hlines!(ax_score, [phase_threshold], color=:crimson, linestyle=:dash, linewidth=1.5, label="Threshold ($phase_threshold)")
    band!(ax_score, res.centers, zeros(length(res.centers)), max.(0.0, replace(res.phasing_score, NaN => 0.0)), color=(:royalblue, 0.25))
    lines!(ax_score, res.centers, res.phasing_score, color=:royalblue, linewidth=2.0, label="Local Phasing Score")
    xlims!(ax_score, x_min, x_max)
    rowsize!(f.layout, row_idx, Fixed(110))
    row_idx += 1

    # 4. 2D Local Autocorrelation Heatmap (Position x Lag)
    ax_ac = nothing
    if show_autocor
        ax_ac = Axis(
            f[row_idx, 1],
            title = "Local Detrended Autocorrelation Ã(x, lag)",
            xlabel = "Genomic Position ($(res.chrom))",
            ylabel = "Lag τ (bp)",
            xgridvisible = false,
            ygridvisible = false
        )
        hm_ac = heatmap!(ax_ac, res.centers, res.lags_bp, res.autocor_matrix', colormap=:balance, colorrange=(-0.25, 0.25))
        Colorbar(f[row_idx, 2], hm_ac, label = "Detrended AC")
        lines!(ax_ac, phased_x, phased_y, color=:black, linewidth=2.0, linestyle=:dot)
        lines!(ax_ac, phased_x, 2.0 .* phased_y, color=(:black, 0.6), linewidth=1.5, linestyle=:dot)
        xlims!(ax_ac, x_min, x_max)
        ylims!(ax_ac, first(res.lags_bp), last(res.lags_bp))
    end

    axes_to_link = Axis[ax_spec, ax_score]
    !isnothing(ax_track) && pushfirst!(axes_to_link, ax_track)
    !isnothing(ax_ac) && push!(axes_to_link, ax_ac)
    linkxaxes!(axes_to_link...)

    return f
end

function plot_read_local_nrl(
    cohort::LocalNRLCohortResult,
    read_idx::Int,
    reads::Union{Nothing, AbstractVector{ModAllFireRead}} = nothing;
    kwargs...
)
    rd = (!isnothing(reads) && 1 <= read_idx <= length(reads)) ? reads[read_idx] : nothing
    return plot_read_local_nrl(cohort.read_results[read_idx], rd; kwargs...)
end

"""
    plot_local_nrl_overview(cohort::LocalNRLCohortResult; order=1:length(cohort.read_results), f=Figure(size=(1250, 820)))

Visualizes the local NRL landscape across a cohort of reads spanning `cohort.loc`:
1. **Top Panel**: Population-mean 2D NRL spectrogram across all reads, plus the fraction of phased reads at each genomic position.
2. **Middle Panel**: Read-by-position heatmap of local phased NRLs (`nrl` in bp; unphased segments shown in light grey).
3. **Bottom Panel**: Read-by-position heatmap of local spectral phasing scores.
"""
function plot_local_nrl_overview(
    cohort::LocalNRLCohortResult;
    order::AbstractVector{<:Integer} = 1:length(cohort.read_results),
    nrl_colorrange::Tuple{Real, Real} = (155, 215),
    f = Figure(size = (1250, 840))
)
    n_windows, n_reads = size(cohort.phasing_score)
    (n_windows == 0 || n_reads == 0) && return f

    x_min = first(cohort.centers)
    x_max = last(cohort.centers)

    # 1. Population Mean Spectrogram
    ax_spec = Axis(
        f[1, 1],
        title = "Population Mean Local NRL Spectrogram (Feature :$(cohort.feature), window=$(cohort.window_bp) bp, step=$(cohort.step_bp) bp, dx=$(cohort.dx) bp)",
        ylabel = "Candidate NRL (bp)",
        xticklabelsvisible = false,
        xgridvisible = false,
        ygridvisible = false
    )
    hm_spec = heatmap!(ax_spec, cohort.centers, cohort.periods, cohort.mean_spectrogram', colormap=:viridis)
    Colorbar(f[1, 2], hm_spec, label = "Mean Power")
    xlims!(ax_spec, x_min, x_max)
    ylims!(ax_spec, first(cohort.periods), last(cohort.periods))
    rowsize!(f.layout, 1, Relative(0.28))

    # 2. Single-Molecule Local NRL Heatmap (Phased Windows)
    nrl_float = Matrix{Float64}(undef, n_windows, length(order))
    score_ord = Matrix{Float64}(undef, n_windows, length(order))
    @inbounds for (new_r, old_r) in enumerate(order)
        for w in 1:n_windows
            v = cohort.nrl[w, old_r]
            nrl_float[w, new_r] = ismissing(v) ? NaN : Float64(v)
            score_ord[w, new_r] = cohort.phasing_score[w, old_r]
        end
    end

    ax_nrl = Axis(
        f[2, 1],
        title = "Single-Molecule Local Phased NRL Regimes (bp, unphased = grey)",
        ylabel = "Reads (n=$(length(order)))",
        xticklabelsvisible = false,
        xgridvisible = false,
        ygridvisible = false,
        backgroundcolor = :grey92
    )
    hm_nrl = heatmap!(
        ax_nrl, cohort.centers, 1:length(order), nrl_float,
        colormap = :turbo, colorrange = nrl_colorrange, nan_color = (:grey90, 0.5)
    )
    Colorbar(f[2, 2], hm_nrl, label = "Local NRL (bp)")
    ylims!(ax_nrl, length(order) + 0.5, 0.5)
    xlims!(ax_nrl, x_min, x_max)

    # 3. Single-Molecule Local Phasing Score Heatmap
    ax_score = Axis(
        f[3, 1],
        title = "Single-Molecule Local Phasing Score",
        xlabel = "Genomic Coordinate",
        ylabel = "Reads (n=$(length(order)))",
        xgridvisible = false,
        ygridvisible = false
    )
    hm_score = heatmap!(
        ax_score, cohort.centers, 1:length(order), score_ord,
        colormap = :magma, colorrange = (0.0, 0.12), nan_color = :grey90
    )
    Colorbar(f[3, 2], hm_score, label = "Phasing Score")
    ylims!(ax_score, length(order) + 0.5, 0.5)
    xlims!(ax_score, x_min, x_max)

    linkxaxes!(ax_spec, ax_nrl, ax_score)
    return f
end

