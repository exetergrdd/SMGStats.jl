


####### FEATURE MATRIX GENERATION ########

using ImageFiltering

@inline inc_pos!(M, i, col, weight, region_length) = (1 <= i <= region_length) && @inbounds M[i, col] += weight

@inline function inc_positions!(M, positions, col, start_coord, region_length; weight=1.0, incfun=inc_pos!, strand=true, decrement_neg=false)
    for p in positions
        if !strand && decrement_neg
            p -= 1
        end
        incfun(M, p - start_coord + 1, col, weight, region_length)
    end
    M
end

function mean_spatial_bin(H::AbstractMatrix{T}, dx::Int) where T
    dx == 1 && return H
    n, m = size(H)
    w = cld(n, dx)
    AH = zeros(T, w, m)
    te = zeros(Int, w)
    
    @inbounds for i = 1:m
        for j = 1:n
            wj = cld(j, dx)
            (i == 1) && (te[wj] += 1)
            AH[wj, i] += H[j, i]
        end
    end

    # Stabilize tiny remainder bins at the right edge to prevent noise
    if te[end] < 0.75 * dx && w > 1
        @inbounds for i = 1:m
            AH[end, i] += AH[end-1, i]
        end
        te[end] += te[end-1]
    end

    # Calculate true mean per bin
    @inbounds for i = 1:m
        for j = 1:w
            AH[j, i] /= te[j]
        end
    end
    
    return AH
end


"""
    featurematrix(loc, reads::Vector{ModAllFireRead}; kwargs...)

Construct a feature matrix for a specific genomic location (`loc`) across a set of `reads`.
Returns a named tuple containing the combined feature matrix `M` constructed at resolution `dx`, 
as well as the individual matrices for FIRE elements, linker, nucleosomes, and modifications.
"""
function featurematrix(loc, reads::Vector{ModAllFireRead}; features=(:fire, :mod_6mA, :mod_5mC, :mod_5hmC), dx = 15, 
                       smooth_sigma = 0.0, 
                       mod_inc_fun=inc_pos!,
                       weight_fire=1.0, weight_linker=0.75, weight_nuc=0.25, weight_mod=1.0, unify_gc_strand=true)
    
    n = length(loc)
    m = length(reads)
    ndx = cld(n, dx)
    loc_first = first(loc)
    levels = 1:m
    
    do_fire       = :fire       ∈ features
    do_fire_elem  = :fire_elem  ∈ features
    do_linker     = :linker     ∈ features
    do_nuc        = :nuc        ∈ features
    do_nuc_dyad   = :nuc_dyad   ∈ features
    do_mod_6mA    = :mod_6mA    ∈ features
    do_mod_5mC    = :mod_5mC    ∈ features
    do_mod_5hmC   = :mod_5hmC   ∈ features
    do_unmod_6mA  = :unmod_6mA  ∈ features
    do_unmod_5mC  = :unmod_5mC  ∈ features
    do_unmod_5hmC = :unmod_5hmC ∈ features
    
    # --- 1. Matrix Allocation ---
    # Everything is built natively at 1bp resolution
    mat_fire      = do_fire      ? zeros(Float64, n, m) : nothing
    mat_fire_elem = do_fire_elem ? zeros(Float64, n, m) : nothing
    mat_linker    = do_linker    ? zeros(Float64, n, m) : nothing
    mat_nuc       = do_nuc       ? zeros(Float64, n, m) : nothing
    mat_nuc_dyad  = do_nuc_dyad  ? zeros(Float64, n, m) : nothing

    mat_6mA       = do_mod_6mA    ? zeros(Float64, n, m) : nothing
    mat_5mC       = do_mod_5mC    ? zeros(Float64, n, m) : nothing
    mat_5hmC      = do_mod_5hmC   ? zeros(Float64, n, m) : nothing
    umat_6mA      = do_unmod_6mA  ? zeros(Float64, n, m) : nothing
    umat_5mC      = do_unmod_5mC  ? zeros(Float64, n, m) : nothing
    umat_5hmC     = do_unmod_5hmC ? zeros(Float64, n, m) : nothing

    # --- 2. Population Loop (1bp Resolution) ---
    for (r, l) in zip(reads, levels)
        if do_fire || do_fire_elem || do_linker || do_nuc
            if do_fire || do_fire_elem || do_linker
                for (s, e, q) in r.msps
                    is_fire_elem = q > 0.9 * 255
                    weight_elem  = is_fire_elem ? weight_fire : weight_linker
                    iv = s:e
                    do_fire      &&                inc_positions!(mat_fire,      iv, l, loc_first, n, weight=weight_elem)
                    do_fire_elem &&  is_fire_elem && inc_positions!(mat_fire_elem, iv, l, loc_first, n, weight=weight_elem)
                    do_linker    && !is_fire_elem && inc_positions!(mat_linker,    iv, l, loc_first, n, weight=weight_elem)
                end
            end
            
            if do_fire || do_nuc
                for (s, e) in r.nucs
                    iv = s:e
                    do_fire && inc_positions!(mat_fire, iv, l, loc_first, n, weight=weight_nuc)
                    do_nuc  && inc_positions!(mat_nuc,  iv, l, loc_first, n, weight=weight_nuc)
                end
            end
        end
        
        do_mod_6mA    && inc_positions!(mat_6mA,   r.mods_6mA,     l, loc_first, n, weight=weight_mod, incfun=mod_inc_fun)
        do_unmod_6mA  && inc_positions!(umat_6mA,  r.un_mods_6mA,  l, loc_first, n, weight=weight_mod, incfun=mod_inc_fun)
        
        do_mod_5mC    && inc_positions!(mat_5mC,   r.mods_5mC,     l, loc_first, n, weight=weight_mod, incfun=mod_inc_fun, strand=r.strand, decrement_neg=unify_gc_strand)
        do_unmod_5mC  && inc_positions!(umat_5mC,  r.un_mods_5mC,  l, loc_first, n, weight=weight_mod, incfun=mod_inc_fun, strand=r.strand, decrement_neg=unify_gc_strand)
        do_mod_5hmC   && inc_positions!(mat_5hmC,  r.mods_5hmC,    l, loc_first, n, weight=weight_mod, incfun=mod_inc_fun, strand=r.strand, decrement_neg=unify_gc_strand)
        do_unmod_5hmC && inc_positions!(umat_5hmC, r.un_mods_5hmC, l, loc_first, n, weight=weight_mod, incfun=mod_inc_fun, strand=r.strand, decrement_neg=unify_gc_strand)
    end
    
    # --- 3. Smoothing and Binning Phase ---
    ker = Kernel.gaussian((smooth_sigma, 0.0))
    mat_fire      = smooth_matrix(mat_fire, smooth_sigma, dx, ker)
    mat_fire_elem = smooth_matrix(mat_fire_elem, smooth_sigma, dx, ker)
    mat_linker    = smooth_matrix(mat_linker, smooth_sigma, dx, ker)
    mat_nuc       = smooth_matrix(mat_nuc, smooth_sigma, dx, ker)
    mat_nuc_dyad  = smooth_matrix(mat_nuc_dyad, smooth_sigma, dx, ker)
    mat_6mA       = smooth_matrix(mat_6mA, smooth_sigma, dx, ker)
    mat_5mC       = smooth_matrix(mat_5mC, smooth_sigma, dx, ker)
    mat_5hmC      = smooth_matrix(mat_5hmC, smooth_sigma, dx, ker)
    umat_6mA      = smooth_matrix(umat_6mA, smooth_sigma, dx, ker)
    umat_5mC      = smooth_matrix(umat_5mC, smooth_sigma, dx, ker)
    umat_5hmC     = smooth_matrix(umat_5hmC, smooth_sigma, dx, ker)
    
    # --- 4. Final Aggregation ---
    M = zeros(Float64, ndx*length(features), m)   
    
    for (k, feat) in enumerate(features)
        ind = (1:ndx) .+ ndx*(k-1)
        if     feat == :fire;       M[ind, :] .= mat_fire
        elseif feat == :fire_elem;  M[ind, :] .= mat_fire_elem
        elseif feat == :linker;     M[ind, :] .= mat_linker
        elseif feat == :nuc;        M[ind, :] .= mat_nuc
        elseif feat == :mod_6mA;    M[ind, :] .= mat_6mA
        elseif feat == :mod_5mC;    M[ind, :] .= mat_5mC
        elseif feat == :mod_5hmC;   M[ind, :] .= mat_5hmC
        elseif feat == :unmod_6mA;  M[ind, :] .= umat_6mA
        elseif feat == :unmod_5mC;  M[ind, :] .= umat_5mC
        elseif feat == :unmod_5hmC; M[ind, :] .= umat_5hmC
        else
            @warn "Feature :$feat is not recognized and was skipped."
        end
    end
    
    # Notice: The global M ./= dx has been completely removed.
    
    return (; M, mat_fire, mat_fire_elem, mat_linker, mat_nuc, 
              mat_6mA, mat_5mC, mat_5hmC, umat_6mA, umat_5mC, umat_5hmC, 
              loc, dx, ldx=loc_first:dx:last(loc), blockind=1:ndx, features=features)
end

function smooth_matrix(M_high, smooth_sigma, dx, ker = Kernel.gaussian((smooth_sigma, 0.0)))
    isnothing(M_high) && return nothing
    if smooth_sigma > 0.0
        M_high = imfilter(M_high, ker, "replicate")
    end
    return mean_spatial_bin(M_high, dx)
end
    

function featureblock(FM, entry::Symbol; M = FM.M)
    i = 0
    for (k, f) in enumerate(FM.features)
        if entry == f
            i = k
            break
        end
    end
    (i == 0) && error("$entry not found in $(FM.features)")
    
    ind = FM.blockind .+ (i - 1)*length(FM.blockind)
    view(M, ind, :)
end

############# ASSESS READ DEPTH ##############################################


function collectreadextents(chrom, loc, reader::ModGroupFileReader ; wex=0)
    
    starts = Int[]
    stops = Int[]
    first_feature = Int[]
    last_feature = Int[]
    labels = String[]
    
    wloc = (first(loc) - wex):(last(loc) + wex)

    recorddata = reader.recorddata

    for (reader, label) in zip(reader.readers, reader.labels)
        for (k, record) in enumerate(eachintersection(reader, chrom, wloc))
            validflag(record) || continue
            processread!(record, recorddata) || continue
            
            lp = SMGReader.leftposition(record) + 1
            rp = SMGReader.rightposition(recorddata)
            push!(starts, lp)
            push!(stops, rp)
            
            firegenomecoords = firegenome(record, recorddata)
            nucs = collect(firegenomecoords(firenucs(record, recorddata)))
            msps = collect(firegenomecoords(firemsps(record, recorddata)))
            
            ff = typemax(Int)
            lf = 0
            for (s, e) in nucs
                ff = min(s, ff)
                lf = max(e, lf)
            end
            for (s, e, q) in msps
                ff = min(s, ff)
                lf = max(e, lf)
            end
            push!(first_feature, ff)
            push!(last_feature, lf)

            push!(labels, label)
        end
    
    end
    (; starts, stops, first_feature, last_feature, labels)
end


function read_depth_table(chrom, loc, reader ; window_sizes = 0:50:5000)
    mid_pos = div(first(loc) + last(loc), 2)
    readextents = DataFrame(collectreadextents(chrom, loc, reader))
    ### determine the number of reads fully in which the fire elements fully span the interval

    depths = mapreduce(gdf -> DataFrame(W=window_sizes, Depth=[sum((gdf.first_feature .< (mid_pos - w)) .& (gdf.last_feature .> (mid_pos + w))) for w in window_sizes], Label=gdf.labels[1]), vcat, groupby(readextents, :labels))

    depths
end



################# BAYESIAN INFORMATION CRITERION FOR KMEANS ########################

## Strategy for clustering use kmeans implied gaussian mixture model to calculate BIC to identify cluster number
## its easy to get a decent clustering, but hard to get an optimal BIC as single outlying molecules that should
## group into small clusters are hard to place, therefore we have a complex hill climbing strategy to find the 
## optimal clustering
# 1. Initial clustering hierachical with single linkage, this sets a good starting point
# 2. Iterative local search hill climbing + kmeans lloyd algorithm iterations for cost function
# 3. Strict hillclimbing, with periodic allowing non-optimal moves to jump out of basin, 
# 4. Multiple random restarts of the whole problem.
## This appears to robustly ensure optimal clustering is found, but its computational expensive



#### Data structures 

struct KMeansData
    centers::Matrix{Float64}
    distances::Matrix{Float64}
    assignments::Vector{Int}
    counts::Vector{Int}
    to_update::Vector{Bool}
end

function KMeansData(n_features::Int, n_samples::Int, k::Int)
    return KMeansData(
        zeros(Float64, n_features, k), 
        zeros(Float64, k, n_samples),  
        ones(Int, n_samples), 
        zeros(Int, k), 
        zeros(Bool, k) 
    )
end

struct KMeansOptResult
    centers::Matrix{Float64}
    assignments::Vector{Int}
    costs::Vector{Float64}       # Per-sample costs
    counts::Vector{Int}
    totalcost::Float64
    probabilities::Matrix{Float64}
    h_reads::Vector{Float64}
end

function KMeansOptResult(n_features::Int, n_samples::Int, k::Int)
    return KMeansOptResult(
        zeros(Float64, n_features, k), 
        ones(Int, n_samples), 
        zeros(Float64, n_samples),
        zeros(Int, k), 
        Inf,
        zeros(Float64, k, n_samples),
        zeros(Float64, n_samples)
    ) 
end

Base.copy(kmd::KMeansData) = KMeansData(copy(kmd.centers), copy(kmd.distances), copy(kmd.assignments), copy(kmd.counts), copy(kmd.to_update))

function Base.copyto!(dest::KMeansData, src::KMeansData)
    dest.centers .= src.centers
    dest.distances .= src.distances
    dest.assignments .= src.assignments
    dest.counts .= src.counts
    dest.to_update .= src.to_update
    return dest
end

@inline function sample_two(rng::AbstractRNG, data, num_samples::Int)
    a = rand(rng, data)
    b = a
    if num_samples > 1
        while a == b
            b = rand(rng, data)
        end
    end
    return a, b
end

@inline function update_cluster_center!(w::KMeansData, M::Matrix{Float64}, rng::AbstractRNG; logiter = false, logtrace=[])
    n_samples = size(M, 2)
    k = length(w.counts)
    num_affected = 0
    
    # 1. Zero out affected centers and counts before accumulation
    @inbounds for c in 1:k
        if w.to_update[c] || w.counts[c] == 0
            w.counts[c] = 0
            @views w.centers[:, c] .= 0.0
        end
    end

    # 2. Accumulate new coordinates
    @inbounds for j in 1:n_samples
        c = w.assignments[j]
        if w.to_update[c]
            w.counts[c] += 1
            @views w.centers[:, c] .+= M[:, j]
        end
    end

    # 3. Average or respawn
    
    @inbounds for c in 1:k
        if w.counts[c] == 0
            j = rand(rng, 1:n_samples)
            if logiter
               push!(logtrace, (:randomspawn, c, sum(w.centers[:, c]), sum(M[:, j]), 0.0, j, 0, 0, 0, 0, 0))
            end
            @views w.centers[:, c] .= M[:, j]
            w.to_update[c] = true
            num_affected += 1
        elseif w.to_update[c]
            @views w.centers[:, c] ./= w.counts[c]
            num_affected += 1
        end
    end
    return num_affected
end

function converge_kmeans!(M::Matrix{Float64}, w::KMeansData, rng::AbstractRNG; logtrace=[], logiter=false, maxiter=300)
    n_samples = size(M, 2)
    k = size(w.centers, 2)
    
    total_cost = Inf
    converged = false
    num_affected = sum(w.to_update)
    
        if logiter
           push!(logtrace, (:kmeansstart, 0, 0.0, sum(w.distances), sum(w.centers), sum(w.to_update), 0, 0, 0, 0, 0))
        end
    
    for iter in 1:maxiter
        old_cost = total_cost
        total_cost = 0.0
        changed = false
        
        # 1. Calculate distances on things that need to be updated
        if (num_affected > 0.75 * k )
            pairwise!(w.distances, SqEuclidean(), w.centers, M, dims=2)
        else
            dists_view = view(w.distances, w.to_update, :)
            centers_view = view(w.centers, :, w.to_update)
            pairwise!(dists_view, SqEuclidean(), centers_view, M, dims=2)
        end
        
        # Reset update flags
        fill!(w.to_update, false)
        
        # 2. Assignments
        @inbounds for j in 1:n_samples
            best_dist = Inf
            best_cluster = 1
            for c in 1:k
                dist = w.distances[c, j]
                if dist < best_dist
                    best_dist = dist
                    best_cluster = c
                end
            end
            
            total_cost += best_dist
            prev_cluster = w.assignments[j]
            
            # If a data point switches cluster, update the source and destination
            if prev_cluster != best_cluster
                w.assignments[j] = best_cluster
                if 1 <= prev_cluster <= k
                    w.to_update[prev_cluster] = true
                end
                w.to_update[best_cluster] = true
                changed = true
            end
        end
        if logiter
           push!(logtrace, (:kmeans, iter, old_cost, total_cost, sum(w.centers), sum(w.to_update), 0, 0, 0, 0, 0))
        end
        
        if !changed || (old_cost - total_cost < 1e-7)
            converged = true
            break
        end
        
        # 3. Center Updates 
        num_affected = update_cluster_center!(w, M, rng, logiter=logiter, logtrace=logtrace)
    end
    
    return total_cost, converged
end

function hillclimb_kmeans!(M::Matrix{Float64}, cost::Float64, worker::KMeansData, currentbest::KMeansData, rng::AbstractRNG; max_rejections=150, logtrace=[], logiter=false, epochs=1000, maxiter=300, mutstep=:two_sample, restart_seed=1)
    k = size(worker.centers, 2)
    clusters = 1:k
    n_samples = size(M, 2)
    samples = 1:n_samples
    failures = 0
    
    currentbest_cost = cost
    resetstate = false
    
    consecutive_rejections = 0
    
    for epoch in 1:epochs
        # 1. State Restoration (triggers if previous mutation failed to find a lower cost)
        if resetstate
            copyto!(worker, currentbest) 
        end
         # copyto!(worker, currentbest) 
        # 2. Mutation
        fill!(worker.to_update, false)
#         n_mutations = min((rand(rng) < 0.80 ? 1 : 2), k)
        
#         ca, cb = sample_two(rng, clusters, n_mutations)
#         da, db = sample_two(rng, samples, n_mutations)
    
#         @views worker.centers[:, ca] .= M[:, da]
#         worker.to_update[ca] = true
#         if n_mutations > 1
#             @views worker.centers[:, cb] .= M[:, db]
#             worker.to_update[cb] = true
#         end
        

        n_mutations = min((rand(rng) < 0.8 ? 1 : 2), k) 
        centers_to_replace = sample(rng, 1:k, n_mutations, replace=false)
        new_points = sample(rng, 1:n_samples, n_mutations, replace=false)
                    
        for (c, p) in zip(centers_to_replace, new_points)
            @views worker.centers[:, c] .= M[:, p]
            worker.to_update[c] = true
        end
        if n_mutations == 1
            ca = centers_to_replace[1]
            cb = ca
            da = new_points[1]
            db = da
        else
            ca, cb = centers_to_replace
            da, db = new_points
        end
        
        # 3. Convergence on new clusters
        new_cost, converged = converge_kmeans!(M, worker, rng, maxiter=maxiter, logtrace=logtrace, logiter=logiter)
        
        
        if !converged
            failures += 1
        end
        
        if logiter
           push!(logtrace, (:hillclimb, epoch,currentbest_cost, new_cost, sum(worker.centers), restart_seed, n_mutations, ca, cb, da, db))
        end
        # 4. Take better costs
        if new_cost < currentbest_cost
            copyto!(currentbest, worker)
            currentbest_cost = new_cost
            consecutive_rejections = 0
        else
            consecutive_rejections += 1
        end
         
        
        resetstate = new_cost >= currentbest_cost
        (consecutive_rejections .> max_rejections) && break
    end
    
    return currentbest_cost, failures
end



function extract_single_molecule_metrics(D::AbstractMatrix{Float64}, totalcost::Float64)
    k, n_samples = size(D) # Assuming D is the K x N distance matrix
    
    
    T = (totalcost / n_samples) + 1e-12
    
    # Pre-allocate outputs for the optimal model
    probabilities = zeros(Float64, k, n_samples)
    h_reads = zeros(Float64, n_samples)
    
    @inbounds for j in 1:n_samples
        dj = view(D, :, j)
        min_d = minimum(dj)
        
        Z = 0.0
        for i in 1:k
            Z += exp(-(dj[i] - min_d) / T)
        end
        
        H_read = 0.0
        for i in 1:k
            p = exp(-(dj[i] - min_d) / T) / Z
            probabilities[i, j] = p
            
            if p > 1e-12
                H_read -= p * log(p)
            end
        end
        h_reads[j] = H_read
    end
    
    return (; probabilities, h_reads)
end

"""
    select_k_kmeans_bic(M::Matrix{Float64}, ks::AbstractVector{Int}; kwargs...)

Identifies the optimal number of clusters `K` from the candidate vector `ks` using 
Bayesian Information Criterion (BIC). It explores clustering assignments with Iterated 
Local Search (ILS) and multiple restarts to avoid local minima.
"""
function select_k_kmeans_bic(M::Matrix{Float64}, ks::AbstractVector{Int}; n_ils_restarts::Int=1, n_restarts::Int=20, random_seed_offset=0, epochs::Int=1000, rejection_neighbourhood_fraction=0.0, max_rejections=150, maxiter::Int=300, linkage=:single, rng=Xoshiro(), logiter=false)
    n_features, n_samples = size(M)
    

    best_bics = zeros(Float64, length(ks))
    
    best_bic = Inf
    best_model = KMeansOptResult(0, 0, 0)
    
    # best_models = Vector{KmeansResult}(undef, length(ks))
    best_seeds = zeros(Int, length(ks)) 
    convergence_failures = zeros(Int, length(ks))
    mean_assignment_entropy = zeros(Float64, length(ks))
    max_read_assignment_entropy = zeros(Float64, n_samples)
    max_assignment_entropy = 0.0
    
    D_dist = pairwise(SqEuclidean(), M, dims=2)
    hc = hclust(D_dist, linkage=linkage)
    
    logtrace = []
    if logiter
        logtrace = Vector{Tuple{Symbol, Int, Float64, Float64, Float64, Int, Int, Int, Int, Int, Int}}(undef, 0)
    end
    
    for (i, k) in enumerate(ks)
        if rejection_neighbourhood_fraction > 0
            max_rejections = max(100, Int(round(n_samples*k*rejection_neighbourhood_fraction)))
        end

        init_clustering = KMeansData(n_features, n_samples, k)
        init_clustering.assignments .= cutree(hc, k=k)
        fill!(init_clustering.to_update, true)
        
        update_cluster_center!(init_clustering, M, rng)
        initial_cost, initial_converged = converge_kmeans!(M, init_clustering, rng, maxiter=maxiter, logiter=logiter, logtrace=logtrace)
        
        currentbest = copy(init_clustering)
        bestoverall = copy(init_clustering)
        worker      = copy(init_clustering)
        globalbest  = copy(init_clustering)
        
        globalbest_cost = initial_cost
        bestoverall_cost = initial_cost
        best_restart = 0
        total_failures = initial_converged ? 0 : 1
        
        
        for seed_offset in 1:n_ils_restarts
            copyto!(worker, init_clustering)
            copyto!(bestoverall, init_clustering)
            bestoverall_cost = initial_cost
            for restart_seed in 1:n_restarts
                current_seed = restart_seed + (seed_offset - 1)*n_restarts + random_seed_offset
                Random.seed!(rng, current_seed)

                # Reset strictly to the single-linkage base state for every restart
                # copyto!(worker, init_clustering)
                # copyto!(currentbest, init_clustering)
                copyto!(worker, bestoverall)
                currentbest_cost, failures = hillclimb_kmeans!(M, initial_cost, worker, currentbest, rng, epochs=epochs, max_rejections=max_rejections, maxiter=maxiter, logiter=logiter, logtrace=logtrace, restart_seed=restart_seed)
                total_failures += failures

                if currentbest_cost < bestoverall_cost
                    copyto!(bestoverall, currentbest)
                    bestoverall_cost = currentbest_cost
                    # best_restart = restart_seed
                end
            end
            
            if bestoverall_cost < globalbest_cost
                copyto!(globalbest, bestoverall)
                globalbest_cost = bestoverall_cost
                best_restart = seed_offset
            end
        end
        
        # final_model = kmeans!(M, globalbest.centers, maxiter=1)
        # @assert final_model.totalcost ≈ globalbest_cost
        
        σ2 = (globalbest_cost / (n_samples * n_features)) + 1e-12
        log_L = -0.5 * n_samples * n_features * (log(2 * π * σ2) + 1.0)
        num_params = k * n_features
        
        best_bics[i] = num_params * log(n_samples) - 2.0 * log_L
        best_seeds[i] = best_restart
        convergence_failures[i] = total_failures
        

        
        info_metrics = extract_single_molecule_metrics(globalbest.distances, globalbest_cost)
        mean_assignment_entropy[i] = mean(info_metrics.h_reads)
        
        if max_assignment_entropy < mean_assignment_entropy[i]
            max_assignment_entropy = mean_assignment_entropy[i]
            max_read_assignment_entropy .= info_metrics.h_reads
            
        end
        
        if best_bics[i] < best_bic
            point_costs = Vector{Float64}(undef, n_samples)
            @inbounds for j in 1:n_samples
                c = globalbest.assignments[j]
                point_costs[j] = globalbest.distances[c, j]
            end
            
            

            best_model = KMeansOptResult(
                copy(globalbest.centers),     
                copy(globalbest.assignments),
                point_costs,
                copy(globalbest.counts),
                globalbest_cost,
                info_metrics.probabilities,
                info_metrics.h_reads,
            )
            best_bic = best_bics[i]
        end
    end
    
    best_idx = argmin(best_bics)
    
    return (; 
        optimal_k = ks[best_idx], 
        optimal_model = best_model,
        optimal_seed = best_seeds[best_idx],
        optimal_bic = best_bic,
        bics = best_bics,
        seeds = best_seeds,
        failures = convergence_failures,
        ks,
        mean_assignment_entropy,
        max_read_assignment_entropy,
        logtrace
    )
end


###### Results tables
function clustercenters(M, model)
    centers = zeros(Float64, size(M, 1), length(model.counts))
    @inbounds for (i, a) in enumerate(model.assignments)
       @views centers[:, a] .+= M[:, i]
    end
    centers./model.counts'
end


function characterise_nucleosomes_spectral!(ac_signal::AbstractVector{Float64}, lags::AbstractVector; 
                                           min_nrl=140.0, max_nrl=240.0, 
                                           phase_threshold=0.04, n = length(ac_signal), detrended = zeros(Float64, n)) #, f = Figure(), k = 1) # Note: Scale is different, needs tuning
    
    n = length(ac_signal)
    
    # 1. Detrend the signal to remove the global exponential decay.
    # We use a moving average roughly the size of one nucleosome (~160bp).
    # This centers the true structural oscillations around 0.0.
    detrended .= 0
    
    # Estimate how many bins make up ~160bp dynamically based on the lags array
    dx_approx = lags[2] - lags[1]
    half_w = max(1, round(Int, 160 / dx_approx) ÷ 2)
    
    @inbounds for i in 1:n
        min_idx = max(1, i - half_w)
        max_idx = min(n, i + half_w)
        
        m = 0.0
        for j in min_idx:max_idx
            m += ac_signal[j]
        end
        m /= (max_idx - min_idx + 1)
        
        detrended[i] = ac_signal[i] - m
    end
    
   
    
    # 2. Targeted Discrete-Time Fourier Transform (Matched Filter)
    best_score = -Inf
    best_nrl = 0.0
    
    # We can step by 1.0bp because we are generating the wave analytically, 
    # completely bypassing FFT resolution limits.
    for test_nrl in min_nrl:1.0:max_nrl
        score = 0.0
        
        @inbounds for i in 1:n
            # A pure cosine perfectly models an autocorrelation signal:
            # Peak at 0, 1*nrl, 2*nrl... Trough at 0.5*nrl, 1.5*nrl...
            template_wave = cos(2.0 * π * lags[i] / test_nrl)
            score += detrended[i] * template_wave
        end
        
        # Normalize by length to get an interpretable power score
        score /= n
        
        if score > best_score
            best_score = score
            best_nrl = test_nrl
        end
    end
    # f = Figure(size=(400, 300))
    # ax = Axis(f[cld(k, 5), mod(k - 1, 5) + 1], title="Cluster=$k\nbest_nrl = $best_nrl\nphasing_score= $best_score")
    #  lines!(lags, ac_signal, label="ac_signal")
    # lines!(lags, detrended, label="detrend")
    # vlines!(best_nrl, label="$best_nrl")
    # axislegend()
    # # display(f)
    
    is_closed = best_score > phase_threshold
    
    # 3. Optional: Determine if the 2N harmonic is dominant (like your previous peak_type)
    peak_type = :none
    if is_closed
        # Fast linear interpolation to check the detrended heights at 1N vs 2N
        idx_1N = searchsortedfirst(lags, best_nrl)
        idx_2N = searchsortedfirst(lags, best_nrl * 2.0)
        
        val_1N = (1 <= idx_1N <= n) ? detrended[idx_1N] : 0.0
        val_2N = (1 <= idx_2N <= n) ? detrended[idx_2N] : 0.0
        
        peak_type = val_2N > val_1N ? :dinucleosome : :mononucleosome
    end
    
    return (; nrl=is_closed ? round(Int, best_nrl) : missing, 
              phasing_score=best_score, 
              is_closed=is_closed,
              peak_type=peak_type)
end

function interval_count(data, loc)
    start = max(data[1], first(loc))
    stop = min(data[2], last(loc))
    
    max(stop - start + 1, 0)
end



"""
    readtable(chrom, loc, reads, clusterresult, featmat; kwargs...)

Builds a DataFrame summarizing statistics per read (molecule) over the specified region.
Includes metrics such as read length, modification rates, feature densities, and 
spectral nucleosome repeat lengths (NRLs).
"""
function readtable(chrom, loc, reads, clusterresult, featmat; nrl_feature=:fire, unifystrand=true, labelfun = x -> replace(basename(x), r".bam|.cram|.haplotagged" => ""))
    n = length(reads)
    cluster = clusterresult.optimal_model.assignments
    cost   = clusterresult.optimal_model.costs
    max_assignment_entropy = clusterresult.max_read_assignment_entropy
    strand = [r.strand for r in reads]
    
    chrom = [chrom for r in reads]
    sample = [labelfun(r.label) for r in reads]
    start = [max(r.lp, first(loc)) for r in reads]
    stop  = [min(r.rp, last(loc)) for r in reads]
    hap   = [r.haplotype for r in reads]
    
    readlength = stop .- start .+ 1
    
    total_mod_6mA   = [sum((p ∈ loc for p in r.mods_6mA); init=0) for r in reads]
    total_unmod_6mA = [sum((p ∈ loc for p in r.unmods_6mA); init=0) for r in reads]
    rate_mod_6mA = total_mod_6mA./(total_mod_6mA .+ total_unmod_6mA)
    
    total_mod_5mC   = [sum((p ∈ loc for p in r.mods_5mC); init=0) for r in reads]
    total_unmod_5mC = [sum((p ∈ loc for p in r.unmods_5mC); init=0) for r in reads]
    rate_mod_5mC = total_mod_5mC./(total_mod_5mC .+ total_unmod_5mC)
    
    total_mod_5hmC   = [sum((p ∈ loc for p in r.mods_5hmC); init=0) for r in reads]
    total_unmod_5hmC = [sum((p ∈ loc for p in r.unmods_5hmC); init=0) for r in reads]
    rate_mod_5hmC = total_mod_5hmC./(total_mod_5hmC .+ total_unmod_5hmC)
    
    total_msp  = [sum((interval_count(m, loc) for m in r.msps if m[3] <= 0.9*255); init=0) for r in reads]
    total_fire = [sum((interval_count(m, loc) for m in r.msps if m[3]  > 0.9*255); init=0) for r in reads]
    total_nuc  = [sum(interval_count(m, loc) for m in r.nucs) for r in reads]
    
    rate_msp = total_msp./readlength
    rate_fire = total_fire./readlength
    rate_nuc = total_nuc./readlength
    
    ### read level nrls
    
    FM = featureblock(featmat, nrl_feature)
    A = autocor(FM, 1:(size(FM, 1)-1))
    detrended = Vector{Float64}(undef, size(FM, 1) - 1)
    
    nrls = [characterise_nucleosomes_spectral!(a, featmat.ldx[1:end-1], detrended=detrended) for a in eachcol(A)]
    

    index = 1:n
    (; index, cluster, cost, max_assignment_entropy,
        sample, chrom, start, stop, strand, hap, readlength,
        total_mod_6mA,
        total_unmod_6mA,
        rate_mod_6mA,
        total_mod_5mC,
        total_unmod_5mC,
        rate_mod_5mC,
        total_mod_5hmC,
        total_unmod_5hmC,
        rate_mod_5hmC,
        total_msp,
        total_fire,
        total_nuc,
        rate_msp,
        rate_fire,
        rate_nuc,
        has_fire = rate_fire .> 0,
        nrl=[n.nrl for n in nrls],
        phasing_score=[n.phasing_score for n in nrls],
        phasing_implies_closed=[n.is_closed for n in nrls],
        phasing_peak_type=[n.peak_type for n in nrls]) |> DataFrame
end


function filtmean(x, filtfun=ismissing)
    total = sum(filtfun, x)
    (total == length(x)) && return missing
    
    mean(v for v in x if !filtfun(v))
end

function countstring(v::AbstractVector{T}) where {T}
    counts = Dict{T, Int}()
    
    for a in v
        c = get(counts, a, 0)
        counts[a] = c + 1
    end
    
    ks = sort(collect(keys(counts)))
    
    iob = IOBuffer()
    for (i, k) in enumerate(ks)
       print(iob, k, ":", counts[k], (i < length(ks)) ? "," : "") 
    end
    String(take!(iob))
end


"""
    clustertable(readtable, chrom, loc, reads, clusterresult, featmat; kwargs...)

Aggregates read-level statistics from `readtable` into a cluster-level summary DataFrame.
Provides metrics such as mean modification rates, phasing scores, and feature depths 
for each discovered cluster in the locus.
"""
function clustertable(readtable, chrom, loc, reads, clusterresult, featmat; nrl_feature=:fire)
    gdf = combine(groupby(readtable, [:chrom, :start, :stop, :cluster]),
        :cost => mean => :meancost,
        :max_assignment_entropy => mean => :mean_max_assignment_entropy,
        :sample => countstring => :counts_sample,
        :strand => mean => :proportion_positive_strand,
        :hap => countstring => :counts_hap,
        :readlength => mean => :mean_readlength,
        :rate_mod_6mA => mean => :mean_rate_mod_6mA,
        :rate_mod_5mC => mean => :mean_rate_mod_5mC,
        :rate_mod_5hmC => mean => :mean_rate_mod_5hmC,
        :rate_nuc => mean => :mean_rate_nuc_bp,
        :rate_msp => mean => :mean_rate_msp_bp,
        :rate_fire => mean => :mean_rate_fire_bp,
        :has_fire => mean => :proportion_has_fire,
        :nrl => (n -> sum(!ismissing, n)) => :total_valid_nrls,
        :nrl => (n -> filtmean(n, ismissing)) => :mean_read_nrl,
        :phasing_score => (p -> filtmean(p, isinf)) => :mean_read_phasing_score,
        )
    
    ### Build cluster-level metrics DataFrame
    k = length(clusterresult.optimal_model.counts)
    cluster_metrics = DataFrame(cluster = 1:k)

    ### fire depth
    fire_pile = clustercenters(featurematrix(loc, reads; features=(:fire_elem, ), weight_fire=1, dx=1).M, clusterresult.optimal_model)

    cluster_metrics.max_fire_depth_rate = replace(vec(maximum(fire_pile, dims=1)), NaN => 0.0)
    cluster_metrics.max_fire_depth = Int.(round.(clusterresult.optimal_model.counts .* cluster_metrics.max_fire_depth_rate))
    
    #### nrl_stats
    FM = featureblock(featmat, nrl_feature, M=clusterresult.optimal_model.centers)

    A = autocor(FM, 1:(size(FM, 1)-1))
    detrended = Vector{Float64}(undef, size(FM, 1) - 1)
    
    nrls = [characterise_nucleosomes_spectral!(a, featmat.ldx[1:end-1], detrended=detrended) for a in eachcol(A)]
    
    cluster_metrics = hcat(cluster_metrics, DataFrame(nrls))
    
    # Safely merge cluster metrics back to the read-populated clusters in gdf
    leftjoin(gdf, cluster_metrics, on=:cluster)
end

"""
    modclustertable_5mC(chrom, loc, reads, clusterresult; kwargs...)

Generates a DataFrame summarizing 5mC modification calls and rates per position 
across the given clusters, filtered by the minimum number of modification calls.
"""
function modclustertable_5mC(chrom, loc, reads, clusterresult; unifystrand=true, min_mods=5)
    readid = Int[]
    cluster = Int[]
    pos = Int[]
    mod = Bool[]

    
    
    for (k, (r, a)) in enumerate(zip(reads, clusterresult.optimal_model.assignments))
        
        for p in r.mods_5mC
            if p in loc
                v = (r.strand || !unifystrand) ? p : p - 1
                push!(readid, k)
                push!(cluster, a)
                push!(pos, v)
                push!(mod, true)
            end
        end
        for p in r.unmods_5mC
            if p in loc
                v = (r.strand || !unifystrand) ? p : p - 1
                push!(readid, k)
                push!(cluster, a)
                push!(pos, v)
                push!(mod, false)
            end
        end
    end
    
    mct = DataFrame(readid=readid, cluster=cluster, pos=pos, mod=mod)
    mct.chrom = fill(chrom, nrow(mct))
    mct.start = fill(first(loc), nrow(mct))
    mct.stop  = fill(last(loc), nrow(mct))
    
    mcounts = combine(groupby(mct, [:chrom, :start, :stop, :cluster, :pos]), 
        :mod => (x -> count(identity, x)) => :mod,
        :mod => (x -> count(!, x)) => :unmod
    )
    
    if nrow(mcounts) == 0
        mcounts.total = Int[]
        mcounts.rate = Float64[]
        mcounts.total_mod = Int[]
        mcounts.total_unmod = Int[]
        mcounts.total_rate = Float64[]
        mcounts.pvalue = Float64[]
        return mcounts
    end

    mcounts.total = mcounts.mod .+ mcounts.unmod
    mcounts.rate = mcounts.mod ./ mcounts.total
    mcounts = mcounts[mcounts.total .> min_mods, :]
    
    if nrow(mcounts) == 0
        mcounts.total_mod = Int[]
        mcounts.total_unmod = Int[]
        mcounts.total_rate = Float64[]
        mcounts.pvalue = Float64[]
        return mcounts
    end

    maxmods = combine(groupby(mcounts, :pos), :mod => maximum => :max)
    pos_above_5 = maxmods[maxmods.max .> min_mods, :pos]
    mcounts = mcounts[mcounts.pos .∈ Ref(Set(pos_above_5)), :]
    
    if nrow(mcounts) == 0
        mcounts.total_mod = Int[]
        mcounts.total_unmod = Int[]
        mcounts.total_rate = Float64[]
        mcounts.pvalue = Float64[]
        return mcounts
    end

    transform!(groupby(mcounts, :pos), :mod => sum => :total_mod, :unmod => sum => :total_unmod)
    mcounts.total_rate = mcounts.total_mod ./ (mcounts.total_mod .+ mcounts.total_unmod)
    bts = BinomialTest.(mcounts.mod, mcounts.total, mcounts.total_rate)
    mcounts.pvalue = pvalue.(bts)
    return mcounts
end

###### Main loop


"""
    cluster_peaks(peakfile, bams, outdir; kwargs...)

Run stencil clustering on all peaks defined in a BED file. This is the main entry point 
for the clustering workflow. Computes clustering for each region and writes results to `outdir`.
It executes in parallel across regions if Julia is started with multiple threads.
"""
function cluster_peaks(peakfile, bams, outdir;  dw=600, features=(:fire, :mod_6mA), dx=15, smooth_sigma=0.0, strictfirefeatures=1, unify_gc_strand=true, ks=[2:20 ; 32:2:60], n_ils_restarts=2, n_restarts=20, random_seed_offset = 0,
                epochs=1000, rejection_neighbourhood_fraction=0.08, max_rejections=150, maxiter=300, linkage=:single, logiter=false, label="")


    smglog("Stencil Clustering  : ")
    smglog("Peakfile            : ", peakfile)
    smglog("BAMs                : ", join(bams, ", "))
    smglog("Outdir              : ", outdir)
    ### log for parameters?
    peaks = loadbed(peakfile)
    
    smglog("Loaded              : ", size(peaks, 1), " regions")

    chroms = peaks.chrom
    starts = peaks.start
    stops  = peaks.stop

    starttime = time()
    Threads.@threads for i in 1:size(peaks, 1)
        
        chrom = String(chroms[i])
        peakloc = starts[i]:stops[i]
        smglog("Clustering              : ", chrom, ":", first(peakloc), "-", last(peakloc)) ### smglog is thread safe

        rng = Xoshiro()
        clustering = cluster_region(chrom, peakloc, bams; dw=dw, features=features, dx=dx, smooth_sigma=smooth_sigma, strictfirefeatures=strictfirefeatures, unify_gc_strand=unify_gc_strand, ks=ks,
                 n_ils_restarts=n_ils_restarts, n_restarts=n_restarts, random_seed_offset = random_seed_offset,
                epochs=epochs, rejection_neighbourhood_fraction=rejection_neighbourhood_fraction, max_rejections=max_rejections, maxiter=maxiter, 
                linkage=linkage, rng=rng, logiter=logiter)

        saveclustering(outdir, clustering, label)
    end

    smglog("Clustering complete in :", time() - starttime, " seconds")

end


to_portable(r::ModAllFireRead) = NamedTuple{propertynames(r)}(Tuple(getproperty(r, p) for p in propertynames(r)))
to_portable(o::KMeansOptResult) = NamedTuple{propertynames(o)}(Tuple(getproperty(o, p) for p in propertynames(o)))

# 2. Re-pack the payload before saving
function save_portable_clustering(cluster_path, payload)
    
    # Surgically convert only the optimal_model field inside the NamedTuple
    portable_clustering = merge(
        payload.optimal_clustering, 
        (; optimal_model = to_portable(payload.optimal_clustering.optimal_model))
    )

    # Re-pack the main payload
    portable_payload = (; 
        chrom = payload.chrom, 
        peakloc = payload.peakloc, 
        loc = payload.loc, 
        reads = map(to_portable, payload.reads), 
        CM = payload.CM, 
        optimal_clustering = portable_clustering, # Use the patched NamedTuple
        readtable = payload.readtable, 
        clustertable = payload.clustertable, 
        modclustertable = payload.modclustertable
    )
    
    JLD2.save_object(cluster_path, portable_payload)
end

"""
    saveclustering(outdir, clustering, label)

Serializes the clustering results and writes the summary tables (readtable, 
clustertable, modclustertable) to compressed TSV files.
"""
function saveclustering(outdir, clustering, label)

    regionlabel = string("smgstats_cluster_", 
        isempty(label) ? join(sort(unique(clustering.readtable.sample)), ".") : label, "_",
        clustering.chrom, "_",
        first(clustering.peakloc), "_",
        last(clustering.peakloc), "_",
        first(clustering.loc), "_",
        last(clustering.loc), "_",
        join(clustering.CM.features, "."))
    smglog("Saving             :  ", regionlabel)

    mkpath(outdir)
    
    cluster_path         = joinpath(outdir, regionlabel * ".jld2")
    readtable_path       = joinpath(outdir, regionlabel * ".readtable.tsv.gz")
    clustertable_path    = joinpath(outdir, regionlabel * ".clustertable.tsv.gz")
    modclustertable_path = joinpath(outdir, regionlabel * ".modclustertable.5mC.tsv.gz")

    save_portable_clustering(cluster_path, clustering)
    CSV.write(readtable_path, clustering.readtable, delim='\t', compress=true)
    CSV.write(clustertable_path, clustering.clustertable, delim='\t', compress=true)
    CSV.write(modclustertable_path, clustering.modclustertable, delim='\t', compress=true)


    smglog("Written        : ", cluster_path)
    smglog("Written        : ", readtable_path)
    smglog("Written        : ", clustertable_path)
    smglog("Written        : ", modclustertable_path)



end

"""
    cluster_region(chrom, peakloc, bams; kwargs...)

Performs single-molecule clustering on a specific genomic region.
Extracts reads, computes the feature matrix, determines the optimal K using BIC 
and iterated local search k-means, and returns the models and metric tables.
"""
function cluster_region(chrom, peakloc, bams; dw=600, features=(:fire, :mod_6mA), dx=10, smooth_sigma = 0.0, strictfirefeatures=1,
              unify_gc_strand=true, ks=2:2:60, n_ils_restarts=1, n_restarts=20, random_seed_offset = 0,
               max_rejections=150, maxiter=300, linkage=linkage, rng=Xoshiro(), logiter=false,
                epochs=1000, rejection_neighbourhood_fraction=0.08)

    loc = (first(peakloc) - dw):(last(peakloc) + dw)


    reads = collectstencilreads(bams, chrom, loc, strictfirefeatures=strictfirefeatures)
    CM = featurematrix(loc, reads, dx=dx, features=features, unify_gc_strand=unify_gc_strand, smooth_sigma=smooth_sigma)
    
    optimal_clustering  = select_k_kmeans_bic(CM.M, ks, n_ils_restarts=n_ils_restarts, n_restarts=n_restarts, epochs=epochs, linkage=linkage, logiter=logiter, 
                                rejection_neighbourhood_fraction=rejection_neighbourhood_fraction, max_rejections=max_rejections, maxiter=maxiter, random_seed_offset=random_seed_offset, rng=rng);
    ### result tables

    rt  = readtable(chrom, loc, reads, optimal_clustering, CM)
    
    ct  = clustertable(rt, chrom, loc, reads, optimal_clustering, CM)
    mct = modclustertable_5mC(chrom, loc, reads, optimal_clustering) 


    ( ; chrom, peakloc, loc, reads, CM, optimal_clustering, readtable=rt, clustertable=ct, modclustertable=mct)                        

end



