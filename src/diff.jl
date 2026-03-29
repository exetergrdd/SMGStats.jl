
    # Helper to sum counts across a group of BAMs
    function sum_counts(files, interval_set, label_prefix)
        # We need to accumulate generic counts.
        # countfire returns (fire_counts, read_counts) as Vectors.
        total_fire = zeros(Int, length(interval_set))
        total_read = zeros(Int, length(interval_set))
        
        for file in files
            smglog("Processing $label_prefix: $file")
            fc, rc = countfire(file, interval_set; total_reads=-1)
            total_fire .+= fc
            total_read .+= rc
        end
        return total_fire, total_read
    end



function fisher(a, b, c, d)
    try
        f = FisherExactTest(a, b, c, d)
        or = f.ω
        logor_se = sqrt.(1.0./a .+ 1.0./b .+ 1.0./c .+ 1.0./d)
        fisher_lpv = pvalue(f, tail=:left)
        fisher_bpv = pvalue(f, tail=:both)
        fisher_rpv = pvalue(f, tail=:right)
        return (;or, logor_se, fisher_lpv, fisher_bpv, fisher_rpv)
    catch
        smglog("FisherExactTest failed for $a, $b, $c, $d")
        return (; or=NaN, logor_se=NaN, fisher_lpv=NaN, fisher_bpv=NaN, fisher_rpv=NaN)
    end

end

function binomialtest(x, n, p)
   try
        bt = BinomialTest(x, n, p)
        return (;lpv=pvalue(bt, tail=:left), bpv=pvalue(bt, tail=:both), rpv=pvalue(bt, tail=:right))
    catch
        smglog("BinomialTest failed for $x, $n, $p")
        return (; lpv=NaN, bpv=NaN, rpv=NaN)
    end
end

"""
    differentialaccessibility(bedfile, bamsA, bamsB, output; method=:fisher, chromlabel=:chrom, startlabel=:start, stoplabel=:stop)

Performs differential accessibility analysis between two groups of BAM files (A and B).
Counts fires and reads for all files, pools them by group, and calculates p-values.

# Arguments
- `bedfile`: Path to peaks BED file.
- `bamsA`: Vector of BAM files for Group A (Control).
- `bamsB`: Vector of BAM files for Group B (Treatment).
- `output`: Output filename.
- `method`: Statistical method (`:fisher`, `:binomial_pooled`, `:binomial_ref`, `:alll`).
"""
function differentialaccessibility(bedfile::String, bamsA, bamsB, output::String; method=:fisher, chromlabel=:chrom, startlabel=:start, stoplabel=:stop)

    smglog("               Differential Accessibility Analysis")
    smglog("Group A:      ", bamsA)
    smglog("Group B:      ", bamsB)
    smglog("Output:       ", output)
    smglog("Method:       ", method)
    # Load peaks
    bdf = loadbed(bedfile)
    intervals = bedintervals(bdf; chromlabel=Symbol(chromlabel), startlabel=Symbol(startlabel), stoplabel=Symbol(stoplabel))
    smglog("Loaded $(length(intervals)) intervals for differential analysis")



    # Process Group A
    smglog("Counting fires for Group A ($(length(bamsA)) files)...")
    fireA, readA = sum_counts(bamsA, intervals, "Group A")

    # Process Group B
    smglog("Counting fires for Group B ($(length(bamsB)) files)...")
    fireB, readB = sum_counts(bamsB, intervals, "Group B")
    
    # Construct Results DataFrame
    
    bdf.FireA = fireA
    bdf.ReadA = readA
    bdf.FireB = fireB
    bdf.ReadB = readB
    
    # Perform Stats
    smglog("Calculating statistics using method: $method")
    
    
    
    if method == :fisher || method == :all
        # Fisher's Exact Test
        # Table:
        #           Fire    NoFire
        # Group A   Fa      Ra - Fa
        # Group B   Fb      Rb - Fb

        # fts = [FisherExactTest(fireA[i], readA[i] - fireA[i], fireB[i], readB[i] - fireB[i]) for i in 1:length(intervals)]
        fts = [fisher(fireA[i], readA[i] - fireA[i], fireB[i], readB[i] - fireB[i]) for i in 1:length(intervals)]
        bdf.a = fireA
        bdf.b = readA .- fireA
        bdf.c = fireB
        bdf.d = readB .- fireB
        bdf.OR = [f.or  for f in fts]
        bdf.logor_se = [f.logor_se  for f in fts]
        bdf.fisher_lpv = [f.fisher_lpv  for f in fts]
        bdf.fisher_bpv = [f.fisher_bpv  for f in fts]
        bdf.fisher_rpv = [f.fisher_rpv  for f in fts]

    end
        
    if method == :binomial_pooled || method == :all
        # Shared rate p = (Fa + Fb) / (Ra + Rb)
        bdf.TotalReads = readA .+ readB
        bdf.TotalFires = fireA .+ fireB
        bdf.p_pool = bdf.TotalFires ./ bdf.TotalReads
        bts = [binomialtest(fireB[i], readB[i], bdf.p_pool[i]) for i in 1:length(intervals)]
        bdf.binpool_lpv = [bt.lpv for bt in bts]
        bdf.binpool_bpv = [bt.bpv for bt in bts]
        bdf.binpool_rpv = [bt.rpv for bt in bts]

    end
    if method == :binomial_ref || method == :all
        # Reference rate derived from A: p = Fa / Ra
        bdf.p_ref = fireA ./ readA
        bts = [binomialtest(fireB[i], readB[i], bdf.p_ref[i]) for i in 1:length(intervals)]
        bdf.binref_lpv = [bt.lpv for bt in bts]
        bdf.binref_bpv = [bt.bpv for bt in bts]
        bdf.binref_rpv = [bt.rpv for bt in bts]

    end

    
    # Write output
    smglog("Writing differential analysis results to: $output")
    if endswith(output, ".gz")
        CSV.write(output, bdf, delim='\t', compress=true)
    else
        CSV.write(output, bdf, delim='\t')
    end
end
