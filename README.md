# SMGStats


**SMGStats** is a Julia package for Single Molecule Genomics Statistics. It provides tools for analyzing single-molecule sequencing data, including standard statistics reporting, coaccessibility analysis, fire counting, and differential accessibility testing.

## Installation

If you are a julia user, you can install SMGStats.jl using the following command in the Julia REPL:

```julia
] # to enter pkg mode
add "https://github.com/owensnick/SMGStats.jl"
```

If you predominantly want to run at the command line, you can install the bash wrapper script by running

``` bash
cd ~/bin
git clone https://github.com/owensnick/SMGStats.jl.git
cd SMGStats.jl/
julia
] # to enter pkg mode
add activate .
instantiate
Ctrl + d # to exit 
cd ..
ln -s SMGStats.jl/bin/smgstats smgstats
```

Now you can run `smgstats --help` to see the available commands.


### 1. Standard Statistics (`stats`)

Generates comprehensive statistics and HTML reports for a BAM/CRAM file.

```bash
smgstats stats -c config.yaml [-o output_dir] input.bam
```

### 2. Coaccessibility (`coaccess`)

Analyzes coaccessibility of fibers across defined peaks.

```bash
smgstats coaccess --peaks peaks.bed --output coaccess_results.tsv.gz sample1.bam sample2.bam ...
```

### 3. Fire Counting (`firecount`)

Counts "fires" (MSP interactions) within defined peaks.

```bash
smgstats firecount  --peaks peaks.bed --output fire_counts.tsv.gz  sample1.bam sample2.bam ...
```

### 4. Differential Accessibility (`diff`)

Performs differential accessibility analysis between two groups of samples (Control vs Treatment).

```bash
smgstats diff --peaks peaks.bed --groupA control1.bam,control2.bam --groupB treatment1.bam,treatment2.bam --output diff_results.tsv.gz --method all
```
*   **Methods**: `fisher` (default), `binomial_pooled`, `binomial_ref`, `all`.

---

## Interactive Usage (Julia REPL)

You can also use the package interactively in Julia.

### Standard Statistics

```julia
using SMGStats

# Run stats on a single file
runstats("input.bam", "config.yaml"; outputdir="stats_output")
```

### Coaccessibility

```julia
using SMGStats

# Define files and peaks
files = ["sample1.bam", "sample2.bam"]
peaks = "peaks.bed"

# Run wrapper
coaccess_wrapper(files, peaks, "coaccess_results.tsv.gz")
```

### Fire Counting

```julia
using SMGStats

files = ["sample1.bam", "sample2.bam"]
countfire(files, "peaks.bed", "fire_counts.tsv.gz")
```

### Differential Accessibility

```julia
using SMGStats

groupA = ["control1.bam", "control2.bam"]
groupB = ["treat1.bam", "treat2.bam"]

# Run differential analysis
differentialaccessibility("peaks.bed", groupA, groupB, "diff_results.tsv.gz"; method=:fisher)
```
