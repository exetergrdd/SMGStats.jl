# Comprehensive Guide to `smgstats cluster`

`smgstats cluster` is a single-molecule chromatin fiber clustering pipeline designed for long-read sequencing technologies (such as PacBio Fiber-seq and Oxford Nanopore single-molecule foot-printing / stencilling). 

Rather than averaging accessibility and epigenetic modifications across bulk cell populations, `smgstats cluster` evaluates individual sequenced DNA molecules spanning genomic regions of interest (e.g., regulatory elements, promoters, enhancers, or CTCF sites). It decomposes complex, heterogeneous chromatin mixtures into distinct molecular sub-states—revealing transcription factor binding footprints, nucleosome positions, heterochromatin closure, and allele-specific epigenetic patterns.

---

## Table of Contents

1. [Quick Start & Command Synopsis](#1-quick-start--command-synopsis)
2. [Input Arguments & CLI Options](#2-input-arguments--cli-options)
3. [Algorithmic Principles & Methodology](#3-algorithmic-principles--methodology)
   - [3.1 Stencil Matrix Construction (Spatial Binning & Smoothing)](#31-stencil-matrix-construction-spatial-binning--smoothing)
   - [3.2 The Single-Molecule Clustering Challenge](#32-the-single-molecule-clustering-challenge)
   - [3.3 Initial Hierarchical Clustering](#33-initial-hierarchical-clustering)
   - [3.4 Iterated Local Search (ILS) with Stochastic Hill Climbing](#34-iterated-local-search-ils-with-stochastic-hill-climbing)
   - [3.5 Model Selection via Bayesian Information Criterion (BIC)](#35-model-selection-via-bayesian-information-criterion-bic)
   - [3.6 Assignment Probabilities & Read Entropy](#36-assignment-probabilities--read-entropy)
   - [3.7 Spectral Nucleosome Repeat Length (NRL) & Phasing](#37-spectral-nucleosome-repeat-length-nrl--phasing)
   - [3.8 Spatial Cluster Ordering & Post-Hoc Relabeling](#38-spatial-cluster-ordering--post-hoc-relabeling)
   - [3.9 Cluster Composition & Fisher Skew Tests (Strand, Haplotype, Sample)](#39-cluster-composition--fisher-skew-tests-strand-haplotype-sample)
4. [Output Files & Directory Structure](#4-output-files--directory-structure)
5. [Output Table Schemas](#5-output-table-schemas)
   - [5.1 Read Table (`*.readtable.tsv.gz`)](#51-read-table-readtabletsvgz)
   - [5.2 Cluster Table (`*.clustertable.tsv.gz`)](#52-cluster-table-clustertabletsvgz)
   - [5.3 5mC Modification Table (`*.modclustertable.5mC.tsv.gz`)](#53-5mc-modification-table-modclustertable5mctsvgz)
   - [5.4 Serialized Portable Data (`*.jld2`)](#54-serialized-portable-data-jld2)
6. [Generating Interactive HTML Reports & Master Catalog (`smgstats report`)](#6-generating-interactive-html-reports--master-catalog-smgstats-report)
   - [6.1 Report CLI Options](#61-report-cli-options)
   - [6.2 Multi-Process Parallel Report Rendering (`--workers` & `--chunk`)](#62-multi-process-parallel-report-rendering---workers---chunk)
7. [Glossary of Key Terms](#7-glossary-of-key-terms)
8. [Recommended Parameter Presets](#8-recommended-parameter-presets)

---

## 1. Quick Start & Command Synopsis

`smgstats cluster` can be invoked using the `smgstats` command-line executable. Multi-threading is supported and configured via the `-t` or `--threads` option.

```bash
# Basic invocation with 8 CPU threads
smgstats -t 8 cluster \
    --peaks peaks.bed \
    --output ./cluster_results \
    sample1.bam sample2.bam
```

### Typical Output Files Per Region
For each region in `peaks.bed`, four files are generated:
```
cluster_results/
├── smgstats_cluster_sample1_chr10_69347828_69350020_69347228_69350620_fire.mod_6mA.readtable.tsv.gz
├── smgstats_cluster_sample1_chr10_69347828_69350020_69347228_69350620_fire.mod_6mA.clustertable.tsv.gz
├── smgstats_cluster_sample1_chr10_69347828_69350020_69347228_69350620_fire.mod_6mA.modclustertable.5mC.tsv.gz
└── smgstats_cluster_sample1_chr10_69347828_69350020_69347228_69350620_fire.mod_6mA.jld2
```

---

## 2. Input Arguments & CLI Options

The CLI accepts both required positional/named parameters and optional algorithmic tuning flags.

### 2.1 Core Inputs

| Argument | Flag | Type | Default | Description |
| :--- | :--- | :--- | :--- | :--- |
| `files` | *(positional)* | Strings (space-separated) | **Required** | One or more indexed BAM or CRAM files containing single-molecule methylation and footprint calls (e.g., tags `MM`/`ML` for 6mA/5mC, or FIRE tags). |
| `--peaks` | `-p` | File Path | **Required** | BED or TSV file defining genomic regions to cluster. Must contain at least `chrom`, `start`, and `stop` (0-based or 1-based, with or without header). Optional `PeakName` (e.g. `ctcf_1_fp_1`) and `Category` (e.g. `AccessibleTSS`, `FirePeak`) columns are automatically preserved and propagated to all output tables and `.jld2` payloads; if omitted, they default to `peak_1, peak_2, ...` and `ClusterPeak`. |
| `--output` | `-o` | Directory Path | **Required** | Directory where output TSV tables and JLD2 files will be written. Directory is created automatically if it does not exist. |
| `--samples` | `-s` | Comma-separated Strings | *Auto-inferred* | Optional custom sample labels for each input BAM/CRAM file (e.g. `endo_as1,endo_as2,endo_wgs` or `050526=endo_as1,200126=endo_as2`). When omitted, sample names are automatically extracted from standard ONT/Dorado directory hierarchies (`*/sample_id/sample_id_rep/run_ind/dorado*/...`) and guaranteed to be unique per BAM file. |
| `--threads` | `-t` | Integer / `auto` | `auto` | Number of CPU worker threads used to cluster genomic regions in parallel. Set in the `smgstats` bash script. |

---

### 2.2 Genomic Window & Feature Extraction Options

| Argument | Flag | Type | Default | Description |
| :--- | :--- | :--- | :--- | :--- |
| `--dw` | | Integer | `600` | Flanking window extension (bp) added to both sides of each peak: $\text{Locus} = [\text{start} - dw, \text{stop} + dw]$. Allows contextual chromatin (flanking nucleosome arrays) to inform the clustering. |
| `--features` | `-f` | Comma-delimited | `fire,mod_6mA` | Features extracted into the clustering matrix. Choices: `fire`, `fire_elem`, `linker`, `nuc`, `nuc_dyad`, `mod_6mA`, `mod_5mC`, `mod_5hmC`, `unmod_6mA`, `unmod_5mC`, `unmod_5hmC`. |
| `--dx` | | Integer | `15` | Spatial binning resolution in base pairs. Converts 1 bp resolution footprints into mean feature rates over bins of size $\Delta x$. |
| `--smooth_sigma` | | Float | `0.0` | Standard deviation for Gaussian smoothing across genomic positions prior to binning. Set to `0.0` for no smoothing. |
| `--strictfirefeatures` | | Integer | `1` | Minimum span requirement (bp). Only reads whose annotated features extend at least this many bp beyond both boundaries of the locus window are kept. Prevents truncated read artifacts. |
| `--unify_gc_strand` | | Boolean (`true`/`false`) | `true` | For CpG modifications (`5mC`, `5hmC`), subtracts 1 from negative-strand coordinates so that symmetrical CpGs on opposite strands map to identical coordinate bins. |

---

### 2.3 Clustering & Optimization Parameters

| Argument | Flag | Type | Default | Description |
| :--- | :--- | :--- | :--- | :--- |
| `--ks` | | String range(s) | `1:20,22:2:60` | Candidate cluster numbers $k$ evaluated by BIC. Comma-separated integers or Julia-style ranges (e.g. `1:20`, `2:2:30`). Automatically constrained by `--max_k_fraction`. |
| `--max_k_fraction` | | Float | `0.65` | Maximum candidate $k$ as a fraction of total spanning reads $N$, constraining evaluated cluster counts to $1 \le k \le \max(1, \lfloor N \cdot \text{max\_k\_fraction} \rfloor)$. |
| `--n_ils_restarts` | | Integer | `2` | Number of outer Iterated Local Search (ILS) restarts. Higher values reduce the chance of remaining trapped in suboptimal basins. |
| `--n_restarts` | | Integer | `20` | Number of inner random hill-climbing restart passes per ILS cycle. |
| `--random_seed_offset` | | Integer | `0` | Offset added to internal PRNG seeds for reproducible stochastic runs. |
| `--epochs` | | Integer | `1000` | Maximum hill-climbing mutation attempts per restart cycle. |
| `--rejection_neighbourhood_fraction` | | Float | `0.08` | Adaptive stopping threshold fraction. Maximum allowable consecutive rejections during hill climbing is set to $\text{round}(N \cdot k \cdot f)$. |
| `--max_rejections` | | Integer | `150` | Fallback maximum consecutive rejections when `--rejection_neighbourhood_fraction` is set to `0.0`. |
| `--maxiter` | | Integer | `300` | Maximum Lloyd iterations permitted for each individual k-means convergence step. |
| `--linkage` | | String | `single` | Linkage criterion used to build the initial hierarchical clustering tree (`single`, `complete`, `average`). |
| `--logiter` | | Flag | `false` | Verbose logging flag. Records detailed trace tuples of every mutation, distance update, and cost change. |
| `--label` | | String | `""` | Optional prefix string added to output file names. If omitted, sample names from the BAM headers are used. |

---

## 3. Algorithmic Principles & Methodology

Understanding the algorithms inside `smgstats cluster` helps when tuning parameters or interpreting unexpected cluster numbers.

```
                      +-----------------------------+
                      | Single-Molecule Reads (BAM) |
                      +-----------------------------+
                                     |
                                     v
             [ Feature Extraction & Spatial Binning (dx=15bp) ]
                                     |
                                     v
                        +--------------------------+
                        |  Feature Matrix M (D x N) |
                        +--------------------------+
                                     |
        +----------------------------+----------------------------+
        |                                                         |
        v                                                         v
 [ Pairwise Distance Matrix ]                              For each candidate k in ks:
        |                                                         |
        v                                                         v
 [ Hierarchical Clustering ]                               [ cutree(hc, k) Initial State ]
 (Single-Linkage Tree)                                            |
                                                                  v
                                                           [ Lloyd's K-Means Convergence ]
                                                                  |
                                                                  v
                                                           [ Iterated Local Search (ILS) ]
                                                           * Random center teleportation
                                                           * Greedy hill climbing
                                                           * Adaptive rejection stopping
                                                                  |
                                                                  v
                                                           [ Calculate BIC(k) & Entropy ]
                                                                  |
                                                                  v
                                            [ Select Optimal k* = argmin BIC(k) ]
                                                                  |
                                                                  v
                                    [ Spatial Cluster Ordering (default_cluster_order) ]
                                    * Active clusters (> mu + 1.0*sigma) ordered L->R
                                    * Closed/nucleosomal clusters ordered L->R
                                                                  |
                                     +----------------------------+----------------------------+
                                     |                                                         |
                                     v                                                         v
                         [ Dual Spectral NRL Phasing ]                              [ Generate Output Tables ]
                         * FIRE & 6mA autocorrelation                               * readtable.tsv.gz
                         * Moving-average detrending                                * clustertable.tsv.gz
                         * Cosine matched filtering                                 * modclustertable.5mC.tsv.gz
```

---

### 3.1 Stencil Matrix Construction (Spatial Binning & Smoothing)

For a genomic window of length $L = \text{stop} - \text{start} + 1$, each sequenced read covering the region is represented across the specified `--features`.

1. **Native 1 bp Representation**: Each read is encoded on a 1 bp genomic coordinate grid. 
   - Methylation calls (`mod_6mA`, `mod_5mC`, `mod_5hmC`) are thresholded at a probability $\ge 0.9 \times 255 \approx 230$.
   - FIRE accessible elements are weighted by confidence:
     $$\text{Weight} = \begin{cases} 1.0 & \text{if } q > 0.9 \times 255 \text{ (confident FIRE element)} \\ 0.75 & \text{if } q \le 0.9 \times 255 \text{ (linker patch)} \end{cases}$$
   - Nucleosome footprints are assigned a base weight of $0.25$.
2. **Optional Gaussian Smoothing**: A 1D Gaussian filter with standard deviation $\sigma$ (`--smooth_sigma`) can blur discrete base-level events along the genomic axis.
3. **Spatial Binning (`mean_spatial_bin`)**: The 1 bp signal is divided into non-overlapping bins of width $\Delta x = \text{dx}$ (default 15 bp). The value in each bin is the exact arithmetic mean of its constituent base pairs:
   $$\bar{x}_{\text{bin}} = \frac{1}{\Delta x} \sum_{i \in \text{bin}} x_i$$
   Edge bins smaller than $0.75 \times \Delta x$ are merged into the preceding bin to prevent boundary noise.
4. **Stacked Feature Matrix ($M$)**: The resulting binned vectors for each feature channel are stacked vertically into a combined feature matrix $M$ with dimensions:
   $$D \times N$$
   where $N$ is the number of reads, and $D = (\text{number of features}) \times \lceil L / \Delta x \rceil$.

---

### 3.2 The Single-Molecule Clustering Challenge

Chromatin fiber datasets exhibit high dimensionality ($D \sim 100\text{--}400$ bins) and biological sparsity. Crucially, chromatin states are often **asymmetric in population size**:
- A prominent state might contain 80% of fibers (e.g. fully closed nucleosome array).
- A critical regulatory state (e.g., transcription factor binding, open promoter) might represent only 2–5% of molecules.

Standard k-means clustering (Lloyd's algorithm) initialized with random centers or k-means++ tends to:
1. Split large, homogeneous clusters into arbitrary subsets.
2. Miss small, rare biological sub-populations by absorbing them into adjacent large clusters.
3. Become trapped in local energy minima, producing inconsistent cluster counts.

To overcome this, `smgstats cluster` pairs **hierarchical initialization** with **Iterated Local Search (ILS) stochastic hill climbing**.

---

### 3.3 Initial Hierarchical Clustering

Before testing candidate cluster counts, `smgstats cluster` computes the full $N \times N$ pairwise squared Euclidean distance matrix between all reads:
$$d(r_a, r_b) = \|M_{:, a} - M_{:, b}\|^2$$
It constructs a hierarchical tree using **single-linkage** clustering (by default). 

Single linkage groups reads based on minimal point-to-point distances, making it sensitive to continuous density manifolds and distinct, compact outlier groups. For each evaluated $k$, the tree is sliced at $k$ branches (`cutree(hc, k)`), and initial cluster centroids are computed directly from these initial groups. This provides a structured starting configuration superior to random centroid guesses.

---

### 3.4 Iterated Local Search (ILS) with Stochastic Hill Climbing

Starting from the hierarchical baseline, `smgstats cluster` performs iterative local search optimization:

```
                  +-----------------------------------+
                  | Baseline Clustering (cost = C_0)  |
                  +-----------------------------------+
                                    |
                                    v
                       [ Perturbation / Mutation ]
            Select 1 or 2 centers; replace them with randomly
                      sampled read feature vectors
                                    |
                                    v
                     [ Fast Lloyd's K-Means Steps ]
                   Converge centers to local minimum (C_new)
                                    |
                                    v
                     Is C_new < Current Best Cost?
                                   / \
                            YES   /   \   NO
                                 /     \
                                v       v
               [ Accept State ]          [ Reject State ]
               * Update current best     * Revert to previous best
               * Reset rejections = 0    * Increment rejections
                                \       /
                                 \     /
                                  v   v
                 Has consecutive rejections exceeded limit?
                                   / \
                            YES   /   \   NO
                                 /     \
                                v       v
                         [ Terminate ]   [ Next Epoch ]
```

1. **Local Convergence (`converge_kmeans!`)**: Lloyd's algorithm iteratively alternates between assigning each read to the nearest centroid and updating each centroid to the mean of its assigned reads. Center updates and distance recalculations are selectively tracked to skip unaffected clusters.
2. **Perturbation (Mutation)**: With 80% probability, 1 cluster center is selected; with 20% probability, 2 cluster centers are selected. These centers are **relocated directly to the coordinates of randomly chosen reads** in $M$. This acts as a targeted structural disruption that forces centroids out of shallow local basins.
3. **Local Re-convergence**: The perturbed state is run through Lloyd's algorithm until convergence, yielding a new total squared Euclidean error $\text{Cost}_{\text{new}}$.
4. **Greedy Hill Climbing**:
   - If $\text{Cost}_{\text{new}} < \text{Cost}_{\text{best}}$, the move is accepted, the best state is overwritten, and the rejection counter is reset to 0.
   - If $\text{Cost}_{\text{new}} \ge \text{Cost}_{\text{best}}$, the mutation is discarded, the worker state is restored from $\text{Cost}_{\text{best}}$, and the consecutive rejection counter increments.
5. **Adaptive Stopping Criterion**: Optimization halts when consecutive rejections exceed:
   $$R_{\max} = \max\left(100, \, \text{round}(N \cdot k \cdot f_{\text{rejection}})\right)$$
   where $f_{\text{rejection}} = \text{rejection\_neighbourhood\_fraction}$ (default 0.08). If $N=100$ and $k=5$, hill climbing allows up to $100 \times 5 \times 0.08 = 40$ (floored at 100) unsuccessful mutations before concluding the search space has been adequately sampled.
6. **Nested Restarts**:
   - **Inner Restarts (`n_restarts = 20`)**: The hill climb is re-initiated from different random seeds.
   - **Outer Restarts (`n_ils_restarts = 2`)**: The full search cycle is repeated across independent outer passes to verify global convergence.

---

### 3.5 Model Selection via Bayesian Information Criterion (BIC)

To identify the optimal number of clusters $k^*$, `smgstats cluster` evaluates the Bayesian Information Criterion (BIC) across all candidate values in `--ks`.

#### Theoretical Foundation
K-means minimization of squared Euclidean distance is mathematically equivalent to maximum likelihood estimation of a **Gaussian Mixture Model (GMM)** where all clusters share an isotropic (spherical) covariance matrix $\Sigma = \sigma^2 I$.

For a chosen $k$, given total squared error $J = \sum_{j=1}^N \|x_j - c_{a(j)}\|^2$, the estimated observational variance is:
$$\sigma^2 = \frac{J}{N \cdot D}$$

The log-likelihood $\ln L$ of this model is:
$$\ln L = -\frac{N \cdot D}{2} \left[ \ln(2\pi \sigma^2) + 1 \right]$$

The number of free parameters $p$ in the model corresponds to the coordinates of the $k$ centroids in $D$-dimensional space:
$$p = k \cdot D$$

The Bayesian Information Criterion is calculated as:
$$\text{BIC}(k) = p \ln(N) - 2 \ln L = \left(k \cdot D \ln N\right) - 2 \ln L$$

#### Interpreting the Tradeoff
- **Data Fit Term ($-2 \ln L$)**: As $k$ increases, centroids move closer to individual data points, $J$ decreases, and $-2 \ln L$ drops.
- **Complexity Penalty ($k \cdot D \ln N$)**: Each additional cluster adds $D$ dimensions of free parameters, penalizing overly complex models proportionally to sample size $\ln N$.
- **Optimal Model ($k^*$)**:
  $$k^* = \arg\min_{k \in \text{ks}} \text{BIC}(k)$$
  The minimum of the BIC curve represents the optimal balance between capturing true biological sub-states and avoiding overfitting experimental noise.

---

### 3.6 Assignment Probabilities & Read Entropy

Because single-molecule reads can occasionally exhibit intermediate or transitioning chromatin configurations, hard k-means assignments can obscure uncertainty. `smgstats cluster` derives a soft, probabilistic assignment for every read using a **Boltzmann-Gibbs statistical mechanics formulation**:

1. **Effective Temperature ($T$)**:
   The average error per read serves as a natural system temperature:
   $$T = \frac{J}{N}$$
2. **Assignment Probability ($P_{ij}$)**:
   The probability that read $j$ belongs to cluster $i$ with squared distance $d_{ij} = \|x_j - c_i\|^2$ is:
   $$P(i \mid j) = \frac{\exp\left(-\frac{d_{ij} - \min_c d_{cj}}{T}\right)}{\sum_{m=1}^k \exp\left(-\frac{d_{mj} - \min_c d_{cj}}{T}\right)}$$
   *(Subtracting $\min_c d_{cj}$ prevents numerical underflow in floating-point operations.)*
3. **Read Assignment Entropy ($H_j$)**:
   The Shannon entropy of the probability distribution for read $j$ is computed as:
   $$H_j = -\sum_{i=1}^k P(i \mid j) \ln P(i \mid j)$$

#### Biological Meaning of Read Entropy
- **Low Entropy ($H_j \approx 0$)**: The read sits near the centroid of its assigned cluster and far from all other clusters. Its chromatin architecture matches a distinct, well-defined biological state.
- **High Entropy ($H_j \to \ln k$)**: The read is equidistant between two or more cluster centroids. This indicates a molecule caught in an intermediate transition (e.g. partial transcription factor eviction, nucleosome sliding, or technical ambiguity).
- **`max_assignment_entropy`**: The column in `readtable` records the read's entropy at the specific candidate $k$ where the population had highest overall uncertainty, providing a conservative metric of assignment robustness.

---

### 3.7 Spectral Nucleosome Repeat Length (NRL) & Phasing

To characterize whether a cluster or single read represents an organized, regular array of nucleosomes (such as heterochromatin or flanking positioned arrays), `smgstats cluster` performs **spectral matched-filter autocorrelation** independently across **both** the `:fire` (footprint stencil) and `:mod_6mA` (raw adenine methylation) channels:

1. **Autocorrelation Function**: The spatial autocorrelation $A(\tau)$ of the feature signal (for both single reads and cluster centroids) is calculated across relative genomic lags $\tau \in [0, L-\Delta x]$. Evaluating both `:fire` and `:mod_6mA` captures complimentary regimes: `:fire` reflects called footprint regularity, whereas `:mod_6mA` retains unthresholded long-range nucleosome linker oscillations even in noisy or partially accessible fibers.
2. **Moving-Average Detrending**: Real single-molecule chromatin signals exhibit a steep exponential drop at short lags due to self-overlap. To isolate true structural oscillations, a moving average with window width $\approx 160$ bp (the approximate size of one nucleosome) is subtracted:
   $$\tilde{A}(\tau) = A(\tau) - \text{MovingAverage}_{160}(\tau)$$
3. **Discrete Matched Filtering (Wave Decomposition)**: Rather than being constrained by the coarse frequency bins of a standard Fast Fourier Transform (FFT), `smgstats cluster` projects $\tilde{A}(\tau)$ onto analytic cosine template waves across continuous candidate period lengths $L_{\text{test}} \in [140\text{ bp}, 240\text{ bp}]$ at 1.0 bp increments:
   $$\text{Score}(L_{\text{test}}) = \frac{1}{n_{\text{lags}}} \sum_{\tau} \tilde{A}(\tau) \cos\left(\frac{2\pi \cdot \tau}{L_{\text{test}}}\right)$$
   The optimal Nucleosome Repeat Length (NRL) is the period achieving maximum power:
   $$\text{NRL} = \arg\max_{L_{\text{test}}} \text{Score}(L_{\text{test}})$$
4. **Phasing Classification**:
   - **Phasing Score**: The power score of the optimal wave (`phasing_score_fire` and `phasing_score_6mA`).
   - **`is_closed` / `phasing_implies_closed`**: If the phasing score exceeds `phase_threshold = 0.035`, the cluster or read is classified as exhibiting structured, phased nucleosomes (`is_closed_fire` / `is_closed_6mA`).
   - **`peak_type`**: Evaluates harmonic heights at $1 \times \text{NRL}$ vs $2 \times \text{NRL}$ (`peak_type_fire` / `peak_type_6mA`):
     - `:mononucleosome`: Standard nucleosome ladder ($1\text{N} > 2\text{N}$).
     - `:dinucleosome`: Di-nucleosome periodicity or alternating packaging ($2\text{N} > 1\text{N}$).
     - `:none`: Unphased chromatin.

---

### 3.8 Spatial Cluster Ordering & Post-Hoc Relabeling

Because k-means and hierarchical tree cutting assign arbitrary integer labels ($1 \dots k^*$), `cluster_region` automatically applies a deterministic, biologically interpretable **spatial left-to-right ordering** (`default_cluster_order` + `reorder_clusters!`) immediately after BIC model selection:

1. **Feature Selection & Coarse Spatial Binning**:
   - Extracts the `:mod_6mA` centroid sub-matrix from `CM` if `:mod_6mA` is present in `CM.features` (otherwise falls back to `first(CM.features)`). Note that all stencil weights are strictly non-negative (nucleosomes $= 0.25$, linkers $= 0.75$, FIRE $= 1.0$).
   - Applies coarse spatial binning via `CB = mean_spatial_bin(CF, 10)` (grouping 10 stencil bins, e.g. $10 \times 15\text{ bp} = 150\text{ bp}$ macro-bins) to suppress high-frequency single-bin noise and isolate focal regulatory accessibility versus phased nucleosome arrays.
2. **Activity Thresholding & Left-to-Right Diagonalization**:
   - Computes an activity threshold across the coarse-binned centroid matrix:
     $$t = \mu(\text{CB}) + s \cdot \sigma(\text{CB}) \quad (\text{default } s = 1.0)$$
   - For each cluster $c \in \{1 \dots k^*\}$, identifies its peak magnitude and coarse bin index $(m_c, p_c) = \text{findmax}(\text{CB}_{:, c})$.
   - With `active_first = true` (default), clusters whose peak exceeds the activity threshold ($m_c > t$, representing focal open/regulatory states) are placed **first**, ordered left-to-right ($5' \to 3'$) by their peak bin position $p_c$.
   - Remaining clusters ($m_c \le t$, representing closed or nucleosomal arrays) follow, also ordered left-to-right by their peak bin position $p_c$.
3. **Propagation & Post-Hoc Relabeling**:
   - `reorder_clusters!` permutes `centers`, `counts`, `probabilities`, and read `assignments` in-place, and propagates the new labels ($1 \dots k^*$) consistently across `readtable`, `clustertable`, and `modclustertable`.
   - When working with saved `.jld2` files in Julia, clusters can be reordered at any time by arbitrary criteria (e.g., haplotype skew, strand bias, 5mC rate, or custom functions) using `cluster_order` and `reorder_clusters!`.

---

### 3.9 Cluster Composition & Fisher Skew Tests (Strand, Haplotype, Sample)

To distinguish true biological chromatin sub-states from technical sequencing artifacts or to discover allele- and sample-specific regulation, `clustertable` automatically computes $2 \times 2$ **Fisher's Exact Tests** for every cluster $c \in \{1 \dots k^*\}$ (when $k^* > 1$), comparing reads in cluster $c$ against all other reads in the locus ($c' \ne c$):

1. **Technical Strand Skew (`strand_or`, `strand_logor_se`, `strand_pvalue`)**:
   - Tests enrichment of positive-strand (`+`) versus negative-strand (`-`) reads in cluster $c$. Significant strand bias ($p < 0.01$) flags potential strand-specific basecalling or modification-calling artifacts.
2. **Allelic / Haplotype Skew (`hap_or`, `hap_logor_se`, `hap_pvalue`)**:
   - Evaluated across phased reads (`hap == 1` vs `hap == 2`) whenever at least 10 phased reads span the locus. Significant allelic skew identifies haplotype-specific chromatin accessibility, transcription factor binding, or genomic imprinting.
3. **Sample / Condition Skew (`sample_or`, `sample_logor_se`, `sample_pvalue`, `sample_min_pvalue`)**:
   - When multiple BAM files/samples are clustered together, tests each cluster for enrichment in specific samples (for $>2$ samples, records the most strongly enriched sample's odds ratio, standard error, and minimum Fisher $p$-value).

---

## 4. Output Files & Directory Structure

Results are written to the directory specified by `--output`. Filenames are formatted according to the following convention:

```
smgstats_cluster_{label}_{chrom}_{peak_start}_{peak_stop}_{loc_start}_{loc_stop}_{features}.{ext}
```

- `{label}`: User `--label`, or concatenated sample names derived from the input BAM files.
- `{chrom}`, `{peak_start}`, `{peak_stop}`: Coordinates of the input peak interval from `peaks.bed`.
- `{loc_start}`, `{loc_stop}`: Actual analyzed window coordinates, including the `--dw` flanking regions.
- `{features}`: Dot-separated list of clustered features (e.g., `fire.mod_6mA`).

### Generated Artifacts Per Locus

| Extension | File Format | Description |
| :--- | :--- | :--- |
| `*.readtable.tsv.gz` | Gzipped Tab-Separated Values | Single-molecule level assignments, modification counts, feature percentages, and individual read phasing metrics. |
| `*.clustertable.tsv.gz` | Gzipped Tab-Separated Values | Cluster-level summaries: mean costs, sample/haplotype composition, feature occupancy, and spectral NRL statistics. |
| `*.modclustertable.5mC.tsv.gz` | Gzipped Tab-Separated Values | Base-pair resolution 5mC CpG methylation rates per cluster, along with statistical binomial tests for cluster-specific differential methylation. |
| `*.jld2` | Portable HDF5-based JLD2 | Serialized Julia named tuple containing all raw reads, feature matrices, optimal model centroids, assignments, and summary DataFrames. |

---

## 5. Output Table Schemas

### 5.1 Read Table (`*.readtable.tsv.gz`)

Contains one row per sequenced molecule spanning the locus window.

| Column Name | Data Type | Description |
| :--- | :--- | :--- |
| `index` | Integer | 1-based sequential read index within this locus. |
| `cluster` | Integer | Assigned cluster index ($1 \dots k^*$) for this molecule under the optimal model (after spatial ordering). |
| `cost` | Float64 | Squared Euclidean distance $\|x_j - c_{a(j)}\|^2$ from this read to its assigned centroid. |
| `max_assignment_entropy` | Float64 | Read assignment entropy $H_j$ evaluated at the candidate $k$ with highest population uncertainty. Higher values indicate borderline/ambiguous assignments. |
| `sample` | String | Sample identifier derived from the source BAM/CRAM file name. |
| `chrom` | String | Chromosome name. |
| `start` | Integer | 1-based start coordinate of the read's overlap with the locus. |
| `stop` | Integer | 1-based stop coordinate of the read's overlap with the locus. |
| `PeakName` | String | Target peak identifier from `--peaks` (or `peak_1, peak_2, ...` if absent). |
| `Category` | String | Target peak category from `--peaks` (or `ClusterPeak` if absent). |
| `strand` | Boolean | Read alignment strand (`true` = positive/forward strand, `false` = negative/reverse strand). |
| `hap` | UInt8 | Phased haplotype assignment from the BAM tag (`0` = unassigned, `1` = Haplotype 1, `2` = Haplotype 2). |
| `readlength` | Integer | Overlapping length (bp) of the molecule across the evaluated locus: $\text{stop} - \text{start} + 1$. |
| `total_mod_6mA` | Integer | Number of 6mA modification calls observed on this molecule inside the locus. |
| `total_unmod_6mA` | Integer | Number of unmethylated adenine calls observed on this molecule inside the locus. |
| `rate_mod_6mA` | Float64 | Molecule 6mA modification rate: $\frac{\text{total\_mod\_6mA}}{\text{total\_mod\_6mA} + \text{total\_unmod\_6mA}}$. |
| `total_mod_5mC` | Integer | Number of 5mC CpG modification calls on this molecule inside the locus. |
| `total_unmod_5mC` | Integer | Number of unmethylated CpG calls on this molecule inside the locus. |
| `rate_mod_5mC` | Float64 | Molecule 5mC modification rate: $\frac{\text{total\_mod\_5mC}}{\text{total\_mod\_5mC} + \text{total\_unmod\_5mC}}$. |
| `total_mod_5hmC` | Integer | Number of 5hmC modification calls on this molecule inside the locus. |
| `total_unmod_5hmC` | Integer | Number of unmethylated calls corresponding to 5hmC positions. |
| `rate_mod_5hmC` | Float64 | Molecule 5hmC modification rate: $\frac{\text{total\_mod\_5hmC}}{\text{total\_mod\_5hmC} + \text{total\_unmod\_5hmC}}$. |
| `total_msp` | Integer | Total base pairs on this read covered by sub-nucleosomal Methyltransferase-Sensitive Patches (MSPs, $q \le 0.9 \times 255$). |
| `total_fire` | Integer | Total base pairs on this read covered by confident accessible FIRE elements ($q > 0.9 \times 255$). |
| `total_nuc` | Integer | Total base pairs on this read covered by nucleosome footprints. |
| `rate_msp` | Float64 | Fraction of read length occupied by MSPs: $\frac{\text{total\_msp}}{\text{readlength}}$. |
| `rate_fire` | Float64 | Fraction of read length occupied by confident FIRE elements: $\frac{\text{total\_fire}}{\text{readlength}}$. |
| `rate_nuc` | Float64 | Fraction of read length occupied by nucleosomes: $\frac{\text{total\_nuc}}{\text{readlength}}$. |
| `has_fire` | Boolean | True if the molecule contains at least one confident FIRE accessible element (`rate_fire > 0`). |
| `nrl` / `nrl_fire` | Int / Missing | Single-molecule FIRE-channel spectral Nucleosome Repeat Length (bp). `missing` if unphased. |
| `phasing_score` / `phasing_score_fire` | Float64 | Single-molecule FIRE-channel spectral power score. |
| `phasing_implies_closed` / `phasing_implies_closed_fire` | Boolean | True if single-molecule FIRE `phasing_score > 0.035`, indicating structured periodic chromatin. |
| `phasing_peak_type` / `phasing_peak_type_fire` | String/Symbol | Single-molecule FIRE spectral harmonic structure: `:mononucleosome`, `:dinucleosome`, or `:none`. |
| `nrl_6mA` | Int / Missing | Single-molecule `:mod_6mA`-channel spectral Nucleosome Repeat Length (bp). `missing` if unphased. |
| `phasing_score_6mA` | Float64 | Single-molecule `:mod_6mA`-channel spectral power score. |
| `phasing_implies_closed_6mA` | Boolean | True if single-molecule `:mod_6mA` `phasing_score_6mA > 0.035`. |
| `phasing_peak_type_6mA` | String/Symbol | Single-molecule `:mod_6mA` spectral harmonic structure: `:mononucleosome`, `:dinucleosome`, or `:none`. |

---

### 5.2 Cluster Table (`*.clustertable.tsv.gz`)

Contains one row per cluster discovered in the locus ($k^*$ rows total), ordered by `cluster` ($1 \dots k^*$).

| Column Name | Data Type | Description |
| :--- | :--- | :--- |
| `cluster` | Integer | Cluster identifier ($1 \dots k^*$), ordered spatially left-to-right by default. |
| `chrom` | String | Chromosome name. |
| `start` | Integer | Locus start coordinate (including flank). |
| `stop` | Integer | Locus stop coordinate (including flank). |
| `PeakName` | String | Target peak identifier from `--peaks` (or `peak_1, peak_2, ...` if absent). |
| `Category` | String | Target peak category from `--peaks` (or `ClusterPeak` if absent). |
| `totalreads` / `count` | Integer | Total number of reads assigned to this cluster. |
| `meancost` | Float64 | Mean squared distance of reads in this cluster to the centroid. Lower values reflect tighter, more homogeneous clusters. |
| `mean_max_assignment_entropy` | Float64 | Mean assignment entropy among reads in this cluster. Identifies whether this cluster contains distinct states or transitional molecules. |
| `counts_sample` | String | Distribution of reads across input samples, formatted as `sample1:count,sample2:count`. |
| `proportion_positive_strand` | Float64 | Fraction of reads in the cluster derived from the $(+)$ strand. Deviations from $\approx 0.5$ indicate technical strand bias. |
| `counts_hap` | String | Distribution of reads across haplotypes, formatted as `0:unassigned,1:hap1,2:hap2`. Crucial for identifying **allele-specific chromatin states**. |
| `mean_readlength` | Float64 | Mean overlapping length (bp) of reads belonging to this cluster. |
| `mean_rate_mod_6mA` | Float64 | Average 6mA accessibility rate across all reads in this cluster. |
| `mean_rate_mod_5mC` | Float64 | Average endogenous 5mC CpG methylation rate across all reads in this cluster. |
| `mean_rate_mod_5hmC` | Float64 | Average 5hmC rate across all reads in this cluster. |
| `mean_rate_nuc_bp` | Float64 | Average fraction of base pairs occupied by nucleosome footprints. |
| `mean_rate_msp_bp` | Float64 | Average fraction of base pairs occupied by linker/sub-nucleosomal patches. |
| `mean_rate_fire_bp` | Float64 | Average fraction of base pairs occupied by confident FIRE elements. High values indicate open/accessible chromatin states. |
| `proportion_has_fire` | Float64 | Proportion of reads in the cluster possessing at least one FIRE element. |
| `total_valid_nrls` | Integer | Number of reads in the cluster with an inferred single-molecule FIRE NRL. |
| `mean_read_nrl` | Float64 / Missing | Mean of individual read FIRE NRLs in this cluster. |
| `mean_read_phasing_score` | Float64 | Mean of individual read FIRE phasing scores in this cluster. |
| `total_valid_nrls_6mA` | Integer | Number of reads in the cluster with an inferred single-molecule `:mod_6mA` NRL. |
| `mean_read_nrl_6mA` | Float64 / Missing | Mean of individual read `:mod_6mA` NRLs in this cluster. |
| `mean_read_phasing_score_6mA` | Float64 | Mean of individual read `:mod_6mA` phasing scores in this cluster. |
| `max_fire_depth_rate` | Float64 | Peak 1 bp accessibility fraction across the cluster centroid. |
| `max_fire_depth` | Integer | Absolute estimated read count at the peak 1 bp accessibility summit within the cluster: $\text{round}(\text{count} \times \text{max\_fire\_depth\_rate})$. |
| `strand_or`, `strand_logor_se`, `strand_pvalue` | Float64 | Fisher's Exact Test odds ratio, log-OR standard error, and two-sided $p$-value for strand enrichment (`+` vs `-`) in this cluster. |
| `hap_or`, `hap_logor_se`, `hap_pvalue` | Float64 | Fisher's Exact Test odds ratio, log-OR standard error, and two-sided $p$-value for haplotype enrichment (`hap == 1` vs `hap == 2`) when $\ge 10$ phased reads span the locus. |
| `sample_or`, `sample_logor_se`, `sample_pvalue`, `sample_min_pvalue` | Float64 | Fisher's Exact Test odds ratio, log-OR standard error, and minimum two-sided $p$-value for sample enrichment in this cluster. |
| `nrl` / `nrl_fire` | Int / Missing | **Cluster-level FIRE Nucleosome Repeat Length (bp)**, computed from the autocorrelation of the `:fire` cluster centroid. |
| `phasing_score` / `phasing_score_fire` | Float64 | Cluster-level `:fire` spectral power score. |
| `is_closed` / `is_closed_fire` | Boolean | True if the `:fire` cluster centroid exhibits phased, closed chromatin structure (`phasing_score > 0.035`). |
| `peak_type` / `peak_type_fire` | String/Symbol | Dominant `:fire` spectral architecture for the cluster centroid: `:mononucleosome`, `:dinucleosome`, or `:none`. |
| `nrl_6mA` | Int / Missing | **Cluster-level 6mA Nucleosome Repeat Length (bp)**, computed from the autocorrelation of the `:mod_6mA` cluster centroid. |
| `phasing_score_6mA` | Float64 | Cluster-level `:mod_6mA` spectral power score. |
| `is_closed_6mA` | Boolean | True if the `:mod_6mA` cluster centroid exhibits phased chromatin structure (`phasing_score_6mA > 0.035`). |
| `peak_type_6mA` | String/Symbol | Dominant `:mod_6mA` spectral architecture for the cluster centroid: `:mononucleosome`, `:dinucleosome`, or `:none`. |

---

### 5.3 5mC Modification Table (`*.modclustertable.5mC.tsv.gz`)

Provides single-CpG resolution quantification across clusters. Only CpG sites with $> 5$ total calls (`min_mods=5`) across the locus are retained.

| Column Name | Data Type | Description |
| :--- | :--- | :--- |
| `chrom` | String | Chromosome name. |
| `start` | Integer | Locus start coordinate. |
| `stop` | Integer | Locus stop coordinate. |
| `cluster` | Integer | Cluster identifier ($1 \dots k^*$). |
| `pos` | Integer | Genomic 1-based coordinate of the cytosine (or unified CpG coordinate if `--unify_gc_strand` is enabled). |
| `mod` | Integer | Number of methylated reads observed at this position in this cluster. |
| `unmod` | Integer | Number of unmethylated reads observed at this position in this cluster. |
| `total` | Integer | Total coverage depth at this position in this cluster: $\text{mod} + \text{unmod}$. |
| `rate` | Float64 | Cluster-specific methylation rate: $\frac{\text{mod}}{\text{total}}$. |
| `total_mod` | Integer | Total methylated reads across **all** clusters combined at this position. |
| `total_unmod` | Integer | Total unmethylated reads across **all** clusters combined at this position. |
| `total_rate` | Float64 | Background locus-wide methylation rate across all reads: $\frac{\text{total\_mod}}{\text{total\_mod} + \text{total\_unmod}}$. |
| `log2FC` | Float64 | Regularized log2 fold-change of cluster methylation relative to locus background: $\log_2(\text{rate} + 0.05) - \log_2(\text{total\_rate} + 0.05)$. |
| `pvalue` | Float64 | Two-sided Binomial Test $p$-value comparing cluster `mod` counts against the background `total_rate`. Tests whether this position is differentially methylated in this specific cluster. |

---

### 5.4 Serialized Portable Data (`*.jld2`)

The `.jld2` file contains the complete execution state saved as a standard Julia `NamedTuple`. It can be inspected directly with `JLD2.load_object` or via `SMGStats.load_clustering`:

```julia
using JLD2, SMGStats

# Load the clustering payload (automatically backfills any missing NRL fields)
data = load_clustering("smgstats_cluster_...jld2")

# Fields available in data:
# data.chrom               -> Chromosome (String)
# data.peakloc             -> Peak interval range (UnitRange{Int})
# data.loc                 -> Window interval range including flank (UnitRange{Int})
# data.PeakName            -> Peak identifier string (e.g. "ctcf_1_fp_1")
# data.Category            -> Peak category string (e.g. "FirePeak", "AccessibleTSS")
# data.reads               -> Vector of read named tuples
# data.CM                  -> Stencil feature matrix struct
# data.optimal_clustering  -> Optimization result struct:
#                             - optimal_clustering.optimal_k
#                             - optimal_clustering.optimal_bic
#                             - optimal_clustering.ks (all evaluated k values)
#                             - optimal_clustering.bics (BIC values for each k)
#                             - optimal_clustering.optimal_model.centers (centroids matrix)
#                             - optimal_clustering.optimal_model.assignments (read cluster vector)
#                             - optimal_clustering.optimal_model.probabilities (K x N matrix)
# data.readtable           -> DataFrame of per-read statistics
# data.clustertable        -> DataFrame of cluster summary statistics
# data.modclustertable     -> DataFrame of CpG modification statistics

# Reorder clusters in-place across all tables and model matrices:
reorder_clusters!(data)                              # Default left-to-right spatial order
reorder_clusters!(data; by=:hap_pvalue)              # Order by haplotype skew significance
reorder_clusters!(data; by=:mean_rate_mod_5mC, rev=true) # Order by 5mC methylation (high -> low)
```

---

## 6. Generating Interactive HTML Reports & Master Catalog (`smgstats report`)

Once `smgstats cluster` has generated `.jld2` files, `smgstats report` builds a standalone interactive web catalog (`index.html`, `catalog.json`, and `guide.html`) along with detailed per-locus HTML diagnostic reports (`regions/<slug>/report.html`).

### 6.1 Report CLI Options

| Argument | Flag | Type | Default | Description |
| :--- | :--- | :--- | :--- | :--- |
| `--input` | `-i` | Path | Optional (if rebuilding UI) | Input directory containing `.jld2` files, or a single `.jld2` file. If omitted when `<output>/catalog.json` already exists, rebuilds `index.html` and `guide.html` instantaneously. |
| `--output` | `-o` | Directory Path | **Required** | Output directory for the master catalog dashboard (`index.html`) and `regions/` reports. |
| `--exondata` | `-e` | File Path | `nothing` | Optional GTF (`.gtf.gz`) or Arrow file for gene model tracks, nearest-gene annotation, and TSS distance calculation. |
| `--peaks` | `-p` | File Path | `nothing` | Optional target peaks BED/TSV file to supply or override `PeakName` and `Category` labels when reporting on older `.jld2` files. |
| `--samples` | `-s` | Comma-separated Strings | *Auto-inferred* | Optional sample label overrides when reporting on `.jld2` files (e.g. `050526=endo_as1,200126=endo_as2,181225=endo_wgs` or `endo_as1,endo_as2,endo_wgs`). Even without `--samples`, `load_clustering` automatically recovers Dorado/ONT sample IDs from the full BAM paths stored in `clustering.reads`. |
| `--refresh_samples` | | Flag | `false` | Fast partial update for existing region reports: regenerates **only** `assets/sample.svg` and `report.html` (skipping the other 10 figures per region) after sample label resolution or `--samples` overrides. |
| `--longest_only` | | Boolean (`true`/`false`) | `true` | Plot only the longest transcript isoform per gene on the gene model track. |
| `--workers` | `-w` | Integer | `1` | Number of parallel OS worker processes to spawn for rendering per-region CairoMakie figures and HTML reports. |
| `--chunk` | | Integer | `1` | 1-based worker chunk index ($1 \le \text{chunk} \le \text{nchunks}$) when splitting report generation externally across cluster jobs. |
| `--nchunks` | | Integer | `1` | Total number of chunks when splitting report generation externally. |
| `--no_reports` | | Flag | `false` | Build only the master catalog (`index.html`, `catalog.json`, `guide.html`) and skip rendering individual region reports. |
| `--no_index` | | Flag | `false` | Skip building the master catalog index and only render the individual region reports for the specified chunk. |
| `--force` | | Flag | `false` | Overwrite and regenerate existing `report.html` files (by default, existing region reports are skipped). |
| `--limit` | | Integer | `-1` | Limit the number of regions processed (useful for quick testing). |
| `--title` | | String | `"SMGStats Chromatin Cluster Catalog"` | Custom header title for the master catalog dashboard. |

---

### 6.2 Multi-Process Parallel Report Rendering (`--workers` & `--chunk`)

`smgstats report` uses a two-phase architecture designed to scale cleanly to thousands of clustered loci:
1. **Phase 1 (Instant Catalog Indexing)**: Extracts `LocusSummary` metadata across all `.jld2` files in parallel using threads (`-t`) and writes `index.html`, `catalog.json`, and `guide.html` **immediately** before plotting begins. You can open `index.html` right away while individual reports render.
2. **Phase 2 (Parallel Process Rendering)**: Because `CairoMakie` / `FreeType` font rendering is process-bound, passing `-w N` (`--workers N`) spawns $N$ independent OS worker processes, each rendering a disjoint slice (`w:N:n_total`) of region reports in parallel with zero lock contention.

```bash
# Single-command parallel report generation using 16 threads for indexing and 24 OS worker processes for plotting:
smgstats -t 16 report \
    -i ./cluster_results \
    -o ./cluster_catalog \
    -e gencode.v46.annotation.gtf.gz \
    -p peaks.tsv \
    --workers 24
```

For HPC / SLURM job arrays where you want to distribute report rendering across separate compute nodes:
```bash
# Step 1: Build master catalog index.html in seconds
smgstats -t 8 report -i ./cluster_results -o ./cluster_catalog -e gencode.gtf.gz --no_reports

# Step 2: Run array workers 1..24 across nodes (each renders only its deterministic slice i:24:n_total)
smgstats report -i ./cluster_results -o ./cluster_catalog -e gencode.gtf.gz --no_index --chunk ${SLURM_ARRAY_TASK_ID} --nchunks 24
```

---

## 7. Glossary of Key Terms

### Biological & Genomic Terms
- **Single-Molecule Genomics / SMG**: Sequencing technologies that read intact, native DNA molecules without PCR amplification, preserving long-range multi-kilobase spatial correlations.
- **Fiber-seq / SMF (Single-Molecule Footprinting)**: Assays employing non-specific DNA methyltransferases (e.g. Hia5 for 6mA, or M.CviPI for GpC) to mark accessible chromatin, leaving unmethylated "footprints" where proteins (transcription factors, nucleosomes) were bound.
- **FIRE (Fiber-seq Inferred Regulatory Element)**: Regions of dense, continuous adenine methylation indicating statistically significant, FDR-controlled chromatin accessibility.
- **MSP (Methyltransferase Sensitive Patch)**: Any contiguous stretch of methylated bases. Large MSPs correspond to regulatory elements, while short MSPs correspond to inter-nucleosomal linkers.
- **Nucleosome Footprint**: An unmethylated protection zone approximately $140\text{--}160$ bp in length, created by an intact nucleosome core particle obstructing methyltransferase enzymes.
- **NRL (Nucleosome Repeat Length)**: The average center-to-center distance (in bp) between adjacent nucleosomes along a chromatin fiber. Typically ranges between $160\text{ bp}$ and $210\text{ bp}$ depending on cell state and chromatin compaction.
- **Closed vs. Open Chromatin**: In `smgstats`, closed chromatin is characterized by low overall accessibility (`rate_fire ≈ 0`) and high periodic nucleosome phasing (`is_closed = true`). Open chromatin features prominent FIRE elements and disrupted nucleosome periodicity.
- **Haplotype (Allele-Specific Chromatin)**: Genetic variants (SNPs) phased onto parental chromosomes (Haplotype 1 vs 2). When a cluster is enriched for `hap:1` or `hap:2`, it indicates allele-specific accessibility or genomic imprinting.

### Algorithmic & Statistical Terms
- **Stencil / Feature Matrix**: The multi-channel numerical representation $M$ of single-molecule chromatin data, spatially binned at resolution $\Delta x$ and optionally smoothed with a Gaussian kernel.
- **K-Means Clustering**: Partitioning $N$ observations into $k$ clusters, where each observation belongs to the cluster with the nearest mean centroid.
- **Lloyd's Algorithm**: The standard iterative expectation-maximization heuristic for k-means: alternating between hard nearest-centroid assignments and recalculating centroid means.
- **Iterated Local Search (ILS)**: A metaheuristic optimization algorithm that explores multiple local minima by repeatedly perturbing (mutating) a known local minimum and re-applying local search.
- **Hill Climbing**: An optimization strategy that systematically accepts moves that improve the objective function while rejecting moves that degrade it.
- **Perturbation / Center Teleportation**: Relocating selected cluster centroids directly to randomly sampled data points, enabling the algorithm to jump across energy barriers between distinct clustering configurations.
- **BIC (Bayesian Information Criterion)**: A model selection criterion that balances goodness-of-fit against model complexity:
  $$\text{BIC} = k D \ln N - 2 \ln L$$
  Lower values denote more parsimonious models.
- **Assignment Entropy ($H_j$)**: A Shannon entropy metric quantifying the ambiguity of assigning read $j$ to a cluster. A read with $H_j = 0$ is unambiguously committed to one cluster; high entropy denotes transitional molecules.
- **Spectral Matched Filter**: A signal processing technique that computes cross-correlation between a signal and an analytic template (in this case, pure cosine waves across candidate NRL frequencies) to extract periodic signals without standard FFT binning constraints.

---

## 8. Recommended Parameter Presets

### A. High-Accuracy Production Run (Default)
Best for finalizing results across publication regions. Thoroughly explores the clustering landscape to guarantee identification of rare chromatin states.

```bash
smgstats -t 8 cluster \
    --peaks targets.bed \
    --output ./results_production \
    --dw 600 \
    --dx 15 \
    --features fire,mod_6mA \
    --ks 1:20,22:2:60 \
    --max_k_fraction 0.65 \
    --n_ils_restarts 2 \
    --n_restarts 20 \
    --rejection_neighbourhood_fraction 0.08 \
    --epochs 1000 \
    sample.bam
```

### B. Fast Exploration / Parameter Tuning Run
Best for testing hundreds to thousands of candidate loci quickly. Uses fewer restarts and an early rejection threshold.

```bash
smgstats -t 8 cluster \
    --peaks targets.bed \
    --output ./results_fast \
    --dw 500 \
    --dx 20 \
    --features fire,mod_6mA \
    --ks 2:2:20 \
    --n_ils_restarts 1 \
    --n_restarts 5 \
    --rejection_neighbourhood_fraction 0.02 \
    --epochs 200 \
    sample.bam
```

### C. Methylation-Focused Clustering (CpG Epigenetics)
Best when studying regions where endogenous DNA methylation (5mC/CpG) defines sub-states (e.g., imprinted regions, bivalent promoters, or methylated repeat elements).

```bash
smgstats -t 8 cluster \
    --peaks targets.bed \
    --output ./results_5mC \
    --dw 600 \
    --dx 15 \
    --features fire,mod_6mA,mod_5mC \
    --unify_gc_strand true \
    --ks 2:15 \
    sample.bam
```
