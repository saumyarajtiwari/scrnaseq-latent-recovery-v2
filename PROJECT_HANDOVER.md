# Project Handover — Mapping Latent Biological Subspace Recovery in scRNA-seq Preprocessing under Sparsity, Depth, and Dropout Stress

**Repository:** `saumyarajtiwari/scrnaseq-latent-recovery` (public, branch: `main`)
**Last updated:** 2026-08-30 (Step 6.9 complete — full failure-mode detection and characterization framework closed out)

---

## 1. Project Overview & Scientific Motivation

This is a latent-recovery study. It is deliberately not a broad benchmark
paper, and that distinction matters for how every later step should be
read. The question we're actually answering is: *which scRNA-seq
preprocessing methods preserve the true latent biological subspace, and
at what technical thresholds do they fail?*

The mechanics of the study are simple to state even though the execution
took months: raw count matrices get pushed through six different
preprocessing methods, and we check whether the resulting low-dimensional
embeddings still contain the real biological structure, or whether that
structure has been distorted, hidden, or replaced by technical noise. The
interesting part is not "which method wins" — it's mapping out *where*
each method starts to fail as sparsity goes up, depth goes down, and
dropout gets worse.

Three simulators are used — Splatter, scDesign3, SymSim — specifically so
no single simulator's quirks can drive the conclusions. Six real datasets
back this up. Six preprocessing methods are compared. Seven failure modes
are formally defined and tested for. All of that is described in detail
below, step by step, roughly in the order it actually happened.

---

## 2. Environment & Machine Setup

- **Machine:** Acer Predator Helios Neo 16 (Intel i5-13500HX, 20 threads, 16GB RAM)
- **OS:** Ubuntu 24.04.4 LTS, username `aayush`
- **R:** 4.3.3, with `renv` (environment locked via `renv.lock`)
- **Storage:** simulation and processed outputs live outside the git repo entirely — `/mnt/extra` and `/mnt/extra2` (internal NVMe partitions) plus `/mnt/archive` (an external USB SSD, ext4, must be mounted manually before touching `pca_libnorm`, `pca_log`, or `pca_shiftedlog`, which symlink into it). A `.renvignore` file excludes these large data directories from `renv`'s dependency scan.
- **Key packages:** Splatter, scDesign3, SymSim (`YosefLab/SymSim`), scran, scuttle, bluster, scater, BiocSingular, BiocParallel, TENxPBMCData, irlba, aricode, BiocNeighbors, RSpectra, glmpca, sctransform.
- **Workflow discipline, held to for the entire project:** one command at a time, output pasted back before the next is proposed; no assumptions without terminal evidence; `set.seed(42)` before every stochastic operation; atomic writes (temp file + rename) for anything that modifies existing files; commit after every completed substep; nothing invalid or superseded gets deleted outright — it gets preserved under a timestamped backup or an explicit `_INVALID` suffix instead.

---

## 3. Repository Structure

The repo carries two different code-organization conventions, from two
different phases of the project, and we've decided to leave that as-is
rather than retroactively rename everything — both halves contain
validated, already-relied-upon code, and moving files around now would
introduce risk for a purely cosmetic gain.

```
R/01_simulation/        — Steps 1.1-1.7: simulator calibration and generation scripts
R/02_real_data/         — Step 2: real-dataset download, QC, harmonization
R/03_eda_checkpoint2/   — EDA Checkpoint 2 scripts
scripts/                — Steps 3 through 6: flat step{N}_{substep}_{description}.R convention
docs/                   — structured findings documents, written for direct reuse in the manuscript
data/                   — real/, simulated/, processed/ — raw and processed .rds files excluded
                          from git (see .gitignore), living on external storage instead;
                          a handful of small, critical artifacts are tracked as explicit exceptions
results/tables/, results/figures/ — tracked CSV/PNG outputs from Step 1's EDA checkpoints
figures/                — Step 5 and EDA Checkpoint 3 rendered figures
logs/                   — tracked per-file progress logs for the largest batch runs
renv.lock               — pinned package versions
PROJECT_HANDOVER.md     — this document
PATCH_NOTES.md          — chronological log of fixes applied across the project
```

Files preserved under an `_INVALID` suffix (e.g.
`simulate_scdesign3_v1_INVALID_sparsity_inert.R`) are not dead code to be
cleaned up — they're kept on purpose, as documented evidence of bugs that
were found and fixed, in case anyone ever needs to see exactly what the
broken version did.

---

## 4. Parameter Grid & Design

`data/simulated/param_grid.csv` (10,940 rows) is the backbone of the
whole simulated side of the study. `run_id` is a sequential integer,
1 through 10,940, and it doubles as the seed-recovery mechanism — any
file, from any simulator, can be regenerated deterministically from its
`run_id` alone, given the matching script version.

| Axis | Levels |
|---|---|
| sparsity | 0.70, 0.80, 0.90, 0.95, 0.98 — ordinal severity ranks, not absolute targets |
| depth | 500, 2000, 10000 |
| dropout | none, low, high |
| separability | null, low, medium, high |
| n_cells | 200, 1000, 5000 (capped per simulator/separability where the reference pool is smaller) |
| batch | none, simple, complex |
| gene_strategy | all, hvg500, hvg2000 |
| clipping | none, clip99, log_stabilized |

`n_groups` is fixed at 5 for the main grid, with one exception: scDesign3,
where the separability tier itself determines how many real cell types
get pulled from the PBMC3k reference (null=1, low=2, medium=4, high=5).

Null-controls live in a separate manifest,
`data/simulated/null_control_grid.csv` — 45 rows in total: 15 of them
just point back into the main grid (`is_new=FALSE`, one per simulator per
sparsity level, already generated as part of the main run), and 30 are
genuinely new replicate files (`is_new=TRUE`, two additional replicates
per condition), added specifically so a negative claim ("no biological
structure here") wouldn't rest on a single stochastic draw.

---

## Step 1.1 — Environment & Tooling Setup

The very first task was getting all three simulators actually working
side by side, which turned out to be less trivial than it sounds. R 4.3.3
and `renv` were set up first, then the GitHub remote, then the three
simulator packages themselves: Splatter and scDesign3 installed cleanly,
but SymSim (`YosefLab/SymSim`) hit a real compatibility problem — a
conflict between its own `rank()` function and `BiocGenerics::rank()`
under the current Bioconductor version. A runtime patch didn't stick,
because SymSim's functions are byte-compiled, so the patch had to be
applied at the source level instead, before the package was built. Once
that was sorted out, all three simulators were confirmed working with a
basic smoke test each.

---

## Step 1.2 — Parameter Grid Construction & Calibration

Building the 10,940-row grid meant first figuring out, for each
simulator, which actual parameter genuinely controls sparsity — and this
is where the first real surprise of the project showed up. The obvious
candidate in Splatter, `mean.rate`, turned out to have essentially no
effect on sparsity at all, because Splatter internally renormalizes gene
means in a way that cancels it out. The parameter that actually does the
job is `bcv.common` (biological coefficient of variation), and that only
came out through direct calibration sweeps, not from reading the
documentation.

The second important design decision made here, which shaped a lot of
what came later: sparsity levels (0.70 through 0.98) are ordinal severity
labels, not literal targets the generator is supposed to hit exactly.
This was documented explicitly at this stage so nobody downstream would
mistake "sparsity=0.90" for a promise that the actual achieved sparsity
would land at 0.90 — it's a rank, not a set-point. That distinction ended
up mattering again in EDA Checkpoint 1, much later.

---

## Step 1.3 — Splatter Generation

The original build used `bcv.common` for sparsity (per the fix above) and
`lib.loc`/`lib.scale` for depth. All 10,940 files were generated and
verified against the target grid at the time.

### Corrections made during Step 1.7 validation

**Finding 1 — the dropout `low`/`high` labels were backwards.**
`param_dict.R` had mapped `"low"` to `dropout.mid=3.0` and `"high"` to
`dropout.mid=1.0`. A direct controlled sweep of `dropout.mid` from 0.5 to
5.0, holding everything else fixed, showed the relationship is
unambiguous: higher `dropout.mid` means *more* dropout and less depth
(sparsity 0.9706 at mid=1.0 vs. 0.9925 at mid=3.0). So every file labeled
`dropout="low"` actually had more real dropout than every file labeled
`"high"` — exactly inverted. This tracked perfectly with something
already visible in the validation numbers: depth deviation was 70.8% for
`"low"` and only 40.0% for `"high"`, which made no sense until the
inversion was found.

Since the underlying count matrices were valid simulations of specific,
real `dropout.mid` values — they just had the wrong label attached — the
fix was to relabel rather than regenerate. `param_grid.csv` was backed up
first, then the `dropout` values for the 7,290 affected rows (3,645 each
direction) were swapped, and every one of those files'
`run_params$dropout` field was patched to match (7,290 files, zero
errors, about half an hour). `param_dict.R` was corrected for future
runs. After the fix, the ordering came out right: `none` (0.74% mean
deviation) < `low` (40.0%) < `high` (70.8%).

**Finding 2 — the depth calibration table never accounted for dropout.**
`calibrate_splatter.R`'s own comments admit sparsity was calibrated only
at `dropout.type="none"`, with dropout "layering on top" — a known but
never-quantified coupling. That meant Step 1.7's depth check had been
comparing `dropout=low/high` files against a baseline that didn't apply
to them. A new calibration script,
`calibrate_splatter_depth_dropout.R`, extended the original `lib.loc`
sweep across all three dropout levels using the corrected `dropout.mid`
values, producing a proper 30-row lookup table. After switching the
validation script to use it, flagged files dropped from 62% down to
37.3%, and mean deviation fell from 37.2% to 17.3%.

What's left after both fixes is concentrated almost entirely at
`target_depth < 100`, where percentage-deviation is just naturally
unstable because the denominator is tiny — not a generation problem.

A small number of near-tied rank-order violations (563 out of 2,187
groups) turned out to be dropout-severity saturation at the fine end of
the sparsity gradient, with the worst single violation across the entire
dataset amounting to 0.22 percentage points. That got resolved with a
tolerance in the validation logic, not treated as a defect.

*(Side note, kept for the record: the script's header comment says
`N_WORKERS = 10`, but the configured value is actually `2L`. Harmless,
left alone.)*

---

## Step 1.4 — scDesign3 Generation

The original approach fit one NB-GLM model per unique combination of
depth, dropout, separability, n_cells, and batch — 244 fits in total —
using a PBMC 3k reference, and then generated many files off each fit. All
10,940 files were produced and looked fine at the time.

### The biggest single finding of Step 1.7

Neither `sparsity_label` nor `dropout` actually did anything. The
original script generated one count matrix per fit and saved it
unchanged across every sparsity label and every dropout value sharing
that fit — confirmed directly by comparing `sum(counts)` across 20
randomly sampled fit groups, all of which came back byte-identical. The
root cause: `dropout` fed into `family_use` (switching between negative
binomial and zero-inflated negative binomial), but the actual calibrated
zero-inflation values were never passed through to the generation call.
`sparsity_label` had no generation-time hook at all.

Fitting with `family_use="zinb"` was tried and dropped — it took about 29
minutes just for the fitting and parameter-extraction steps at the
smallest test size, which projected to roughly five days for the full
244-fit grid, and on top of that the extracted zero-inflation matrix came
back entirely `NA` on this package version, which crashed the downstream
simulation call outright.

What actually worked: keep `family_use` fixed at `"nb"` for every fit
(fast — about 13 minutes for the full generation at `n_cells=1000`,
matching the original timing), and implement both `dropout` and
`sparsity_label` as calibrated post-hoc stochastic zero-masking instead —
Bernoulli draws, seeded per row, composed as
`p_combined = 1 - (1-p_dropout)(1-p_sparsity)`. The sparsity ladder ended
up at `{0.7:0.00, 0.8:0.15, 0.9:0.35, 0.95:0.55, 0.98:0.75}`, dropout at
`{none:0.00, low:0.10, high:0.40}` (the pi values that had already been
defined but never used). Both extremes were tested directly — even the
harshest combination (dropout=high, sparsity=0.98, combined mask
probability 0.85) produced zero degenerate empty cells.

One nice side effect: since dropout no longer needed its own model fit,
the number of unique fits dropped from 244 to 82.

The old script was preserved as
`simulate_scdesign3_v1_INVALID_sparsity_inert.R`, and the old data was
backed up (6.0GB) rather than deleted. Full regeneration — 82 fits,
10,940 files, zero errors — took 1 hour 16 minutes, much faster than the
15-to-25-hour estimate that had been floating around, simply because
fewer fits were needed this time.

After the fix, every rank-order check passed cleanly (0 of 2,187 broken,
down from all of them), and depth deviation (14.2% flagged, mean 9.4%) is
fully explained by the masking removing real count mass — the validation
script now compares against a masking-adjusted baseline rather than the
raw pre-masking one.

A smaller bug turned up alongside this: the 30 new null-control files
were missing `gene_strategy` and `clipping` in their metadata — a
generation-script gap, not a data problem, and it was fixed at the
script level. The null-control generation script had the identical
sparsity-inert bug and got the identical masking fix; 10 files were
regenerated and confirmed monotonic across the sparsity ladder.

---

## Step 1.5 — SymSim Generation

The original design used a simulate-once-save-many pattern across 244
unique parameter combinations (separability, n_cells, batch, dropout,
depth), with `Sigma` fixed at 0.4 and 2,000 genes. Depth is severely
sublinear at low `alpha_mean` — the 10,000 UMI target simply isn't
reachable at high dropout, topping out around 3,967–4,000 regardless. All
10,940 files were generated and verified at the time.

### Correction made during Step 1.7

The original script's own comments already admitted `sparsity_label` had
no effect, since `Sigma` was fixed — this is the same underlying mistake
as scDesign3's bug: a parameter that seemed like it should control
sparsity didn't, and rather than searching for the real lever, the axis
got written off as uncontrollable. Dropout, by contrast, was already
confirmed working correctly and needed no change.

The fix was the same post-hoc masking mechanism used for scDesign3,
applied only to `sparsity_label` this time. A worst-case test (depth=500,
dropout=high, separability=high, smallest population size) confirmed zero
empty cells even under the harshest setting. Since sparsity was never
part of the fit key to begin with, the fit count stayed at 244.

The old script was preserved
(`simulate_symsim_v1_INVALID_sparsity_inert.R`), old data backed up
(22GB). Full regeneration completed with zero errors, though runtime
varied wildly by combination — some of the larger, low-dropout fits (5,000
cells) took over an hour and a half each, since `True2ObservedCounts` has
much less zero-inflation to short-circuit through when there's little
technical dropout to begin with.

The run got interrupted once, mid-way through (216 of 244 done) — the
terminal application itself closed, which killed the process despite
`nohup` (which only protects against the shell hanging up, not the whole
terminal app quitting). It resumed cleanly from the checkpoint with zero
data loss. Rank-order came back clean (0/2,187 broken) and depth
deviation barely moved (0.82%), which is reassuring — it means the
masking fix didn't disturb the part that was already correct.

The null-control generation script for SymSim had the same bug and got
the same fix. Its regeneration was accidentally missed during recovery
from the terminal interruption above, caught only when the validation
script threw file-not-found warnings, and run separately afterward.

---

## Step 1.6 — Null-Control Generation

45 files total: 15 reused directly from the main grid (5 per simulator,
flagged `is_null_control=TRUE`), plus 30 genuinely new replicate files
(10 per simulator). The reasoning for the extra replicates is worth
restating since it comes up again much later, in Step 6: a claim of "no
real structure here" that rests on one random draw is weaker than one
confirmed across independent seeds, and a reviewer could reasonably
attribute a lack-of-structure finding to a lucky or unlucky draw
otherwise.

All three simulators' null-control generation scripts turned out to have
bugs, found during Step 1.7: scDesign3 had the sparsity-inert bug plus
the missing metadata fields; Splatter was missing only the metadata
fields (dropout is fixed at `"none"` for null-controls by design, so the
label-inversion bug never touched it); SymSim had the sparsity-inert bug
only, with metadata already correct. All were patched.

---

## Step 1.7 — Output Validation and Inventory

This was the step where most of the above findings actually surfaced.
The validation script checked all 32,865 files (32,820 main-grid plus 45
null-control): sparsity checked descriptively against raw counts plus
monotonic rank-order per simulator (aggregated across replicates first,
to avoid flagging ordinary sampling noise); depth checked against each
simulator's own calibration table with a 20% deviation threshold, chosen
to absorb legitimate variance without missing real failures.

Two infrastructure problems had to be cleared before the real findings
could even be seen. `mclapply` turned out to be unstable when forked
alongside the `Matrix`/`cholmod` C library — the validation logic worked
fine in serial and crashed under forking, so the script just ran serially
instead (three minutes for the whole thing, so parallelism wasn't buying
much anyway). Separately, `data.frame()` construction crashed on `NULL`
`gene_strategy`/`clipping` fields from the 30 new null-control files —
which is exactly what led to discovering that metadata gap.

The actual discovery order, roughly: the scDesign3 sparsity bug came
first and was the biggest single finding; while re-validating that fix,
the null-control metadata gaps in Splatter and SymSim turned up; then the
SymSim sparsity bug, spotted because it showed the exact same
100%-rank-order-broken signature as scDesign3's had; then the missed
SymSim null-control regeneration; then the Splatter dropout-label
inversion, found while chasing persistently high depth deviation; then
the depth-calibration dropout-blindness, found immediately after; and
finally the rank-order tolerance issue, resolved once across all three
simulators with a single empirically-derived tolerance value.

Final validated state:

| Simulator | Files | Flagged | Mean depth deviation | Rank-order broken |
|---|---|---|---|---|
| Splatter | 10,935 | 4,077 (37.3%) | 17.3% | 0/2,187 |
| scDesign3 | 10,935 | 1,557 (14.2%) | 9.4% | 0/2,187 |
| SymSim | 10,935 | 90 (0.8%) | 2.9% | 0/2,187 |
| Null-control | 60 | 0 | — | — |

Every remaining flag has an explanation behind it — real depth reduction
from masking or dropout, small-sample noise at the extreme low-depth
corner — nothing here is an unexplained anomaly. Zero unreadable files,
zero cell-count or group-count mismatches anywhere.

---

## Step 1.8 — Ground-Truth Extraction

Step 1.8's original brief called for storing ground-truth cell labels,
loading vectors, and a true subspace basis alongside each simulated
dataset once everything moved into a unified SingleCellExperiment format.
None of the Step 1.3–1.7 production files had this — the original
lightweight-save decision (counts, cell metadata, run params only)
deliberately deferred it until it was actually needed.

The definition settled on: the true subspace is the linear span of
noise-free group-mean expression vectors, a (G−1)-dimensional signal
subspace per file — standard in the signal-subspace-recovery literature.
Where that ground truth actually comes from differs by simulator, since
each one exposes a different native pre-noise layer: Splatter's
`TrueCounts` assay, SymSim's `SimulateTrueCounts()` output (before
`True2ObservedCounts()` adds technical noise), and for scDesign3 — which
doesn't expose its fitted model directly — the real PBMC3k reference
group means, using the same subsampling and depth-scaling seed as
production. Each extraction was validated on its own terms: Splatter's
against known differentially-expressed genes; SymSim's by confirming
label recovery beat random permutations; scDesign3's against known PBMC
marker genes localizing to the right cell types.

**Ground-truth-linearity caveat.** This ground truth is, by construction,
a *linear* subspace (the span of group-mean vectors). Grassmannian
distance, subspace recovery score, and spectral recovery score are
computed against it for all six preprocessing methods, including the two
whose own transforms are nonlinear (SCTransform v2, GLM-PCA). For those
two methods, recovery scores partly reflect representational mismatch
between a nonlinear embedding and a linear reference, not purely how much
biological signal survived preprocessing. This is the same reason Step
6.5c excludes both methods from its ground-truth-ratio test outright.
Rather than exclude them here too, Step 5.8's threshold table (and any
script consuming it) now carries an explicit `method_family` column and
a `ground_truth_caveat` field for nonlinear-method rows, so cross-method
comparisons of these scores can be reported stratified (linear vs.
nonlinear) rather than silently pooled. This should be stated explicitly
in the manuscript wherever simulated-data recovery scores are compared
across methods.

**Step 5.8 vs. Step 6.8 resolution (added post-review).** Three
(axis, simulator, method) combinations initially appeared to show Step
5.8's geometric recovery boundary ("never fails") contradicting three
independent Step 6 failure-mode tests ("fails 70-93% of the time"). Root
cause, confirmed with direct file-level evidence: Step 5.8's per-axis
thresholds are computed from a one-axis-at-a-time (OFAT) sweep held at a
single baseline point in the 8-dimensional grid, while Step 6.8's
cross-reference averages flag rates across the full marginal grid (e.g.
3,643 rows vs. 1 baseline row, for one such combination). These are two
different statistical claims about the same parameter space, not two
measurements of the same thing -- disagreement between a point estimate
and a marginal average is expected, not contradictory. Full writeup in
docs/step6_9_failure_mode_review.md.

**Gene-panel-size confound (Item 2), tested and found immaterial.** SCTransform
v2 and GLM-PCA cap at 3,000 genes (top mean-expression + all-zero-cell
rescue) while the four linear methods use the full gene set -- raising the
concern that recovery-score differences partly reflect gene-panel size,
not method quality. Tested directly: reran all four linear methods
restricted to the identical 3,000-gene cap, on a 135-file reduced grid
(sparsity x depth x separability, 3 simulators, baseline elsewhere),
comparing same-genes vs. full-genes subspace_recovery_score for each of
540 (file, method) combinations, all succeeding with no errors.

Result: scDesign3 and SymSim's native gene counts (2,000 each, fixed
since Step 1.4/1.5) fall below the 3,000-gene cap, so it never triggers
for them -- the confound is structurally moot there, not merely untested.
Splatter (10,000 genes natively) is the only simulator where the cap
actually applies, and is therefore the only valid test of this concern;
n=180 (45 files x 4 methods). There, the same-genes vs. full-genes
recovery-score delta has median 0.001-0.004 and mean 0.002-0.009 across
the four methods (max 0.06 for one outlier case) -- two orders of
magnitude smaller than the linear-vs-nonlinear recovery-score gaps of
0.3-0.8 documented in the ground-truth-linearity caveat above. Gene-panel
size is confirmed real as a structural difference between methods but
not a material driver of this study's observed recovery-score
differences. Full data in
results/step4_metrics/step4_item2_gene_panel_sensitivity.csv, analysis
script at code/04_metrics/step4_item2_gene_panel_sensitivity.R.

**Multiple-comparisons null calibration (Item 6), Technical Separation
(Step 6.1).** The literal spec threshold (AMI>0.5 against batch_id or
UMI-quartile) had no check against what pure chance would produce.
Calibrated via label-permutation nulls (200 draws, aricode::AMI),
stratified by (n_cells, k_used) for UMI-quartile and additionally by
batch cardinality for the batch covariate -- 16 UMI buckets, 25 batch
buckets (20 simulated + 5 real), covering all ~196,830 rows. Stability
confirmed: 200 vs. 2,000 draws give 95th-percentile estimates within
~15% of each other, same order of magnitude.

Result, in the opposite direction from what a "how much of this is
noise" framing anticipates: the calibrated 95th-percentile null
thresholds are tiny (0.00003 at n=61,292 cells to 0.03 at n=200 cells),
because AMI's null distribution under permutation concentrates sharply
around zero and tightens further as n grows -- standard behavior for
information-theoretic association measures at this scale. Consequently
the calibrated flag count (155,711 of 196,830) is six times LARGER than
the literal AMI>0.5 flag count (25,192), not smaller. This means the
literal 0.5 threshold is not really a statistical-significance
criterion -- it is a de facto effect-size cutoff. Essentially none of
the 25,192 literal flags are attributable to chance; they are a
conservative, high-effect-size subset of a much larger set of
statistically real (if often practically small) batch/UMI associations.
Both the literal flag (`technical_separation_flag`) and the new
calibrated flag (`technical_separation_flag_relative`) are retained
side by side, per Step 6.3's established precedent -- neither replaces
the other. Full calibration tables in
data/processed/step6_1_umi_ami_null_calibration.csv and
data/processed/step6_1_batch_ami_null_calibration.csv; per-row output in
data/processed/step6_1_technical_separation_calibrated.csv; script at
code/06_failure_modes/step6_1_null_calibration.R.

**GLM-PCA chunked-projection flagging (Item 5).** GLM-PCA's chunked
fixed-loadings projection (used for pbmc68k and tabula_sapiens_lung,
since native GLM-PCA fitting exceeds available memory at that scale) was
validated only on Baron -- a dataset that never needed chunking in
production, showing a real ARI degradation (0.652 native vs. 0.421
projected). Left unflagged, these two datasets' GLM-PCA rows could be
silently pooled with natively-fit GLM-PCA results in aggregate tables or
figures. Fixed at the source: a full manifest rebuild
(embedding_manifest.csv, 197,100 files) adds `internal_method` (each
.rds file's own "glmpca_chunked_projection" vs. "glmpca_nb" tag, verbatim)
and `has_validation_note` (boolean), propagated through
step4_9_compile_results.R's explicit column selection into
step4_master_results_table.csv. Exactly 2 of 196,830 rows carry the
chunked-projection flag, matching the known scope precisely. Any
manuscript table or figure aggregating GLM-PCA results should filter or
visibly mark on `internal_method`, the same way `method_family` handles
the linear/nonlinear distinction from Item 1.

**Multiple-comparisons null calibration (Item 6), Cluster Collapse (Step
6.2).** Same motivation as Step 6.1's calibration, but silhouette is a
geometry-dependent statistic, not a label-agreement one, so it is
calibrated via Step 6.3's noise-generation template rather than Step
6.1's label-permutation approach: random Gaussian embeddings at the real
(n_labeled, n_groups_used, nv_used) structure, roughly-equal synthetic
group sizes (same simplifying assumption as Step 6.1; real group sizes
are often unequal, a disclosed limitation not expected to change the
qualitative picture), 200 draws, using the identical simplified_
silhouette() function from step4_4to6_secondary_metrics.R. 29 buckets
covering all ~196,830 rows.

Result, in the OPPOSITE direction from Step 6.1's finding: the
calibrated flag (silhouette below the bucket's null 5th percentile) is
52,293 -- roughly one-third of the literal silhouette<0.2 flag count
(163,844). Null-distribution medians cluster tightly around zero (-0.018
to +0.010 across buckets) while real data's median silhouette (0.022)
sits meaningfully above most of them. This means roughly 111,551 of the
163,844 literally-flagged rows show silhouette statistically
distinguishable from pure noise -- real, detectable separation, just at
a numerically low value, consistent with silhouette's known tendency to
compress toward zero in higher-dimensional embeddings even for
genuinely-separated clusters (a curse-of-dimensionality effect, not
evidence of absent structure).

**This refines, but does not reverse or contradict, docs/step6_9_
failure_mode_review.md's "RETAINED, strongly supported" verdict for
Cluster Collapse**, nor its 52/52 Step 5.8 alignment finding (which
concerns correlation between independently-derived signals, unaffected
by whether underlying silhouette values are modest-but-real or robustly
separated). The precise, calibration-informed claim: of rows flagged
under the literal threshold, roughly two-thirds show statistically real
(if numerically modest) separation, and roughly one-third are genuinely
indistinguishable from random noise. Both flags
(`silhouette_trigger` literal, `silhouette_trigger_relative` calibrated)
retained side by side per established precedent. Full calibration table
in data/processed/step6_2_silhouette_null_calibration.csv; per-row
output in data/processed/step6_2_cluster_collapse_calibrated.csv; script
at code/06_failure_modes/step6_2_null_calibration.R.

**Multiple-comparisons null calibration (Item 6), Variance Hijacking
(Step 6.4).** This test already carried two threshold variants (literal:
Spearman-only for both covariates; relative: Spearman for UMI, eta for
batch, per the ordinality methodological note) -- neither was
null-calibrated; both still used the fixed 0.7 spec threshold. Added a
third, genuinely calibrated variant. Both statistics' null distributions
depend only on (n_cells) for Spearman and (n_cells, n_categories) for
eta -- not on embedding/group structure -- calibrated via permutation
(200 draws), critically taking max(|.|) across the top 3 PCs in the null
exactly as the flag logic does, so the calibration itself absorbs the
top-3 multiple-comparisons structure. 12 Spearman buckets, 17 eta buckets
(12 simulated + 5 real).

Result, same direction as Step 6.1 (not Step 6.2): calibrated null 95th
percentiles are far below 0.7 (0.01 at n=61,292-65,690 to ~0.26 at
n=200), shrinking with n exactly as AMI's did. Calibrated flag (157,010,
79.8%) is more than double the literal flag (70,743, 35.9%) and the
methodology-relative flag (67,341, 34.2%). Same interpretation as Step
6.1: the literal 0.7 threshold is an effect-size cutoff, not a
statistical-significance one -- most of the additional calibrated flags
represent real, detectable (if modest) technical association that the
spec's high bar doesn't capture, not noise inflating the literal count.
All three flags (`hijack_flag_literal`, `hijack_flag_relative`,
`hijack_flag_calibrated`) retained side by side. Full calibration tables
in data/processed/step6_4_spearman_umi_null_calibration.csv and
data/processed/step6_4_eta_batch_null_calibration.csv; per-row output in
data/processed/step6_4_variance_hijacking_calibrated.csv; script at
code/06_failure_modes/step6_4_null_calibration.R.

**Multiple-comparisons null calibration (Item 6), Neighborhood Collapse
(Step 6.6).** Mean k-NN overlap < 0.5 flags collapse; unlike prior
calibrations, this statistic's null distribution under "no relationship
between reference and method embeddings" is purely combinatorial (mean
intersection size of two independent random k-subsets from n-1 items),
depending only on (n_cells, k=15) -- no matrix/embedding simulation
needed. 11 buckets (a 12th, Tabula Sapiens Lung at n=61,292, is entirely
absent from this test's own output -- see the addendum added to
docs/step6_9_failure_mode_review.md, an undisclosed pre-existing gap
found during this work, not something we introduced).

Result, same direction as Step 6.2: null 95th percentiles are tiny
(0.0013 at n=65,690 to 0.093 at n=200) against a real median overlap of
0.208 -- comfortably above chance almost everywhere. Calibrated flag
(indistinguishable from/below chance: 8,741, 4.4%) is a small fraction of
the literal overlap<0.5 flag (183,900, 93.5%). This gives docs/step6_9's
existing qualitative caveat ("should not be read as neighborhood
structure is destroyed almost everywhere... partly reflects the
strictness of a binary criterion") a precise number: roughly 95.6% of
literally-flagged rows show real, statistically-detectable neighborhood
structure, just below an aggressive 0.5 bar; only about 4.4% are
genuinely indistinguishable from random chance. Both flags
(`overlap_flag` literal, `overlap_flag_calibrated`) retained side by
side. Full calibration table in
data/processed/step6_6_overlap_null_calibration.csv; per-row output in
data/processed/step6_6_neighborhood_collapse_calibrated.csv (196,824
rows, excludes the missing Tabula Sapiens Lung rows noted above); script
at code/06_failure_modes/step6_6_null_calibration.R.

**Multiple-comparisons null calibration (Item 6), Subspace Rotation
Slippage (Step 6.7).** Unlike the prior four tests, this one has an exact
closed-form null rather than requiring simulation: for a uniformly random
unit vector in R^p (ambient gene space), cos^2(angle to a fixed rank-r
true subspace) follows a Beta(r/2, (p-r)/2) distribution exactly --
standard directional-statistics result. r comes directly from each row's
own rank_true (already computed by Step 6.7; notably 1 for all scDesign3
rows, reflecting the true group-means matrix's actual numerical rank, not
an assumed n_groups-1=4). p (ambient gene dimension) is deterministic by
source for the 4 linear methods (2,000/10,000, no file reads needed) but
varies substantially per file for SCTransform v2/GLM-PCA -- found, while
building this calibration, to range from 121 to 3,000 genes (median
2,000), well below the nominal 3,000-gene cap discussed in Item 2,
apparently due to additional zero-count-gene removal at high-sparsity
parameter combinations. Read directly from all 65,598 affected files
(parallelized, ~2.6 min) rather than approximated.

Result: null 5th-percentile angle thresholds are enormous (73.9-88.2
degrees, median 86.3) -- expected, since a random vector's angle to a
low-rank subspace concentrates near 90 degrees in high-dimensional gene
space. Calibrated flag (angle statistically indistinguishable from
chance: 14,285, 7.3%) is a small fraction of the literal angle>30 flag
(172,934, 87.9%). This means roughly 92% of literally-flagged
"slippage" rows show real, statistically-significant alignment with the
true subspace -- the 30-degree threshold is a strict effect-size cutoff
in an extremely high-dimensional space, not a significance boundary.
Same qualitative pattern as Steps 6.2 and 6.6, now with an exact rather
than simulated null. Both flags (`flag_pc1` literal,
`flag_pc1_calibrated`) retained side by side. Per-row output (including
p_ambient and the per-row null threshold) in
data/processed/step6_7_subspace_rotation_slippage_calibrated.csv; script
at code/06_failure_modes/step6_7_null_calibration.R.

**Multiple-comparisons null calibration (Item 6), Over-Smoothing --
real-data DPT branch (Step 6.5a).** Only 18 rows total (3 pancreas
datasets x 6 methods), each calibrated exactly per-row (matching that
row's actual progenitor/mature cell-count split, not just total n --
correlation with a binary reference's null variance depends on class
balance) rather than bucketed, at 2,000 draws given the trivial
computational cost at this scale.

Result: the most decisive of all six calibrations. Null 95th percentiles
are tiny (0.019-0.039) against real Spearman correlations of 0.18-0.70 --
even the single weakest real correlation (0.177, GLM-PCA/Baron) sits
~9x above its own null threshold. Calibrated flag: 0/18. Every one of
the 10 literally-flagged (rho<0.6) rows reflects overwhelmingly real,
non-random biological signal between diffusion pseudotime and the
progenitor/mature reference -- none are remotely close to chance. This
confirms and sharpens docs/step6_9's existing framing: the 0.6 threshold
is a demanding effect-size bar requiring strong trajectory preservation,
not a test of whether any real structure survived preprocessing at all.
Both flags (`over_smoothing_flag` literal,
`over_smoothing_flag_calibrated`) retained side by side. Per-row output
in data/processed/step6_5a_pancreas_dpt_calibrated.csv; script at
code/06_failure_modes/step6_5a_null_calibration.R.

**Multiple-comparisons null calibration (Item 6), Over-Smoothing --
simulated variance-ratio branch (Step 6.5c).** This branch found 0/31,640
flagged (ratio<0.20) under the literal spec -- the relevant calibration
question here inverts from the other tests: not "how many flags are
noise" (there are none), but whether the 0.20 threshold is so lenient
the test could never detect genuine over-smoothing even if present, i.e.
whether the null result reflects a real absence of the failure mode or
an unfalsifiable test. Tested via label-permutation (200 draws): keep
the real gene-expression matrix and real embedding fixed, permute
true_group, recompute both pre-PCA and post-PCA signal fractions and
their ratio. Full 31,640-row grid was too expensive to permute in full;
calibrated on a stratified sample (25 files x 4 methods = 100 cases)
spanning all 3 simulators and both target log2FC bands (0.25, 0.50),
matching Item 3's cross-simulator triangulation standard rather than a
single-simulator sample (Splatter has zero eligible files at the 0.50
band -- a genuine scope gap, not a sampling bug, consistent with this
branch's already-narrow documented scope).

Result: null ratios cluster tightly near 1.0 (median 0.98-0.99 across
methods), with the lowest 5th percentile across any method still at 0.75
(Log-PCA) -- far above the 0.20 flag threshold. This confirms the test
has genuine statistical power: were real over-smoothing present (ratio
dropping toward the ~0.75-0.99 range that even a total ABSENCE of signal
produces), the test would detect it. The 0/31,640 result is therefore a
meaningful, well-powered negative finding, not an artifact of an
unfalsifiable threshold. (Null ratios cluster near 1 rather than near 0
because signal_fraction is an R^2-type statistic computed identically
pre- and post-PCA -- even under permuted labels both numerator and
denominator land near the same small random-label baseline, keeping
their ratio near 1 regardless of dimensionality change -- a reassuring
structural property of the ratio construction itself.) Full calibration
table in data/processed/step6_5c_variance_ratio_null_calibration.csv;
script at code/06_failure_modes/step6_5c_null_calibration.R.

scDesign3 and SymSim extracted cleanly — 82/82 and 244/244 fit-keys, zero
failures. Splatter needed one call per row rather than per fit-key, since
its seeding is unique per row, and turned up two separate, real data-
quality findings along the way, both caught because the extraction
script checks every single row rather than sampling.

The first: four files (run_ids 4, 5, 9, 3289) traced back to the very
first invocation of the Splatter generation script, run under
BiocParallel before it was swapped for `mclapply` the same day —
BiocParallel manages its own per-worker RNG substreams, which plausibly
interfered with the seed. These four were regenerated under the current
logic and the fix was verified with a third independent simulation call.

The second is smaller and was left open: six files (run_ids 3284, 3285,
6565, 6569, 9849, 9850) — all sharing sparsity in {0.95, 0.98} and
n_cells=200, though those two conditions together aren't sufficient on
their own to explain it (1,458 rows share both and only six fail) —
couldn't be reproduced from their own recorded parameters. Every
plausible explanation was checked and ruled out directly: metadata
inconsistency, dropout value, file corruption, `gene_strategy` as an
actual causal input (it isn't — Splatter never sees it), RNG carryover
from a prior session, and silent parameter drift in either
`param_dict.R` or the generation script itself (checked against git
history, unchanged since the first commit). These six were simply
excluded from the ground-truth set rather than shipped with a mismatched
label — logged, not silently dropped. In total, 10,929 files were
written and 6 excluded, accounting for all 10,935.

---

## EDA Checkpoint 1 (1.1–1.4)

**Sparsity verification** confirmed what Step 1.2 had already flagged:
`sparsity_label` is an ordinal rank, not a literal target, across all
three simulators. 68.4% of rows deviate from the literal label by more
than 0.03, which is entirely expected under that design and not a
problem.

**Depth verification** showed Splatter and SymSim calibrated accurately
at baseline. scDesign3 showed a genuine, isolated roughly 18%
under-delivery specifically at `depth_label="500"`, not explained by
masking — flagged and deliberately deferred rather than chased down
immediately, since it wasn't blocking anything downstream.

**DE-magnitude vs. separability** computed pairwise |log2FC| between
true group means per file. SymSim showed the cleanest separation between
adjacent separability tiers; scDesign3 and Splatter both showed
substantial overlap between tiers.

**The summary pass/fail table** used peer-group z-scoring (within
matched simulator × depth × sparsity × dropout groups) rather than a
literal-label tolerance, since the latter would have falsely flagged 68%
of the grid given the ordinal-label finding above. Separability showed a
fairly high rank-inconsistency rate (~23%, concentrated in scDesign3 and
Splatter), which was investigated and traced to something specific: no
calibration script for any simulator had ever empirically validated
separability against actual DE-magnitude output, unlike sparsity and
depth, which both went through proper calibration sweeps. Checking back
against the original project design document confirmed this was
intentional — separability was always meant to be characterized properly
later, via the subspace-recovery metrics this whole project exists to
compute, not fixed at the EDA stage with a rough log-fold-change proxy.
No recalibration was done; it was logged as a Step 3–4 watch item
instead.

---

## Step 2 — Real-Data Curation

Six datasets curated, each QC'd with MAD-based thresholds:

| Dataset | Cells | Genes | Cell types | Batches | Fully labeled |
|---|---|---|---|---|---|
| PBMC 68k | 65,690 | 20,387 | 5 | 1 | Yes, after Step 4.7's re-annotation |
| Muraro | 2,403 | 18,197 | 10 | 4 | No — 277 unlabeled |
| Baron | 8,569 | 17,499 | 14 | 4 | Yes |
| Segerstolpe | 2,929 | 23,238 | 14 | 10 | No — 720 unlabeled |
| Tabula Sapiens Lung | 61,292 | 56,139 | 34 | 4 | Yes |
| Tasic 2018 | 20,562 | 42,865 | 24 | 341 | Yes |

These weren't chosen to rank methods on real data in the abstract — the
point was to pick datasets that would independently stress the same kind
of failure modes the simulated grid was built to probe, and see whether
the simulated regimes actually show up in real biology too.

---

## EDA Checkpoint 2 — Real-Dataset Technical Characterization

Phase-space coverage between the simulated grid and the real datasets was
confirmed adequate. A more open-ended investigation went into whether
batch structure was visible in the pancreas datasets (Muraro, Baron,
Segerstolpe) — a series of escalating naive PCA analyses couldn't
visually confirm it. Both obvious explanations were checked and ruled
out (batch confounded with true cell type; batch simply having no effect
at all), and the investigation was deliberately stopped there rather than
pushed further, with a note that it should be re-checked once Step 3's
proper normalization pipeline existed. That re-check ended up happening
much later, formally, as part of Step 6.4 — see below.

---

## Step 3 — Preprocessing & Dimensionality Reduction (3.1–3.10) + EDA Checkpoint 3

Six preprocessing methods were implemented: Raw PCA, Library-Normalized
PCA, Log-PCA, Shifted-Log-PCA, a Pearson-residual (SCTransform v2 style)
PCA, and GLM-PCA. `irlba` was used throughout for memory-safe truncated
PCA, seeded consistently, with `nv=50` for Raw PCA and `nv=30` for
everything else.

GLM-PCA couldn't be fit natively at PBMC68k's or Tabula Sapiens Lung's
scale, so a chunked, fixed-loadings projection approach was built
instead. A GPU/PyTorch reimplementation was also attempted, in the hope
of running full GLM-PCA at scale, but it was abandoned after producing
roughly 49 degrees of subspace misalignment against the CPU-fitted
reference — too far off to trust.

Two real data-corruption incidents happened during this step, both
eventually fully resolved. The bigger one: 13,091 files corrupted during
the loadings backfill for the four linear methods, traced to `/mnt/extra`
filling to 100% capacity mid-write combined with an in-place-overwrite
pattern that left torn files behind when the disk ran out. Everything was
regenerated and independently re-verified by reading back all 32,814
raw-PCA files. A later full-grid integrity scan (Step 3.10) caught two
more files from this same incident that had been missed the first time
around — both real-data files, both fixed. That same integrity scan also
turned up and fixed a smaller, unrelated Raw PCA bug on real data (wrong
`nv` value, missing loadings).

SCTransform v2 needed two separate rescue passes to get everything
working; after both, 16 files still had one genuinely zero-total-count
cell that had to be explicitly excluded from that file's input, with the
exclusion logged directly in the file's own metadata rather than silently
dropped.

**EDA Checkpoint 3** did a visual review of representative embeddings —
best case, worst case, and null-control, across every method and
simulator. One thing flagged here, explicitly as unconfirmed: the
null-control panels for SCTransform v2 and GLM-PCA looked visually
non-unimodal, which was noted as a candidate signal worth a formal test
later, not treated as an established finding at the time. That formal
test came in Step 6.3, and it told a more nuanced story than the visual
flag suggested — more on that below.

---

## Step 4 — Subspace-Recovery Metrics (4.1–4.9)

The primary metrics — Grassmannian distance, principal angles, Subspace
Recovery Score, Spectral Recovery Score — were all computed against each
simulator's per-file ground truth from Step 1.8. Principal angles came
from an SVD of the cross-product between the estimated and true loading
subspaces, with the true subspace's rank determined numerically via an
SVD tolerance rather than assumed to always be one less than the number
of groups.

Secondary metrics — trustworthiness, continuity, ARI, silhouette — needed
a real scalability compromise, since exact versions of these metrics are
quadratic in the number of cells and some of the real datasets have tens
of thousands of them. Trustworthiness and continuity were computed
against a common reference space per file (log-normalized counts, top
500 highly variable genes, reduced to 15 PCs — built once per file, not
once per method), using k=15 neighbors and a capped-rank approximation
for tractability. Silhouette used a centroid-based simplified formulation
from the clustering-validity literature rather than the full pairwise
version. ARI and silhouette both came from k-means with `nstart=25` and
seed 42 — a detail that turned out to matter a great deal later, when
Step 6.1 had to reverse-engineer this exact configuration from the stored
results.

PBMC68k needed real intervention here: it turned out to have zero usable
ground-truth cell-type labels at all going into this step. A full
re-annotation pipeline — Louvain clustering plus marker-based scoring,
with a CD3D-gating fix for a CD8⁺ T-cell/NK marker overlap — resolved it,
landing at 65,690 out of 65,690 cells labeled across five canonical
types. The fix replaced the production file and updated the dataset
inventory, and `embedding_manifest.csv` was re-joined to reflect it — but
that re-join only touched the manifest's own metadata columns, not the
individual embedding files' internally-stored labels. That distinction
turned out to matter quite a bit in Step 6.1, below.

---

## Step 5 — Phase-Space & Failure-Boundary Analysis (5.1–5.9)

Six single-axis recovery plots were produced — sparsity, depth,
separability, batch, dropout, clipping — each with its own round of
visual review that caught real bugs. A scDesign3 dropout-label swap was
found and fixed here. Two other findings were investigated thoroughly and
then deliberately left as documented limitations rather than fixed:
clipping turned out to be generatively inert for both scDesign3 and
SymSim (it simply doesn't do anything in either simulator), and SymSim's
dropout axis showed genuine non-monotonic behavior that couldn't be
explained by a simple label mix-up. SymSim's batch axis separately showed
an unexplained positive relationship going the wrong direction — reported
honestly, flagged as unreliable, and not chased further into the
simulator's internals, since that would have been a much larger
undertaking than the finding warranted.

A 108-row critical-thresholds table was built — six axes times three
simulators times six methods — giving the interpolated stressor value at
which Subspace Recovery Score crosses 0.5 for each combination. 24 of
those 108 rows carry an explicit "unreliable" flag, tracing back to the
clipping-inertness and SymSim anomalies above.

**Cross-simulator triangulation caveat.** This project's core epistemic
argument -- agreement across 3 independent simulators rules out
simulator-specific artifacts -- is degraded for 3 of the 6 stressor axes,
not just the 2 originally called out inline above: clipping is reliable
in only 1/3 simulators (Splatter), while both batch and dropout are
reliable in 2/3 (SymSim excluded from each, for different confirmed
reasons). Full detail and evidence in docs/step5_axis_findings.md's
summary section. This should be stated explicitly in the manuscript
wherever cross-simulator triangulation is invoked as supporting evidence.

Two families of phase-space heatmaps were rendered — sparsity-by-depth
and dropout-by-depth, one per method — as the final visual deliverable of
this step.

---

## Step 6 — Failure Mode Detection & Characterization (6.1–6.9)

This is where the project moved from mapping where recovery fails
geometrically to actually naming and testing seven specific mechanisms
that could explain why: technical separation, cluster collapse, phantom
clustering, variance hijacking, over-smoothing, neighborhood collapse,
and subspace rotation slippage. Every one of these got its own detection
script, and — almost without exception — writing that script surfaced
something that needed fixing before the result could be trusted. That
turned out to be the real character of this step: not "run seven tests,"
but "build seven tests properly, catch what breaks, fix it, and only then
believe the number."

Two loose ends from earlier steps got tied up along the way, not because
they were originally scoped into Step 6, but because Step 6 couldn't
proceed correctly without them. The 30 newer null-control replicate files
from Step 1.6 had never actually been run through Step 3's preprocessing
— only the original 15 had — which meant Step 6.3 would have run on a
third of its intended data. That got fixed by writing a small conversion
script and simply re-running the existing, unmodified Step 3 scripts,
which skip anything already done. And EDA Checkpoint 2's long-deferred
recommendation to re-check pancreas batch visibility once real
normalization was in place — it finally got a proper answer, as part of
Step 6.4.

### 6.1 — Technical Separation

The idea here is simple: does batch identity or sequencing depth end up
correlating with the clusters a method produces, more than real biology
does? Measured via Adjusted Mutual Information between k-means clusters
and two technical covariates — batch and UMI-count quartile.

Getting the k-means configuration right took real work. The stored ARI
values from Step 4.5 had to be reproduced exactly before any new
computation could be trusted, and that meant testing `nstart=1` against
`nstart=25` and checking which one actually reproduced the stored numbers
— `nstart=25` won clearly, matching 6 of 10 test cases exactly against
only 3 for `nstart=1`.

Along the way, a genuine bug in the ARI/AMI package itself turned up: it
silently mis-handles character labels that aren't purely numeric strings,
converting them to integers in a way that quietly produces `NA` instead
of failing loudly. That had to be worked around everywhere `aricode` got
used for the rest of Step 6. A second, separate bug — an `%in% FALSE`
filter that silently dropped every real-data row, since the relevant
field is `NA` rather than `FALSE` for real data — went unnoticed through
several dry runs before it was caught.

The most consequential find, though, was that all six of PBMC68k's
embedding files had completely blank ground-truth labels — every single
cell. This traced straight back to Step 4.7's re-annotation fix: it
updated the manifest, but never touched the embedding files' own internal
copy of the labels. Once caught, the fix was straightforward — patch the
labels back in from the corrected source file — but confirming the fix
was actually correct took some care, since neither file had an
independent cell-identifier to check row alignment against directly. That
got resolved indirectly, by recomputing ARI after the patch and checking
it matched Step 4's stored values almost exactly.

Final result: 196,813 valid rows, 12.8% flagged for technical separation.

### 6.2 — Cluster Collapse

Whether two supposedly distinct biological populations have effectively
merged — measured by comparing the distance between their centroids
against how spread out each population is internally. This one ran
cleanly, no real bugs, and lined up well with expectations: 83.3% of rows
triggered the silhouette-based warning sign, and 1.4 million individual
population pairs were flagged as collapsed across the whole grid. One
specific pair — CD4 and CD8 T-cells, under scDesign3 with GLM-PCA — sat
right at the exact boundary of collapse across many different parameter
settings, which makes biological sense given how similar those two cell
types are transcriptomically.

### 6.3 — Phantom Clustering

Run on the null-control data specifically — files with only one real
biological population — to see whether any method invents apparent
cluster structure where none should exist. Two tests were used: ARI
against a randomly-assigned set of labels, and the Calinski-Harabasz
index.

The ARI-against-random test behaved exactly as its underlying statistics
predict it should — since ARI is chance-corrected by construction, it
came back essentially at zero everywhere, which isn't actually
informative one way or the other. The Calinski-Harabasz test needed more
care: taken at face value, every single file flagged, which seemed too
good (or bad) to be true. A quick simulation confirmed why — pure random
noise, with no structure at all, reliably produces a CH score around 27
just from the geometry of forcing k-means into three groups regardless.
So the real threshold had to be calibrated against that noise floor
rather than trusted at face value.

Once corrected, the picture matched the EDA Checkpoint 3 hunch only
partly. GLM-PCA really is a strong phantom-clustering offender — nearly
universal, and by a wide margin. SCTransform v2, which had been flagged
alongside GLM-PCA on visual inspection, turned out on the formal test to
be the *best*-performing method of all six by this metric. The visual
flag was half right.

### 6.4 — Variance Hijacking

Checking whether the top principal components end up correlating more
with technical covariates than with real biology — the formal version of
EDA Checkpoint 2's deferred batch-visibility question, finally answered
properly against real, finished embeddings rather than an early
exploratory PCA.

One methodological wrinkle needed sorting out: Spearman correlation
assumes some kind of ordering, which is fine for the simulated batch axis
(it only ever has up to three levels, and they're explicitly ordered by
complexity), but doesn't really make sense for real datasets, where batch
is just an arbitrary donor or sample ID with no natural order — Tasic
2018 alone has 341 of them. A different statistic, the correlation ratio,
was used for that case instead.

34.2% of rows flagged overall. One genuinely reassuring number came out
of this too: in nearly half the grid, real biology was already the
dominant signal at the very first principal component, ahead of any
technical covariate.

### 6.5 — Over-Smoothing

This one split cleanly into two very different tests, because the
original specification described something that didn't map directly onto
either our real or simulated data.

For real data, the question was whether a method preserves known
biological trajectories — tested via diffusion pseudotime on the three
pancreas datasets, since pancreatic endocrine differentiation is a
well-established continuous process in the literature (PBMC68k was left
out here, since its five cell types are discrete, terminally
differentiated lineages with no comparable single trajectory). The
standard R package for this, `destiny`, wouldn't install — a broken
dependency chain going several packages deep — so diffusion pseudotime
was implemented directly from the published method instead, using only
packages that were already available. Ten of the eighteen
dataset-method combinations tested (55.6%) flagged as over-smoothed, and
GLM-PCA was the worst performer across all three datasets, including the single
lowest correlation score seen anywhere in this whole step.

For simulated data, the first design that got built turned out to be
mathematically circular — it compared a method's own fitted loadings
against ground truth derived by projecting through those same loadings,
which is guaranteed by basic linear algebra to come out looking almost
perfect no matter what, regardless of whether the method actually did
anything meaningful. That was caught before the result got reported
anywhere, and the whole test was rebuilt from scratch to compare two
genuinely independent things instead — how much of the real signal
survives in full gene space versus in the compressed embedding. Under the
corrected version, restricted to the four simpler linear preprocessing
methods (SCTransform v2 and GLM-PCA don't have a simple enough transform
to build this comparison from), not a single one of 31,640 tested cases
came back flagged. That's a real, useful negative result for that
specific scope — it should not be read as "over-smoothing doesn't happen"
more broadly, since it only covers four methods and two narrow separation
regimes.

### 6.6 — Neighborhood Collapse

Whether a cell's local neighbors in the real embedding still resemble its
neighbors in an independent reference space. This one had the roughest
path to a clean result of any substep in Step 6.

Building the reference space needed computing gene-level variance on a
sparse matrix, and the way that got written initially — using base R's
`apply()` — turned out to force expensive, memory-heavy conversions under
the hood. Combined with four worker processes running at once, it caused
outright crashes on the two largest, gene-richest real datasets, crashes
that didn't even show up as ordinary errors since the whole worker
process was dying, not just throwing an R-level exception. Rewriting that
one calculation to use a sparse-friendly formula, and additionally
pre-building the real datasets' reference spaces one at a time before
letting the parallel run loose, cleared it up.

A second, subtler problem showed up only at full scale: a caching
function was writing its output directly to disk without the
temp-file-then-rename pattern used everywhere else in this project, which
meant two worker processes could occasionally race to write the same
cache file at once and leave a corrupted, half-written file behind for
everyone downstream to trip over. Once that pattern was fixed to match
the rest of the project's convention, the handful of corrupted cache
files were found, deleted, and the affected rows recomputed cleanly.

Final result: 93.5% of rows flagged. Worth being upfront about what that
number really means, though — this is a strict, binary yes/no test on
whether specific neighbors match exactly, which is a much harsher
standard than the graded continuity metric computed back in Step 4. A
neighbor that's just one rank off still counts as a complete miss here,
even though it would barely register on a graded scale. So this high
number reflects the strictness of the test as much as it reflects actual
structural damage — both things are true at once, and the writeup needs
to say so.

### 6.7 — Subspace Rotation Slippage

Whether individual principal components rotate away from the true
biological subspace under stress, and if so, which one goes first.

The existing subspace-recovery metric from Step 4 compared whole
subspaces to each other and had no way to identify which specific,
numbered principal component was responsible — so a new angle calculation
was built, measuring each individual PC's own angle to the true subspace
directly. It was checked against Step 4's stored numbers in the one case
where the two calculations should mathematically agree exactly (a
rank-one true subspace), and they did, down to two decimal places.

A third of an early test batch failed outright with a matrix-dimension
error, and it traced to something specific: two of the six methods —
SCTransform v2 and GLM-PCA — cap the number of genes they use internally,
independent of the general gene-selection setting elsewhere in the grid,
so their loadings didn't line up gene-for-gene with the full ground-truth
matrix for the larger simulated datasets. Matching genes by name instead
of by position fixed it.

The main finding here is one of the cleaner, more publishable results in
the whole step: PC1 rotates past the 30-degree threshold in 88% of rows,
and — more importantly — it does so already at the mildest stress level
tested, across every single axis and method combination. In other words,
within the range of conditions this study actually tested, rotation
isn't a threshold you cross as things get worse — it's already happened
by the time you get to the easiest condition, and what the stressors
actually control is how much worse it gets from there. Separability turned
out to be by far the strongest driver of that severity, for the four
methods that respond to stress at all; Raw PCA and GLM-PCA barely move
regardless of condition, because they're already substantially rotated
from the very start.

### 6.8 — Cross-Referencing Against Step 5's Boundaries

With all seven failure modes computed, the natural next question was
whether they actually line up with where Step 5 had already identified
geometric recovery failure. Three of them — cluster collapse,
neighborhood collapse, and subspace rotation slippage — line up almost
perfectly: everywhere Step 5 found recovery always failing, all three of
these independently agree. Two others — technical separation and
variance hijacking — turned out to be largely decoupled from geometric
failure. That's a real, interesting result on its own: a method can lose
real biological signal without technical covariates ever actually taking
over the leading components. Those are, apparently, two mostly separate
ways preprocessing can go wrong, not the same failure showing up twice.

One loose thread got flagged rather than resolved here: three specific
parameter combinations where Step 5 says recovery never fails, but three
independently-computed failure modes all disagree with that and with each
other's independent computation agree with each other. That's a genuine
open question for the manuscript's limitations section, not something
that got quietly explained away.

### 6.9 — Final Retention Review

The last step was going back over all seven failure modes with the full
weight of evidence in hand and deciding, honestly, what deserves to stay
in the final characterization and what needs a caveat attached. All seven
survived — none of them turned out to be unsupported once actually
tested — but three needed a real interpretive caveat (technical
separation and variance hijacking's decoupling from geometric failure;
neighborhood collapse's test-strictness issue), and over-smoothing needed
to be reported as two separate results rather than blended into one,
since its real-data and simulated-data branches told genuinely different
stories. The full writeup, including a section explaining exactly why
each retained finding can be trusted — every one of them survived a real
bug being found and fixed before it was reported — lives in
`docs/step6_9_failure_mode_review.md`.

---

## Key Learnings & Principles

Some of these go all the way back to Step 1; some are new from Step 6.
Worth keeping all of them together, since a few repeated themselves in
slightly different form more than once.

- An axis having zero real effect, despite being logged as if it does, is
  best caught by direct byte-for-byte comparison across labels that
  share an underlying fit — this single check found both scDesign3's and
  SymSim's sparsity bugs.
- A fit-once-save-many architecture is efficient, but anything not
  included in the fit key needs its own, separately verified downstream
  mechanism — never assume a label implies a real effect just because
  it's recorded.
- A label pointing the wrong direction is worse than a label that does
  nothing at all, since it actively produces a misleading trend rather
  than just an uninformative one. Always verify direction with a
  controlled sweep, never from documentation or a variable's name alone.
- When two things are genuinely coupled (dropout reducing effective
  depth, say), the fix is usually not to decouple them — that coupling
  might be realistic and worth keeping — but to make whatever's comparing
  them aware of it.
- Percentage-deviation metrics get unstable at small denominators. Always
  check whether "high deviation" is really concentrated at some tiny
  extreme corner before treating it as a systemic problem.
- Rank-order and monotonicity checks on single-draw stochastic data need
  an explicit, empirically-derived tolerance — strict inequality will
  flag ordinary sampling noise as if it were a real violation. Derive the
  tolerance from the actual observed ceiling, not a round number picked
  in advance.
- When data is wrong only in its label, relabeling beats regenerating —
  just as rigorous, cheaper, and doesn't introduce fresh risk from a new
  generation run.
- `mclapply` can behave unpredictably alongside certain C-library-backed
  packages even when the same code is fine in serial — always test both
  paths when a fork failure looks strange, and don't assume the bug is in
  your own logic just because it only shows up under parallelism.
- `nohup` protects against the shell hanging up, not against the whole
  terminal application being closed. Checkpointing is the real safety
  net, and it's paid off more than once on this project.
- Packages that quietly coerce or silently drop things instead of
  failing loudly (`aricode`'s character-to-integer coercion, `%in%`'s
  handling of `NA`) are a recurring, genuinely dangerous class of bug —
  worth testing explicitly rather than trusting.
- Atomic writes — temp file, then rename — aren't optional decoration.
  Skipping that pattern even once, in Step 6.6's cache-writing function,
  produced a real race condition and real file corruption under
  concurrent access.
- A test that looks mathematically elegant can still be circular. If a
  new metric can be reduced to an identity that holds regardless of the
  data, it isn't measuring what it appears to measure — this cost a full
  day of work in Step 6.5 before being caught.
- `param_dict.R` remains the single source of truth for every parameter
  mapping, and `run_id` remains the universal seed-recovery mechanism,
  all the way through to the end of Step 6.

---

## Current State & Next Steps

Steps 1 through 6, in full — 1.1 through 6.9 — are complete as of this
writing. Two loose ends from earlier steps got tied off along the way in
Step 6 rather than left hanging: the 30 unprocessed null-control
replicates from Step 1.6, and EDA Checkpoint 2's deferred pancreas
batch-visibility question.

**Genuinely open items, carried forward rather than resolved:**
- The three-combination discrepancy between Step 5.8's thresholds and
  three independent Step 6 failure modes (see Step 6.8 above) — a real
  candidate for the manuscript's limitations section, or a targeted
  follow-up before submission.
- A stale, superseded preprocessing output directory
  (`data/processed/pca_sctransform/`, without the `_v2` suffix) turned up
  during Step 6.5 and was left in place rather than deleted without
  explicit confirmation — still sitting there as of this writing.
- `renv::status()` has reported the environment as out of sync since a
  package install partway through Step 6.1 — not yet reconciled.
- scDesign3's isolated depth under-delivery at `depth_label="500"`,
  flagged back in EDA Checkpoint 1, was deliberately deferred and never
  revisited.

**Repository housekeeping**, done separately from the scientific work
above but worth recording here: the repo's `.gitignore`, commit history,
and folder structure were all audited folder by folder ahead of
publication. Thirty-two stray raw log files that had been committed
directly at the repository root — leftover terminal output that never
made it into the proper `logs/` directory — were removed, and a rule was
added so it can't happen again. The README was rewritten, since it had
been describing a project structure that hadn't matched reality for a
long time. Every other folder — `scripts/`, `R/`, `data/`, `docs/`,
`results/`, `figures/`, `logs/`, `renv/` — was checked directly and found
to already be organized the way it should be.
