# AGENTS.md — Gourgouthakas cave 16S

Working agreement for this repository. Read alongside `PLAN.md`, which
specifies the downstream analysis; this file records what is **actually true
of this dataset and this machine**, including the places where reality departs
from `PLAN.md`.

## What this is

Full-length 16S rRNA (PacBio HiFi / Kinnex) from **Gourgouthakas cave, Crete**
— a 1100 m deep vertical cave system. Sediment and water sampled down a depth
transect on 2022-08-09. Upstream processing is `HiFi-16S-workflow` (Nextflow,
DADA2 ASVs); downstream is a tidyverse + vegan ecology analysis per `PLAN.md`.

## Everything runs in the container

Both the Nextflow pipeline and the R analysis run inside
`localhost/gourgouthakas-16s`. The host carries **podman only** — no conda, no
docker, no host Nextflow. Do not add host-level tool installs.

```
podman build -f Containerfile -t gourgouthakas-16s .
./scripts/run_pipeline.sh download   # once: SILVA + GTDB + GG2
./scripts/run_pipeline.sh run
```

Image layering, each file building on the one above:

| File | Adds |
|---|---|
| `omics-16s.Containerfile` | base: tidyverse, dada2 1.40.0, vegan, SRS, vsearch, fastp — built from `localhost/r-tidyverse` |
| `Containerfile` | R analysis packages (§5 of PLAN) **+** Nextflow **+** pinned pipeline tools in the `hifi16s` micromamba env |

Two R installations coexist on purpose: the **system R** (dada2 1.40.0) for the
downstream analysis, and `/opt/conda/envs/hifi16s` (dada2 **1.38.0**, matching
`HiFi-16S-workflow/env/dada2.yml`) for the pipeline. The env is kept off the
default PATH; `scripts/nextflow.config` prepends it for task execution
only. Don't "simplify" this by merging them — it silently changes the DADA2
version the ASVs are inferred with.

`HiFi-16S-workflow` ships no podman profile and this host has no conda, so
`scripts/nextflow.config` runs every process **natively inside the
container** with both `enable_conda` and per-process containers disabled. It
also caps resources to 14 cpus / 26 GB, because `conf/base.config` requests up
to 64 cpus and 500 GB for some labels and Nextflow would refuse to schedule
them on this 16-core / 31 GB host.

## The data

`data/PB482_SP/` — 51 libraries, **all 51 md5-verified** against
`data/md5sum-16s.txt` (run `md5sum -c` from inside `data/`, the paths are
relative to it). `data/PB482_SP.tar` is the untouched delivery archive.

51 libraries, with **two nested levels of replication** on the sediment side:

- **Biological replicates** — at each of the 9 depth sites (C1–C9), two
  independent sediment samples were taken: the main depth-transect sample
  (`C1`, …) and a parallel isolate-source sample (`C1_I`, …). These are the
  `transect` / `isolate_source` arms of `sample_set` below.
- **Technical replicates** — each of those (and each water sample) went
  through **two independent DNA extractions**, e.g. `C1_1_A` vs `C1_2_A`. This
  is the `replicate` column (`1`/`2`) that `build_inputs.R` parses out of the
  library name.

| Set | Sites | Libraries |
|---|---|---|
| `transect` sediment | C1–C9 | 18 (9 × 2 technical reps) |
| `isolate_source` sediment | C1_I–C9_I | 18 (9 × 2 technical reps) |
| `transect` water | W1–W5 | 10 (5 × 2 technical reps) |
| `control` | 5 | 5 (no replicate structure) |

Water has no biological-replicate arm (no `W1_I`) — the two-level nesting
applies to sediment only.

Controls: `ctr_zymo_com`, `ctr_zymo_log`, `ctr_msa_3001` (mocks),
`ctr_EB` (extraction blank), `ctr_MM` (mastermix blank).

**Column names, and one gap to know about:** `data/metadata.tsv` has `site`
(the *biological*-replicate-level id — `C1` and `C1_I` are two distinct
`site` values, not one), `sample_set` (`transect` / `isolate_source` /
`control` — the column that actually distinguishes the two biological arms),
and `replicate` (`1`/`2`, the *technical*/extraction-replicate id). There is
**no ready-made column pairing `C1` with `C1_I` as replicates of the same
depth location** — derive it by stripping the `_I` suffix from `site`
(`sub("_I$", "", site)`) if a script needs to group them, e.g. for the
biological-replicate concordance check in `PLAN.md` §7 `02_qc_filter.R`.

`scripts/build_inputs.R` (base R, runs on the host) generates
`data/samplesheet.tsv` and `data/metadata.tsv` from the barcode map
plus the two ENA/MIxS checklist sheets. Regenerate rather than hand-editing —
the pipeline's `inspect_metadata` process diffs the two sample-id columns and
dies on any mismatch.

**One row per library, not per site.** Both replicate levels — technical
(extraction) and biological (`sample_set` arm) — stay separate through the
pipeline so concordance can be measured at each level before anything is
averaged away; pool them downstream, deliberately (`PLAN.md` §7
`02_qc_filter.R` covers both levels explicitly).

## Departures from PLAN.md — read before implementing 02–10

1. **No phylogeny.** The workflow emits no tree, and we opted not to build one.
   **Faith's PD (§04) and UniFrac (§06) are out of scope.** Use Bray–Curtis and
   Aitchison distances. `picante` / `GUniFrac` are deliberately not installed.

2. **Spatial analysis (§09) does not apply as written.** Every sample shares
   one lat/lon (35.33213, 24.08346) — it is a single cave. There is no
   distance–decay and no Moran's I over geography. **The gradient is `depth_m`
   (0 → 1100 m), with `elevation_m` its exact complement.** Rework §09 as a
   depth transect: community turnover vs. depth, distance–decay in *depth*
   space. `sf` / `terra` are not installed.

3. **Controls exist, so `decontam` is on** (§11.3 of PLAN). `ctr_EB` and
   `ctr_MM` are the negatives; prevalence method. The three mocks are a
   separate accuracy check, not decontam input — the `control_type` column
   distinguishes them.

4. **Taxonomy: all three DBs** (SILVA 138.2, GTDB r220, GG2 2024.09) via naive
   Bayes, SILVA also via VSEARCH. `db_to_prioritize` stays at the pipeline
   default `GG2`. Rank parsing must handle whichever the merged table carries.

5. **`condition` in PLAN terms is not one variable here.** The real contrasts
   are `sample_type` (sediment vs. water) and the continuous `depth_m`. The
   `sample_set` column separates the transect from the isolate-source
   sediments — these are the **biological replicates** at each site (see
   "The data" above) — do not pool those two arms without checking they agree
   at the ASV level, and don't stop there: check the **technical**
   (extraction) replicates too, at the `replicate` column level, before that.
   `PLAN.md` §7 `02_qc_filter.R` runs both checks explicitly — the ISD
   reference scripts (`isd_archive_scripts/isd_crete_numerical_ecology.R`)
   only ever checked the biological level, which is why this needed calling
   out separately here.

## Reference databases

Downloaded once with `scripts/run_pipeline.sh download` into
`/mnt/data/databases` (11 GB), verified 2026-07-19:

| DB | naive-Bayes trainset | VSEARCH seqs / taxonomy rows |
|---|---|---|
| SILVA 138.2 | 452,055 seqs | 510,495 / 510,495 |
| GTDB r220 | 61,774 seqs | 863,832 / 863,832 |
| GG2 2024.09 | 337,506 seqs | 23,467,470 / 23,450,268 ⚠ |

**GG2's VSEARCH pair is inconsistent** — 17,202 more sequences than taxonomy
rows. Harmless as configured, because `vsearch_databases = ['silva']` and GG2 is
used for naive-Bayes only (which reads the trainset, not these files). Reconcile
before pointing VSEARCH at GG2.

## Environment traps that cost real time

Each produces a silent hang or a misleading error; each is already fixed
(in `Containerfile` / `scripts/run_pipeline.sh` for #1-2, `.Rprofile` for #3,
`scripts/nextflow.config` for #4), documented here so none is re-introduced.

1. **SELinux is enforcing.** Podman bind mounts are inaccessible without a
   label option — every file under `/work` returns "Permission denied", and
   Nextflow reports it as a *remote pipeline not found* error, which points
   nowhere near the real cause. `run_pipeline.sh` passes
   `--security-opt label=disable` rather than the `:z` mount flag, since `:z`
   recursively relabels and would rewrite SELinux contexts across the whole
   4.5 GB raw dataset on every run.

2. **Fedora ships wget2 as `/usr/bin/wget`.** The download processes call
   `wget -O`. zenodo.org resolves to IPv6 first, the rootless container has no
   IPv6 route, and neither wget does Happy Eyeballs — so it connects to the v6
   address and blocks until timeout, retrying forever with zero bytes written
   and nothing on stdout. It looks hung rather than failed; the real error
   (`errno=110`) only appears in `work/*/*/.command.log`. `curl` succeeds
   throughout and masks the problem. Fixed by putting GNU wget 1.25 in the
   pipeline env with `prefer-family = IPv4` in its system `wgetrc`.

Also: Nextflow 26.04's strict config parser rejects the bare `$HOME` in
upstream's `nextflow.config`, so `run_pipeline.sh` sets `NXF_SYNTAX_PARSER=v1`
rather than patching a checkout we do not own. Drop it once upstream quotes
that as `env('HOME')`.

3. **`renv`'s autoloader crashes every plain `Rscript` invocation.** `.Rprofile`
   sources `renv/activate.R`, which tries to bootstrap `BiocManager` into
   `~/.cache/R/renv` (renv's default cache root) to resolve the lockfile's
   Bioconductor entries — and that install fails outright in this image,
   halting before any script code runs (`Error: failed to install
   "BiocManager"`). A `.renv-cache/` with real cached packages sits at the
   project root, but that's not the path renv actually resolves to, and the
   project-local `renv/library/` only ever held `renv` itself — the 357
   packages in `renv.lock` were never installed *through* renv's private
   library; they're the system R library (`/usr/lib64/R/library`), confirmed
   version-for-version identical to the lockfile. So `renv.lock` here is a
   manifest, not a live dependency source. Fixed by setting
   `Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")` at the top of
   `.Rprofile`, before it sources `activate.R` — every script still gets
   exactly the same packages, just without the broken bootstrap attempt.

4. **`conf/base.config`'s `highparallel` label caps `time` at 8h.** The four
   deep sediment libraries (below) all exceed that. Hitting the limit doesn't
   fail the task cleanly — Nextflow's local-executor timeout-kill path races
   with its own exit-status check (`java.lang.IllegalThreadStateException:
   process hasn't exited` in `LocalTaskHandler.checkIfCompleted()`) and
   **aborts the entire session** instead, observed twice at exactly the 8h
   mark (`C1_2_A`, then `C1I_1_A`). `errorStrategy`'s retry-on-exit-143 never
   gets a chance to run because it's an executor-level exception, not a normal
   task failure. Fixed in `scripts/nextflow.config` by overriding
   `withLabel: highparallel { time = 96.h }` and raising `resourceLimits.time`
   to match — do not remove this override, the upstream 8h default is too low
   for this dataset's biggest libraries.

## Runtime: denoising is diversity-bound, not depth-bound

`dada2_denoise_independent` costs scale with **unique sequences**, not reads.
Measured on this dataset:

| Library | Filtered reads | Unique seqs | Denoise time |
|---|---|---|---|
| `ctr_msa_3001` (mock, ~8 species) | 41,281 | 16,142 | ~1m 8s |
| `C1_1_A` (cave sediment) | 38,742 | 27,601 | 49m 9s |
| `C6I_1_A` | 583,374 | 222,653 | 16h 57m |
| `C1I_1_A` | 487,320 | 252,801 | 1d 15h 40m |
| `C1I_2_A` | 588,618 | 308,920 | 2d 6h 27m |
| `C6I_2_A` | 778,394 | 291,931 | 1d 4h 4m |

Same order of magnitude of reads, wildly different times — do not estimate
runtime from read counts *or* from unique-sequence counts alone. `C6I_1_A` and
`C1I_1_A` have near-identical uniqueness (222k vs 253k) but a 2.3× time gap,
and `C6I_2_A` finished faster than `C1I_2_A` despite having more reads and
comparable uniqueness — whatever drives cost also depends on the specific
error/quality structure per site, not just diversity. A full run of all 51
libraries took **~9 days** end to end (resume-to-resume, across the 8h-timeout
interruption above), almost entirely spent on these four libraries.

The stage is also **serialised**: the process carries label `highparallel`
(16 cpus). Even with `resourceLimits` raised to 20 cpus (this dev host has
24 cores / 251 GB, not the 16-core/31 GB target host this file otherwise
assumes — check `nproc`/`free -h` before reusing these limits elsewhere),
16-cpu tasks still can't run two at a time, so exactly one denoise task is
schedulable at once regardless. `denoise_independent.R` also hardcodes
`multithread = TRUE`, so it takes every allocated core regardless of
`task.cpus` — the CPU is saturated, not idle. Reducing cpus/task to raise
concurrency is plausible but would oversubscribe threads; not attempted.

**Always run with `DETACH=1`** (`scripts/run_pipeline.sh`) when launching from
the host. Attached, the run is a child of the calling shell and dies with the
session — that happened once here, mid-denoise. Detached it is owned by
conmon/systemd and survives. (When already inside the container, `DETACH` is a
no-op — `run_pipeline.sh` just execs `nextflow` directly — so back the job with
a session-independent background mechanism instead.)

**Pipeline run status:** completed successfully for all 51 libraries
(0 failures, `[SUCCESS] completed=38 failed=0 cached=283`). Final ASV +
taxonomy tables are in `results/hifi/final/` (`best_tax_merged_freq_tax.tsv`
is the main one — 30,532 ASVs × 51 samples + taxonomy); per-database taxonomy
detail is in `results/hifi/nb_tax/` and `results/hifi/vsearch_tax/`.

## Known data-quality issue

**FIXED: W5's conductivity and temperature were transposed** in
`data/metadata.tsv` — `conductivity_ms = 8.4`, `temperature_c = 195.9`,
impossible in a cave stream (every other water sample sits near 178–188
mS/cm and 6.6–7.9 °C). Corrected directly in `data/metadata.tsv`
(`temperature_c = 8.4`, `conductivity_ms = 195.9` — matches the depth trend:
6.6 → 7.3 → 7.6 → 7.9 → 8.4 °C across W1–W5. Note: the source MIxS sheet's
own value is `196.9`, not `195.9` — a ~1 mS/cm rounding-level discrepancy
from a hand edit rather than a `build_inputs.R` regeneration; immaterial to
any result, but regenerate via `build_inputs.R` rather than hand-editing
again if exactness matters).

Left uncorrected, this wasn't just a cosmetic problem: `conductivity_ms` is
**only measured for water**, so W5's swapped pair was 2 of the only 10
complete-case observations `07_environmental_drivers.R`'s correlation/VIF/
db-RDA/envfit/varpart/Mantel analyses had to work with — it was forcing a
mechanically-exact **-1.00** temperature~conductivity correlation that
wasn't real, and roughly halved the db-RDA R² (0.076 → 0.143 once fixed).
Caught by auditing `07_env_correlation.png` for exactly this kind of
suspicious exact value, the same way the elevation_m/depth_m -1.00 was
legitimate but this one wasn't.

`build_inputs.R` now asserts water `temperature_c` is in a plausible 0–40 °C
range and refuses to write `metadata.tsv` otherwise; `01_import.R` carries
the same guard for whatever `metadata.tsv` is actually on disk. Neither
script silently "fixes" a transposition if one recurs — both stop loudly.

Sediment temperature is missing for C7, C8, C9 (and their `_I` counterparts) —
genuinely absent, not a parsing failure.

## Conventions

Follow `PLAN.md` §6 (snake_case, `here::here()`, samples-as-rows for vegan,
scripts read `data/` + earlier `results/` and write only `results/` + `plots/`).
Beyond that:

- **No phyloseq, no qiime2R, no `.qza`** anywhere. This is non-negotiable in
  `PLAN.md` and the container does not carry them.
- `data/PB482_SP/` is **immutable**.
- `renv` manages the project library, hydrated from the container's system
  library (`renv::hydrate()`) rather than recompiled — see `scripts/00_setup.R`.
- Generated inputs (`data/samplesheet.tsv`, `data/metadata.tsv`) are
  build products of `build_inputs.R`; edit the script, not the outputs.
- Reference databases live at `/mnt/data/databases` on the host and are mounted
  at `/databases` in the container. Override with `DB_DIR=... scripts/run_pipeline.sh`.

## Downstream analysis notes

Statistical choices made while implementing `PLAN.md` §7, and why. Read
before changing any of these -- each was picked for a reason specific to
this dataset, not a generic default.

- **Rarefaction depth (499) comes from the upstream pipeline's own
  `results/hifi/final/rarefaction_depth_suggested.txt`**, not a re-derived
  quantile — reused rather than duplicated. Its sibling
  `alpha_depth_suggested.txt` (116877) is a much higher *saturation* depth
  (how deep the curve needs to go to plateau), not a depth to rarefy to;
  using it would drop all but a handful of the highest-biomass libraries.
  **This value was computed over all 51 libraries including the 5 controls**
  (`HiFi-16S-workflow/bin/final_stats.R` sorts every sample's final read
  count with no `sample_type` awareness and picks rank `floor(0.8*N)`) —
  the near-empty blanks (`ctr_EB`=41, `ctr_MM`=70 raw reads) sit in that
  sorted list and pull the chosen rank down. `03_normalize.R`'s T8 audit
  reports the alternative: recomputed over the 46 real libraries alone, the
  same rank-`floor(0.8*46)` value is **409**, not 499 — lower, not higher
  (removing 5 shallow libraries from the bottom of the sort shifts every
  rank above them down by 5 positions, which moves the 80th-percentile
  cutoff to a shallower point in the real-sample distribution than where
  it sat in the mixed 51-library list). **Stop-and-ask, not changed**: this
  script still rarefies to 499, the upstream value. 12/46 samples fall
  below 499 reads and are excluded from the rarefied matrix only
  (`03_normalize.R`) — still present in relative-abundance/CLR.

- **CLR needs its own, coarser ASV filter, re-derived on the 46
  control-free libraries.** The general QC prevalence filter (now
  `02b_controls.R`'s final gate, ≥2 samples, controls already dropped)
  leaves 15,457 ASVs — far too sparse for compositional zero-replacement:
  `zCompositions::cmultRepl()`'s default thresholds (`z.warning`/`z.delete`
  = 0.8/`TRUE`) **silently dropped samples** the first time this ran,
  because most ASVs are zero in >80% of samples at that resolution. Fixed
  in `03_normalize.R` with a CLR-specific ≥5-of-46-samples (~10.9%)
  prevalence filter (1,201 ASVs, up from 717 when the same ~10% threshold
  was computed over the control-contaminated 51-library table — the mocks'
  low diversity had been inflating apparent sparsity, as the original
  controls plan suspected) plus `z.delete = FALSE` explicitly, so a still-
  sparse sample (the shallowest sediment libraries) gets a warning, never a
  silent drop. Check `setequal(rownames(counts_clr), rownames(counts_clean))`
  after any change here — it should always be `TRUE`.

- **Controls were leaking into the ecological analysis — fixed by
  reordering the filter pipeline and splitting `02_qc_filter.R` into two
  scripts.** The old order prevalence-filtered (≥2 samples) while the 5
  control columns (`ctr_EB`, `ctr_MM`, `ctr_zymo_com`, `ctr_zymo_log`,
  `ctr_msa_3001`) were still present, so an ASV seen only in `ctr_EB` +
  `ctr_MM` and nowhere real could pass on the controls' strength alone, and
  `03_normalize.R`'s old CLR filter (≥10% of 51 samples) had the identical
  problem. `02_qc_filter.R` now stops after lineage-filtering + library-size
  QC + flagging (not removing) `decontam` candidates, and hands off
  `counts_postlineage.rds`/`taxonomy_postlineage.rds` (controls still
  present) to the new **`scripts/02b_controls.R`**, which: (T1) normalizes
  GTDB/GG2 polyphyly-split genus suffixes (`_[A-Z]`, `_[A-Z]_[0-9]+` —
  confirmed both databases need this, not just GTDB, from the same ASV
  reading `Akkermansia_muciniphila_A` under GTDB and
  `Akkermansia_muciniphila_D_776786` under GG2) against
  `data/mock_expected.tsv`; (T2) scores all 4 classification methods against
  the 3 mock controls; (T3) reports the ASV-level false-positive rate from
  the mocks; (T4) inventories every ASV present in the blanks; (T5) computes
  the bidirectional PacBio barcode-cross-talk floor; then applies the
  `decontam` blacklist + bleed floor to the **full 51-library table**,
  drops the 5 control libraries, *then* prevalence-filters, with a
  `stopifnot` guard that no control row survives into `counts_clean.rds`.
  `04`–`10` no longer need their own `sample_type != "control"` filters —
  removed throughout, since controls never reach `metadata_clean.rds` in
  the first place. See `data/mock_expected.tsv` for the versioned,
  web-sourced (Zymo/ATCC product pages, cited by URL) expected compositions
  this all scores against — not reconstructed from memory.

- **Database verdict from the mock scorecard (T2): no change indicated.**
  Scored SILVA-NB/GTDB-NB/GG2-NB/SILVA-VSEARCH against all 3 mocks
  (`results/02b_mock_accuracy.tsv`); mean recall by method: `gg2_nb` = 0.93,
  the best of the four. `PLAN.md`'s default `db_to_prioritize = GG2` is
  competitive on this dataset's own mock controls, not just a generic
  default — kept as-is.

- **`decontam` runs `method="prevalence"`, `threshold=0.5`** (up from an
  earlier, undocumented implicit 0.1), against `ctr_EB` + `ctr_MM` (n=2
  negatives) with no `batch` argument — there is no extraction-date/plate
  field recorded anywhere in `data/` to stratify on, a limitation, not an
  oversight. n=2 negatives is thin statistical power at *any* threshold;
  its output (`results/02b_contaminant_candidates.tsv`) is a **manual-review
  candidate list, not automated removal** — a `decision` column, empty on
  first run, that a human fills in and the script round-trips on every
  subsequent run into `results/02b_contaminant_blacklist.tsv`. **Do not read
  an empty blacklist as a clean bill of health**: `ctr_EB` (41 raw reads)
  and `ctr_MM` (70 raw reads) are shallow negatives, so contamination
  screening here is limited by control read depth, not exhaustive.

- **Bidirectional bleed floor (T5) came back at 0.00000 (0%)** — no
  mock-exclusive ASV was detected in any real sample, and no cave/water-
  characteristic genus (prevalence ≥5 in reals) was detected in any mock,
  across the full 51-library table. This is a real, plausible finding for
  PacBio (barcode demultiplexing has much higher specificity than Illumina
  index-hopping — a different mechanism, not comparable to Illumina hop-rate
  expectations), not a bug in the check. Directly answers the question this
  investigation started from: Akkermansia's sediment abundance sits above
  a floor of 0, i.e. this check gives no evidence for or against it being
  cross-talk — its presence has to be evaluated on its own taxonomic/
  ecological merits, not waved away as bleed.

- **`conductivity_ms` is too sparse (19.6% complete) to be a primary
  predictor.** Including it in `07_environmental_drivers.R`'s db-RDA drops
  casewise-complete N from 34/46 to 10/46, and at N=10 forward selection
  correctly retains no term (`RsquareAdj()` returns `numeric(0)` — guarded
  for explicitly, since `sprintf("%.3f", numeric(0))` silently produces an
  empty string instead of erroring). The primary model uses `depth_m` +
  `temperature_c` (N=34); conductivity gets a clearly-labeled secondary,
  reduced-N sensitivity run alongside it, not silently dropped.

- **`elevation_m` is `depth_m`'s exact complement** (`depth_m + elevation_m
  = 1535` constant, confirmed numerically, cor = -1) — excluded from every
  model, not just deprioritized in VIF, to avoid perfect collinearity.

- **`vegan::adonis2(..., by = "margin")` on a formula with an interaction
  silently drops the main effects.** `06_beta_diversity.R` originally fit
  `~ sample_type * depth_m` with `by = "margin"` and got back a table with
  *only* the interaction row (R²=0.027, marginal p) — not a bug, that's
  documented `adonis2` behavior (main effects inside an interaction aren't
  marginally estimable), but it silently hid a much stronger result: fit
  additively (`~ sample_type + depth_m`, `by = "margin"`), both main effects
  are highly significant on their own (R²=0.053/0.056, p=0.001/0.001).
  `results/06_permanova.tsv` now reports both: the additive model as
  primary (`model == "additive_marginal"`), the full interactive model
  sequentially (`by = "terms"`) alongside it so the interaction term
  (sample_type's effect depending on depth) is still visible, just not at
  the cost of hiding the main effects.

- **`vegan::plot.varpart()` clips long `Xnames`** (`"depth_m"`,
  `"temperature_c"` both ran off the device edge at default margins) —
  shortened to `"depth"`/`"temp"` for the plot only; full names stay in
  `results/07_varpart.tsv`. It also doesn't compose with a ggplot in the
  same device (base `layout()`/`print(<ggplot>)` leave it blank, and
  `grid::grid.grabExpr()` came back blank too on this specific plot) —
  `10_figures.R`'s Fig 3 embeds the already-rendered `plots/07_varpart.png`
  as a raster (`png::readPNG` + `grid::rasterGrob`) instead of recapturing it.

- **ALDEx2 and Maaslin2 substantially disagree** on sediment-vs-water
  differential abundance (q<0.05): ALDEx2 finds 0 of 15,460 ASVs
  significant, Maaslin2 finds 29, intersection 0. This is a real result, not
  a bug — ALDEx2's Monte-Carlo compositional-uncertainty approach is
  markedly more conservative, and N is unbalanced (36 sediment vs. 10
  water). Report both numbers; don't quietly pick the one with hits.

- **`11_faprotax.R` predicts functional guilds (nitrification, sulfate
  respiration, methanotrophy, fermentation, ...) from ASV taxonomy using
  FAPROTAX (Louca et al. 2016)** — the same tool
  `isd_archive_scripts/isd_crete_workflow.sh` used, run the same way (the
  official `collapse_table.py` against the official `FAPROTAX.txt`), but
  new to this project's own `00`-`10` pipeline. Needs Python + `numpy`,
  which the container's system `python3` has neither of (no `pip` either);
  the script builds and reuses a local venv + downloads FAPROTAX 1.2.12
  itself, both under `tools/faprotax/` (`.gitignore`'d — vendored tool, not
  source). `-n none` (raw summed reads per group) is used instead of
  FAPROTAX's own normalization, then relative abundance is computed the
  same way as `05_taxonomic_composition.R` (divide by `library_size`), so
  every script in this pipeline normalizes on the same convention.

  **FAPROTAX's database is written against pre-GTDB (SILVA/NCBI-style)
  names** (`Proteobacteria`, `Firmicutes`, ...); this project's taxonomy is
  GTDB-style (`Pseudomonadota`, `Bacillota`, ...). Measured, not assumed:
  **64.2% of reads (71.6% of ASVs) matched no functional group at all**
  (`results/11_faprotax_unassigned.tsv`, corroborated by
  `results/11_faprotax_report.txt`'s own tally). Spot-checked directly
  against `Nitrospira` (a real, non-trivial genus in the sediment samples,
  ~0.3–0.5% relative abundance per `results/05_genus_relabund.tsv`) — it
  lands in **no** FAPROTAX group, for two compounding reasons confirmed by
  reading `FAPROTAX.txt` itself: (1) FAPROTAX's nitrification entries name
  specific comammox species (`Nitrospira nitrosa/nitrificans/inopinata`),
  not the bare genus — DADA2 ASVs unresolved to species (`s__NA`, common
  here) can never match those regardless of naming convention; (2) the one
  broader rule is class-level (`*Nitrospirae*Nitrospiria*`) and uses the
  pre-GTDB phylum spelling `Nitrospirae`, which doesn't match this
  project's GTDB `Nitrospirota` string — a direct, confirmed instance of
  the naming-mismatch concern, not a hypothetical one. Read
  `results/11_faprotax_relabund.tsv` with this substantial, disclosed
  under-matching in mind — it undercounts guild membership, particularly
  for species-unresolved ASVs and any group defined above the genus level.

- **Replicate concordance found real disagreement, not just noise:** 14/18
  biological-replicate pairs (`sample_set` transect vs. isolate_source, same
  site + tech_rep) are discordant by the `02_qc_filter.R` criterion, vs. only
  4/23 technical (extraction) pairs. The two `sample_set` arms are not
  interchangeable at this site — treat that as a finding to report, not
  a QC problem to pool away. See `results/replicate_concordance_decision.tsv`
  (both replicate levels are **not pooled by default**) and `PLAN.md` §11.7.

## Layout

Everything executable lives in `scripts/`; there is no separate `workflow/`.

```
Containerfile              the image used for everything (layers on omics-16s.Containerfile)
fetch-16s-databases.sh     container-side DB helper (see caveat below)
scripts/                   build_inputs.R, run_pipeline.sh, nextflow.config,
                           then the 00–10 analysis pipeline of PLAN.md §4
data/                      PB482_SP (raw fastq), MIxS sheets, md5sums,
                           samplesheet.tsv + metadata.tsv (generated), processed/
HiFi-16S-workflow/         upstream Nextflow pipeline (separate git repo — do not edit)
isd_archive_scripts/       reference implementation from ISD Crete; read for idiom, don't run
results/ plots/            outputs
```

`fetch-16s-databases.sh` writes a **flat** `silva/ gtdb/ ncbi/
greengenes/` layout for the container's own `/opt/databases`. It is *not* what
the Nextflow pipeline consumes — that needs `<db>/nb/` + `<db>/vsearch/` and is
produced by `run_pipeline.sh download`. Keep the two straight.
