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
  15/51 samples fall below 499 reads and are excluded from the rarefied
  matrix only (`03_normalize.R`) — still present in relative-abundance/CLR.

- **CLR needs its own, coarser ASV filter.** The general QC prevalence
  filter (`02_qc_filter.R`, ≥2 samples) leaves 15,460 ASVs — far too sparse
  for compositional zero-replacement: `zCompositions::cmultRepl()`'s default
  thresholds (`z.warning`/`z.delete` = 0.8/`TRUE`) **silently dropped 33/51
  samples** the first time this ran, because most ASVs are zero in >80% of
  samples at that resolution. Fixed in `03_normalize.R` with a CLR-specific
  ≥10%-of-samples prevalence filter (717 ASVs) plus `z.delete = FALSE`
  explicitly, so a still-sparse sample (the low-diversity mocks, by design;
  the shallowest sediment libraries) gets a warning, never a silent drop.
  Check `setequal(rownames(counts_clr), rownames(counts_clean))` after any
  change here — it should always be `TRUE`.

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
