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

51 libraries = 23 sites × 2 technical replicates + 5 controls:

| Set | Sites | Libraries |
|---|---|---|
| `transect` sediment | C1–C9 | 18 |
| `isolate_source` sediment | C1_I–C9_I | 18 |
| `transect` water | W1–W5 | 10 |
| `control` | 5 | 5 |

Controls: `ctr_zymo_com`, `ctr_zymo_log`, `ctr_msa_3001` (mocks),
`ctr_EB` (extraction blank), `ctr_MM` (mastermix blank).

`scripts/build_inputs.R` (base R, runs on the host) generates
`data/samplesheet.tsv` and `data/metadata.tsv` from the barcode map
plus the two ENA/MIxS checklist sheets. Regenerate rather than hand-editing —
the pipeline's `inspect_metadata` process diffs the two sample-id columns and
dies on any mismatch.

**One row per library, not per site.** Technical replicates stay separate
through the pipeline so replicate concordance can be measured before it is
averaged away; pool them downstream, deliberately.

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
   sediments — do not pool those two arms without checking they agree.

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

## Two environment traps that cost real time

Both produce silent hangs or misleading errors; both are already fixed in
`Containerfile` / `scripts/run_pipeline.sh`, documented here so they are not
re-introduced.

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

3. **`conf/base.config`'s `highparallel` label caps `time` at 8h.** The four
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

**W5's conductivity and temperature are transposed in the source sheet**:
`conductivity_ms = 8.4`, `temperature_c = 195.9`. 195.9 °C is impossible in a
cave stream, and every other water sample sits near 188 mS/cm and 6.6–7.9 °C.
`build_inputs.R` carries the values through **verbatim** rather than silently
correcting them. Handle it explicitly in `02_qc_filter.R` and say so in the
output; do not quietly swap them upstream.

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
