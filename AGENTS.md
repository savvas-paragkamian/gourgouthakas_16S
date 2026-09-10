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

**Read/amplicon length**: full-length 16S rRNA (V1-V9), not a short-read
amplicon like V4. `HiFi-16S-workflow/nextflow.config` filters DADA2 input to
`min_len = 1000` / `max_len = 1600` (`main.nf --min_len`/`--max_len`
defaults). Realized post-filter ASV length distribution (30,477 ASVs,
`taxonomy_postlineage.rds`): median 1,453 bp, IQR 1,437-1,472 bp, 99% range
1,389-1,532 bp — centered on the ~1,500 bp expected for a full bacterial 16S
gene, as expected for PacBio HiFi/Kinnex circular consensus reads (not the
raw subread length, which is longer per pass but irrelevant post-CCS).

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

- **Library size is reported as a covariate (or an above-floor
  restriction) in every community-level model, not just visualized.**
  Motivation is a confirmed, not hypothetical, confound: the
  replicate-concordance check (`02b_controls.R`) found technical-pair
  Bray-Curtis distance strongly anti-correlated with the pair's minimum
  library size (Spearman r = -0.86 — shallow libraries look artificially
  dissimilar), and sediment/water differ systematically in sequencing
  depth too, so a `sample_type`/`depth_m` effect could in principle be a
  library-size effect in disguise. `metadata_clean.rds` carries
  `library_size`/`log_library_size` (log10) as first-class columns,
  computed once in `02b_controls.R` from the final matrix, name-joined —
  every script below reads the same numbers rather than re-deriving them.
  "Above-floor" restriction always means the *same* threshold
  `03_normalize.R` already uses for rarefaction eligibility (currently 499
  reads, `data/processed/rarefaction_excluded_samples.rds`), not a second,
  competing floor. Per-script treatment:
  - `06_beta_diversity.R` PERMANOVA: two extra models beyond the existing
    additive/interaction pair — `additive_marginal_libsize` (full N=46,
    `+ log_library_size` as its own covariate) and
    `additive_marginal_restricted` (above-floor N=34, no covariate needed —
    the restriction itself removes the confound). Result: `log_library_size`
    is itself significant (R²=0.057, p=0.001), but `sample_type`/`depth_m`
    barely move under adjustment (0.052→0.053, 0.056→0.055) and get
    *stronger* under restriction (0.067/0.081) — the main effects are
    robust, not artifacts.
  - `07_environmental_drivers.R`: a `libsize_adjusted` db-RDA variant
    (`Condition(log_library_size)`, vegan's standard covariate-partialling
    mechanism for constrained ordination — R²=0.131 vs. primary's 0.141,
    barely moved) and a library-size partial Mantel (`community ~ depth |
    log_library_size`: r=0.280 vs. unpartialled r=0.272, if anything
    slightly stronger).
  - `08_differential_abundance.R`: Maaslin2 gets `log_library_size` added
    to `fixed_effects` (natively supported, standard practice for this
    tool) — `results/da_maaslin2_all_terms.tsv` has both terms' full
    results, `da_maaslin2_sample_type.tsv` stays filtered to the
    `sample_type` contrast only (mixing the two terms' significant hits
    together would have been a real bug, caught before it shipped). ALDEx2
    does **not** get the covariate: its CLR transform divides each
    sample's counts by that same sample's own geometric mean before
    logging, so it's compositionally invariant to total library size by
    construction — adding one would mean switching the primary DA method
    to `aldex.glm()` with a full design matrix, a materially bigger
    rewrite for a confound this method doesn't actually have.
  - `04_alpha_diversity.R`: no change. Alpha diversity runs on rarefied
    counts (every included sample equalized to exactly 499 reads), so
    library size has zero variance in that set by construction and can't
    be a covariate there — rarefaction *is* the library-size control for
    this metric, not a gap.

- **Mock cross-talk index (`02b_controls.R` T6,
  `results/02b_crosstalk_index.tsv` + `02b_crosstalk_index.png`, plus
  `02b_mock_recovery_species.tsv`/`02b_mock_crosstalk.tsv`) — the measured
  index-hopping/barcode-cross-talk rate, and what it does and doesn't
  license claiming about species-level resolution.** Built from the three
  mocks (`data/mock_expected.tsv` gives near-mutually-exclusive membership
  across them, confirmed: 6 genera exclusive to `mock_env`, 3 exclusive to
  `mock_even`+`mock_log`), not asserted from general PacBio literature.
  - **Mock-to-mock cross-talk (the direct index-hopping estimate)**: reads
    in one mock's well assigned to a genus expected only in a *different*
    mock. Measured **0.09-0.21%** of a mock's reads (mock_even 0.09%,
    mock_log 0.13%, mock_env 0.21%) — this is the concrete number behind
    "PacBio barcode cross-talk is normally low": here it *is* low, and now
    it's a citable rate rather than an assumed one. The Listeria-in-mock_env
    false positive visible in `02b_mock_recovery.png` is exactly this
    signal; T6 is its formalization across all three mocks, not a new
    finding.
  - **Species-level recovery is 19% (5/26 expected genus×mock pairs
    resolved to the correct species), vs. 92% at genus level (T2).** This
    is the number that actually bears on the "full-length reads give
    species-level resolution" claim — genus-level recovery is strong,
    species-level recovery with this classifier (`gg2_nb`, naive Bayes
    against GTDB-backed Greengenes2) is not, on this dataset. Of the 26
    expected pairs, 21 are "not detected" (the ASV classifies confidently to
    genus but `gg2_nb` doesn't commit to a species) and the *other* mock
    ASVs contribute 54 "false positive" species rows — nearly all
    same-genus/different-species misclassification (e.g. `Escherichia_coli`
    reads landing on `Escherichia_albertii`/`boydii`/`sonnei`,
    `Listeria_monocytogenes` on `Listeria_A_marthii`), **not** cross-mock
    contamination — that signal is already isolated separately in the
    mock-to-mock cross-talk number above. **Practical implication: claim
    species-level taxonomic *reads* (full-length 16S genuinely resolves
    more than V4 short reads) but not species-level *classifier accuracy*
    with this NB/GTDB pipeline on this dataset** — any species-level claim
    in the writeup should cite the 19% figure, not assume full-length =
    accurate species calls.
  - **Chimera rate**: mocks 2.5-6.1% (mock_even 2.53%, mock_env 5.74%,
    mock_log 6.07%) vs. real samples 2.55% and extraction/mastermix blanks
    3.20% — mocks are in line with the rest of the run, not unusually
    clean or dirty; chimera removal isn't hiding anything mock-specific.
  - **ASV excess rate** (T3's excess-ASV count as a % of ASVs observed):
    mock_env 5.7%, mock_even 34.4%, **mock_log 93.1%** (578 ASVs observed
    vs. ~16-40 expected). This is the number that directly motivates the
    genus-level-primary decision below — it is *not* cross-talk (mock-to-
    mock cross-talk for mock_log is only 0.13%): it's singletons and rRNA-
    operon copy-number variants of the correctly-identified genera, exactly
    as expected for a 16S-only ASV table with no copy-number correction at
    the ASV level.
  - **`read_method_tax()` bug fixed while building this**: its genus/species
    extraction regex (`sub("^.*g__(...)...")`) silently fell through to the
    *raw, full taxon string* — not `NA` — whenever a rank prefix was absent
    entirely (common: `silva_vsearch`'s Taxon field has no `g__`/`s__`
    prefixes at all; `gg2_nb` genus-classifies but doesn't always
    species-classify an ASV). This corrupted the species-level recovery
    table (raw taxon strings appearing as "species") and, downstream,
    crashed `score_mock_method()` for `silva_vsearch` (`stats::aggregate()`
    returns a `count` column typed `list()`, not `numeric`, when its
    grouping vector is *entirely* `NA` — `sum()` then errors on that list).
    Fixed at the source (`stringr::str_extract()`, which returns `NA` on no
    match, instead of base `sub()`'s silent no-op) plus an explicit
    NA-drop-before-aggregate guard everywhere `methods$*$genus`/`$species`
    feeds `stats::aggregate()`.

- **Depth-honest replicate concordance (`02b_controls.R`,
  `02_replicate_concordance.png` + `results/02_replicate_vs_null_band.tsv`,
  `02_replicate_null_band_summary.tsv`, `02_replicate_floor_by_level.tsv`)
  — an empirical, measured read-count floor, not the ad hoc 500/499 already
  in use. STOP-AND-ASK finding: read before touching either existing
  threshold.** Two changes to the original (now superseded)
  `02_replicate_concordance.png`:
  1. **Pairwise-rarefied Bray-Curtis** (`bray_rarefied`, added alongside the
     existing full-depth `bray` in `results/replicate_concordance_*.tsv`):
     both samples in a pair are repeatedly (30x) subsampled down to their
     own shared `min(lib_a, lib_b)` via `vegan::rrarefy()` before computing
     Bray-Curtis, and the result averaged — so a 24-read vs. 400,000-read
     "pair" is compared on equal footing at the depth it can actually
     achieve, not penalized for the shallow side's sampling noise alone.
  2. **Subsampling null band**: the single deepest real sediment sample and
     deepest real water sample are each repeatedly (50x per depth,
     log-spaced from 20 to 400,000 reads) split into two INDEPENDENT
     `rrarefy()` subsamples, and Bray-Curtis computed between them — a null
     distribution of "how dissimilar do two subsamples of ONE identical
     community look, from pure sampling noise alone, at depth N." Plotted
     as a median+10th-90th-percentile ribbon under the real within-/
     between-pair points, now positioned on that same depth axis by each
     pair's own `min_lib` (the categorical technical/biological x-axis is
     gone).
  - **Result — genuinely different stories for technical vs. biological
    pairs, not one number**:
    - **Technical pairs (`02_replicate_floor_by_level.tsv`: 4/23
      indistinguishable from noise, all shallow) split cleanly by depth.**
      Below ~8,000-8,300 reads, a technical pair's Bray-Curtis distance
      cannot be told apart from what two subsamples of the *same* community
      would produce by chance alone — "this pair looks discordant" and
      "this pair is just shallow" are not distinguishable claims down here.
      Above that depth, pairs separate from the (rapidly shrinking) null
      band, but their own absolute `bray_rarefied` keeps *falling* as depth
      increases (e.g. one technical pair: 0.44 at 8,122 reads → 0.10 at
      467,624 reads) — the "reproduce above it" half of the story: real
      reproducibility, visible once sampling noise is no longer big enough
      to hide it, not a failure mode that shallow sampling was masking.
    - **Biological pairs (0/18 indistinguishable from noise, at ANY depth
      tested, including the shallowest — 24 reads) are a different story
      entirely, not a depth-threshold one.** Their `bray_rarefied` (0.7-1.0)
      already clearly exceeds the null band at the very shallowest depths
      in the dataset, where the band is at its widest and most forgiving —
      there's no depth at which more sequencing would have changed this
      verdict. This is *expected*, not a data-quality problem: a
      "biological pair" here is same `location`+`tech_rep` but a different
      arm (transect vs. isolate_source) — deliberately different physical
      samples, not a duplicate extraction of the same one. Real, immediate
      spatial heterogeneity between arms, not evidence replication failed
      or a case the read floor below is relevant to.
  - **Combined empirical floor (`results/02_replicate_floor_by_level.tsv`
    has the technical/biological split) — SUPERSEDED, see below.** Originally
    measured at ~8,000-8,300 reads, computed on `counts_filtered` (the
    post-prevalence-filter matrix). Once Part 3a (below) moved this whole
    computation earlier — onto `counts_noctrl`, pre-prevalence-filter, so
    its own result could gate that filter without circularity — the same
    method measured **in the tens of reads (observed 24-78 across repeated
    runs, Monte Carlo noise in which 1-2 borderline pairs happen to fall
    inside the shrinking null band — not a precise fixed number), not
    ~8,272.** This drop wasn't noise in the "ignore it" sense, though: the
    original ~8,272 figure was itself an artifact of running on an
    already-prevalence-filtered matrix (the exact ordering problem Part 3a
    exists to fix) — on the full, unfiltered ASV set, the deepest
    previously-"indistinguishable" technical pair (`C9`, 8,272 reads) turned
    out to show real discordance once every ASV was counted, not just the
    ones that happened to survive prevalence filtering (`bray_rarefied`
    rose from 0.183 to 0.28-0.29 at nearly the same depth across repeated
    runs, while the null band's ceiling at that depth barely moved). **Per
    direct instruction, the prevalence filter's eligibility floor (below) is
    set manually to 8,000 reads rather than tracking this recomputed,
    unstable few-dozen-reads figure** — that figure is the more
    methodologically honest self-vs-self noise floor, but (a) a filter gated
    there excludes almost no libraries (2/36 sediment, 0/10 water) and (b)
    its own run-to-run instability is a second, independent reason not to
    build a hard filter threshold on it directly. Both numbers are real and
    both are reported (`empirical_read_floor` = the noisy few-dozen-reads
    figure, exact value varies by run, in the code and
    `results/02_replicate_floor_by_level.tsv`; `prevalence_eligibility_floor`
    = 8,000, the manual value actually used to gate the prevalence filter)
    — they measure genuinely different things (self-vs-self noise floor vs.
    a chosen conservative cutoff) and shouldn't be conflated.
  - Null-band mechanics sanity-checked with `stopifnot()`: the sediment
    reference's null median must fall as depth increases and exceed 0.05 at
    the shallow end — both hold.

- **Prevalence filter reworked (`02b_controls.R`, Part 3a of
  `golden-napping-breeze.md`) — within-type, read-count-aware, floor-gated,
  replacing a pooled `prevalence >= 2 of 46 samples, >0 reads` rule with
  three independent, real problems.** Direct objection, not a style
  preference: the old rule let a taxon pass on e.g. "1 sediment read + 1
  water read" (sediment and water share almost no taxa — a pooled count is
  not ecologically meaningful), treated a single read as "present" (Part 1
  measured mock-to-mock cross-talk directly in this dataset at 0.09-0.21% of
  reads — a lone read is indistinguishable from that), and let libraries
  with 24-200 reads vote on prevalence before any read floor had been
  established. New rule: an ASV is kept if it has ≥2 reads in ≥2 samples of
  the SAME type, counting only samples at/above the 8,000-read eligibility
  floor (`prevalence_eligibility_floor` — manual, see the empirical-floor
  note above for why this isn't the noisy, recomputed few-dozen-reads
  figure). Result: **14,453/30,430
  ASVs kept, down from 15,455/30,430 under the old rule** — a real,
  reported reduction, not a rounding difference. Eligible-sample counts per
  type are asymmetric by design/data: sediment 17/36, water only 2/10 (an
  honest reflection of water's much shallower sequencing in this dataset,
  not softened for that reason — see the eligibility-floor decision above).
  The entire replicate-concordance diagnostic block (pairwise-rarefied
  Bray-Curtis, subsampling null band) was relocated earlier in the script to
  run on `counts_noctrl` instead of `counts_filtered`, specifically so its
  own `empirical_read_floor` output could feed this filter without
  circularity (compute the floor on a matrix, then use that floor to build
  the SAME matrix, would have been backwards).

- **Intragenomic-variant test (`02b_controls.R`, Part 3b of
  `golden-napping-breeze.md`, `results/02b_intragenomic_variant_test.tsv` +
  `_summary.tsv` + `plots/02b_intragenomic_variant_test.png`) — the direct
  evidence for genus-level aggregation (Part 3c), not an assumption.**
  Question: do ASVs sharing an identical species-level GTDB label
  (`species_norm`) co-occur across samples with correlated abundance (the
  signature of intragenomic rRNA-operon copies of ONE genome, present in a
  fixed ratio) more than random cross-species-label ASV pairs do? Run on the
  Part-3a-corrected `counts_filtered`/`taxonomy_filtered` (not the earlier
  reconnaissance numbers, which were computed on the uncorrected table and
  are superseded). **Result: yes, clearly.** Same-species-label pairs
  (44,063 pairs, 1,144 groups sharing a label, capped at 20 ASVs/group):
  median Jaccard co-occurrence 0.40 vs. a cross-species null's 0.00; median
  Spearman abundance correlation 0.60 vs. the null's -0.05. Both Wilcoxon
  p≈0, rank-biserial effect ≈0.40-0.41 (positive = same-species pairs
  stochastically exceed the null, the direction the hypothesis predicts).
  **Not a clean unimodal signal, though** — the figure shows both
  distributions are bimodal, and same-species-label pairs still have real
  mass down near 0 on both metrics: a meaningful fraction of same-label ASV
  pairs behave like independent organisms (a coarse nearest-reference label
  shared by genuinely distinct taxa), not like intragenomic copies of one
  genome. Case study `GMQP-bins7_sp004366385` (the label pointed at
  directly): 180 ASVs pre-Part-3a collapsed to 20 post-3a (capped for the
  pairwise computation) — a concrete illustration of how much the
  uncorrected prevalence filter had been letting through. **Conclusion for
  Part 3c: the evidence supports genus-level aggregation as reducing real
  ASV-level redundancy, not manufacturing an assumption — but "mostly
  intragenomic variants" would overstate it; "a substantial, statistically
  clear excess of same-species co-occurrence over what distinct organisms
  sharing a label would produce, alongside a real population of
  genuinely-distinct co-labeled organisms" is the accurate summary.**

- **Genus-level is the primary unit for core community structure
  (`02b_controls.R`/`03_normalize.R`/`04_alpha_diversity.R`/
  `06_beta_diversity.R`/`07_environmental_drivers.R`, Part 3c of
  `golden-napping-breeze.md`) — ASV-level kept as supplementary, not
  dropped.** Built from the Part-3a-corrected ASV table
  (`counts_filtered`/`taxonomy_filtered`, not the earlier uncorrected one)
  and justified by Part 3b's measured evidence, not assumed.
  - **`02b_controls.R`**: `counts_genus_clean.rds`/`taxonomy_genus_clean.rds`
    — ASVs summed by `genus_norm`; genus-NA ASVs excluded (821/14,453 ASVs,
    4.587% of reads — reported, not silently dropped). **46 samples x 1,406
    genera, ASV:genus ratio 10.3:1.**
  - **`03_normalize.R`**: genus-level rarefied/relabund/CLR
    (`counts_genus_rarefied.rds`/`_relabund.rds`/`_clr.rds`), same 499-read
    depth and same re-derived (>=10% of samples) CLR-prevalence logic as
    ASV-level, but re-applied to the genus matrix's own column count and
    sparsity, not reusing the ASV-level threshold's absolute number.
    `rarefaction_excluded_samples_genus.rds` is its OWN list, recomputed
    from the genus matrix's own row sums, not reused from the ASV-level
    exclusion list — genus-level totals are lower (genus-NA reads excluded),
    so the same depth excludes a different set (13/46 samples, vs. 12-13/46
    at ASV-level depending on the exact run — one additional sample,
    `C3_2_A`, dips below the floor only once genus-NA reads are removed).
  - **`04_alpha_diversity.R`/`06_beta_diversity.R`**: refactored into one
    function each (`run_alpha_diversity()`/`run_beta_diversity_level()`),
    called twice — genus-level first, writing the PRIMARY, unsuffixed
    result files (`alpha_diversity.tsv`, `06_permanova.tsv`,
    `dist_bray.rds`, etc. — the same names these scripts have always used);
    ASV-level second, writing `*_asv`-suffixed supplementary files. This
    means every downstream consumer of those unsuffixed names (`10_figures.R`,
    `09_spatial_analysis.R`, `07_environmental_drivers.R`) automatically
    started reading genus-level data with **zero code changes in those
    three scripts** — confirmed by re-running each standalone.
  - **`07_environmental_drivers.R`, `09_spatial_analysis.R`,
    `10_figures.R`: genuinely unchanged**, per the point above.
  - **`08_differential_abundance.R`, `11_faprotax.R`, `05_taxonomic_composition.R`:
    unchanged, deliberately** — 08 stays ASV-level (scope decision,
    `golden-napping-breeze.md` Context section: doubling ALDEx2's ~4-5 min
    runtime for a second full model wasn't judged worth it against the
    payoff); 05 and 11 already operate at genus level in practice.
  - **06's dissimilarity heatmap (`06_dissimilarity_heatmap.png`) uses the
    genus-level (primary) Bray distance only** — not duplicated at ASV
    level, judged modest incremental value for a fairly heavy
    pheatmap+PNG-composition diagnostic. Its library-size color annotation
    stays ASV-level raw read counts (`counts_clean.rds`) regardless — that's
    the actual sequencing depth per sample, not something that changes
    meaning at genus level.
  - **Result direction, for context** (not the headline — Part 3b's
    evidence is): genus-level R² is consistently a bit higher than the
    ASV-level supplementary result across every model checked (PERMANOVA
    Bray additive: sample_type 0.076 vs. 0.057, depth_m 0.074 vs. 0.056;
    db-RDA primary: 0.183 vs. 0.145) — consistent with genus collapse
    reducing ASV-level noise rather than just discarding signal.

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
  candidate list, not automated removal** — a `decision` column that a human
  fills in and the script round-trips on every subsequent run into
  `results/02b_contaminant_blacklist.tsv`. **Do not read an empty blacklist
  as a clean bill of health**: `ctr_EB` (41 raw reads) and `ctr_MM` (70 raw
  reads) are shallow negatives, so contamination screening here is limited
  by control read depth, not exhaustive.

  **Reviewed** (all 72 candidates; `results/02b_contaminant_candidates.tsv`
  and `_blacklist.tsv`): cross-referenced against `data/mock_expected.tsv`
  first — 40/72 hit a Zymo/ATCC mock genus, confirming most of this list is
  PacBio barcode bleed from the mock wells into `ctr_EB`/`ctr_MM`, not a
  reagent contaminant reaching real samples. Decision rule, applied
  uniformly rather than judged organism-by-organism (not confident enough
  in every genus's soil/cave ecology to override the data): **remove** if
  absent from every real sample (`prevalence_real_samples == 0`, 32/72) or
  present at only trace level (`max_relabund_real_samples < 0.001`, 4 more
  — Staphylococcus 0.014%, Parafrigoribacterium 0.04%, Lactobacillus
  0.002%, Enterenecus 0.09%), **1** other trace-level mock-genus row also
  removed for the same reason; **keep** everything else (25/72) — the
  observed multi-sample real abundance outweighs a statistically thin n=2
  test. This explicitly keeps `Akkermansia` (33/46 real samples, up to 24%
  relative abundance — the organism this whole controls investigation
  started from) and `Parabacteroides` (21/46, up to 17%): both flagged by
  decontam on blank-prevalence pattern alone, both overwhelmingly real by
  direct observation. 47 ASVs blacklisted; final gate removed
  15457→15455/30430 ASVs (most of the 47 were already going to fail the
  ≥2-sample prevalence filter regardless, being absent from real samples).

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

- **Replicate concordance found real disagreement, not just noise:** 16/18
  biological-replicate pairs (`sample_set` transect vs. isolate_source, same
  site + tech_rep — sediment only, water has no second `sample_set` arm)
  are discordant by `02b_controls.R`'s criterion (a within-pair Bray
  distance at or above the 5th percentile of the between-pair reference
  distribution), vs. 5/23 technical (extraction) pairs. The two
  `sample_set` arms are not interchangeable at this site — treat that as a
  finding to report, not a QC problem to pool away. The between-pair
  reference is restricted to pairs of the **same `sample_type`** (a
  sediment-vs-water comparison is trivially, definitionally dissimilar and
  doesn't belong in the noise floor a same-type replicate pair gets judged
  against) **and the discordance threshold itself is computed per
  `sample_type`**, not pooled across sediment+water — both refinements
  moved the counts from an original 3/23 and 14/18, to 5/23 and 15/18 with
  same-type-only between-pairs, to the current 5/23 and 16/18 once the
  threshold was also split by type.
  `plots/02_replicate_concordance.png` is two side-by-side panels
  (sediment | water, patchwork, not `facet_wrap`) with an explicit
  within-pair/between-pair definition in the caption; water's panel
  correctly has no "biological" x-tick at all (water was never eligible
  for biological pairing, so there's no within-pair to judge a between-pair
  distribution against there — an earlier version of this figure showed an
  orphaned between-only "biological" box for water, which was wrong, not
  just unlabeled). See `results/replicate_concordance_decision.tsv` (both
  replicate levels are **not pooled by default**) and `PLAN.md` §11.7.

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
