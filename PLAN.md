# Analysis Plan — Downstream R Analysis of HiFi 16S Soil Microbiome
### tidyverse + vegan native — no phyloseq, no qiime2R

A plan for Claude Code to scaffold and implement an R analysis that consumes the
**HiFi-16S-workflow** outputs (PacBio HiFi full-length 16S → DADA2 ASVs) and runs
a soil-microbiome ecology analysis in the style of **crete-soil-health** —
describing microbial diversity and identifying the environmental drivers behind
it. The data model is plain tidy tables + numeric matrices + an `ape` tree.
No `phyloseq` S4 objects and no `qiime2R`/`.qza` imports anywhere.

---

## 1. Objective

A reproducible, script-driven R pipeline that:

1. Imports flat ASV outputs into a small set of canonical tidy objects.
2. Cleans, filters, and normalizes the ASV count matrix.
3. Characterizes **alpha diversity**, **taxonomic composition**, and **beta diversity**.
4. Identifies **environmental drivers** of community structure (constrained ordination, PERMANOVA, variance partitioning) — the analytical core.
5. Runs **differential abundance** across conditions/gradients.
6. Adds a **spatial** layer (sample maps, distance–decay) for island / landscape sampling.
7. Emits publication-quality figures to `plots/` and tables/objects to `results/`.

Each script runs stand-alone and in sequence, reads only from `data/` and prior
`results/`, and writes only to `results/` and `plots/`. Relative paths only; no
manual steps.

---

## 2. Inputs (flat files — qiime-free)

Land these in `data/raw/`. They are the exported/native products of the
workflow (or of DADA2 directly). Read them in R with `readr`, `biomformat`, and
`ape` — **not** `qiime2R`, and never as `.qza`.

| File | What it is | Reader |
|------|------------|--------|
| `feature-table.tsv` **or** `feature-table.biom` | ASV × sample counts | `readr::read_tsv` / `biomformat::read_biom` |
| `taxonomy.tsv` | Per-ASV lineage (SILVA/GTDB/GG2) + confidence | `readr::read_tsv` |
| `tree.nwk` | Rooted phylogeny (UniFrac / Faith's PD) | `ape::read.tree` |
| `dna-sequences.fasta` | ASV representative sequences (optional) | `Biostrings::readDNAStringSet` |
| `metadata.tsv` | Sample metadata (`sample_id`, `condition`, …) | `readr::read_tsv` |
| `soil_chemistry.tsv` | pH, moisture, C/N, texture, nutrients, coordinates | `readr::read_tsv` |

> If the workflow only gives you `.qza`, export **once at the shell** with
> `qiime tools export` (or unzip and grab `data/`), so the R code stays
> qiime-free. A `.biom` is read directly in R via `biomformat` with no QIIME.
> If you run DADA2 yourself, the `seqtab` matrix and `assignTaxonomy` output are
> already flat R objects — feed those in directly.

---

## 3. Canonical data model (replaces the phyloseq object)

`01_import` produces these and saves each to `data/processed/*.rds`. Everything
downstream operates on them:

- `counts` — **integer matrix, samples as rows, ASVs as columns** (vegan's expected orientation; flag this explicitly, it's the classic bug).
- `counts_long` — tidy tibble `sample_id, asv_id, count` for ggplot/dplyr work.
- `taxonomy` — tibble: `asv_id` + `domain, phylum, class, order, family, genus, species` (split from the lineage string) + `confidence`.
- `metadata` — tibble: `sample_id`, `condition`, soil-chemistry columns, `lon`/`lat`, plus the replicate-structure columns needed for §7's replicate-concordance step: `site` (the depth/location identifier shared across replicates, e.g. `C1`), `bio_rep` (biological replicate arm — e.g. `C1` vs `C1I`; parsed from `sample_set` in the upstream metadata), and `tech_rep` (technical/extraction replicate — the `_1`/`_2` suffix on `sample_id`). Parse these once in `01_import.R` rather than re-deriving them from `sample_id` strings in later scripts.
- `tree` — `ape::phylo`, pruned to the ASVs present in `counts`.
- (optional) `refseqs` — `DNAStringSet`.

Keep `counts`, `taxonomy`, `metadata`, and `tree` **aligned** on the same ASV
and sample ids at every step; write a tiny `assert_aligned()` helper and call it
after each filtering operation.

---

## 4. Directory structure

```
crete-soil-health-analysis/
├── README.md
├── AGENTS.md                 # working agreement + conventions for Claude Code
├── renv.lock                 # pinned package versions
├── run_all.R                 # sources scripts 00–10 in order
├── Makefile                  # optional: per-target incremental rebuilds
├── data/
│   ├── raw/                  # flat inputs (read-only)
│   └── processed/            # counts.rds, taxonomy.rds, metadata.rds, tree.rds, dists, clr…
├── scripts/                  # (== src/) numbered, ordered steps
│   ├── 00_setup.R
│   ├── 01_import.R
│   ├── 02_qc_filter.R
│   ├── 03_normalize.R
│   ├── 04_alpha_diversity.R
│   ├── 05_taxonomic_composition.R
│   ├── 06_beta_diversity.R
│   ├── 07_environmental_drivers.R
│   ├── 08_differential_abundance.R
│   ├── 09_spatial_analysis.R
│   ├── 10_figures.R
│   └── functions.R           # shared helpers, sourced by 00_setup.R
├── plots/                    # figures (.pdf vector + .png raster)
└── results/                  # tables (.tsv/.csv), stats objects (.rds), session_info
```

Use `scripts/` (treat `src/` as a synonym; keep one). `data/raw/` is immutable.

---

## 5. Environment & reproducibility

- R ≥ 4.3. Write `sessionInfo()` to `results/session_info.txt` at the end of `run_all.R`.
- `renv` lockfile. Bioconductor packages (`biomformat`, `ALDEx2`, `DESeq2`, `decontam`, `Biostrings`) via `BiocManager`.
- Global `set.seed(42)` in `00_setup.R` (rarefaction, NMDS, permutations).
- Every script begins with `source("scripts/00_setup.R")` (paths, libraries, seed, ggplot theme, palettes, helpers).
- **Packages:** `tidyverse`, `vegan`, `ape`, `picante` (Faith's PD), `GUniFrac` (UniFrac), `compositions`/`zCompositions` (CLR), `decontam`, `Maaslin2` + `ALDEx2` (differential abundance), `DESeq2` (optional cross-check, matrix interface), `sf`, `terra`, `patchwork`, `ggpubr`, `ggrepel`, `scales`, `here`, `biomformat`, `Biostrings`.

None of these requires `phyloseq`. `decontam`, `ALDEx2`, and `DESeq2` all accept
plain matrices; `Maaslin2` takes plain data frames.

---

## 6. Conventions

- **Naming:** `snake_case`; scripts zero-padded and ordered.
- **Paths:** relative via `here::here()` / constants in `00_setup.R`. No absolute paths.
- **Matrix orientation:** samples = rows for all vegan calls; assert it.
- **I/O contract:** a script reads `data/` and earlier `results/`; writes only `results/` + `plots/`. Heavy intermediates → `data/processed/*.rds`.
- **Figures:** save `.pdf` (vector) + `.png` (preview) with explicit `width`/`height`/`dpi`; one shared `theme_set()`.
- **Tables:** tidy, one `write_tsv()` per result.
- **Style:** tidyverse style guide; comment the *why*.

---

## 7. Pipeline — script by script

### `00_setup.R`
Libraries, seed, path constants, shared ggplot theme, palettes (taxa/conditions),
`source("scripts/functions.R")`. Fail loudly on any missing package.

### `01_import.R`
Read the flat files (§2) → build the canonical objects (§3). Split the taxonomy
lineage string into ranks; strip rank prefixes (`d__`, `p__`, …). Join soil
chemistry into `metadata` by `sample_id`, reporting key mismatches. Prune `tree`
to observed ASVs (`ape::keep.tip`). Save all objects to `data/processed/`; write
a sample × variable completeness table to `results/`.

### `02_qc_filter.R`
- Drop non-target lineages: mitochondria, chloroplast, Eukaryota, ASVs unassigned at phylum (string-filter on `taxonomy`, then subset `counts`).
- Library sizes: ordered depth plot; flag/remove very shallow samples.
- If negative/extraction controls exist, run **`decontam::isContaminant`** on the count matrix (prevalence and/or frequency) and remove contaminants.
- Prevalence filter (ASV in ≥ 2 samples) + optional low-count filter.
- **Replicate concordance, at both levels present in this design — do this before any pooling decision, and before any other script pools or averages replicates:**
  - Sample IDs carry **two nested levels of replication**: a **biological** pair per site (`bio_rep`, e.g. `C1` vs `C1I` — two independent sediment samples at the same site) and, within each of those, a **technical** pair of DNA extractions (`tech_rep`, e.g. `C1_1` vs `C1_2`). The ISD reference scripts (`isd_archive_scripts/isd_crete_numerical_ecology.R`) only ever check the biological level (`loc_1`/`loc_2` pairwise dissimilarity) — deliberately extend that here to also check the technical level, since it's the finer-grained and more diagnostic of the two (a bad extraction shows up here first).
  - For each level, compute a **within-pair Bray–Curtis (and Jaccard, for presence/absence) distance** on relative abundance — same distances used downstream in `06_beta_diversity.R`, just computed early enough to gate pooling — and compare it against the distribution of **between-pair** (different site) distances at the same level. Report both as a table (`results/replicate_concordance_technical.tsv`, `results/replicate_concordance_biological.tsv`: `site`, `pair_id`, `bray`, `jaccard`, `level`) and a paired box/strip plot (within-pair vs. between-pair distance) to `plots/`.
  - Also report simpler concordance signals per pair: correlation of ASV relative abundances (Spearman), and ΔASV richness / Δlibrary size — cheap sanity checks that catch a failed extraction independent of the distance metric.
  - Flag (don't silently drop) any pair whose within-pair distance falls inside the between-pair distribution — that's a replicate that didn't reproduce, and `03_normalize.R` onward should know about it via a `flag_discordant_pair` column on `metadata` rather than have it disappear into an average.
  - This step's output is a **decision, not just a diagnostic**: record in `results/replicate_concordance_decision.tsv` whether technical replicates are pooled (mean/sum counts per `bio_rep × site`) before `03_normalize.R`, and whether the two `bio_rep` arms are pooled or kept as a deliberate contrast downstream — either way, keep both the pooled and unpooled count matrices in `data/processed/` so the decision is reversible.
- Re-align all objects; `assert_aligned()`. Save `*_clean.rds`; write a read/ASV retention table.

### `03_normalize.R`
- **Rarefaction curves** (`vegan::rarecurve`) → `plots/`.
- Produce three matrices, explicit about where each is used:
  - **rarefied** (`vegan::rrarefy` at a documented depth) → alpha diversity,
  - **relative abundance** (`sweep`/`decostand(method="total")`) → composition & unconstrained ordination,
  - **CLR** (`compositions::clr` after a zero-treatment via `zCompositions::cmultRepl`) → Aitchison beta diversity & DA input.
- Note the rarefaction-depth rationale in `results/`.

### `04_alpha_diversity.R`
- `vegan::specnumber` (Observed), `vegan::diversity` (Shannon, Simpson), Pielou evenness (derive), **Faith's PD** via `picante::pd(counts, tree)`.
- Test across `condition` (Kruskal–Wallis / Wilcoxon or LMs); regress metrics on continuous soil variables.
- Outputs: `results/alpha_diversity.tsv`, `results/alpha_stats.tsv`; boxplots + scatter-vs-gradient to `plots/`.

### `05_taxonomic_composition.R`
- Aggregate counts to Phylum and Genus with dplyr (`counts_long` ⨝ `taxonomy` → `group_by` → `summarise`); relative abundances.
- Stacked barplots of top-N taxa by `condition` (remainder → "Other").
- **Core taxa**: prevalence/abundance thresholds via dplyr → membership table.
- Outputs: composition tables to `results/`; barplots + core plot to `plots/`.

### `06_beta_diversity.R`
- Uses whichever count matrix `02_qc_filter.R` decided on (pooled or unpooled technical/biological replicates — see `results/replicate_concordance_decision.tsv`); if kept unpooled, colour ordinations by `tech_rep`/`bio_rep` as a visual re-check that the concordance call was right before interpreting any other grouping.
- Distances: **Bray–Curtis** (`vegan::vegdist`), **weighted & unweighted UniFrac** (`GUniFrac::GUniFrac` with `counts` + `tree`), **Aitchison** (`vegdist(clr, "euclidean")`).
- Ordinations: **PCoA** (`ape::pcoa` / `cmdscale`) and **NMDS** (`vegan::metaMDS`, report stress); color by `condition` and gradients.
- **PERMANOVA** (`vegan::adonis2`, `by="margin"`) + **`betadisper`**/`permutest` for dispersion homogeneity.
- Save distance matrices to `data/processed/`; coords + stats to `results/`; plots to `plots/`.

### `07_environmental_drivers.R` — analytical core
- Build & scale the environmental matrix; screen collinearity (correlation heatmap, VIF via `usdm`/`car`).
- **db-RDA** (`vegan::capscale` / `dbrda` on Bray or Aitchison) and/or **CCA**; forward-select with `ordiR2step`; permutation-test axes/terms.
- **`envfit`** of soil variables (and key taxa) onto the ordination.
- **Variance partitioning** (`vegan::varpart`): chemistry vs. spatial vs. climate.
- **Mantel / partial Mantel** (`vegan::mantel`): community vs. environmental vs. geographic distance.
- Outputs: fitted models (`results/*.rds`), drivers summary table (R², p, % variance), triplot(s) to `plots/`.

### `08_differential_abundance.R`
- Primary: **ALDEx2** (`aldex` on the count matrix + `condition`) — compositional, matrix-native, no phyloseq.
- Cross-check: **Maaslin2** (feature data frame + metadata data frame) and optionally **DESeq2** (`DESeqDataSetFromMatrix`). Report the intersection.
- Volcano / effect-size plots; significant-taxa table with effect sizes and adjusted p.
- Outputs: `results/da_<method>_<contrast>.tsv`; plots to `plots/`.

### `09_spatial_analysis.R`
- `sf` points from `lon`/`lat`; basemap via `terra` (DEM/land cover) or a coastline outline.
- Maps colored by richness / dominant taxon / key gradient.
- **Distance–decay** of community similarity vs. geographic distance; **Moran's I** (`ape::Moran.I`) on alpha diversity and leading ordination axes.
- Outputs: maps + distance–decay plot to `plots/`; spatial stats to `results/`.

### `10_figures.R`
Assemble multi-panel publication figures with `patchwork` (Fig 1 map+alpha;
Fig 2 ordination+PERMANOVA; Fig 3 db-RDA triplot+varpart; Fig 4 composition+DA).
Consistent theme, labeled panels, journal dimensions. Write
`results/figure_index.tsv`.

---

## 8. Shared helpers (`functions.R`)

- `read_feature_table(path)` — TSV or biom → integer matrix, samples as rows.
- `split_taxonomy(tbl)` — lineage string → rank columns, prefixes stripped.
- `assert_aligned(counts, taxonomy, metadata, tree)` — stops on id mismatch.
- `save_plot(plot, name, w, h)` — matched `.pdf` + `.png`.
- `write_result(x, name)` — standardized `write_tsv` into `results/`.
- `theme_soil()` + palette getters; tidy wrappers for `adonis2`/`envfit` → data.frame.

---

## 9. Orchestration

`run_all.R` sources `00`→`10` with logging + runtime summary, then writes
`session_info.txt`. Optional `Makefile` with a target per script keyed on the
`.rds`/`.tsv` each consumes, so `make all` reruns only what changed.

---

## 10. Deliverables checklist

- [x] Repo scaffolded per §4 with filled `README.md` + `AGENTS.md`.
- [x] `renv.lock` + one-command setup (renv's autoloader needed a fix -- see AGENTS.md environment trap #3).
- [x] Scripts `00`–`10` + `functions.R`, each independently runnable, **phyloseq- and qiime-free**.
- [x] `run_all.R` reproducing everything end-to-end (no `Makefile` -- optional per §9, not built).
- [x] Figures in `plots/` (pdf+png), tables/objects in `results/`, `session_info.txt`.
- [x] Every statistical choice (rarefaction depth, filters, model terms) documented -- see AGENTS.md "Downstream analysis notes".

---

## 11. Decisions to confirm before coding

1. **Taxonomy DB** upstream (SILVA / GTDB / Greengenes2) — sets rank parsing + lineage-filter strings.
2. Feature table delivered as **TSV or biom** — picks the reader.
3. Presence of **negative controls** — enables/disables `decontam`.
4. **Grouping variable(s)** (`condition`) + available **continuous soil variables** — defines contrasts + the environmental matrix.
5. **Coordinates** present — enables/disables `09_spatial_analysis.R`.
6. Preferred **DA method** if not defaulting to ALDEx2.
7. **Replicate pooling threshold** — how discordant a technical (extraction) or biological pair has to be, on the §7 `02_qc_filter.R` metrics, before it's flagged rather than pooled; and whether discordant pairs are dropped, kept unpooled and carried through as-is, or pooled anyway with a caveat.

---

## Kickoff prompt for Claude Code

> Read `PLAN.md`. Scaffold the repository exactly as in §4, initialize `renv`,
> and create `00_setup.R`, `functions.R`, and `01_import.R` first. Build the
> canonical objects from §3 (no phyloseq, no qiime2R) and stop after `01` so we
> can verify the imported matrix/tree/metadata against the real files in
> `data/raw/` before continuing to `02`–`10`.
