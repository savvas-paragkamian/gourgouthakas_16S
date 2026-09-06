# 01_import.R — read the actual HiFi-16S-workflow outputs, build the
# canonical objects (PLAN.md §3), save to data/processed/*.rds.
#
# Adapted from PLAN.md §2's generic qiime-style file list: this project's
# real inputs are results/hifi/final/best_tax_merged_freq_tax.tsv (counts +
# taxonomy combined in one wide file, not separate feature-table.tsv /
# taxonomy.tsv) and data/metadata.tsv (soil chemistry already joined in by
# build_inputs.R, so there's no separate soil_chemistry.tsv). No tree.nwk:
# this dataset has no phylogeny (AGENTS.md departure #1), so `tree` stays
# NULL throughout and Faith's PD / UniFrac are out of scope everywhere below.

source("scripts/00_setup.R")

# --- read --------------------------------------------------------------

stopifnot(file.exists(path_asv_table), file.exists(path_metadata))

ft <- read_feature_table(path_asv_table)
counts <- ft$counts # samples x ASVs, integer

taxonomy <- ft$taxonomy |>
  split_taxonomy(lineage_col = "lineage") |>
  dplyr::select(asv_id, domain, phylum, class, order, family, genus, species,
                confidence, lineage, sequence)

metadata_raw <- readr::read_tsv(path_metadata, show_col_types = FALSE)

# --- metadata: canonical names + the replicate-structure columns PLAN.md §3
# now expects (added alongside the technical/biological replicate-concordance
# work in §7 02_qc_filter.R) ---------------------------------------------
#
# `site` in the raw metadata is already the *biological*-replicate-level id
# (C1 and C1_I are distinct `site` values, per AGENTS.md "The data"). There is
# no ready-made column pairing C1 with C1_I as the same depth location, so
# `location` is derived here by stripping the "_I" suffix — this is the join
# key 02_qc_filter.R uses for the biological-replicate concordance check.
metadata <- metadata_raw |>
  dplyr::rename(sample_id = sample_name) |>
  dplyr::mutate(
    location = sub("_I$", "", site), # shared key across the two bio_rep arms
    bio_rep = sample_set,            # "transect" / "isolate_source" / "control"
    tech_rep = replicate,            # 1 / 2 (DNA extraction); NA for controls
    condition = sample_type          # nearest equivalent to PLAN.md's `condition`
  ) |>
  dplyr::relocate(sample_id, site, location, bio_rep, tech_rep, sample_type,
                   sample_set, condition)

# W5's conductivity_ms/temperature_c used to be transposed in data/
# metadata.tsv (195.9 degC, impossible in a cave stream) -- fixed at the
# source now, not handled here. This is a standing guard against a future
# data refresh reintroducing it silently, not the fix itself.
temp_bad <- metadata$sample_id[!is.na(metadata$temperature_c) & metadata$sample_type == "water" &
                                  (metadata$temperature_c < 0 | metadata$temperature_c > 40)]
if (length(temp_bad) > 0) {
  stop(sprintf(
    "[01_import] Implausible temperature_c for water sample(s) %s -- looks transposed with conductivity_ms again (see AGENTS.md 'Known data-quality issue'). Fix data/metadata.tsv (or its source sheet), don't carry it through.",
    paste(temp_bad, collapse = ", ")
  ))
}

# --- read tracking (read counts through cutadapt -> filter -> denoise ->
# non-chimeric), for the completeness table below and later QC reporting ---
tracking <- if (file.exists(path_tracking)) {
  readr::read_tsv(path_tracking, show_col_types = FALSE) |>
    dplyr::rename(sample_id = sample)
} else {
  NULL
}

# --- representative sequences (optional; kept for provenance / potential
# future primer-trim or length sanity checks, not otherwise used downstream)
refseqs <- Biostrings::DNAStringSet(stats::setNames(taxonomy$sequence, taxonomy$asv_id))

# tree stays NULL: no phylogeny in this dataset (AGENTS.md departure #1).
tree <- NULL

# --- long form for ggplot/dplyr work ------------------------------------
counts_long <- counts |>
  as.data.frame() |>
  tibble::rownames_to_column("sample_id") |>
  tidyr::pivot_longer(-sample_id, names_to = "asv_id", values_to = "count")

# --- align + save --------------------------------------------------------
assert_aligned(counts, taxonomy |> dplyr::rename(asv_id = asv_id), metadata, tree = NULL)

saveRDS(counts, file.path(path_processed, "counts.rds"))
saveRDS(counts_long, file.path(path_processed, "counts_long.rds"))
saveRDS(taxonomy, file.path(path_processed, "taxonomy.rds"))
saveRDS(metadata, file.path(path_processed, "metadata.rds"))
saveRDS(refseqs, file.path(path_processed, "refseqs.rds"))
saveRDS(tree, file.path(path_processed, "tree.rds"))
if (!is.null(tracking)) saveRDS(tracking, file.path(path_processed, "tracking.rds"))

# --- completeness table ---------------------------------------------------
completeness <- metadata |>
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(!is.na(.x)))) |>
  tidyr::pivot_longer(dplyr::everything(), names_to = "variable", values_to = "fraction_complete") |>
  dplyr::arrange(fraction_complete)

write_result(completeness, "metadata_completeness")

message(sprintf(
  "[01_import] %d samples x %d ASVs; %d taxonomy ranks parsed; tree = NULL (no phylogeny)",
  nrow(counts), ncol(counts), sum(!is.na(taxonomy$phylum))
))
