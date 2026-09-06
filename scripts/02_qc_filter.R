# 02_qc_filter.R — lineage filtering, library-size QC, decontam, prevalence
# filter, and (PLAN.md §7 update) two-level replicate concordance.

source("scripts/00_setup.R")

counts <- readRDS(file.path(path_processed, "counts.rds"))
taxonomy <- readRDS(file.path(path_processed, "taxonomy.rds"))
metadata <- readRDS(file.path(path_processed, "metadata.rds"))

# W5's conductivity_ms/temperature_c transposition (AGENTS.md "Known
# data-quality issue") is fixed at the source (data/metadata.tsv) now -- no
# runtime correction needed here any more. It mattered well beyond looking
# wrong: conductivity_ms is only measured for water, so W5's swapped pair
# used to be 2 of the only 10 complete-case observations 07's correlation/
# VIF/db-RDA/envfit/varpart/Mantel analyses had to work with, forcing a
# mechanically-exact -1.00 temperature~conductivity correlation that wasn't
# real. 01_import.R now asserts this can't recur silently.

n_asv_start <- ncol(counts)
n_reads_start <- sum(counts)

# --- 1. drop non-target lineages -----------------------------------------
# Mitochondria/chloroplast (keyword search across all ranks -- GG2/SILVA/GTDB
# don't agree on which rank carries it), Eukaryota, and ASVs unassigned at
# phylum.
is_offtarget <- taxonomy$domain == "Eukaryota" |
  is.na(taxonomy$phylum) |
  apply(taxonomy[, c("domain", "phylum", "class", "order", "family", "genus", "species")],
        1, function(r) any(grepl("mitochond|chloroplast", r, ignore.case = TRUE), na.rm = TRUE))
is_offtarget[is.na(is_offtarget)] <- FALSE

taxonomy_kept <- taxonomy[!is_offtarget, ]
counts_kept <- counts[, taxonomy_kept$asv_id, drop = FALSE]

message(sprintf(
  "[02_qc_filter] dropped %d/%d off-target ASVs (mito/chloroplast/Eukaryota/unassigned-phylum), %d reads",
  sum(is_offtarget), n_asv_start, sum(counts) - sum(counts_kept)
))

# --- 2. library sizes ------------------------------------------------------
lib_sizes <- tibble::tibble(
  sample_id = rownames(counts_kept),
  library_size = rowSums(counts_kept)
) |>
  dplyr::left_join(metadata, by = "sample_id") |>
  dplyr::arrange(library_size)

p_libsize <- ggplot2::ggplot(lib_sizes, ggplot2::aes(
  x = forcats::fct_reorder(sample_id, library_size), y = library_size, fill = sample_type
)) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_manual(values = palette_sample_type()) +
  ggplot2::scale_y_log10() +
  ggplot2::coord_flip() +
  ggplot2::labs(x = NULL, y = "Library size (log10)", fill = "Sample type",
                title = "Library sizes after lineage filtering") +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6))
save_plot(p_libsize, "02_library_sizes", w = 7, h = 10)

# Controls (blanks especially) are *expected* to be shallow -- don't flag
# those on the same footing as a failed biological library. Floor is
# documented, not derived, because there's no natural break in this data;
# revisit if a future run shows one.
shallow_floor <- 500L
shallow_biological <- lib_sizes |>
  dplyr::filter(sample_type != "control", library_size < shallow_floor)

if (nrow(shallow_biological) > 0) {
  message(sprintf(
    "[02_qc_filter] %d biological sample(s) below the %d-read shallow floor (flagged, not auto-dropped): %s",
    nrow(shallow_biological), shallow_floor, paste(shallow_biological$sample_id, collapse = ", ")
  ))
}
write_result(lib_sizes |> dplyr::select(sample_id, library_size, sample_type, sample_set, site),
             "02_library_sizes")

# --- 3. decontam (prevalence method) --------------------------------------
# ctr_EB (extraction blank) and ctr_MM (mastermix blank) are the negatives
# (AGENTS.md departure #3); the three mocks are an accuracy check, not
# decontam input.
is_neg <- metadata$control_type[match(rownames(counts_kept), metadata$sample_id)] %in%
  c("extraction_blank", "mastermix_blank")
stopifnot(sum(is_neg) == 2L)

contam <- decontam::isContaminant(counts_kept, method = "prevalence", neg = is_neg)
is_contaminant <- contam$contaminant & !is.na(contam$contaminant)

message(sprintf(
  "[02_qc_filter] decontam (prevalence, neg = ctr_EB + ctr_MM): %d/%d ASVs flagged as contaminants",
  sum(is_contaminant), ncol(counts_kept)
))

counts_decontam <- counts_kept[, !is_contaminant, drop = FALSE]
taxonomy_decontam <- taxonomy_kept[!is_contaminant, ]

write_result(
  tibble::tibble(asv_id = rownames(contam), contam) |> dplyr::filter(is_contaminant),
  "02_decontam_flagged_asvs"
)

# Standard decontam diagnostic (its own vignette's PA plot): prevalence of
# each ASV among the negative controls vs. among every other sample, colored
# by the contaminant call. A real contaminant should sit high on the
# negative-control axis and/or low on the true-sample axis; anything flagged
# that instead lands in the bottom-right (common in real samples, rare in
# blanks) is worth a second look at the prevalence-method threshold.
pa_neg <- colSums(counts_kept[is_neg, , drop = FALSE] > 0)
pa_pos <- colSums(counts_kept[!is_neg, , drop = FALSE] > 0)
decontam_pa <- tibble::tibble(
  asv_id = colnames(counts_kept),
  prevalence_controls = pa_neg,
  prevalence_samples = pa_pos,
  contaminant = is_contaminant
)
write_result(decontam_pa, "02_decontam_prevalence")

p_decontam <- ggplot2::ggplot(decontam_pa, ggplot2::aes(
  x = prevalence_controls, y = prevalence_samples, color = contaminant
)) +
  ggplot2::geom_jitter(width = 0.08, height = 0.4, alpha = 0.4, size = 1) +
  ggplot2::scale_color_manual(values = c(`TRUE` = "#E15759", `FALSE` = "grey60")) +
  ggplot2::scale_x_continuous(breaks = 0:sum(is_neg)) +
  ggplot2::labs(x = "Prevalence (negative controls: ctr_EB, ctr_MM)",
                y = "Prevalence (all other samples)", color = "Contaminant\n(decontam)",
                title = "Decontam diagnostic: ASV prevalence, controls vs. samples",
                subtitle = sprintf("%d/%d ASVs flagged, prevalence method", sum(is_contaminant), ncol(counts_kept)))
save_plot(p_decontam, "02_decontam_prevalence", w = 6, h = 5)

# --- 4. prevalence filter (+ optional low-count filter) --------------------
prevalence <- colSums(counts_decontam > 0)
keep_prevalence <- prevalence >= 2L

counts_filtered <- counts_decontam[, keep_prevalence, drop = FALSE]
taxonomy_filtered <- taxonomy_decontam[keep_prevalence, ]

message(sprintf(
  "[02_qc_filter] prevalence filter (>=2 samples): %d/%d ASVs kept",
  sum(keep_prevalence), length(keep_prevalence)
))

assert_aligned(counts_filtered, taxonomy_filtered, metadata)

# ===========================================================================
# 5. Replicate concordance -- both levels, before any pooling decision
#
# Sample IDs carry two nested levels of replication (AGENTS.md "The data"):
#   - biological: `bio_rep` (sample_set: transect / isolate_source), paired
#     via `location` (site with the "_I" suffix stripped)
#   - technical:  `tech_rep` (the 1/2 DNA-extraction replicate), paired
#     within a single `site`
# The ISD reference scripts (isd_archive_scripts/isd_crete_numerical_ecology.R)
# only ever checked the biological level (loc_1/loc_2 pairwise dissimilarity).
# This extends that to the technical level too -- it's the finer-grained
# check and the more diagnostic of a bad extraction.
# ===========================================================================

rel_ab <- counts_filtered / rowSums(counts_filtered)

# concordance_at(): given a grouping key that identifies replicate PAIRS
# (e.g. `site` for technical reps, `location` for biological reps) and a
# column that distinguishes the two members of a pair (e.g. `tech_rep`,
# `bio_rep`), compute within-pair vs. between-pair Bray-Curtis/Jaccard plus
# the cheaper sanity signals (Spearman correlation, delta-richness,
# delta-library-size). Sites/locations without exactly 2 members are skipped.
concordance_at <- function(rel_ab, counts, metadata, pair_key, level_label) {
  bray_d <- as.matrix(vegan::vegdist(rel_ab, method = "bray"))
  jac_d <- as.matrix(vegan::vegdist(counts > 0, method = "jaccard"))

  ids <- rownames(rel_ab)
  meta_i <- metadata[match(ids, metadata$sample_id), ]
  richness <- rowSums(counts > 0)
  libsize <- rowSums(counts)

  groups <- split(ids, meta_i[[pair_key]])
  groups <- groups[lengths(groups) == 2L] # only real pairs

  pair_rows <- lapply(names(groups), function(g) {
    a <- groups[[g]][1]; b <- groups[[g]][2]
    tibble::tibble(
      level = level_label,
      pair_id = g,
      sample_a = a, sample_b = b,
      bray = bray_d[a, b],
      jaccard = jac_d[a, b],
      spearman_cor = suppressWarnings(stats::cor(rel_ab[a, ], rel_ab[b, ], method = "spearman")),
      delta_richness = richness[a] - richness[b],
      delta_library_size = libsize[a] - libsize[b]
    )
  })
  within_pair <- dplyr::bind_rows(pair_rows)

  # Between-pair background: same distance metric, all cross-group sample
  # pairs that are NOT one of the within-pair comparisons above.
  all_pairs <- t(utils::combn(ids, 2))
  within_pair_set <- paste(within_pair$sample_a, within_pair$sample_b)
  is_within <- paste(all_pairs[, 1], all_pairs[, 2]) %in% within_pair_set
  between_pair <- tibble::tibble(
    level = level_label,
    sample_a = all_pairs[!is_within, 1],
    sample_b = all_pairs[!is_within, 2],
    bray = bray_d[all_pairs[!is_within, , drop = FALSE]],
    jaccard = jac_d[all_pairs[!is_within, , drop = FALSE]]
  )

  list(within = within_pair, between = between_pair)
}

tech_conc <- concordance_at(rel_ab, counts_filtered, metadata, pair_key = "site", level_label = "technical")

# Biological pairing is NOT `location` alone -- a sediment location has 4
# samples (2 bio_rep x 2 tech_rep), not 2, so grouping on `location` either
# skips every sediment site (group size 4) or -- the bug this comment is
# replacing -- silently picks up water's (location) groups of size 2 instead,
# which are actually *technical* pairs, not biological ones (water has no
# bio_rep arm at all). The real biological pair is the same location AND the
# same tech_rep, differing only in bio_rep: C1_1_A vs C1I_1_A, C1_2_A vs
# C1I_2_A. Water/control rows get a unique key each (no bio_rep counterpart)
# so they can't accidentally form a pair here.
metadata_biopair <- metadata |>
  dplyr::mutate(bio_pair_key = paste(location, tech_rep, sep = "_"))
bio_conc <- concordance_at(rel_ab, counts_filtered, metadata_biopair, pair_key = "bio_pair_key", level_label = "biological")

write_result(tech_conc$within, "replicate_concordance_technical")
write_result(bio_conc$within, "replicate_concordance_biological")

# Flag (don't drop) any within-pair distance that falls inside the
# between-pair distribution -- i.e. no more similar than two unrelated
# samples. That's a replicate that didn't reproduce.
flag_discordant <- function(within, between) {
  thresh <- stats::quantile(between$bray, 0.05, na.rm = TRUE) # 5th pct of *unrelated*-pair distances
  within |> dplyr::mutate(discordant = bray >= thresh)
}
tech_flagged <- flag_discordant(tech_conc$within, tech_conc$between)
bio_flagged <- flag_discordant(bio_conc$within, bio_conc$between)

p_concordance <- dplyr::bind_rows(
  tech_conc$within |> dplyr::mutate(pair_type = "within-pair"),
  tech_conc$between |> dplyr::mutate(pair_type = "between-pair"),
  bio_conc$within |> dplyr::mutate(pair_type = "within-pair"),
  bio_conc$between |> dplyr::mutate(pair_type = "between-pair")
) |>
  ggplot2::ggplot(ggplot2::aes(x = level, y = bray, fill = pair_type)) +
  ggplot2::geom_boxplot(outlier.size = 0.5) +
  ggplot2::geom_jitter(ggplot2::aes(color = pair_type), width = 0.15, alpha = 0.4, size = 0.8) +
  ggplot2::labs(x = NULL, y = "Bray-Curtis distance",
                title = "Replicate concordance: within-pair vs. between-pair",
                subtitle = "Lower is better -- a within-pair distance near the between-pair distribution means that replicate didn't reproduce")
save_plot(p_concordance, "02_replicate_concordance", w = 6, h = 5)

n_discordant_tech <- sum(tech_flagged$discordant, na.rm = TRUE)
n_discordant_bio <- sum(bio_flagged$discordant, na.rm = TRUE)
message(sprintf(
  "[02_qc_filter] replicate concordance: %d/%d technical pairs discordant, %d/%d biological pairs discordant",
  n_discordant_tech, nrow(tech_flagged), n_discordant_bio, nrow(bio_flagged)
))

metadata_flagged <- metadata |>
  dplyr::mutate(
    flag_discordant_pair = sample_id %in% c(
      unlist(tech_flagged |> dplyr::filter(discordant) |> dplyr::select(sample_a, sample_b)),
      unlist(bio_flagged |> dplyr::filter(discordant) |> dplyr::select(sample_a, sample_b))
    )
  )

# --- pooling decision -----------------------------------------------------
# Recorded, not silently acted on: both pooled and unpooled matrices are kept
# in data/processed/ so this is reversible. Default here is to NOT pool
# either level automatically -- 03_normalize.R onward runs on the unpooled,
# per-library matrix, coloring by tech_rep/bio_rep, and pooling is left as an
# explicit decision to confirm (PLAN.md §11.7) once a human has looked at
# 02_replicate_concordance.png.
pooling_decision <- tibble::tibble(
  level = c("technical", "biological"),
  n_pairs = c(nrow(tech_flagged), nrow(bio_flagged)),
  n_discordant = c(n_discordant_tech, n_discordant_bio),
  pooled_by_default = c(FALSE, FALSE),
  note = "not pooled by default -- see plots/02_replicate_concordance.png and confirm PLAN.md §11.7 before changing this"
)
write_result(pooling_decision, "replicate_concordance_decision")

# --- 6. re-align + save ---------------------------------------------------
assert_aligned(counts_filtered, taxonomy_filtered, metadata_flagged)

saveRDS(counts_filtered, file.path(path_processed, "counts_clean.rds"))
saveRDS(taxonomy_filtered, file.path(path_processed, "taxonomy_clean.rds"))
saveRDS(metadata_flagged, file.path(path_processed, "metadata_clean.rds"))

retention <- tibble::tibble(
  stage = c("import", "lineage_filter", "decontam", "prevalence_filter"),
  n_asv = c(n_asv_start, ncol(counts_kept), ncol(counts_decontam), ncol(counts_filtered)),
  n_reads = c(n_reads_start, sum(counts_kept), sum(counts_decontam), sum(counts_filtered))
)
write_result(retention, "02_retention")

message(sprintf(
  "[02_qc_filter] final: %d samples x %d ASVs, %d reads",
  nrow(counts_filtered), ncol(counts_filtered), sum(counts_filtered)
))
