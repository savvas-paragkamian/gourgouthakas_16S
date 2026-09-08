# 02b_controls.R -- the control-sample accuracy check AGENTS.md has always
# said should exist but never did (mocks as "a separate accuracy check, not
# decontam input"). Second half of QC (see scripts/02_qc_filter.R's header
# for why it's split): T1-T5 use the control libraries while they're still
# present, then this script applies the reviewed contaminant blacklist + the
# bleed-through floor, drops the 5 control libraries, runs the prevalence
# filter, runs replicate concordance, and saves *_clean.rds -- the objects
# 03_normalize.R onward actually read.
#
# Every table here that reports an abundance also reports the underlying
# read count next to it (T4's point, generalized): a genus at 44% of a
# 300-read library and the same genus at 44% of a 300,000-read library are
# different findings, and this script never collapses that distinction away.

source("scripts/00_setup.R")

counts <- readRDS(file.path(path_processed, "counts_postlineage.rds"))
taxonomy <- readRDS(file.path(path_processed, "taxonomy_postlineage.rds"))
metadata <- readRDS(file.path(path_processed, "metadata.rds")) # all 51 rows still
mock_expected <- readr::read_tsv("data/mock_expected.tsv", show_col_types = FALSE)

is_mock <- metadata$sample_type == "control" & metadata$control_type == "mock"
is_neg <- metadata$control_type %in% c("extraction_blank", "mastermix_blank")
is_real <- metadata$sample_type != "control"

mock_ids <- stats::setNames(metadata$sample_id[is_mock], metadata$control_type[is_mock])
# metadata$control_type for mocks is just "mock" (build_inputs.R) -- the finer
# mock_even/mock_log/mock_env label lives in mock_expected.tsv, not metadata.
mock_sample_for <- c(mock_even = "ctr_zymo_com", mock_log = "ctr_zymo_log", mock_env = "ctr_msa_3001")

# ============================================================================
# T1 -- name reconciliation (GTDB *and* GG2 both split/rename genera)
# ============================================================================
# Pattern confirmed this session on a real ASV: GTDB calls it
# "Akkermansia_muciniphila_A", GG2 calls it "Akkermansia_muciniphila_D_776786"
# -- both append a polyphyly-split suffix a vendor's NCBI-style name won't
# match without stripping. Order matters: strip the longer "_LETTER_DIGITS"
# pattern before the bare "_LETTER" one.
normalize_epithet <- function(x) {
  x <- sub("_[A-Z]_[0-9]+$", "", x)
  x <- sub("_[A-Z]$", "", x)
  x
}

taxonomy <- taxonomy |>
  dplyr::mutate(genus_norm = normalize_epithet(genus), species_norm = normalize_epithet(species))

db_name_map <- taxonomy |>
  dplyr::filter(!is.na(genus), genus != genus_norm) |>
  dplyr::distinct(genus, genus_norm) |>
  dplyr::arrange(genus_norm, genus)
write_result(db_name_map, "02b_db_name_map")
message(sprintf("[02b_controls] T1: %d genus name(s) normalized (GTDB/GG2 polyphyly suffixes stripped)", nrow(db_name_map)))

# ============================================================================
# T2 -- database selection from the mocks
# ============================================================================
# Per-method taxonomy tables are keyed by full ASV SEQUENCE, not the md5
# asv_id (confirmed this session) -- join through taxonomy$sequence.
read_method_tax <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE) |>
    dplyr::rename(sequence = 1, taxon = Taxon) |>
    dplyr::mutate(genus = normalize_epithet(sub("^.*g__([^;]*).*$", "\\1", taxon)))
}
methods <- list(
  silva_nb = read_method_tax("results/hifi/nb_tax/silva_nb.tsv"),
  gtdb_nb = read_method_tax("results/hifi/nb_tax/gtdb_nb.tsv"),
  gg2_nb = read_method_tax("results/hifi/nb_tax/gg2_nb.tsv"),
  silva_vsearch = read_method_tax("results/hifi/vsearch_tax/silva_vsearch.tsv")
)

aitchison_dist_2 <- function(observed, expected) {
  # Aitchison distance between two composition vectors = Euclidean distance
  # of their CLR transforms. Zero-replace with a small pseudocount (no
  # cmultRepl here -- this is a 2-row, often-sparse comparison, not a
  # dataset-wide compositional matrix).
  m <- rbind(observed, expected)
  m[m == 0] <- min(m[m > 0], na.rm = TRUE) * 0.65
  m <- m / rowSums(m)
  clr_m <- t(apply(m, 1, function(r) log(r) - mean(log(r))))
  sqrt(sum((clr_m[1, ] - clr_m[2, ])^2))
}

score_mock_method <- function(mock_type, method_name, method_tax) {
  sample_id <- mock_sample_for[[mock_type]]
  expected <- mock_expected |> dplyr::filter(control_type == mock_type)
  expected_genus <- unique(expected$expected_genus)

  asv_seq <- taxonomy$sequence
  counts_row <- counts[sample_id, ]
  present <- counts_row[counts_row > 0]
  seq_for_present <- asv_seq[match(names(present), taxonomy$asv_id)]

  observed_genus <- method_tax$genus[match(seq_for_present, method_tax$sequence)]
  genus_counts <- stats::aggregate(as.numeric(present), by = list(genus = observed_genus), FUN = sum)
  names(genus_counts) <- c("genus", "count")
  genus_counts$rel_abund <- genus_counts$count / sum(genus_counts$count)

  recovered <- intersect(expected_genus, genus_counts$genus[genus_counts$rel_abund > 0])
  false_positive <- genus_counts |> dplyr::filter(rel_abund > 0.001, !genus %in% expected_genus)

  aitchison <- NA_real_
  lowest_input_detected <- NA_character_
  if (mock_type %in% c("mock_even", "mock_env")) {
    ref_col <- if (all(is.na(expected$pct_16s_corrected))) "genomic_dna_pct" else "pct_16s_corrected"
    exp_vec <- stats::setNames(expected[[ref_col]], expected$expected_genus)
    exp_vec <- exp_vec[!is.na(exp_vec)]
    obs_vec <- stats::setNames(rep(0, length(exp_vec)), names(exp_vec))
    common <- intersect(names(exp_vec), genus_counts$genus)
    obs_vec[common] <- genus_counts$rel_abund[match(common, genus_counts$genus)]
    if (length(exp_vec) >= 2) aitchison <- aitchison_dist_2(obs_vec, exp_vec)
  } else {
    ranked <- expected |> dplyr::arrange(genomic_dna_pct)
    detected_rank <- ranked |> dplyr::filter(expected_genus %in% recovered)
    if (nrow(detected_rank) > 0) lowest_input_detected <- detected_rank$expected_genus[1]
  }

  tibble::tibble(
    mock = mock_type, sample_id = sample_id, method = method_name,
    library_size = sum(present),
    n_expected_genera = length(expected_genus),
    n_recovered = length(recovered),
    n_false_positive_genera = nrow(false_positive),
    false_positive_genera = paste(false_positive$genus, collapse = ";"),
    aitchison_distance = aitchison,
    lowest_input_detected = lowest_input_detected
  )
}

mock_accuracy <- purrr::map_dfr(names(mock_sample_for), function(mt) {
  purrr::map_dfr(names(methods), function(m) score_mock_method(mt, m, methods[[m]]))
})
write_result(mock_accuracy, "02b_mock_accuracy")

p_mock_accuracy <- ggplot2::ggplot(mock_accuracy, ggplot2::aes(method, n_recovered / n_expected_genera, fill = method)) +
  ggplot2::geom_col() +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%d/%d\n(%d FP)", n_recovered, n_expected_genera, n_false_positive_genera)),
                      vjust = -0.2, size = 3) +
  ggplot2::facet_wrap(~mock) +
  ggplot2::scale_y_continuous(limits = c(0, 1.15), labels = scales::percent) +
  ggplot2::labs(x = NULL, y = "Expected genera recovered", fill = "Method",
                title = "Mock accuracy by taxonomy method",
                subtitle = "FP = false-positive genera at >0.1% relative abundance") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
save_plot(p_mock_accuracy, "02b_mock_accuracy", w = 8, h = 5)

# Informational only -- does NOT change db_to_prioritize (that's upstream
# HiFi-16S-workflow config). If this table suggests a different DB, that's a
# stop-and-ask finding for AGENTS.md, not a silent change (PLAN.md/AGENTS.md
# per this session's plan).
gg2_recall <- mock_accuracy |> dplyr::filter(method == "gg2_nb") |> dplyr::summarise(mean(n_recovered / n_expected_genera)) |> dplyr::pull()
best_method <- mock_accuracy |> dplyr::group_by(method) |> dplyr::summarise(mean_recall = mean(n_recovered / n_expected_genera), .groups = "drop") |> dplyr::slice_max(mean_recall, n = 1)
message(sprintf(
  "[02b_controls] T2: mean recall by method -- gg2_nb (pipeline default priority) = %.2f; best = %s (%.2f). %s",
  gg2_recall, best_method$method[1], best_method$mean_recall[1],
  if (best_method$method[1] != "gg2_nb" && best_method$mean_recall[1] > gg2_recall + 0.05)
    "NOTE: another method scores meaningfully higher -- see results/02b_mock_accuracy.tsv and consider whether db_to_prioritize should change (do not change it silently)."
  else "gg2_nb is competitive; no change indicated."
))

# ============================================================================
# T3 -- ASV-level false-positive rate (DADA2, whole-run property)
# ============================================================================
# Same off-target-lineage filter as 02_qc_filter.R §1 already applied
# (counts/taxonomy here are counts_postlineage). 10 genomes with multiple
# rRNA operons should denoise to roughly 20-50 true ASVs.
asv_fpr <- purrr::map_dfr(names(mock_sample_for), function(mt) {
  sample_id <- mock_sample_for[[mt]]
  n_asv <- sum(counts[sample_id, ] > 0)
  n_genomes <- dplyr::n_distinct((mock_expected |> dplyr::filter(control_type == mt))$expected_genus)
  expected_lo <- n_genomes * 2L
  expected_hi <- n_genomes * 5L
  tibble::tibble(
    mock = mt, sample_id = sample_id, library_size = sum(counts[sample_id, ]),
    n_genomes_expected = n_genomes, n_asv_observed = n_asv,
    expected_asv_range_lo = expected_lo, expected_asv_range_hi = expected_hi,
    excess_asvs = pmax(0, n_asv - expected_hi)
  )
})
write_result(asv_fpr, "02b_mock_asv_fpr")
message(sprintf("[02b_controls] T3: %s", paste(sprintf("%s: %d ASVs observed (expected ~%d-%d), excess=%d",
                asv_fpr$mock, asv_fpr$n_asv_observed, asv_fpr$expected_asv_range_lo, asv_fpr$expected_asv_range_hi, asv_fpr$excess_asvs),
                collapse = "; ")))

# ============================================================================
# T4 -- negative-control ASV inventory (EB vs. MM, two-level diagnostic)
# ============================================================================
neg_ids <- metadata$sample_id[is_neg]
real_ids <- metadata$sample_id[is_real]
neg_asvs <- colnames(counts)[colSums(counts[neg_ids, , drop = FALSE] > 0) > 0]

real_relabund <- counts[real_ids, , drop = FALSE] / rowSums(counts[real_ids, , drop = FALSE])

neg_inventory <- tibble::tibble(
  asv_id = neg_asvs,
  reads_ctr_EB = as.numeric(counts[neg_ids[metadata$control_type[match(neg_ids, metadata$sample_id)] == "extraction_blank"], neg_asvs]),
  reads_ctr_MM = as.numeric(counts[neg_ids[metadata$control_type[match(neg_ids, metadata$sample_id)] == "mastermix_blank"], neg_asvs]),
  in_both_blanks = reads_ctr_EB > 0 & reads_ctr_MM > 0,
  prevalence_real_samples = colSums(counts[real_ids, neg_asvs, drop = FALSE] > 0),
  max_relabund_real_samples = apply(real_relabund[, neg_asvs, drop = FALSE], 2, max)
) |>
  dplyr::left_join(
    taxonomy |> dplyr::select(asv_id, domain, phylum, genus_silva_nb = genus_norm),
    by = "asv_id"
  ) |>
  dplyr::left_join(methods$silva_nb |> dplyr::select(sequence, silva_nb = taxon), by = c("asv_id" = "sequence")) |>
  dplyr::arrange(dplyr::desc(reads_ctr_EB + reads_ctr_MM))
# note: the join above on asv_id==sequence for silva_nb is intentionally
# skipped (methods tables are keyed by sequence, not asv_id) -- taxonomy
# strings per method are read via the sequence key properly, done below.
neg_inventory$silva_nb <- NULL
seq_lookup <- taxonomy$sequence[match(neg_inventory$asv_id, taxonomy$asv_id)]
for (m in names(methods)) {
  neg_inventory[[m]] <- methods[[m]]$taxon[match(seq_lookup, methods[[m]]$sequence)]
}
write_result(neg_inventory, "02b_negative_control_asvs")

message(sprintf(
  "[02b_controls] T4: %d ASV(s) present in ctr_EB and/or ctr_MM; %d in both (reagent/PCR/prep signature), %d in one only (extraction-kit signature). ctr_EB=%d reads, ctr_MM=%d reads total -- shallow, screening power is limited (AGENTS.md).",
  length(neg_asvs), sum(neg_inventory$in_both_blanks), sum(!neg_inventory$in_both_blanks),
  sum(counts[neg_ids[metadata$control_type[match(neg_ids, metadata$sample_id)] == "extraction_blank"], ]),
  sum(counts[neg_ids[metadata$control_type[match(neg_ids, metadata$sample_id)] == "mastermix_blank"], ])
))

p_control_depth <- tibble::tibble(sample_id = metadata$sample_id[is_mock | is_neg],
                                    library_size = rowSums(counts[metadata$sample_id[is_mock | is_neg], , drop = FALSE]),
                                    control_type = metadata$control_type[match(metadata$sample_id[is_mock | is_neg], metadata$sample_id)]) |>
  ggplot2::ggplot(ggplot2::aes(forcats::fct_reorder(sample_id, library_size), library_size, fill = control_type)) +
  ggplot2::geom_col() +
  ggplot2::scale_y_log10() +
  ggplot2::coord_flip() +
  ggplot2::labs(x = NULL, y = "Library size (log10)", fill = "Control type",
                title = "Control read depth", subtitle = "Read counts, not just relative abundance -- see T4")
save_plot(p_control_depth, "02b_control_read_depth", w = 6, h = 3)

# ============================================================================
# T5 -- bidirectional bleed floor (PacBio barcode cross-talk, not Illumina
# index-hopping -- described that way deliberately)
# ============================================================================
mock_ids_vec <- metadata$sample_id[is_mock]
mock_exclusive <- colnames(counts)[colSums(counts[mock_ids_vec, , drop = FALSE] > 0) > 0 &
                                      colSums(counts[real_ids, , drop = FALSE] > 0) == 0]
bleed_in <- tibble::tibble(
  sample_id = real_ids,
  library_size = rowSums(counts[real_ids, , drop = FALSE]),
  mock_exclusive_reads = rowSums(counts[real_ids, mock_exclusive, drop = FALSE]),
  bleed_in_frac = mock_exclusive_reads / library_size
)

# Cave/water-characteristic genera: present in >=5 real samples, absent from
# every mock -- a conservative "clearly not from the mocks" set, not an
# exhaustive one.
real_prevalence <- colSums(counts[real_ids, , drop = FALSE] > 0)
mock_prevalence <- colSums(counts[mock_ids_vec, , drop = FALSE] > 0)
cave_characteristic <- colnames(counts)[real_prevalence >= 5 & mock_prevalence == 0]
bleed_out <- tibble::tibble(
  sample_id = mock_ids_vec,
  control_type = metadata$control_type[match(mock_ids_vec, metadata$sample_id)],
  library_size = rowSums(counts[mock_ids_vec, , drop = FALSE]),
  cave_reads = rowSums(counts[mock_ids_vec, cave_characteristic, drop = FALSE]),
  bleed_out_frac = cave_reads / library_size
)

write_result(dplyr::bind_rows(
  bleed_in |> dplyr::mutate(direction = "bleed_in_to_real_sample"),
  bleed_out |> dplyr::rename(mock_exclusive_reads = cave_reads, bleed_in_frac = bleed_out_frac) |>
    dplyr::mutate(direction = "bleed_out_to_mock")
), "02b_bleed_estimate")

bleed_floor <- max(c(bleed_in$bleed_in_frac, bleed_out$bleed_out_frac), na.rm = TRUE)

# The actual question this investigation started from: does the floor
# clear Akkermansia's sediment abundance?
akk_asv <- taxonomy$asv_id[!is.na(taxonomy$genus_norm) & taxonomy$genus_norm == "Akkermansia"]
akk_sediment_ids <- metadata$sample_id[metadata$sample_type == "sediment"]
akk_relabund <- if (length(akk_asv) > 0) {
  rs <- rowSums(counts[akk_sediment_ids, akk_asv, drop = FALSE])
  rs / rowSums(counts[akk_sediment_ids, , drop = FALSE])
} else {
  numeric(0)
}
akk_below_floor <- akk_relabund[akk_relabund > 0 & akk_relabund < bleed_floor]

decisions <- tibble::tibble(
  decision = c("bleed_floor", "akkermansia_min_sediment_relabund", "akkermansia_n_sediment_samples_below_floor"),
  value = c(bleed_floor, if (length(akk_relabund) > 0) min(akk_relabund[akk_relabund > 0]) else NA_real_, length(akk_below_floor)),
  note = c(
    sprintf("max(bleed_in, bleed_out) across all real samples and mocks; applied as a per-cell relative-abundance floor below in this script"),
    "lowest nonzero Akkermansia relative abundance among sediment samples",
    "sediment samples where Akkermansia is present but below the bleed floor -- these would read as bleed-through, not real presence, at that abundance"
  )
)
write_result(decisions, "02b_decisions")

message(sprintf(
  "[02b_controls] T5: bleed floor = %.5f (%.3f%%). Akkermansia sediment relabund range does %s the floor for %d/%d nonzero sediment samples.",
  bleed_floor, bleed_floor * 100,
  if (length(akk_below_floor) > 0) "fall below" else "stay above",
  length(akk_below_floor), sum(akk_relabund > 0)
))

# ============================================================================
# Apply: blacklist + bleed floor -> drop controls -> prevalence filter
# ============================================================================
candidates <- readr::read_tsv("results/02b_contaminant_candidates.tsv", show_col_types = FALSE)
blacklist <- candidates$asv_id[!is.na(candidates$decision) & candidates$decision == "remove"]
write_result(candidates |> dplyr::filter(asv_id %in% blacklist), "02b_contaminant_blacklist")
message(sprintf(
  "[02b_controls] applying blacklist: %d/%d reviewed candidate(s) marked \"remove\" (%d unreviewed/kept, treated as not-removed)",
  length(blacklist), nrow(candidates), sum(is.na(candidates$decision) | candidates$decision != "remove")
))

counts_deblacklisted <- counts[, !(colnames(counts) %in% blacklist), drop = FALSE]
taxonomy_deblacklisted <- taxonomy[!(taxonomy$asv_id %in% blacklist), ]

# Bleed floor applied to REAL samples only (mocks/blanks have already done
# their job feeding T1-T5; the floor describes cross-talk *into* real
# samples, so it's real-sample cells that get zeroed below it).
rel_full <- counts_deblacklisted / rowSums(counts_deblacklisted)
below_floor_real <- rel_full[real_ids, , drop = FALSE] < bleed_floor & rel_full[real_ids, , drop = FALSE] > 0
n_zeroed <- sum(below_floor_real)
counts_floored <- counts_deblacklisted
counts_floored[real_ids, ][below_floor_real] <- 0L
message(sprintf("[02b_controls] bleed floor zeroed %d real-sample x ASV cell(s) below %.5f relative abundance", n_zeroed, bleed_floor))

counts_noctrl <- counts_floored[real_ids, , drop = FALSE] # by name, not position -- counts' row order need not match metadata's
metadata_noctrl <- metadata[is_real, ]

prevalence <- colSums(counts_noctrl > 0)
keep_prevalence <- prevalence >= 2L
counts_filtered <- counts_noctrl[, keep_prevalence, drop = FALSE]
taxonomy_filtered <- taxonomy_deblacklisted[keep_prevalence, ]

message(sprintf(
  "[02b_controls] final gate: %d controls dropped, blacklist -%d ASVs, prevalence filter (>=2 samples) %d/%d ASVs kept -> %d samples x %d ASVs, %d reads",
  sum(is_mock | is_neg), length(blacklist), sum(keep_prevalence), length(keep_prevalence),
  nrow(counts_filtered), ncol(counts_filtered), sum(counts_filtered)
))

stopifnot("control row(s) survived to counts_clean/metadata_clean" = !any(metadata_noctrl$sample_type == "control"))

# ============================================================================
# Replicate concordance -- both levels (moved here from 02_qc_filter.R: needs
# the final post-blacklist/floor/prevalence matrix, not the raw one)
# ============================================================================
rel_ab <- counts_filtered / rowSums(counts_filtered)

# One row per pair of samples (not built separately for within vs. between --
# a single table, filtered two ways, is easier to follow and to check).
# "Within" = the two samples share `pair_key` (the replicate-pair identifier
# -- `site` for technical, `location`+`tech_rep` for biological). "Between"
# is the reference/noise-floor distribution a within-pair distance gets
# judged against, restricted to pairs of the SAME sample_type: a
# sediment-vs-water pair is trivially, definitionally dissimilar and says
# nothing about whether a same-type replicate pair reproduced, so it doesn't
# belong in that reference distribution.
concordance_at <- function(rel_ab, counts, metadata, pair_key, level_label) {
  bray_d <- as.matrix(vegan::vegdist(rel_ab, method = "bray"))
  jac_d <- as.matrix(vegan::vegdist(counts > 0, method = "jaccard"))
  richness <- rowSums(counts > 0)
  libsize <- rowSums(counts)

  ids <- rownames(rel_ab)
  meta_i <- metadata[match(ids, metadata$sample_id), ]
  i <- utils::combn(seq_along(ids), 2)[1, ]
  j <- utils::combn(seq_along(ids), 2)[2, ]

  pairs <- tibble::tibble(
    sample_a = ids[i], sample_b = ids[j],
    key_a = meta_i[[pair_key]][i], key_b = meta_i[[pair_key]][j],
    type_a = meta_i$sample_type[i], type_b = meta_i$sample_type[j]
  ) |>
    dplyr::mutate(
      level = level_label,
      bray = bray_d[cbind(sample_a, sample_b)],
      jaccard = jac_d[cbind(sample_a, sample_b)],
      spearman_cor = purrr::map2_dbl(sample_a, sample_b, \(a, b) suppressWarnings(stats::cor(rel_ab[a, ], rel_ab[b, ], method = "spearman"))),
      delta_richness = richness[sample_a] - richness[sample_b],
      delta_library_size = libsize[sample_a] - libsize[sample_b],
      is_within = key_a == key_b,
      same_type = type_a == type_b
    )

  within_pair <- pairs |>
    dplyr::filter(is_within) |>
    dplyr::transmute(level, pair_id = key_a, sample_a, sample_b, bray, jaccard,
                      spearman_cor, delta_richness, delta_library_size)
  between_pair <- pairs |>
    dplyr::filter(!is_within, same_type) |>
    dplyr::select(level, sample_a, sample_b, bray, jaccard)

  list(within = within_pair, between = between_pair)
}

tech_conc <- concordance_at(rel_ab, counts_filtered, metadata_noctrl, pair_key = "site", level_label = "technical")
metadata_biopair <- metadata_noctrl |> dplyr::mutate(bio_pair_key = paste(location, tech_rep, sep = "_"))
bio_conc <- concordance_at(rel_ab, counts_filtered, metadata_biopair, pair_key = "bio_pair_key", level_label = "biological")

write_result(tech_conc$within, "replicate_concordance_technical")
write_result(bio_conc$within, "replicate_concordance_biological")

flag_discordant <- function(within, between) {
  thresh <- stats::quantile(between$bray, 0.05, na.rm = TRUE)
  within |> dplyr::mutate(discordant = bray >= thresh)
}
tech_flagged <- flag_discordant(tech_conc$within, tech_conc$between)
bio_flagged <- flag_discordant(bio_conc$within, bio_conc$between)

p_concordance <- dplyr::bind_rows(
  tech_conc$within |> dplyr::mutate(pair_type = "within-pair"), tech_conc$between |> dplyr::mutate(pair_type = "between-pair"),
  bio_conc$within |> dplyr::mutate(pair_type = "within-pair"), bio_conc$between |> dplyr::mutate(pair_type = "between-pair")
) |>
  ggplot2::ggplot(ggplot2::aes(x = level, y = bray, fill = pair_type)) +
  # geom_boxplot dodges its two pair_type boxes apart per level (default
  # position_dodge2); geom_jitter doesn't know about that grouping on its
  # own and was jittering around the plain, un-dodged x position, scattering
  # points across both boxes instead of onto their own. position_dodge2 and
  # position_jitterdodge don't share the same width semantics, so both
  # layers are pinned to an explicit, identical position_dodge(width=0.75)
  # (jitterdodge's own jitter is layered on top via jitter.width/height).
  ggplot2::geom_boxplot(outlier.size = 0.5, position = ggplot2::position_dodge(width = 0.75)) +
  ggplot2::geom_point(ggplot2::aes(color = pair_type),
                       position = ggplot2::position_jitterdodge(dodge.width = 0.75, jitter.width = 0.15),
                       alpha = 0.4, size = 0.8) +
  ggplot2::labs(x = NULL, y = "Bray-Curtis distance", title = "Replicate concordance: within-pair vs. between-pair",
                subtitle = "Lower is better -- a within-pair distance near the between-pair distribution means that replicate didn't reproduce")
save_plot(p_concordance, "02_replicate_concordance", w = 6, h = 5)

n_discordant_tech <- sum(tech_flagged$discordant, na.rm = TRUE)
n_discordant_bio <- sum(bio_flagged$discordant, na.rm = TRUE)
message(sprintf("[02b_controls] replicate concordance: %d/%d technical pairs discordant, %d/%d biological pairs discordant",
                n_discordant_tech, nrow(tech_flagged), n_discordant_bio, nrow(bio_flagged)))

metadata_flagged <- metadata_noctrl |>
  dplyr::mutate(flag_discordant_pair = sample_id %in% c(
    unlist(tech_flagged |> dplyr::filter(discordant) |> dplyr::select(sample_a, sample_b)),
    unlist(bio_flagged |> dplyr::filter(discordant) |> dplyr::select(sample_a, sample_b))
  ))

pooling_decision <- tibble::tibble(
  level = c("technical", "biological"), n_pairs = c(nrow(tech_flagged), nrow(bio_flagged)),
  n_discordant = c(n_discordant_tech, n_discordant_bio), pooled_by_default = c(FALSE, FALSE),
  note = "not pooled by default -- see plots/02_replicate_concordance.png and confirm PLAN.md before changing this"
)
write_result(pooling_decision, "replicate_concordance_decision")

# --- save the objects 03_normalize.R onward actually read -------------------
assert_aligned(counts_filtered, taxonomy_filtered, metadata_flagged)
saveRDS(counts_filtered, file.path(path_processed, "counts_clean.rds"))
saveRDS(taxonomy_filtered, file.path(path_processed, "taxonomy_clean.rds"))
saveRDS(metadata_flagged, file.path(path_processed, "metadata_clean.rds"))

retention <- tibble::tibble(
  stage = c("postlineage", "blacklist_removed", "controls_dropped", "prevalence_filter"),
  n_asv = c(ncol(counts), ncol(counts_deblacklisted), ncol(counts_noctrl), ncol(counts_filtered)),
  n_samples = c(nrow(counts), nrow(counts_deblacklisted), nrow(counts_noctrl), nrow(counts_filtered)),
  n_reads = c(sum(counts), sum(counts_deblacklisted), sum(counts_noctrl), sum(counts_filtered))
)
write_result(retention, "02b_retention")

message(sprintf("[02b_controls] final: %d samples x %d ASVs, %d reads (0 controls)",
                nrow(counts_filtered), ncol(counts_filtered), sum(counts_filtered)))
