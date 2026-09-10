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
  # NOTE: sub("^.*g__([^;]*).*$", ...) silently falls through to the RAW
  # taxon string, unmodified, whenever that rank prefix is absent from the
  # string entirely (no match -> sub() returns x as-is, not NA) -- and GG2
  # commonly classifies an ASV only to genus, leaving no "s__" field at all.
  # Caught via a raw "d__Bacteria;p__...;g__X" string leaking into the
  # species column of 02b_mock_recovery_species.tsv. str_extract() returns
  # NA on no match instead, which is the correct "no assignment at this
  # rank" value.
  extract_rank <- function(taxon, prefix) {
    epithet <- sub(paste0("^", prefix), "", stringr::str_extract(taxon, paste0(prefix, "[^;]*")))
    dplyr::na_if(epithet, "")
  }
  readr::read_tsv(path, show_col_types = FALSE) |>
    dplyr::rename(sequence = 1, taxon = Taxon) |>
    dplyr::mutate(
      genus = normalize_epithet(extract_rank(taxon, "g__")),
      species = normalize_epithet(extract_rank(taxon, "s__"))
    )
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
  # NA genus (unclassified at this rank -- e.g. silva_vsearch's Taxon field
  # has no g__ prefix at all, so every ASV is genus-NA for that method) is
  # dropped before aggregating: it's not a group, it's "no genus signal."
  # aggregate()'s FUN can't infer a return type from zero groups, so an
  # all-NA input makes the `count` column silently become list() instead of
  # numeric, and sum() downstream then errors -- guarded explicitly here
  # rather than relying on aggregate() to do something graceful with none.
  keep <- !is.na(observed_genus)
  if (!any(keep)) {
    genus_counts <- tibble::tibble(genus = character(0), count = numeric(0), rel_abund = numeric(0))
  } else {
    genus_counts <- stats::aggregate(as.numeric(present[keep]), by = list(genus = observed_genus[keep]), FUN = sum)
    names(genus_counts) <- c("genus", "count")
    genus_counts$rel_abund <- genus_counts$count / sum(genus_counts$count)
  }

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

# --- Mock recovery panel: expected vs. observed abundance, per genus -------
# The flagship control-validation figure this session's mock work was
# missing: not just a recall count (p_mock_accuracy above) but *how well*
# each expected organism's abundance was actually recovered. These three
# mocks are well-sequenced (library_size below; ~165K/49K/36K reads --
# results/02b_mock_accuracy.tsv) and deserve it.
#
# gg2_nb only, not all 4 methods -- the pipeline's own db_to_prioritize,
# T2's confirmed-competitive default. silva_vsearch is excluded on purpose:
# its taxonomy file (results/hifi/vsearch_tax/silva_vsearch.tsv) turns out to
# have malformed rows -- checked directly, 95% of lines have more than the
# expected 5 tab-separated fields, because some species-level strain
# descriptions (e.g. "Bacillus subtilis subsp. spizizenii str. W23") contain
# literal embedded tab characters instead of spaces in the reference
# database this pinned submodule shipped. That's the source of the "one or
# more parsing issues" warning that's shown up in every full-pipeline run
# this session -- not a bug in this project's own code, and not something to
# hand-edit in a pinned upstream file, but worth having finally tracked down.
# Vendor-name -> GTDB-genus crosswalk for genuine RECLASSIFICATIONS, not the
# polyphyly-split suffix pattern T1's normalize_epithet() already handles
# (that's "_A"/"_D_776786" tacked onto the *same* name, e.g. Listeria ->
# Listeria_A, and it's already applied to methods$gg2_nb$genus below).
# Checked systematically against every expected_species in mock_expected.tsv
# (search gg2_nb.tsv for each species epithet, compare genus): the only real
# hit is Lactobacillus fermentum, which GTDB places in Limosilactobacillus
# (the 2020 Zheng et al. split of the old Lactobacillus genus) -- confirmed
# directly (`grep fermentum results/hifi/nb_tax/gg2_nb.tsv` ->
# g__Limosilactobacillus), not assumed. Scoped to exactly this dataset's
# expected list (mock_expected.tsv has only Lactobacillus fermentum under
# that genus, no other species that would need a different split-genus
# target), not a general-purpose Lactobacillus crosswalk.
vendor_to_gtdb_genus <- c(Lactobacillus = "Limosilactobacillus")

mock_recovery <- purrr::map_dfr(names(mock_sample_for), function(mt) {
  sample_id <- mock_sample_for[[mt]]
  expected <- mock_expected |>
    dplyr::filter(control_type == mt) |>
    dplyr::mutate(expected_genus = dplyr::coalesce(vendor_to_gtdb_genus[expected_genus], expected_genus))
  ref_col <- if (all(is.na(expected$pct_16s_corrected))) "genomic_dna_pct" else "pct_16s_corrected"

  asv_seq <- taxonomy$sequence
  counts_row <- counts[sample_id, ]
  present <- counts_row[counts_row > 0]
  seq_for_present <- asv_seq[match(names(present), taxonomy$asv_id)]
  observed_genus <- methods$gg2_nb$genus[match(seq_for_present, methods$gg2_nb$sequence)]
  # Genus-NA reads (unclassified at this rank) are a negligible fraction in
  # every mock (<=0.01%, confirmed) but explicitly dropped before
  # aggregating rather than left to aggregate()'s default na.action, both
  # for clarity and because an all-NA input makes aggregate()'s `count`
  # column silently become list() instead of numeric (see score_mock_method
  # above, where this bit for real with silva_vsearch).
  keep <- !is.na(observed_genus)
  genus_reads <- stats::aggregate(as.numeric(present[keep]), by = list(genus = observed_genus[keep]), FUN = sum)
  names(genus_reads) <- c("genus", "reads")
  genus_reads$observed_pct <- 100 * genus_reads$reads / sum(genus_reads$reads)

  dplyr::full_join(
    expected |> dplyr::transmute(genus = expected_genus, expected_pct = .data[[ref_col]], ref_column = ref_col),
    tibble::as_tibble(genus_reads),
    by = "genus"
  ) |>
    dplyr::mutate(
      mock = mt, sample_id = sample_id, library_size = sum(present),
      status = dplyr::case_when(
        !is.na(expected_pct) & !is.na(observed_pct) ~ "recovered",
        !is.na(expected_pct) & is.na(observed_pct) ~ "not detected",
        TRUE ~ "false positive"
      )
    )
})
write_result(mock_recovery, "02b_mock_recovery")

# Undetected/false-positive genera can't be placed on a log axis at their
# true value (0 / NA) -- each mock's panel gets its own floor, a decade
# below that panel's own smallest nonzero recovered value, so "not
# detected" reads as "below what this panel can even show" rather than
# arbitrarily fixed. False positives are floored on the x axis the same way
# and restricted to >0.1% observed abundance (T2's existing FP threshold) --
# without that floor, mock_log's ASV-level noise (T3: 578 observed vs
# ~16-40 expected) would flood the panel with trace-level clutter.
mock_recovery_plot_df <- mock_recovery |>
  dplyr::filter(status != "false positive" | (observed_pct > 0.1)) |>
  dplyr::group_by(mock) |>
  dplyr::mutate(
    floor_pct = min(c(expected_pct, observed_pct), na.rm = TRUE) / 10,
    x = ifelse(is.na(expected_pct), floor_pct, expected_pct),
    y = ifelse(is.na(observed_pct), floor_pct, observed_pct)
  ) |>
  dplyr::ungroup()

p_mock_recovery <- ggplot2::ggplot(mock_recovery_plot_df, ggplot2::aes(x = x, y = y, color = status)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  ggplot2::geom_point(size = 2.2, alpha = 0.85) +
  ggrepel::geom_text_repel(ggplot2::aes(label = genus), size = 3, fontface = "italic",
                            show.legend = FALSE, max.overlaps = Inf, seed = 42) +
  ggplot2::facet_wrap(~mock, scales = "free", nrow = 1) +
  ggplot2::scale_x_log10(labels = scales::label_percent(scale = 1)) +
  ggplot2::scale_y_log10(labels = scales::label_percent(scale = 1)) +
  ggplot2::scale_color_manual(values = c(recovered = "#2E7DA6", "not detected" = "#B0794A", "false positive" = "#E15759")) +
  ggplot2::labs(x = "Expected (16S-corrected genomic %, or genomic % where uncorrected)", y = "Observed (gg2_nb relative abundance)",
                color = NULL, title = "Mock community recovery: expected vs. observed abundance",
                subtitle = "Dashed line = perfect recovery (y=x). Points at each panel's floor = not detected / not expected.") +
  ggplot2::theme(
    strip.text = ggplot2::element_text(size = 13, face = "bold"),
    plot.title = ggplot2::element_text(size = 18, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 11),
    legend.position = "bottom"
  )
save_plot(p_mock_recovery, "02b_mock_recovery", w = 15, h = 6, dpi = 350)

n_not_detected <- sum(mock_recovery$status == "not detected")
n_false_pos_shown <- sum(mock_recovery_plot_df$status == "false positive")
message(sprintf(
  "[02b_controls] mock recovery: %d/%d expected genus-mock pairs not detected at all; %d false-positive genera shown (>0.1%% observed) -- see results/02b_mock_recovery.tsv, plots/02b_mock_recovery.png",
  n_not_detected, sum(mock_recovery$status != "false positive"), n_false_pos_shown
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
# T6 -- mock cross-talk index: species-level recovery, mock-to-mock
# cross-talk, chimera rate, ASV FPR -- a formal index-hopping estimate, not
# just genus-level recall (T2 above) or the mock-to-real bleed floor (T5).
# ============================================================================
# Species-level recovery: same design as the genus-level 02b_mock_recovery
# above, keyed on species instead. The same vendor->GTDB crosswalk applies
# at both ranks (built once above for genus; the species string carries the
# same genus prefix, so it's rebuilt from the already-resolved mapped genus
# rather than duplicating the crosswalk lookup).
mock_recovery_species <- purrr::map_dfr(names(mock_sample_for), function(mt) {
  sample_id <- mock_sample_for[[mt]]
  expected <- mock_expected |>
    dplyr::filter(control_type == mt, !is.na(expected_species)) |>
    dplyr::mutate(
      mapped_genus = dplyr::coalesce(vendor_to_gtdb_genus[expected_genus], expected_genus),
      # NOTE: pattern must NOT depend on expected_genus row-by-row -- base
      # sub()'s `pattern` arg is not vectorized over its input (only the
      # first element is used for the whole column, silently), which briefly
      # corrupted this for every mock with >1 expected genus. expected_
      # species is always "Genus_epithet" with the epithet containing no
      # underscore, so stripping "everything up to the first underscore" is
      # both correct and pattern-row-independent.
      species_epithet = sub("^[^_]+_", "", expected_species),
      expected_species_mapped = paste0(mapped_genus, "_", species_epithet)
    )
  ref_col <- if (all(is.na(expected$pct_16s_corrected))) "genomic_dna_pct" else "pct_16s_corrected"

  asv_seq <- taxonomy$sequence
  counts_row <- counts[sample_id, ]
  present <- counts_row[counts_row > 0]
  seq_for_present <- asv_seq[match(names(present), taxonomy$asv_id)]
  observed_species <- methods$gg2_nb$species[match(seq_for_present, methods$gg2_nb$sequence)]
  keep <- !is.na(observed_species)
  species_reads <- stats::aggregate(as.numeric(present[keep]), by = list(species = observed_species[keep]), FUN = sum)
  names(species_reads) <- c("species", "reads")
  species_reads$observed_pct <- 100 * species_reads$reads / sum(species_reads$reads)

  dplyr::full_join(
    expected |> dplyr::transmute(species = expected_species_mapped, expected_pct = .data[[ref_col]]),
    tibble::as_tibble(species_reads),
    by = "species"
  ) |>
    dplyr::mutate(
      mock = mt, sample_id = sample_id,
      status = dplyr::case_when(
        !is.na(expected_pct) & !is.na(observed_pct) ~ "recovered",
        !is.na(expected_pct) & is.na(observed_pct) ~ "not detected",
        TRUE ~ "false positive"
      )
    )
})
write_result(mock_recovery_species, "02b_mock_recovery_species")

n_species_recovered <- sum(mock_recovery_species$status == "recovered")
n_species_expected <- sum(mock_recovery_species$status != "false positive")
message(sprintf(
  "[02b_controls] T6: species-level recovery %d/%d (%.0f%%) -- see results/02b_mock_recovery_species.tsv (genus-level was T2/02b_mock_recovery.tsv, %.0f%%)",
  n_species_recovered, n_species_expected, 100 * n_species_recovered / n_species_expected,
  100 * gg2_recall
))

# Mock-to-mock cross-talk: the most direct index-hopping signal this dataset
# has. These are synthetic communities with mostly mutually-exclusive
# membership by design (confirmed: 6 genera exclusive to mock_env, 3
# exclusive to mock_even+mock_log -- data/mock_expected.tsv) -- any read
# assigned to a genus expected in a DIFFERENT mock but not this one, inside
# this mock's own well, is barcode cross-talk on this SMRT cell, not a
# biological explanation. This is what the "Listeria in mock_env" false
# positive in 02b_mock_recovery.png actually is; this formalizes it into a
# rate.
mock_genus_lists <- mock_expected |>
  dplyr::distinct(control_type, expected_genus) |>
  dplyr::mutate(expected_genus = dplyr::coalesce(vendor_to_gtdb_genus[expected_genus], expected_genus)) |>
  dplyr::group_by(control_type) |>
  dplyr::summarise(genera = list(unique(expected_genus)), .groups = "drop")
genus_lookup <- stats::setNames(mock_genus_lists$genera, mock_genus_lists$control_type)

mock_crosstalk <- purrr::map_dfr(names(mock_sample_for), function(mt) {
  sample_id <- mock_sample_for[[mt]]
  this_genera <- genus_lookup[[mt]]
  foreign_genera <- setdiff(unique(unlist(genus_lookup[setdiff(names(genus_lookup), mt)])), this_genera)

  asv_seq <- taxonomy$sequence
  counts_row <- counts[sample_id, ]
  present <- counts_row[counts_row > 0]
  seq_for_present <- asv_seq[match(names(present), taxonomy$asv_id)]
  observed_genus <- methods$gg2_nb$genus[match(seq_for_present, methods$gg2_nb$sequence)]
  is_foreign <- observed_genus %in% foreign_genera

  tibble::tibble(
    mock = mt, sample_id = sample_id, library_size = sum(present),
    n_foreign_genera_checked = length(foreign_genera),
    foreign_genera_detected = paste(sort(unique(observed_genus[is_foreign])), collapse = ";"),
    crosstalk_reads = sum(present[is_foreign]),
    crosstalk_frac = sum(present[is_foreign]) / sum(present)
  )
})
write_result(mock_crosstalk, "02b_mock_crosstalk")

# Chimera rate (tracking.rds: denoised_reads -> non_chimeric_reads is DADA2's
# own chimera-removal step) for the mocks specifically, vs. real samples as
# context -- are the mocks unusually clean or dirty relative to the rest of
# the run, or in line with everything else? Loaded here (not just at the
# read-tracking figure further down) since T6 needs it now.
tracking <- readRDS(file.path(path_processed, "tracking.rds"))
chimera_rates <- tracking |>
  dplyr::mutate(chimera_reads = denoised_reads - non_chimeric_reads, chimera_rate = chimera_reads / denoised_reads) |>
  dplyr::mutate(group = dplyr::case_when(
    sample_id %in% mock_sample_for ~ "mock",
    sample_id %in% neg_ids ~ "blank",
    TRUE ~ "real"
  ))
write_result(chimera_rates, "02b_chimera_rates")
chimera_summary <- chimera_rates |> dplyr::group_by(group) |> dplyr::summarise(mean_chimera_rate = mean(chimera_rate), .groups = "drop")

# Synthesis: one table combining every cross-talk-relevant signal per mock --
# mock-to-mock cross-talk (above), the ASV-level excess as a rate (T3, not
# just a raw count), chimera rate (this run's DADA2 output, not mock-
# specific contamination), and T5's already-computed mock-to-real bleed
# floor for completeness.
crosstalk_index <- mock_crosstalk |>
  dplyr::select(mock, sample_id, library_size, crosstalk_reads, crosstalk_frac) |>
  dplyr::left_join(
    asv_fpr |> dplyr::transmute(mock, asv_fpr_pct = 100 * excess_asvs / n_asv_observed),
    by = "mock"
  ) |>
  dplyr::left_join(
    chimera_rates |> dplyr::filter(sample_id %in% mock_sample_for) |> dplyr::select(sample_id, chimera_rate),
    by = "sample_id"
  ) |>
  dplyr::mutate(
    crosstalk_pct = 100 * crosstalk_frac,
    chimera_pct = 100 * chimera_rate,
    mock_to_real_bleed_floor_pct = 100 * bleed_floor
  ) |>
  dplyr::select(mock, sample_id, library_size, crosstalk_pct, asv_fpr_pct, chimera_pct, mock_to_real_bleed_floor_pct)
write_result(crosstalk_index, "02b_crosstalk_index")

p_crosstalk_index <- crosstalk_index |>
  tidyr::pivot_longer(c(crosstalk_pct, asv_fpr_pct, chimera_pct), names_to = "metric", values_to = "value") |>
  dplyr::mutate(metric = factor(metric, levels = c("crosstalk_pct", "chimera_pct", "asv_fpr_pct"),
                                  labels = c("Mock-to-mock\ncross-talk (%)", "Chimera rate (%)", "ASV excess\nrate (%)"))) |>
  ggplot2::ggplot(ggplot2::aes(mock, value, fill = mock)) +
  ggplot2::geom_col() +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f%%", value)), vjust = -0.3, size = 3) +
  ggplot2::facet_wrap(~metric, scales = "free_y") +
  ggplot2::labs(x = NULL, y = NULL, title = "Mock cross-talk index",
                subtitle = "Mock-to-mock cross-talk is the direct index-hopping estimate; chimera / ASV-excess rates are context, not cross-talk") +
  ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
save_plot(p_crosstalk_index, "02b_crosstalk_index", w = 11, h = 4.5)

message(sprintf(
  "[02b_controls] T6 cross-talk index: mock-to-mock = %s; chimera rate mocks=%.2f%% vs real=%.2f%% vs blank=%.2f%%; mock-to-real bleed floor (T5) = %.3f%%. See results/02b_crosstalk_index.tsv.",
  paste(sprintf("%s=%.3f%%", crosstalk_index$mock, crosstalk_index$crosstalk_pct), collapse = ", "),
  chimera_summary$mean_chimera_rate[chimera_summary$group == "mock"] * 100,
  chimera_summary$mean_chimera_rate[chimera_summary$group == "real"] * 100,
  chimera_summary$mean_chimera_rate[chimera_summary$group == "blank"] * 100,
  bleed_floor * 100
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

# ============================================================================
# Replicate concordance -- both levels (moved here from 02_qc_filter.R: needs
# the post-blacklist/bleed-floor/control-drop matrix). Runs on `counts_noctrl`
# (NOT `counts_filtered`) DELIBERATELY, and BEFORE the prevalence filter
# below -- moved here specifically so its output (`empirical_read_floor`)
# is available to gate the prevalence filter itself (see Part 3a of
# golden-napping-breeze.md: "run this filter after the read floor, not
# before" -- a taxon clearing the gate on a detection in a 24-read library
# is exactly what this ordering prevents). This also makes the concordance
# diagnostic itself more correct, incidentally: it's no longer truncated by
# a prevalence filter that, at the point it used to run, hadn't been
# justified by anything yet.
# ============================================================================
rel_ab <- counts_noctrl / rowSums(counts_noctrl)

# One row per pair of samples (not built separately for within vs. between --
# a single table, filtered two ways, is easier to follow and to check).
# "Within" = the two samples share `pair_key` (the replicate-pair identifier
# -- `site` for technical, `location`+`tech_rep` for biological). "Between"
# is the reference/noise-floor distribution a within-pair distance gets
# judged against, restricted to pairs of the SAME sample_type: a
# sediment-vs-water pair is trivially, definitionally dissimilar and says
# nothing about whether a same-type replicate pair reproduced, so it doesn't
# belong in that reference distribution.
concordance_at <- function(rel_ab, counts, metadata, pair_key, level_label, rarefied_lookup) {
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
      same_type = type_a == type_b,
      min_lib = pmin(libsize[sample_a], libsize[sample_b])
    ) |>
    # rarefied_lookup was precomputed once for every same-sample_type pair
    # (cross-type pairs are never within- or between-pairs, so were never
    # computed) -- joined here rather than recomputed per level_label call,
    # since it depends only on which two samples are paired, not which
    # replicate-pair key grouped them (see precompute block above this
    # function).
    dplyr::left_join(rarefied_lookup, by = c("sample_a", "sample_b"))

  within_pair <- pairs |>
    dplyr::filter(is_within) |>
    dplyr::transmute(level, sample_type = type_a, pair_id = key_a, sample_a, sample_b, bray, bray_rarefied, jaccard,
                      spearman_cor, delta_richness, delta_library_size, min_lib)
  between_pair <- pairs |>
    dplyr::filter(!is_within, same_type) |>
    dplyr::transmute(level, sample_type = type_a, sample_a, sample_b, bray, bray_rarefied, jaccard, min_lib)

  list(within = within_pair, between = between_pair)
}

# Pairwise-rarefied Bray-Curtis (Part 2 of the mock cross-talk / depth-honest
# concordance / genus-collapse plan): the full-depth `bray` column above
# conflates real biological/technical discordance with pure sampling-depth
# artifact -- a pair with one 24-read side and one 400,000-read side looks
# "dissimilar" partly just because the shallow side is a noisy subsample,
# not because the underlying communities differ. Fix: repeatedly subsample
# BOTH samples in a pair down to their own shared min(lib_a, lib_b) and
# average Bray-Curtis across reps, so every pair is compared on equal
# footing at its own achievable depth. Computed ONCE for every same-
# sample_type pair (542 of the 1035 total C(46,2) pairs -- cross-type pairs
# are never used as either a within- or between-pair) and shared by both the
# technical and biological concordance_at() calls below, since it depends
# only on which two samples are paired, not on which pair_key grouped them.
pairwise_rarefied_bray <- function(count_a, count_b, n_reps = 30) {
  depth <- min(sum(count_a), sum(count_b))
  if (depth < 1) return(NA_real_)
  both <- rbind(count_a, count_b)
  vals <- vapply(seq_len(n_reps), function(k) {
    sub <- vegan::rrarefy(both, sample = depth)
    as.numeric(vegan::vegdist(sub, method = "bray"))
  }, numeric(1))
  mean(vals)
}

ids_all <- rownames(counts_noctrl)
type_all <- metadata_noctrl$sample_type[match(ids_all, metadata_noctrl$sample_id)]
idx_all <- utils::combn(seq_along(ids_all), 2)
same_type_mask <- type_all[idx_all[1, ]] == type_all[idx_all[2, ]]
idx_st <- idx_all[, same_type_mask, drop = FALSE]

# Each pair is independent, so this is an embarrassingly parallel loop.
# mclapply() forks (Linux-only, fine here; copy-on-write means
# counts_noctrl isn't duplicated per worker) rather than needing a cluster.
# (Runs on counts_noctrl now, not the smaller post-prevalence-filter
# counts_filtered -- more columns than the single-threaded 209s benchmark
# this comment used to cite, hence parallelized rather than left serial.)
mc_cores <- max(1L, parallel::detectCores() - 1L)
t_rarefy_start <- Sys.time()
rarefied_lookup <- tibble::tibble(
  sample_a = ids_all[idx_st[1, ]],
  sample_b = ids_all[idx_st[2, ]]
)
rarefied_lookup$bray_rarefied <- unlist(parallel::mclapply(
  seq_len(nrow(rarefied_lookup)),
  function(i) pairwise_rarefied_bray(counts_noctrl[rarefied_lookup$sample_a[i], ], counts_noctrl[rarefied_lookup$sample_b[i], ], n_reps = 30),
  mc.cores = mc_cores
))
message(sprintf("[02b_controls] pairwise-rarefied Bray-Curtis: %d same-type pairs x 30 reps in %.1fs (%d cores)",
                 nrow(rarefied_lookup), as.numeric(Sys.time() - t_rarefy_start, units = "secs"), mc_cores))

tech_conc <- concordance_at(rel_ab, counts_noctrl, metadata_noctrl, pair_key = "site", level_label = "technical", rarefied_lookup = rarefied_lookup)
metadata_biopair <- metadata_noctrl |> dplyr::mutate(bio_pair_key = paste(location, tech_rep, sep = "_"))
bio_conc <- concordance_at(rel_ab, counts_noctrl, metadata_biopair, pair_key = "bio_pair_key", level_label = "biological", rarefied_lookup = rarefied_lookup)

write_result(tech_conc$within, "replicate_concordance_technical")
write_result(bio_conc$within, "replicate_concordance_biological")

# Threshold is per sample_type now that between-pair is too (see
# concordance_at() above) -- a technical pair's discordance call should be
# judged against ITS OWN type's between-pair noise floor, not a floor
# blended from both sediment's and water's between-distances together.
flag_discordant <- function(within, between) {
  thresh_by_type <- between |>
    dplyr::group_by(sample_type) |>
    dplyr::summarise(thresh = stats::quantile(bray, 0.05, na.rm = TRUE), .groups = "drop")
  within |>
    dplyr::left_join(thresh_by_type, by = "sample_type") |>
    dplyr::mutate(discordant = bray >= thresh) |>
    dplyr::select(-thresh)
}
tech_flagged <- flag_discordant(tech_conc$within, tech_conc$between)
bio_flagged <- flag_discordant(bio_conc$within, bio_conc$between)

# between_pair (see concordance_at() above) comes from ALL pairwise
# comparisons under that level's key, so a (level, sample_type) combo with
# zero real within-pairs -- e.g. "biological" x water, since water has no
# second bio_rep arm to ever be a replicate of -- still produces between-pair
# rows (two water samples are trivially "not a replicate pair" under a key
# neither of them was ever eligible to share). That's a real number but not
# a meaningful one here: there's no "did the replicate reproduce" question
# to judge it against without a within-pair to compare it to, so it's
# dropped from the figure via this semi_join rather than shown as an
# orphaned box.
valid_combos <- dplyr::bind_rows(tech_conc$within, bio_conc$within) |> dplyr::distinct(level, sample_type)

concordance_all <- dplyr::bind_rows(
  tech_conc$within |> dplyr::mutate(pair_type = "within-pair"), tech_conc$between |> dplyr::mutate(pair_type = "between-pair"),
  bio_conc$within |> dplyr::mutate(pair_type = "within-pair"), bio_conc$between |> dplyr::mutate(pair_type = "between-pair")
) |>
  dplyr::semi_join(valid_combos, by = c("level", "sample_type")) |>
  # Explicit factor, not left as character: water's panel only ever has
  # level == "technical" rows (no second bio_rep arm to pair), and
  # scale_shape_manual(..., drop = FALSE) below only preserves an unused
  # FACTOR level -- a character column with a value simply absent has no
  # "unused level" for drop=FALSE to retain, so the two panels' shape
  # scales still silently disagreed without this and patchwork couldn't
  # merge them into one "Level" legend.
  dplyr::mutate(level = factor(level, levels = c("technical", "biological")))
# min_lib is already carried through from concordance_at() above -- the
# smaller of the pair's two library sizes, the bottleneck: a 24-read vs.
# 400,000-read "pair" is limited by the 24-read side, and this is exactly
# what makes many within-pair distances look artificially high (Spearman
# r = -0.86 between min library size and technical within-pair Bray
# distance -- shallow libraries are noisy, not necessarily discordant).

# ----------------------------------------------------------------------
# Subsampling null band: "how dissimilar do two subsamples of ONE real,
# identical community look, from pure sampling noise alone, at depth N."
# Built from the single deepest real sediment sample and the single
# deepest real water sample (the two "true communities" this dataset can
# most confidently subsample down from) -- at each of a log-spaced series
# of depths, draw two INDEPENDENT rrarefy() subsamples from that same
# reference count vector and compute Bray-Curtis between them. This is a
# null distribution, not a model: at low depth it should be high (two
# small subsamples of the same community still look different by chance),
# and it should fall toward 0 as depth approaches the reference's own
# total (a full-depth "subsample" of itself is just itself twice).
# ----------------------------------------------------------------------
ref_sediment_id <- names(which.max(rowSums(counts_noctrl[metadata_noctrl$sample_id[metadata_noctrl$sample_type == "sediment"], , drop = FALSE])))
ref_water_id <- names(which.max(rowSums(counts_noctrl[metadata_noctrl$sample_id[metadata_noctrl$sample_type == "water"], , drop = FALSE])))
ref_samples <- c(sediment = ref_sediment_id, water = ref_water_id)

null_band_at_depth <- function(ref_counts, depth, n_reps = 50) {
  if (depth > sum(ref_counts)) return(rep(NA_real_, n_reps))
  vapply(seq_len(n_reps), function(k) {
    s1 <- vegan::rrarefy(ref_counts, sample = depth)
    s2 <- vegan::rrarefy(ref_counts, sample = depth)
    as.numeric(vegan::vegdist(rbind(s1, s2), method = "bray"))
  }, numeric(1))
}

depth_seq <- unique(round(10^seq(log10(20), log10(400000), length.out = 24)))

# One (sample_type, depth) combination per task -- 48 tasks across mc_cores
# (set above, next to the pairwise-rarefied Bray-Curtis parallelization),
# each still doing its own 50 reps internally.
null_tasks <- tidyr::expand_grid(st = names(ref_samples), d = depth_seq)

t_null_start <- Sys.time()
null_band <- dplyr::bind_rows(parallel::mclapply(seq_len(nrow(null_tasks)), function(i) {
  st <- null_tasks$st[i]; d <- null_tasks$d[i]
  ref_counts <- counts_noctrl[ref_samples[[st]], ]
  vals <- null_band_at_depth(ref_counts, d, n_reps = 50)
  if (all(is.na(vals))) return(NULL) # depth exceeds this reference's own total -- can't subsample past it
  tibble::tibble(sample_type = st, ref_sample_id = ref_samples[[st]], depth = d, bray_null = vals)
}, mc.cores = mc_cores))
message(sprintf("[02b_controls] subsampling null band: %d depths x 50 reps x 2 references in %.1fs (%d cores; refs: sediment=%s @ %d reads, water=%s @ %d reads)",
                 length(depth_seq), as.numeric(Sys.time() - t_null_start, units = "secs"), mc_cores,
                 ref_sediment_id, sum(counts_noctrl[ref_sediment_id, ]), ref_water_id, sum(counts_noctrl[ref_water_id, ])))
write_result(null_band, "02_replicate_null_band_draws")

null_summary <- null_band |>
  dplyr::group_by(sample_type, depth) |>
  dplyr::summarise(median = stats::median(bray_null), q10 = stats::quantile(bray_null, 0.10),
                    q25 = stats::quantile(bray_null, 0.25), q75 = stats::quantile(bray_null, 0.75),
                    q90 = stats::quantile(bray_null, 0.90), .groups = "drop")
write_result(null_summary, "02_replicate_null_band_summary")
stopifnot(
  "null band should fall as depth increases (sanity check on subsampling mechanics)" =
    with(null_summary[null_summary$sample_type == "sediment", ], median[which.max(depth)] < median[which.min(depth)]),
  "null band should rise toward the shallow end" =
    with(null_summary[null_summary$sample_type == "sediment", ], median[which.min(depth)] > 0.05)
)

# Empirical read floor: at each real within-pair's own min_lib, interpolate
# the null band's median/q10/q90 (linear in log10-depth) for that pair's
# sample_type, and ask whether the pair's own pairwise-rarefied Bray-Curtis
# distance is statistically indistinguishable from pure subsampling noise
# (inside [q10, q90]) or clearly separated from it (above q90 -- real
# discordance beyond what depth alone explains). This is what converts the
# figure from "replicates look bad" into an actual depth threshold.
interp_null <- function(min_lib, st) {
  nb <- null_summary[null_summary$sample_type == st, ]
  nb <- nb[order(nb$depth), ]
  list(
    median = stats::approx(log10(nb$depth), nb$median, xout = log10(min_lib), rule = 2)$y,
    q10 = stats::approx(log10(nb$depth), nb$q10, xout = log10(min_lib), rule = 2)$y,
    q90 = stats::approx(log10(nb$depth), nb$q90, xout = log10(min_lib), rule = 2)$y
  )
}
within_vs_null <- dplyr::bind_rows(tech_conc$within, bio_conc$within) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    null_median = interp_null(min_lib, sample_type)$median,
    null_q10 = interp_null(min_lib, sample_type)$q10,
    null_q90 = interp_null(min_lib, sample_type)$q90,
    indistinguishable_from_noise = bray_rarefied <= null_q90
  ) |>
  dplyr::ungroup()
write_result(within_vs_null, "02_replicate_vs_null_band")

# Empirical floor = the largest min_lib among within-pairs that are still
# indistinguishable from pure subsampling noise at their own depth -- below
# it, "this pair didn't reproduce" and "this pair is just shallow" can't be
# told apart; above it, a pair with bray_rarefied clearly above the null
# band is showing real discordance, not a depth artifact. Reported here,
# not silently substituted for the existing 500-read shallow_floor
# (02_qc_filter.R) or 499-read rarefaction depth (03_normalize.R) -- see
# AGENTS.md for the comparison and the stop-and-ask decision.
noisy_pairs <- within_vs_null |> dplyr::filter(indistinguishable_from_noise)
empirical_read_floor <- if (nrow(noisy_pairs) > 0) max(noisy_pairs$min_lib) else NA_real_
n_below_floor_pairs <- sum(within_vs_null$min_lib <= empirical_read_floor, na.rm = TRUE)
message(sprintf(
  "[02b_controls] empirical read floor: %s reads (largest min_lib among within-pairs statistically indistinguishable from subsampling noise; %d/%d within-pairs are at or below it). Existing floors: shallow_floor=500 (02_qc_filter.R), rarefaction_depth=499 (03_normalize.R).",
  format(empirical_read_floor, big.mark = ","), n_below_floor_pairs, nrow(within_vs_null)
))

# Technical and biological pairs tell different stories once separated --
# a single pooled floor hides that. Technical pairs (same site, two
# independent extractions) are expected to converge toward the null band as
# depth increases; biological pairs (same location/tech_rep, different arm
# -- transect vs. isolate_source) are expected to stay elevated above the
# null band at every depth, since they're deliberately different physical
# samples, not duplicate extractions -- that's real spatial heterogeneity,
# not evidence of a reproducibility failure.
floor_by_level <- within_vs_null |>
  dplyr::group_by(level) |>
  dplyr::summarise(
    n_pairs = dplyr::n(),
    n_indistinguishable_from_noise = sum(indistinguishable_from_noise),
    max_min_lib_indistinguishable = if (any(indistinguishable_from_noise)) max(min_lib[indistinguishable_from_noise]) else NA_real_,
    min_min_lib_discordant = if (any(!indistinguishable_from_noise)) min(min_lib[!indistinguishable_from_noise]) else NA_real_,
    .groups = "drop"
  )
write_result(floor_by_level, "02_replicate_floor_by_level")
message(sprintf(
  "[02b_controls] by level: %s",
  paste(sprintf("%s: %d/%d indistinguishable from noise (floor %s reads)",
                floor_by_level$level, floor_by_level$n_indistinguishable_from_noise, floor_by_level$n_pairs,
                format(floor_by_level$max_min_lib_indistinguishable, big.mark = ",")),
        collapse = "; ")
))

size_breaks <- scales::breaks_log(n = 5)(range(concordance_all$min_lib))
size_limits <- range(concordance_all$min_lib)

# ONE ggplot object, faceted by sample_type, rather than two separately-
# constructed panels combined via patchwork -- two attempted fixes
# (declaring `level` as an explicit factor with both levels; forcing one
# literal shared scale object into both panels) still left patchwork's
# guides="collect" rendering two separate "Level" legend blocks instead of
# merging them (a real patchwork guide-merge fragility, not further chased
# here). A single faceted plot has exactly one shape scale and one legend
# by construction, sidestepping the merge question entirely -- and
# size_limits/breaks above are already a single shared range across both
# sample_types, so facet_wrap's default fixed x-scale reproduces the same
# "identical x-axis in both panels" behavior the two-panel version had.
# null_summary's own sample_type column lines the ribbon/line up with the
# matching facet automatically.
p_concordance <- ggplot2::ggplot(concordance_all, ggplot2::aes(x = min_lib, y = bray_rarefied)) +
  ggplot2::geom_ribbon(data = null_summary, ggplot2::aes(x = depth, ymin = q10, ymax = q90), inherit.aes = FALSE,
                        fill = "grey60", alpha = 0.35) +
  ggplot2::geom_line(data = null_summary, ggplot2::aes(x = depth, y = median), inherit.aes = FALSE,
                      color = "grey30", linewidth = 0.8) +
  ggplot2::geom_point(ggplot2::aes(color = pair_type, shape = level), size = 2.4, alpha = 0.75) +
  ggplot2::scale_shape_manual(values = c(technical = 17, biological = 16), drop = FALSE) +
  ggplot2::scale_x_log10(labels = scales::label_comma(), breaks = size_breaks, limits = size_limits) +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  ggplot2::facet_wrap(~sample_type) +
  ggplot2::labs(x = "min(library size) of the pair", y = "Bray-Curtis distance (pairwise-rarefied)",
                color = "Pair type", shape = "Level",
                title = "Depth-honest replicate concordance: pairwise-rarefied Bray-Curtis vs. subsampling null band",
                subtitle = paste(
                  sprintf(
                    "Grey ribbon/line: null dist. (median, 10th-90th pctile) of Bray-Curtis between 2 INDEPENDENT subsamples of ONE deep reference (sediment=%s, water=%s) -- pure sampling noise only.",
                    ref_sediment_id, ref_water_id
                  ),
                  "Points: real within-/between-pairs, pairwise-rarefied to their own min(lib_a,lib_b) and plotted at that depth on the same x-axis.",
                  sprintf("A point inside the ribbon is statistically indistinguishable from subsampling noise at that depth. Empirical read floor (largest min_lib still indistinguishable from noise): %s reads.",
                          format(empirical_read_floor, big.mark = ",")),
                  sep = "\n"
                )) +
  ggplot2::theme(
    axis.text = ggplot2::element_text(size = 12),
    axis.title = ggplot2::element_text(size = 15),
    plot.title = ggplot2::element_text(size = 20, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 10.5),
    strip.text = ggplot2::element_text(size = 16, face = "bold"),
    legend.text = ggplot2::element_text(size = 12),
    legend.title = ggplot2::element_text(size = 13),
    panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1)
  )
save_plot(p_concordance, "02_replicate_concordance", w = 14, h = 8.5, dpi = 400)

n_discordant_tech <- sum(tech_flagged$discordant, na.rm = TRUE)
n_discordant_bio <- sum(bio_flagged$discordant, na.rm = TRUE)
message(sprintf("[02b_controls] replicate concordance: %d/%d technical pairs discordant, %d/%d biological pairs discordant",
                n_discordant_tech, nrow(tech_flagged), n_discordant_bio, nrow(bio_flagged)))

metadata_flagged <- metadata_noctrl |>
  dplyr::mutate(flag_discordant_pair = sample_id %in% c(
    unlist(tech_flagged |> dplyr::filter(discordant) |> dplyr::select(sample_a, sample_b)),
    unlist(bio_flagged |> dplyr::filter(discordant) |> dplyr::select(sample_a, sample_b))
  ))

pooling_decision <- dplyr::bind_rows(tech_flagged, bio_flagged) |>
  dplyr::group_by(level, sample_type) |>
  dplyr::summarise(n_pairs = dplyr::n(), n_discordant = sum(discordant, na.rm = TRUE), .groups = "drop") |>
  dplyr::mutate(
    pooled_by_default = FALSE,
    note = "not pooled by default -- see plots/02_replicate_concordance.png and confirm PLAN.md before changing this"
  )
write_result(pooling_decision, "replicate_concordance_decision")

# ============================================================================
# Prevalence filter -- floor-gated, within-sample_type, read-count-aware
# (Part 3a of golden-napping-breeze.md). Supersedes the old pooled
# "prevalence >= 2 of 46 samples, >0 reads" gate for three independent
# reasons, all real problems with the old rule, not stylistic preferences:
#  1. WITHIN TYPE: sediment and water share almost no taxa (confirmed
#     throughout this session) -- a pooled count let e.g. "1 sediment read +
#     1 water read" pass the >=2-samples bar for a reason that has nothing
#     to do with either community.
#  2. READ-COUNT-AWARE: a single read is indistinguishable from a
#     cross-talk/error event -- Part 1 measured mock-to-mock cross-talk
#     directly in THIS dataset at 0.09-0.21% of reads, and the mock ASV
#     excess (mock_log: 578 observed vs. ~16-40 expected genomes) is exactly
#     what single-read/rRNA-copy-variant noise looks like. >=2 reads is
#     required before a sample "votes" for an ASV's presence.
#  3. FLOOR-GATED: only samples at/above a read-count floor get to vote at
#     all -- a taxon can no longer clear the gate on the strength of a
#     detection in a 24-200-read library that contributes almost no real
#     information.
# Applied to BOTH sample types identically (per-type eligible-sample counts
# differ sharply -- confirmed sediment >> water -- but the floor itself is
# not softened for water; see AGENTS.md).
#
# NOTE on the floor value: `empirical_read_floor` (computed just above, on
# this same counts_noctrl -- correctly ordered, per this exact section's own
# fix) measured at only 78 reads, not the ~8,272 first estimated in Part 2.
# That's not noise: the 8,272 figure was itself computed on the
# already-prevalence-filtered matrix (the exact ordering violation this
# section exists to fix) -- on the full, unfiltered ASV set, the deepest
# previously-"indistinguishable" technical pair (C9, 8,272 reads) turned out
# to show real discordance once every ASV was counted, not just the ones
# that survive prevalence filtering. 78 reads is the more honest number, but
# it makes the floor-gating dimension of this filter nearly a no-op (only
# 2/36 sediment, 0/10 water samples excluded) -- a much weaker gate than
# intended. Per direct instruction, `prevalence_eligibility_floor` is set
# manually to 8,000 (not derived from `empirical_read_floor`) to keep the
# gate as strict as originally intended; `empirical_read_floor` is still
# computed and reported above as its own, separately meaningful finding.
prevalence_eligibility_floor <- 8000L
min_reads_present <- 2L
min_samples_per_type <- 2L
lib_size_noctrl <- rowSums(counts_noctrl)
eligible <- lib_size_noctrl >= prevalence_eligibility_floor
type_noctrl <- metadata_noctrl$sample_type[match(rownames(counts_noctrl), metadata_noctrl$sample_id)]
is_present_strong <- counts_noctrl >= min_reads_present

sediment_elig <- eligible & type_noctrl == "sediment"
water_elig <- eligible & type_noctrl == "water"
prevalence_sediment <- colSums(is_present_strong[sediment_elig, , drop = FALSE])
prevalence_water <- colSums(is_present_strong[water_elig, , drop = FALSE])
keep_prevalence <- (prevalence_sediment >= min_samples_per_type) | (prevalence_water >= min_samples_per_type)

counts_filtered <- counts_noctrl[, keep_prevalence, drop = FALSE]
taxonomy_filtered <- taxonomy_deblacklisted[keep_prevalence, ]

# Old-vs-new comparison, reported directly rather than buried -- this is a
# big, visible change to the final ASV table, not an implementation detail.
keep_prevalence_old <- colSums(counts_noctrl > 0) >= 2L
message(sprintf(
  "[02b_controls] final gate: %d controls dropped, blacklist -%d ASVs, prevalence filter (floor-gated >=%s reads [manual -- measured empirical_read_floor was %s], within-type, >=%d/%d samples): sediment eligible=%d/%d samples, water eligible=%d/%d samples -> %d/%d ASVs kept (was %d/%d under the old pooled >=2-of-%d-samples, >0-read rule) -> %d samples x %d ASVs, %d reads",
  sum(is_mock | is_neg), length(blacklist), format(prevalence_eligibility_floor, big.mark = ","), format(empirical_read_floor, big.mark = ","), min_samples_per_type, min_samples_per_type,
  sum(sediment_elig), sum(type_noctrl == "sediment"), sum(water_elig), sum(type_noctrl == "water"),
  sum(keep_prevalence), length(keep_prevalence), sum(keep_prevalence_old), length(keep_prevalence_old), nrow(metadata_noctrl),
  nrow(counts_filtered), ncol(counts_filtered), sum(counts_filtered)
))

stopifnot("control row(s) survived to counts_clean/metadata_clean" = !any(metadata_noctrl$sample_type == "control"))

# library_size / log_library_size: a first-class metadata column from here
# on, computed once from the final matrix and joined by sample_id (not
# positional), so every downstream script (06/07/08) picks up the exact same
# numbers automatically instead of re-deriving rowSums() independently in
# each one. Library size is a confirmed confound in this dataset -- the
# replicate-concordance check earlier found technical-pair Bray distance
# strongly anti-correlated with the pair's minimum library size (Spearman
# r = -0.86) -- so it's reported as a covariate (or restriction) in every
# community-level model from here on; see AGENTS.md.
lib_size_by_sample <- rowSums(counts_filtered)
metadata_flagged <- metadata_flagged |>
  dplyr::mutate(
    library_size = lib_size_by_sample[sample_id],
    log_library_size = log10(library_size)
  )

# --- save the objects 03_normalize.R onward actually read -------------------
assert_aligned(counts_filtered, taxonomy_filtered, metadata_flagged)
saveRDS(counts_filtered, file.path(path_processed, "counts_clean.rds"))
saveRDS(taxonomy_filtered, file.path(path_processed, "taxonomy_clean.rds"))
saveRDS(metadata_flagged, file.path(path_processed, "metadata_clean.rds"))

# ============================================================================
# Genus-level primary aggregation (Part 3c of golden-napping-breeze.md) --
# built from the Part-3a-corrected counts_filtered/taxonomy_filtered, and
# justified by Part 3b's measured evidence, not an assumption: ASVs sharing
# an identical species-level label show a real, statistically clear excess
# of co-occurrence/correlation over a cross-species null (median Jaccard
# 0.40 vs 0.00, median Spearman-cor 0.60 vs -0.05, both p~0 -- see
# AGENTS.md), alongside a genuine population of distinct organisms sharing
# a coarse label. Genus-level becomes the primary unit for alpha diversity,
# PERMANOVA/ordination/betadisper, and db-RDA/Mantel (04/06/07) -- the
# community-structure claims most exposed to ASV-level noise; ASV-level
# results are kept as supplementary throughout. Differential abundance (08)
# and FAPROTAX (11) stay ASV-level (scope decision, golden-napping-breeze.md
# Context section).
# ============================================================================
genus_na <- is.na(taxonomy_filtered$genus_norm)
genus_na_reads <- sum(counts_filtered[, taxonomy_filtered$asv_id[genus_na], drop = FALSE])
message(sprintf(
  "[02b_controls] genus aggregation: %d/%d ASVs have no genus (genus_norm NA), excluding %d reads (%.3f%% of %d total) from the genus-level table",
  sum(genus_na), length(genus_na), genus_na_reads, 100 * genus_na_reads / sum(counts_filtered), sum(counts_filtered)
))

# rowsum() sums ROWS sharing a group; counts_filtered is samples x ASVs, so
# transpose to ASVs x samples first (group ASVs into genera, sum reads),
# then transpose back to samples x genera to match counts_filtered's own
# orientation.
counts_genus_clean <- t(rowsum(t(counts_filtered[, !genus_na, drop = FALSE]), taxonomy_filtered$genus_norm[!genus_na]))
storage.mode(counts_genus_clean) <- "integer"

taxonomy_genus_clean <- taxonomy_filtered |>
  dplyr::filter(!is.na(genus_norm)) |>
  dplyr::group_by(genus_norm) |>
  dplyr::summarise(
    domain = dplyr::first(domain), phylum = dplyr::first(phylum), class = dplyr::first(class),
    order = dplyr::first(order), family = dplyr::first(family), n_asvs = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::rename(genus = genus_norm)

stopifnot(
  "counts_genus_clean row sums must equal counts_filtered row sums minus only genus-NA reads" =
    isTRUE(all.equal(unname(rowSums(counts_genus_clean)), unname(rowSums(counts_filtered) - rowSums(counts_filtered[, genus_na, drop = FALSE])))),
  "counts_genus_clean columns must match taxonomy_genus_clean rows (as a set)" =
    setequal(colnames(counts_genus_clean), taxonomy_genus_clean$genus)
)

saveRDS(counts_genus_clean, file.path(path_processed, "counts_genus_clean.rds"))
saveRDS(taxonomy_genus_clean, file.path(path_processed, "taxonomy_genus_clean.rds"))
message(sprintf("[02b_controls] genus-level table: %d samples x %d genera (ASV:genus ratio %.1f:1)",
                nrow(counts_genus_clean), ncol(counts_genus_clean), ncol(counts_filtered) / ncol(counts_genus_clean)))

retention <- tibble::tibble(
  stage = c("postlineage", "blacklist_removed", "controls_dropped", "prevalence_filter"),
  n_asv = c(ncol(counts), ncol(counts_deblacklisted), ncol(counts_noctrl), ncol(counts_filtered)),
  n_samples = c(nrow(counts), nrow(counts_deblacklisted), nrow(counts_noctrl), nrow(counts_filtered)),
  n_reads = c(sum(counts), sum(counts_deblacklisted), sum(counts_noctrl), sum(counts_filtered))
)
write_result(retention, "02b_retention")

message(sprintf("[02b_controls] final: %d samples x %d ASVs, %d reads (0 controls)",
                nrow(counts_filtered), ncol(counts_filtered), sum(counts_filtered)))

# ============================================================================
# Intragenomic-variant test (Part 3b of golden-napping-breeze.md) -- the
# direct test for genus-level aggregation's actual premise, run here (on the
# CORRECTED counts_filtered/taxonomy_filtered from Part 3a) rather than
# assumed. Question: do ASVs sharing an identical species-level GTDB
# assignment co-occur across samples with correlated abundance (the
# signature of intragenomic rRNA-operon copies of ONE genome, which exist in
# a fixed ratio) or not (the signature of distinct organisms that happen to
# share a coarse nearest-reference label)? Reconnaissance earlier this
# session ran this same method on the OLD, uncorrected ASV table -- not
# reused here, since that table was itself contaminated by the ordering
# problem Part 3a fixed.
# ============================================================================
sp_filtered <- taxonomy_filtered$species_norm
sp_tab_filtered <- table(sp_filtered[!is.na(sp_filtered)])
shared_species_filtered <- names(sp_tab_filtered[sp_tab_filtered > 1])

# Cap group size at 20 ASVs (random subsample if larger) so one big group
# (GMQP-bins7_sp004366385 alone had 180 ASVs pre-3a) can't dominate the
# pair count -- C(20,2)=190 pairs is still plenty per group.
max_group_size <- 20L
group_asvs <- purrr::map(shared_species_filtered, function(s) {
  asvs <- taxonomy_filtered$asv_id[!is.na(sp_filtered) & sp_filtered == s]
  if (length(asvs) > max_group_size) asvs <- sample(asvs, max_group_size)
  asvs
})
names(group_asvs) <- shared_species_filtered

same_species_pairs <- purrr::imap_dfr(group_asvs, function(asvs, s) {
  if (length(asvs) < 2) return(NULL)
  idx <- utils::combn(length(asvs), 2)
  tibble::tibble(species = s, asv_a = asvs[idx[1, ]], asv_b = asvs[idx[2, ]])
})

# Null: an equal-sized draw of random ASV pairs from DIFFERENT species
# labels (same species-labeled pool, cross-label only) -- not from the full
# ASV set, so the comparison isn't confounded by unclassified-vs-classified
# differences, only by same-label vs. different-label.
species_labeled_asvs <- taxonomy_filtered$asv_id[!is.na(sp_filtered)]
asv_to_species <- stats::setNames(sp_filtered[!is.na(sp_filtered)], species_labeled_asvs)
set.seed(42)
null_pairs <- tibble::tibble(
  asv_a = sample(species_labeled_asvs, nrow(same_species_pairs), replace = TRUE),
  asv_b = sample(species_labeled_asvs, nrow(same_species_pairs), replace = TRUE)
) |>
  dplyr::filter(asv_to_species[asv_a] != asv_to_species[asv_b]) # re-drawn pairs that landed same-species by chance are simply dropped, not resampled -- negligible given >1900 species labels

jaccard_cor_pair <- function(a, b) {
  ca <- counts_filtered[, a]; cb <- counts_filtered[, b]
  pa <- ca > 0; pb <- cb > 0
  inter <- sum(pa & pb); uni <- sum(pa | pb)
  jacc <- if (uni == 0) NA_real_ else inter / uni
  cor_ab <- if (stats::sd(ca) == 0 || stats::sd(cb) == 0) NA_real_ else suppressWarnings(stats::cor(ca, cb, method = "spearman"))
  c(jaccard = jacc, spearman_cor = cor_ab)
}

t_intragenomic_start <- Sys.time()
same_species_result <- dplyr::bind_cols(
  same_species_pairs,
  as.data.frame(t(simplify2array(parallel::mclapply(
    seq_len(nrow(same_species_pairs)),
    function(i) jaccard_cor_pair(same_species_pairs$asv_a[i], same_species_pairs$asv_b[i]),
    mc.cores = mc_cores
  ))))
) |> dplyr::mutate(group = "same_species_label")
null_result <- dplyr::bind_cols(
  null_pairs,
  as.data.frame(t(simplify2array(parallel::mclapply(
    seq_len(nrow(null_pairs)),
    function(i) jaccard_cor_pair(null_pairs$asv_a[i], null_pairs$asv_b[i]),
    mc.cores = mc_cores
  ))))
) |> dplyr::mutate(species = NA_character_, group = "null_cross_species")
message(sprintf("[02b_controls] intragenomic-variant test: %d same-species-label pairs (%d groups) + %d null pairs in %.1fs",
                 nrow(same_species_result), length(shared_species_filtered), nrow(null_result),
                 as.numeric(Sys.time() - t_intragenomic_start, units = "secs")))

intragenomic_test <- dplyr::bind_rows(same_species_result, null_result)
write_result(intragenomic_test, "02b_intragenomic_variant_test")

# Formal test: Wilcoxon rank-sum, same-species-label pairs vs. null, for
# both metrics -- statistic, p-value, and rank-biserial effect size, not
# just "elevated vs. null" by eyeball.
# Standard convention: positive means the FIRST wilcox.test() argument
# (same_species_result) is stochastically greater than the second (null) --
# matches the direction the hypothesis actually predicts, so a positive
# number here means "supports intragenomic variance," not the reverse.
rank_biserial <- function(w_stat, n1, n2) (2 * w_stat) / (n1 * n2) - 1
wilcox_jaccard <- stats::wilcox.test(same_species_result$jaccard, null_result$jaccard)
wilcox_cor <- stats::wilcox.test(same_species_result$spearman_cor, null_result$spearman_cor)
intragenomic_summary <- tibble::tibble(
  metric = c("jaccard", "spearman_cor"),
  median_same_species = c(stats::median(same_species_result$jaccard, na.rm = TRUE), stats::median(same_species_result$spearman_cor, na.rm = TRUE)),
  median_null = c(stats::median(null_result$jaccard, na.rm = TRUE), stats::median(null_result$spearman_cor, na.rm = TRUE)),
  wilcox_statistic = c(wilcox_jaccard$statistic, wilcox_cor$statistic),
  wilcox_p = c(wilcox_jaccard$p.value, wilcox_cor$p.value),
  rank_biserial_effect = c(
    rank_biserial(wilcox_jaccard$statistic, sum(!is.na(same_species_result$jaccard)), sum(!is.na(null_result$jaccard))),
    rank_biserial(wilcox_cor$statistic, sum(!is.na(same_species_result$spearman_cor)), sum(!is.na(null_result$spearman_cor)))
  )
)
write_result(intragenomic_summary, "02b_intragenomic_variant_summary")

# Case study: GMQP-bins7_sp004366385, the label pointed at directly --
# recomputed here on the CORRECTED table, not the contaminated
# reconnaissance numbers from before Part 3a.
gmqp_case <- same_species_result |> dplyr::filter(species == "GMQP-bins7_sp004366385")
gmqp_n_asvs <- length(unique(c(gmqp_case$asv_a, gmqp_case$asv_b)))

p_intragenomic <- intragenomic_test |>
  tidyr::pivot_longer(c(jaccard, spearman_cor), names_to = "metric", values_to = "value") |>
  dplyr::mutate(metric = factor(metric, levels = c("jaccard", "spearman_cor"),
                                  labels = c("Jaccard co-occurrence", "Spearman abundance correlation"))) |>
  ggplot2::ggplot(ggplot2::aes(x = group, y = value, fill = group)) +
  ggplot2::geom_violin(alpha = 0.6, na.rm = TRUE) +
  ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA, na.rm = TRUE) +
  ggplot2::facet_wrap(~metric, scales = "free_y") +
  ggplot2::scale_x_discrete(labels = c(same_species_label = "Same species\nlabel", null_cross_species = "Null (cross-\nspecies)")) +
  ggplot2::labs(x = NULL, y = NULL,
                title = "Intragenomic-variant test: same-species-label ASV pairs vs. cross-species null",
                subtitle = sprintf(
                  "Jaccard p=%.2g, r=%.3f; Spearman-cor p=%.2g, r=%.3f (rank-biserial). Case study GMQP-bins7_sp004366385: %d ASVs (capped at %d).",
                  wilcox_jaccard$p.value, intragenomic_summary$rank_biserial_effect[1],
                  wilcox_cor$p.value, intragenomic_summary$rank_biserial_effect[2],
                  gmqp_n_asvs, max_group_size
                )) +
  ggplot2::theme(legend.position = "none")
save_plot(p_intragenomic, "02b_intragenomic_variant_test", w = 11, h = 5)

message(sprintf(
  "[02b_controls] intragenomic-variant test: Jaccard median same-species=%.3f vs null=%.3f (Wilcoxon p=%.3g, effect=%.3f); Spearman-cor median same-species=%.3f vs null=%.3f (Wilcoxon p=%.3g, effect=%.3f). GMQP-bins7_sp004366385 case: %d ASVs post-3a (was 180 pre-3a).",
  intragenomic_summary$median_same_species[1], intragenomic_summary$median_null[1], intragenomic_summary$wilcox_p[1], intragenomic_summary$rank_biserial_effect[1],
  intragenomic_summary$median_same_species[2], intragenomic_summary$median_null[2], intragenomic_summary$wilcox_p[2], intragenomic_summary$rank_biserial_effect[2],
  gmqp_n_asvs
))

# ============================================================================
# Read-tracking funnel: raw -> QC -> denoise -> taxonomy -> decontamination
# ============================================================================
# Samples on the y axis (with depth brackets, rotated 90 deg from
# 02_qc_filter.R's horizontal version -- 51 samples read better as a list
# than crammed along a shared x axis across 5 facets), faceted by pipeline
# step. Sources:
#   raw            = tracking$input_reads          (upstream cutadapt input)
#   qc             = tracking$filtered_reads        (DADA2 filterAndTrim)
#   denoise        = tracking$non_chimeric_reads     (DADA2 denoise + chimera removal --
#                     the actual final ASV table upstream hands off)
#   taxonomy       = rowSums(counts)                (this script's `counts` --
#                     counts_postlineage.rds, i.e. after both the upstream
#                     "best taxonomy" reconciliation AND 02_qc_filter.R's own
#                     on-target lineage filter; two taxonomy-related steps
#                     folded into one, since the user asked for 5 named
#                     stages, not 6)
#   decontamination = rowSums(counts_filtered)       (this script's final gate --
#                     blacklist + bleed floor + control-drop + prevalence
#                     filter; controls have NO row here at all, not a zero
#                     one -- they're genuinely absent from this stage, and a
#                     missing bar says that more honestly than a 0-length
#                     one would, which log(0) can't render anyway)
tracking <- readRDS(file.path(path_processed, "tracking.rds"))

read_tracking <- dplyr::bind_rows(
  tracking |> dplyr::transmute(sample_id, step = "raw", reads = input_reads),
  tracking |> dplyr::transmute(sample_id, step = "qc", reads = filtered_reads),
  tracking |> dplyr::transmute(sample_id, step = "denoise", reads = non_chimeric_reads),
  tibble::tibble(sample_id = rownames(counts), step = "taxonomy", reads = rowSums(counts)),
  tibble::tibble(sample_id = rownames(counts_filtered), step = "decontamination", reads = rowSums(counts_filtered))
) |>
  dplyr::mutate(step = factor(step, levels = c("raw", "qc", "denoise", "taxonomy", "decontamination"))) |>
  dplyr::left_join(metadata |> dplyr::select(sample_id, sample_type, depth_m), by = "sample_id")

# Samples ordered by sample_type then depth_m (this figure has no
# facet_grid(~sample_type) doing that grouping for free, since the facet
# dimension here is `step` -- ordered explicitly instead), shared across
# every facet and the bracket panel below.
sample_order_rt <- metadata |> dplyr::arrange(sample_type, depth_m) |> dplyr::pull(sample_id)
read_tracking$sample_id_f <- factor(read_tracking$sample_id, levels = rev(sample_order_rt))
metadata_rt <- metadata |> dplyr::mutate(sample_id_f = factor(sample_id, levels = rev(sample_order_rt)))

p_read_tracking <- ggplot2::ggplot(read_tracking, ggplot2::aes(y = sample_id_f, x = reads, fill = sample_type)) +
  ggplot2::geom_col() +
  ggplot2::facet_wrap(~step, nrow = 1) +
  ggplot2::scale_fill_manual(values = palette_sample_type()) +
  ggplot2::scale_x_log10(breaks = 10^(1:6), labels = scales::label_comma()) +
  ggplot2::labs(x = "Reads", y = NULL, fill = "Sample type",
                title = "Read tracking: raw → QC → denoise → taxonomy → decontamination") +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 7),
    axis.text.x = ggplot2::element_text(size = 9, angle = 90, hjust = 1, vjust = 0.5),
    axis.title = ggplot2::element_text(size = 14),
    plot.title = ggplot2::element_text(size = 18, face = "bold"),
    strip.text = ggplot2::element_text(size = 12, face = "bold"),
    legend.text = ggplot2::element_text(size = 11),
    legend.title = ggplot2::element_text(size = 12)
  )

# Depth brackets, rotated: one vertical bracket per group of samples sharing
# a depth_m, spanning that group's rows with end-ticks and a depth label,
# same construction as 02_qc_filter.R's horizontal version with x/y swapped.
depth_groups_rt <- metadata_rt |>
  dplyr::filter(!is.na(depth_m)) |>
  dplyr::group_by(sample_type, depth_m) |>
  dplyr::summarise(
    y_start = dplyr::first(sample_id_f), y_end = dplyr::last(sample_id_f),
    y_mid = sample_id_f[ceiling(dplyr::n() / 2)],
    label = sprintf("%g", dplyr::first(depth_m)),
    .groups = "drop"
  )

p_depth_brackets_v <- ggplot2::ggplot(metadata_rt, ggplot2::aes(y = sample_id_f, x = 1)) +
  ggplot2::geom_blank() +
  ggplot2::geom_segment(data = depth_groups_rt, ggplot2::aes(y = y_start, yend = y_end, x = 1, xend = 1), inherit.aes = FALSE) +
  ggplot2::geom_segment(data = depth_groups_rt, ggplot2::aes(y = y_start, yend = y_start, x = 1, xend = 0.55), inherit.aes = FALSE) +
  ggplot2::geom_segment(data = depth_groups_rt, ggplot2::aes(y = y_end, yend = y_end, x = 1, xend = 0.55), inherit.aes = FALSE) +
  ggplot2::geom_text(data = depth_groups_rt, ggplot2::aes(y = y_mid, x = 0.3, label = label), inherit.aes = FALSE, size = 2.6, hjust = 1) +
  ggplot2::scale_x_continuous(limits = c(-1.8, 1.3)) +
  ggplot2::labs(x = "Depth\n(m)") +
  ggplot2::theme_void() +
  ggplot2::theme(axis.title.x = ggplot2::element_text(size = 10))

# Retention-vs-raw panel: the absolute-count bars above are correct (every
# one of 51 samples genuinely loses reads at every step, confirmed directly
# -- min 5, max 23,258 reads dropped denoise->taxonomy alone) but a 5-20%
# step-to-step loss is close to invisible on a log axis spanning 5 orders of
# magnitude across samples 10,000x apart in depth. A linear 0-100%-of-raw
# line makes that same, already-correct attrition actually visible instead
# of looking like "these are all the same numbers."
retention_all <- read_tracking |>
  tidyr::pivot_wider(id_cols = sample_id, names_from = step, values_from = reads) |>
  dplyr::left_join(metadata |> dplyr::select(sample_id, sample_type), by = "sample_id") |>
  tidyr::pivot_longer(cols = c(raw, qc, denoise, taxonomy, decontamination), names_to = "step", values_to = "reads") |>
  dplyr::mutate(step = factor(step, levels = c("raw", "qc", "denoise", "taxonomy", "decontamination"))) |>
  dplyr::group_by(sample_id) |>
  dplyr::mutate(pct_of_raw = 100 * reads / reads[step == "raw"]) |>
  dplyr::ungroup()

p_retention <- ggplot2::ggplot(retention_all, ggplot2::aes(x = step, y = pct_of_raw, group = sample_id, color = sample_type)) +
  ggplot2::geom_line(alpha = 0.5) +
  ggplot2::geom_point(size = 1) +
  ggplot2::scale_color_manual(values = palette_sample_type()) +
  ggplot2::scale_y_continuous(labels = scales::label_percent(scale = 1), limits = c(0, 100)) +
  ggplot2::labs(x = NULL, y = "% of raw reads", color = "Sample type",
                title = "Retention per step, relative to raw (same data as above, linear scale)") +
  ggplot2::theme(
    axis.text = ggplot2::element_text(size = 11),
    axis.title = ggplot2::element_text(size = 13),
    plot.title = ggplot2::element_text(size = 15, face = "bold"),
    legend.position = "none" # already shown in the panel above
  )

p_read_tracking_full <- (p_depth_brackets_v | p_read_tracking) / p_retention +
  patchwork::plot_layout(heights = c(3, 1))
save_plot(p_read_tracking_full, "02b_read_tracking", w = 22, h = 14, dpi = 350)

# Joined by sample_id, not paired positionally -- the raw/qc/denoise blocks
# above are in tracking.rds's own row order, taxonomy in rownames(counts)'s,
# decontamination in rownames(counts_filtered)'s, three different orderings
# that a plain positional divide would silently mismatch.
retention_pct <- read_tracking |>
  dplyr::filter(step %in% c("raw", "decontamination")) |>
  tidyr::pivot_wider(id_cols = sample_id, names_from = step, values_from = reads) |>
  dplyr::filter(!is.na(decontamination)) |> # controls: no decontamination row at all
  dplyr::mutate(pct = 100 * decontamination / raw)

message(sprintf(
  "[02b_controls] read tracking: mean retention raw->decontamination (real samples only) = %.1f%%",
  mean(retention_pct$pct)
))

# Per-step summary table (raw / qc / denoise / taxonomy / decontamination),
# reads and retention-vs-raw, split by group (real / mock / blank) since
# controls drop out entirely at the decontamination step (no row, not a
# zero) and mixing them into one "all samples" row would be misleading.
# This is the numeric table version of p_read_tracking + p_retention above,
# not a new computation -- same read_tracking / retention_all objects.
read_tracking_summary <- retention_all |>
  dplyr::mutate(group = dplyr::case_when(
    sample_id %in% mock_sample_for ~ "mock",
    sample_id %in% neg_ids ~ "blank",
    TRUE ~ "real"
  )) |>
  dplyr::group_by(step, group) |>
  dplyr::summarise(
    n_samples = dplyr::n(),
    total_reads = sum(reads),
    mean_reads = mean(reads),
    median_reads = stats::median(reads),
    mean_pct_of_raw = mean(pct_of_raw),
    median_pct_of_raw = stats::median(pct_of_raw),
    .groups = "drop"
  ) |>
  dplyr::arrange(step, group)
write_result(read_tracking_summary, "02b_read_tracking_summary")

message("[02b_controls] read-tracking step summary written to results/02b_read_tracking_summary.tsv")
