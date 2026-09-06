# 02_qc_filter.R -- lineage filtering, library-size QC, decontam (candidate-
# flagging only). This is the FIRST HALF of QC; scripts/02b_controls.R is the
# second half and does the actual gating (apply the reviewed contaminant
# blacklist + the bleed-through floor, drop control libraries, run the
# prevalence filter, run replicate concordance, and save *_clean.rds).
#
# Why the split: decontam's negative controls AND 02b's bleed-floor/mock
# scoring both need the control libraries as *input*, but the final analysis
# matrix must not contain them (AGENTS.md; see also PLAN.md §7's replicate-
# concordance addition). A single top-to-bottom script can't do both -- it
# would either filter controls before decontam can use them, or leak them
# past decontam into the final save. So: 02_qc_filter.R stops right after
# decontam produces a *candidate* list (not yet applied), and 02b_controls.R
# picks up from there once the candidates have a chance to be reviewed.
#
# Filter order this enforces (previously backwards -- an ASV seen only in
# ctr_EB + ctr_MM and nowhere real used to pass the >=2-sample prevalence
# filter on the controls' strength alone):
#   decontam candidates + bleed floor (both need controls) -> apply both +
#   drop controls (02b) -> prevalence filter (02b) -> everything in 03-10.

source("scripts/00_setup.R")

counts <- readRDS(file.path(path_processed, "counts.rds"))
taxonomy <- readRDS(file.path(path_processed, "taxonomy.rds"))
metadata <- readRDS(file.path(path_processed, "metadata.rds"))

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
# revisit if a future run shows one. ctr_EB (41 reads) and ctr_MM (70 reads)
# are shallow enough that decontam's power against them is limited
# regardless of threshold -- see the threshold note in §3 and AGENTS.md.
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

# --- 3. decontam (prevalence method) -- candidates only, not applied -------
# ctr_EB (extraction blank) and ctr_MM (mastermix blank) are the negatives
# (AGENTS.md departure #3); the three mocks are a separate accuracy check
# (scripts/02b_controls.R), not decontam input.
#
# threshold = 0.5, not decontam's default 0.1: n=2 negatives has very
# little statistical power at any threshold (Davis et al. 2018, the decontam
# paper, are explicit that prevalence mode wants more than 2 negatives), and
# 0.1 is tuned for larger negative sets. 0.5 ("more prevalent in negatives
# than in samples") is the more defensible cutoff at n=2, but even so this
# is a CANDIDATE list, not an automated removal -- a human reviews it
# (results/02b_contaminant_candidates.tsv's `decision` column) and
# 02b_controls.R applies only what's marked "remove". No `batch` argument:
# no extraction-batch field exists anywhere in data/ to stratify by.
is_neg <- metadata$control_type[match(rownames(counts_kept), metadata$sample_id)] %in%
  c("extraction_blank", "mastermix_blank")
stopifnot(sum(is_neg) == 2L)

decontam_threshold <- 0.5
contam <- decontam::isContaminant(counts_kept, method = "prevalence", neg = is_neg,
                                   threshold = decontam_threshold)
is_candidate <- contam$contaminant & !is.na(contam$contaminant)

message(sprintf(
  "[02_qc_filter] decontam (prevalence, threshold=%.1f, neg = ctr_EB + ctr_MM): %d/%d ASVs flagged as contaminant CANDIDATES (review required -- see results/02b_contaminant_candidates.tsv)",
  decontam_threshold, sum(is_candidate), ncol(counts_kept)
))

# --- candidate list, with review-decision round-tripping --------------------
# First run: no prior file, every candidate's `decision` is NA -- 02b_controls.R
# treats NA the same as "keep" (nothing removed) until a human sets
# `decision = "remove"` for specific rows and reruns. This is deliberate:
# an empty/unset decision column on a fresh run means "not yet reviewed",
# not "reviewed and found clean" -- don't read silence as a clean bill of
# health, especially given how shallow ctr_EB/ctr_MM are.
candidates_path <- "results/02b_contaminant_candidates.tsv"
new_candidates <- tibble::tibble(asv_id = rownames(contam), contam) |>
  dplyr::filter(is_candidate) |>
  dplyr::left_join(taxonomy_kept |> dplyr::select(asv_id, domain, phylum, class, order, family, genus, species),
                    by = "asv_id")

if (file.exists(candidates_path)) {
  prior <- readr::read_tsv(candidates_path, show_col_types = FALSE)
  prior_decisions <- prior |> dplyr::select(asv_id, decision)
  new_candidates <- new_candidates |>
    dplyr::left_join(prior_decisions, by = "asv_id") |>
    dplyr::relocate(decision, .after = asv_id)
  n_carried <- sum(!is.na(new_candidates$decision))
  message(sprintf("[02_qc_filter] carried forward %d prior review decision(s) from %s", n_carried, candidates_path))
} else {
  new_candidates <- new_candidates |> dplyr::mutate(decision = NA_character_) |> dplyr::relocate(decision, .after = asv_id)
  message(sprintf("[02_qc_filter] no prior %s -- every candidate is unreviewed (decision = NA) this run", candidates_path))
}
write_result(new_candidates, "02b_contaminant_candidates")

# Standard decontam diagnostic (its own vignette's PA plot): prevalence of
# each ASV among the negative controls vs. among every other sample, colored
# by the candidate call. A real contaminant should sit high on the
# negative-control axis and/or low on the true-sample axis.
pa_neg <- colSums(counts_kept[is_neg, , drop = FALSE] > 0)
pa_pos <- colSums(counts_kept[!is_neg, , drop = FALSE] > 0)
decontam_pa <- tibble::tibble(
  asv_id = colnames(counts_kept),
  prevalence_controls = pa_neg,
  prevalence_samples = pa_pos,
  candidate = is_candidate
)
write_result(decontam_pa, "02_decontam_prevalence")

p_decontam <- ggplot2::ggplot(decontam_pa, ggplot2::aes(
  x = prevalence_controls, y = prevalence_samples, color = candidate
)) +
  ggplot2::geom_jitter(width = 0.08, height = 0.4, alpha = 0.4, size = 1) +
  ggplot2::scale_color_manual(values = c(`TRUE` = "#E15759", `FALSE` = "grey60")) +
  ggplot2::scale_x_continuous(breaks = 0:sum(is_neg)) +
  ggplot2::labs(x = "Prevalence (negative controls: ctr_EB, ctr_MM)",
                y = "Prevalence (all other samples)", color = "Contaminant\ncandidate",
                title = "Decontam diagnostic: ASV prevalence, controls vs. samples",
                subtitle = sprintf("%d/%d ASVs are candidates (threshold=%.1f) -- review, not yet removed",
                                    sum(is_candidate), ncol(counts_kept), decontam_threshold))
save_plot(p_decontam, "02_decontam_prevalence", w = 6, h = 5)

# --- handoff to 02b_controls.R ---------------------------------------------
# Deliberately NOT dropping controls, NOT prevalence-filtering, NOT saving
# *_clean.rds here -- see header. 02b_controls.R reads these and metadata.rds
# (which still has all 51 rows) directly.
saveRDS(counts_kept, file.path(path_processed, "counts_postlineage.rds"))
saveRDS(taxonomy_kept, file.path(path_processed, "taxonomy_postlineage.rds"))

retention <- tibble::tibble(
  stage = c("import", "lineage_filter"),
  n_asv = c(n_asv_start, ncol(counts_kept)),
  n_reads = c(n_reads_start, sum(counts_kept))
)
write_result(retention, "02_retention")

message(sprintf(
  "[02_qc_filter] handing off to 02b_controls.R: %d samples (incl. 5 controls) x %d ASVs",
  nrow(counts_kept), ncol(counts_kept)
))
