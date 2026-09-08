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
  dplyr::arrange(depth_m) # controls (no depth_m) sort last, dplyr's default NA handling

# Slide-ready version: 16:9, higher dpi, larger text throughout -- x axis
# ordered by depth_m (controls, with no depth concept, sort last) but
# labeled with sample_id, faceted into control/sediment/water panels
# (space="free_x" so each panel width matches its own sample count -- 5
# controls shouldn't get the same width as 36 sediment samples), y axis in
# plain read counts (not "1e+05") at every power of ten for a fast
# order-of-magnitude read. sample_id_f is a single shared factor (levels
# fixed once, in depth order) so the bracket panel below lines up exactly.
lib_sizes$sample_id_f <- forcats::fct_inorder(lib_sizes$sample_id)

p_libsize <- ggplot2::ggplot(lib_sizes, ggplot2::aes(x = sample_id_f, y = library_size, fill = sample_type)) +
  ggplot2::geom_col() +
  ggplot2::facet_grid(~sample_type, scales = "free_x", space = "free_x") +
  ggplot2::scale_fill_manual(values = palette_sample_type()) +
  ggplot2::scale_y_log10(breaks = 10^(1:6), labels = scales::label_comma()) +
  ggplot2::labs(x = NULL, y = "Library size (reads)", fill = "Sample type",
                title = "Library sizes after lineage filtering") +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.text.y = ggplot2::element_text(size = 14),
    axis.title = ggplot2::element_text(size = 18),
    plot.title = ggplot2::element_text(size = 22),
    strip.text = ggplot2::element_text(size = 16, face = "bold"),
    legend.text = ggplot2::element_text(size = 14),
    legend.title = ggplot2::element_text(size = 16)
  )

# Depth brackets below the x axis: one bracket per group of samples sharing
# the same depth_m (typically a site's two tech-rep extractions), spanning
# that group's bars with end-ticks and a centered depth label. Controls have
# no depth_m and get no bracket. Built from the full sample set via
# geom_blank() first (invisible, one row per sample) so every panel gets the
# same x domain -- and so the same per-facet width -- as p_libsize above;
# the visible brackets are a second, depth-only layer on top.
depth_groups <- lib_sizes |>
  dplyr::filter(!is.na(depth_m)) |>
  dplyr::group_by(sample_type, depth_m) |>
  dplyr::summarise(
    x_start = dplyr::first(sample_id_f),
    x_end = dplyr::last(sample_id_f),
    x_mid = sample_id_f[ceiling(dplyr::n() / 2)],
    label = sprintf("%g m", dplyr::first(depth_m)),
    .groups = "drop"
  )

p_depth_brackets <- ggplot2::ggplot(lib_sizes, ggplot2::aes(x = sample_id_f, y = 1)) +
  ggplot2::geom_blank() +
  ggplot2::geom_segment(data = depth_groups, ggplot2::aes(x = x_start, xend = x_end, y = 1, yend = 1), inherit.aes = FALSE) +
  ggplot2::geom_segment(data = depth_groups, ggplot2::aes(x = x_start, xend = x_start, y = 1, yend = 0.55), inherit.aes = FALSE) +
  ggplot2::geom_segment(data = depth_groups, ggplot2::aes(x = x_end, xend = x_end, y = 1, yend = 0.55), inherit.aes = FALSE) +
  ggplot2::geom_text(data = depth_groups, ggplot2::aes(x = x_mid, y = 0.45, label = label),
                      inherit.aes = FALSE, size = 3.2, angle = 90, hjust = 1) +
  ggplot2::facet_grid(~sample_type, scales = "free_x", space = "free_x") +
  ggplot2::scale_y_continuous(limits = c(-1.3, 1.3)) +
  ggplot2::theme_void() +
  ggplot2::theme(strip.text = ggplot2::element_blank())

p_libsize_full <- p_libsize / p_depth_brackets + patchwork::plot_layout(heights = c(6, 1.8))
save_plot(p_libsize_full, "02_library_sizes", w = 13.333, h = 9.5, dpi = 400)

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
