# 08_differential_abundance.R — ALDEx2 (primary) + Maaslin2 (cross-check),
# contrast = sample_type (sediment vs. water; controls already dropped
# upstream in 02b_controls.R). Report the intersection.

source("scripts/00_setup.R")

counts <- readRDS(file.path(path_processed, "counts_clean.rds"))
taxonomy <- readRDS(file.path(path_processed, "taxonomy_clean.rds"))
metadata <- readRDS(file.path(path_processed, "metadata_clean.rds"))

eco_samples <- metadata$sample_id
counts_eco <- counts[eco_samples, , drop = FALSE]
meta_eco <- metadata
conditions <- meta_eco$sample_type
stopifnot(dplyr::n_distinct(conditions) == 2L) # aldex.ttest/effect assume a 2-group contrast

message(sprintf("[08_differential_abundance] contrast = sample_type (%s), N = %d (%s)",
                 paste(unique(conditions), collapse = " vs. "), length(eco_samples),
                 paste(table(conditions), collapse = "/")))

# --- ALDEx2 (primary) -------------------------------------------------------
# ALDEx2 wants features as ROWS, samples as COLUMNS -- the opposite of this
# project's canonical orientation (PLAN.md §3), so transpose here only.
#
# No library_size covariate here, deliberately, not an oversight: ALDEx2's
# CLR transform (aldex.clr(), run internally by aldex()) divides each
# sample's counts by that SAME sample's own geometric mean before taking
# logs, so a sample's total library size cancels out of the transform by
# construction -- it's compositionally invariant to library size already.
# Bolting a library_size covariate onto the t-test framework aldex() uses
# would need switching to aldex.glm() with a full design matrix, a
# materially bigger, riskier rewrite of the primary DA method for a
# confound this method doesn't actually have. (Maaslin2 below isn't CLR-based
# and does get the covariate, for exactly that reason.)
reads_for_aldex <- t(counts_eco)

aldex_fit <- ALDEx2::aldex(reads_for_aldex, conditions, mc.samples = 128,
                            test = "t", effect = TRUE, denom = "all")

aldex_result <- tibble::as_tibble(aldex_fit, rownames = "asv_id") |>
  dplyr::left_join(taxonomy |> dplyr::select(asv_id, phylum, genus, species), by = "asv_id") |>
  dplyr::arrange(wi.eBH)

write_result(aldex_result, "da_aldex2_sample_type")

# --- Maaslin2 (cross-check) --------------------------------------------
# Maaslin2 wants samples as rows, features as columns (this project's
# canonical orientation) -- no transpose needed; it does its own
# normalization (default TSS + LOG) from raw counts, which is NOT
# library-size-invariant the way ALDEx2's CLR is above -- so log_library_size
# is added as its own fixed effect here, standard practice for this tool,
# so a sample_type effect can't just be tracking library size instead.
maaslin_input_data <- as.data.frame(counts_eco)
maaslin_input_metadata <- as.data.frame(meta_eco[, c("sample_id", "sample_type", "log_library_size")])
rownames(maaslin_input_metadata) <- maaslin_input_metadata$sample_id

maaslin_fit <- Maaslin2::Maaslin2(
  input_data = maaslin_input_data,
  input_metadata = maaslin_input_metadata,
  output = "results/maaslin2_sample_type",
  fixed_effects = c("sample_type", "log_library_size"),
  min_prevalence = 0, # already prevalence-filtered in 02_qc_filter.R
  standardize = FALSE,
  plot_heatmap = FALSE,
  plot_scatter = FALSE
)

maaslin_all_terms <- tibble::as_tibble(maaslin_fit$results)
write_result(maaslin_all_terms, "da_maaslin2_all_terms") # both fixed effects, unfiltered -- for inspecting the log_library_size term itself

# fixed_effects now has two terms (sample_type, log_library_size), so
# maaslin_fit$results has a row per (feature, term) pair -- filtered to the
# sample_type contrast specifically here, or sig_maaslin below would silently
# mix in ASVs significant for log_library_size instead of/as well as
# sample_type.
maaslin_result <- maaslin_all_terms |>
  dplyr::filter(metadata == "sample_type") |>
  dplyr::rename(asv_id = feature) |>
  dplyr::left_join(taxonomy |> dplyr::select(asv_id, phylum, genus, species), by = "asv_id") |>
  dplyr::arrange(qval)

n_sig_libsize <- sum((maaslin_all_terms |> dplyr::filter(metadata == "log_library_size"))$qval < 0.05, na.rm = TRUE)
message(sprintf(
  "[08_differential_abundance] Maaslin2: %d/%d ASVs significantly associated with log_library_size itself (q<0.05) -- see results/da_maaslin2_all_terms.tsv",
  n_sig_libsize, nrow(maaslin_all_terms |> dplyr::filter(metadata == "log_library_size"))
))

write_result(maaslin_result, "da_maaslin2_sample_type")

# --- intersection ----------------------------------------------------------
sig_aldex <- aldex_result$asv_id[aldex_result$wi.eBH < 0.05]
sig_maaslin <- maaslin_result$asv_id[maaslin_result$qval < 0.05]
sig_both <- intersect(sig_aldex, sig_maaslin)

intersection_tbl <- taxonomy |>
  dplyr::filter(asv_id %in% sig_both) |>
  dplyr::select(asv_id, phylum, genus, species) |>
  dplyr::left_join(aldex_result |> dplyr::select(asv_id, aldex_effect = effect, aldex_qval = wi.eBH), by = "asv_id") |>
  dplyr::left_join(maaslin_result |> dplyr::select(asv_id, maaslin_coef = coef, maaslin_qval = qval), by = "asv_id") |>
  dplyr::arrange(aldex_qval)

write_result(intersection_tbl, "da_intersection_sample_type")

message(sprintf(
  "[08_differential_abundance] significant (q<0.05): ALDEx2 = %d, Maaslin2 = %d, intersection = %d",
  length(sig_aldex), length(sig_maaslin), length(sig_both)
))

# --- plots -----------------------------------------------------------------
p_volcano <- aldex_result |>
  dplyr::mutate(significant = wi.eBH < 0.05) |>
  ggplot2::ggplot(ggplot2::aes(x = effect, y = -log10(wi.eBH), color = significant)) +
  ggplot2::geom_point(alpha = 0.5, size = 1) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ggplot2::scale_color_manual(values = c(`TRUE` = "#E15759", `FALSE` = "grey70")) +
  ggplot2::labs(x = "ALDEx2 effect size", y = "-log10(BH-adjusted p, Wilcoxon)",
                title = "Differential abundance: sediment vs. water (ALDEx2)", color = "q < 0.05")
save_plot(p_volcano, "08_da_volcano", w = 6, h = 5)

if (nrow(intersection_tbl) > 0) {
  p_effects <- intersection_tbl |>
    dplyr::mutate(label = ifelse(!is.na(genus), genus, asv_id)) |>
    ggplot2::ggplot(ggplot2::aes(x = forcats::fct_reorder(label, aldex_effect), y = aldex_effect, fill = phylum)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = stats::setNames(palette_taxa(dplyr::n_distinct(intersection_tbl$phylum)),
                                                          unique(intersection_tbl$phylum))) +
    ggplot2::labs(x = NULL, y = "ALDEx2 effect size (+ = higher in water, given factor order)",
                  title = "Taxa significant in both ALDEx2 and Maaslin2 (q < 0.05)")
  save_plot(p_effects, "08_da_intersection_effects", w = 7, h = max(4, 0.25 * nrow(intersection_tbl)))
}
