# 04_alpha_diversity.R — Observed richness, Shannon, Simpson, Pielou evenness.
#
# No Faith's PD here: this dataset has no phylogeny (AGENTS.md departure #1),
# so `picante` is deliberately not installed and PD is out of scope.
#
# No library_size covariate here, deliberately (see AGENTS.md's "library
# size as a covariate" note): every sample in the rarefied matrix was
# equalized to exactly the same 499-read depth, so library_size has zero
# variance in this set by construction -- rarefaction already *is* the
# library-size control for these metrics, not a gap to fill.
#
# GENUS-LEVEL IS PRIMARY (Part 3c of golden-napping-breeze.md, evidence from
# Part 3b): metrics/tests/plots below run on `counts_genus_rarefied` first
# and write the primary result files. ASV-level (the original, un-collapsed
# computation) is kept and run second, writing to `*_asv` names -- kept as
# supplementary, not dropped, since the specific-lineage discussion still
# wants ASV resolution per AGENTS.md's genus-collapse scope note.

source("scripts/00_setup.R")

metadata <- readRDS(file.path(path_processed, "metadata_clean.rds"))
metrics <- c("observed", "shannon", "simpson", "pielou")

# One function, run twice (genus primary, ASV supplementary) rather than
# duplicating the metric/test/plot logic -- `suffix`/`label` route output
# filenames and plot titles, everything else is identical.
run_alpha_diversity <- function(counts_rarefied, excluded_samples_file, suffix, label) {
  observed <- vegan::specnumber(counts_rarefied)
  shannon <- vegan::diversity(counts_rarefied, index = "shannon")
  simpson <- vegan::diversity(counts_rarefied, index = "simpson")
  pielou <- shannon / log(observed) # 0 richness would be div/0; none here post-rarefaction

  alpha <- tibble::tibble(
    sample_id = rownames(counts_rarefied),
    observed = observed,
    shannon = shannon,
    simpson = simpson,
    pielou = pielou
  ) |>
    dplyr::left_join(metadata, by = "sample_id")

  write_result(
    alpha |> dplyr::select(sample_id, observed, shannon, simpson, pielou,
                            sample_type, sample_set, site, location, tech_rep, depth_m),
    paste0("alpha_diversity", suffix)
  )

  # Real contrasts per AGENTS.md departure #5: sample_type (sediment vs
  # water) and continuous depth_m. Controls never reach this script --
  # 02b_controls.R drops them before metadata_clean/the rarefied matrices
  # exist -- so `alpha` here is already the ecological subset, no filter
  # needed.
  kw_tests <- purrr::map_dfr(metrics, function(m) {
    ft <- stats::kruskal.test(alpha[[m]] ~ alpha$sample_type)
    tibble::tibble(metric = m, test = "kruskal.test ~ sample_type",
                    statistic = unname(ft$statistic), df = unname(ft$parameter), p_value = ft$p.value)
  })

  depth_tests <- purrr::map_dfr(metrics, function(m) {
    fit <- stats::lm(alpha[[m]] ~ depth_m * sample_type, data = alpha)
    s <- summary(fit)
    tibble::tibble(
      metric = m, test = "lm(metric ~ depth_m * sample_type)",
      term = rownames(s$coefficients),
      estimate = s$coefficients[, "Estimate"],
      p_value = s$coefficients[, "Pr(>|t|)"],
      r_squared = s$r.squared
    )
  })

  write_result(kw_tests, paste0("alpha_stats_sample_type", suffix))
  write_result(depth_tests, paste0("alpha_stats_depth", suffix))

  p_box <- alpha |>
    tidyr::pivot_longer(dplyr::all_of(metrics), names_to = "metric", values_to = "value") |>
    ggplot2::ggplot(ggplot2::aes(x = sample_type, y = value, fill = sample_type)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
    ggplot2::facet_wrap(~metric, scales = "free_y") +
    ggplot2::scale_fill_manual(values = palette_sample_type()) +
    ggplot2::labs(x = NULL, y = NULL, title = sprintf("Alpha diversity by sample type (%s)", label)) +
    ggplot2::theme(legend.position = "none")
  save_plot(p_box, paste0("04_alpha_boxplot", suffix), w = 7, h = 6)

  p_depth <- alpha |>
    tidyr::pivot_longer(dplyr::all_of(metrics), names_to = "metric", values_to = "value") |>
    ggplot2::ggplot(ggplot2::aes(x = depth_m, y = value, color = sample_type)) +
    ggplot2::geom_point(size = 1.5, alpha = 0.7) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, formula = y ~ x) +
    ggplot2::facet_wrap(~metric, scales = "free_y") +
    ggplot2::scale_color_manual(values = palette_sample_type()) +
    ggplot2::labs(x = "Depth (m)", y = NULL, color = "Sample type",
                  title = sprintf("Alpha diversity vs. depth (%s)", label))
  save_plot(p_depth, paste0("04_alpha_vs_depth", suffix), w = 8, h = 6)

  message(sprintf("[04_alpha_diversity] %s: %d samples (rarefied @ %d), %d excluded as too shallow",
                   label, nrow(alpha), unique(rowSums(counts_rarefied))[1],
                   length(readRDS(file.path(path_processed, excluded_samples_file)))))

  invisible(alpha)
}

# --- primary: genus-level ---------------------------------------------------
counts_genus_rarefied <- readRDS(file.path(path_processed, "counts_genus_rarefied.rds"))
run_alpha_diversity(counts_genus_rarefied, "rarefaction_excluded_samples_genus.rds", suffix = "", label = "genus-level, primary")

# --- supplementary: ASV-level (original, un-collapsed) ---------------------
counts_rarefied <- readRDS(file.path(path_processed, "counts_rarefied.rds"))
run_alpha_diversity(counts_rarefied, "rarefaction_excluded_samples.rds", suffix = "_asv", label = "ASV-level, supplementary")
