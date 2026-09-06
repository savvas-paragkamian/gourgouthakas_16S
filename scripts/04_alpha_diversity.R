# 04_alpha_diversity.R — Observed richness, Shannon, Simpson, Pielou evenness.
#
# No Faith's PD here: this dataset has no phylogeny (AGENTS.md departure #1),
# so `picante` is deliberately not installed and PD is out of scope.

source("scripts/00_setup.R")

counts_rarefied <- readRDS(file.path(path_processed, "counts_rarefied.rds"))
metadata <- readRDS(file.path(path_processed, "metadata_clean.rds"))

# --- metrics ---------------------------------------------------------------
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
  "alpha_diversity"
)

# --- tests -------------------------------------------------------------
# Real contrasts per AGENTS.md departure #5: sample_type (sediment vs water)
# and continuous depth_m. Controls are excluded from ecological testing (QC
# artifacts, not part of the environmental gradient) but were reported above.
alpha_eco <- alpha |> dplyr::filter(sample_type != "control")

metrics <- c("observed", "shannon", "simpson", "pielou")

kw_tests <- purrr::map_dfr(metrics, function(m) {
  ft <- stats::kruskal.test(alpha_eco[[m]] ~ alpha_eco$sample_type)
  tibble::tibble(metric = m, test = "kruskal.test ~ sample_type",
                  statistic = unname(ft$statistic), df = unname(ft$parameter), p_value = ft$p.value)
})

depth_tests <- purrr::map_dfr(metrics, function(m) {
  fit <- stats::lm(alpha_eco[[m]] ~ depth_m * sample_type, data = alpha_eco)
  s <- summary(fit)
  tibble::tibble(
    metric = m, test = "lm(metric ~ depth_m * sample_type)",
    term = rownames(s$coefficients),
    estimate = s$coefficients[, "Estimate"],
    p_value = s$coefficients[, "Pr(>|t|)"],
    r_squared = s$r.squared
  )
})

write_result(kw_tests, "alpha_stats_sample_type")
write_result(depth_tests, "alpha_stats_depth")

# --- plots -----------------------------------------------------------------
p_box <- alpha_eco |>
  tidyr::pivot_longer(dplyr::all_of(metrics), names_to = "metric", values_to = "value") |>
  ggplot2::ggplot(ggplot2::aes(x = sample_type, y = value, fill = sample_type)) +
  ggplot2::geom_boxplot(outlier.shape = NA) +
  ggplot2::geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
  ggplot2::facet_wrap(~metric, scales = "free_y") +
  ggplot2::scale_fill_manual(values = palette_sample_type()) +
  ggplot2::labs(x = NULL, y = NULL, title = "Alpha diversity by sample type") +
  ggplot2::theme(legend.position = "none")
save_plot(p_box, "04_alpha_boxplot", w = 7, h = 6)

p_depth <- alpha_eco |>
  tidyr::pivot_longer(dplyr::all_of(metrics), names_to = "metric", values_to = "value") |>
  ggplot2::ggplot(ggplot2::aes(x = depth_m, y = value, color = sample_type)) +
  ggplot2::geom_point(size = 1.5, alpha = 0.7) +
  ggplot2::geom_smooth(method = "lm", se = TRUE, formula = y ~ x) +
  ggplot2::facet_wrap(~metric, scales = "free_y") +
  ggplot2::scale_color_manual(values = palette_sample_type()) +
  ggplot2::labs(x = "Depth (m)", y = NULL, color = "Sample type",
                title = "Alpha diversity vs. depth")
save_plot(p_depth, "04_alpha_vs_depth", w = 8, h = 6)

message(sprintf("[04_alpha_diversity] %d samples (rarefied @ %d), %d excluded as too shallow",
                 nrow(alpha), unique(rowSums(counts_rarefied))[1],
                 length(readRDS(file.path(path_processed, "rarefaction_excluded_samples.rds")))))
