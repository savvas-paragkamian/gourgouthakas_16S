# 09_spatial_analysis.R — reworked as a depth transect (AGENTS.md departure
# #2): every sample shares one lat/lon (single cave), so there is no map, no
# geographic distance-decay, and no geographic Moran's I. `sf`/`terra` are
# deliberately not installed. The gradient here is depth_m; this script
# does the depth-space equivalents PLAN.md's generic §09 asks for.

source("scripts/00_setup.R")

metadata <- readRDS(file.path(path_processed, "metadata_clean.rds"))
dist_bray <- readRDS(file.path(path_processed, "dist_bray.rds"))
alpha <- readr::read_tsv("results/alpha_diversity.tsv", show_col_types = FALSE)
ord_coords <- readr::read_tsv("results/06_ordination_coords.tsv", show_col_types = FALSE) |>
  dplyr::filter(distance == "bray")

eco_samples <- metadata$sample_id[metadata$sample_type != "control" & !is.na(metadata$depth_m)]

# --- distance-decay in depth space -----------------------------------------
# The figure PLAN.md §09 asks for, in depth space rather than geographic
# space; 07_environmental_drivers.R's Mantel test is the numeric companion
# to this plot (same relationship, tested rather than just visualized).
bray_mat <- as.matrix(dist_bray)[eco_samples, eco_samples]
depth_v <- metadata$depth_m[match(eco_samples, metadata$sample_id)]
depth_dist_mat <- as.matrix(stats::dist(depth_v))
dimnames(depth_dist_mat) <- dimnames(bray_mat)

pairs_upper <- which(upper.tri(bray_mat), arr.ind = TRUE)
decay_df <- tibble::tibble(
  sample_a = eco_samples[pairs_upper[, 1]],
  sample_b = eco_samples[pairs_upper[, 2]],
  depth_distance = depth_dist_mat[pairs_upper],
  bray_dissimilarity = bray_mat[pairs_upper]
) |>
  dplyr::mutate(
    type_a = metadata$sample_type[match(sample_a, metadata$sample_id)],
    type_b = metadata$sample_type[match(sample_b, metadata$sample_id)],
    pair_type = ifelse(type_a == type_b, type_a, "sediment-water")
  )

write_result(decay_df, "09_distance_decay")

p_decay <- ggplot2::ggplot(decay_df, ggplot2::aes(x = depth_distance, y = bray_dissimilarity)) +
  ggplot2::geom_point(ggplot2::aes(color = pair_type), alpha = 0.4, size = 1) +
  ggplot2::geom_smooth(method = "loess", formula = y ~ x, color = "black", se = TRUE) +
  ggplot2::labs(x = "Depth distance (m)", y = "Bray-Curtis dissimilarity",
                color = "Pair type",
                title = "Distance-decay in depth space",
                subtitle = "Community turnover vs. depth separation (Mantel test: see results/07_mantel.tsv)")
save_plot(p_decay, "09_distance_decay", w = 7, h = 5)

# --- Moran's I, substituting depth distance for geographic distance --------
# Classic Moran's I needs a spatial weights matrix; there is no geography
# here (one lat/lon for the whole cave), so depth is used as the 1D
# coordinate instead -- inverse-depth-distance weights, diagonal zero,
# +1 in the denominator so exact depth ties (e.g. C1 vs C1_I, both 0 m)
# get a large but finite weight instead of Inf.
depth_weights <- 1 / (depth_dist_mat + 1)
diag(depth_weights) <- 0

moran_targets <- list(
  observed = alpha$observed[match(eco_samples, alpha$sample_id)],
  shannon = alpha$shannon[match(eco_samples, alpha$sample_id)],
  pcoa1 = ord_coords$PCoA1[match(eco_samples, ord_coords$sample_id)]
)

moran_results <- purrr::imap_dfr(moran_targets, function(x, name) {
  keep <- !is.na(x)
  if (sum(keep) < 3) {
    return(tibble::tibble(variable = name, observed = NA, expected = NA, sd = NA, p_value = NA,
                            note = "too few non-NA samples"))
  }
  fit <- ape::Moran.I(x[keep], depth_weights[keep, keep])
  tibble::tibble(variable = name, observed = fit$observed, expected = fit$expected,
                  sd = fit$sd, p_value = fit$p.value, note = NA_character_)
})
write_result(moran_results, "09_morans_i_depth")

message(sprintf(
  "[09_spatial_analysis] Moran's I (depth-weighted): %s",
  paste(sprintf("%s: I=%.3f, p=%.3f", moran_results$variable, moran_results$observed, moran_results$p_value), collapse = "; ")
))

message(sprintf("[09_spatial_analysis] distance-decay pairs: %d", nrow(decay_df)))
