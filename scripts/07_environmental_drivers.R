# 07_environmental_drivers.R — environmental matrix + collinearity, db-RDA,
# envfit, variance partitioning, Mantel. Analytical core (PLAN.md §7).
#
# Adapted per AGENTS.md departure #2: there is no spatial/climate layer (one
# cave, one lat/lon) -- the gradient is depth_m. Variance partitioning and
# Mantel below are reframed around depth_m vs. the other chemistry variables,
# not "chemistry vs. spatial vs. climate" as PLAN.md's generic template has it.

source("scripts/00_setup.R")

metadata <- readRDS(file.path(path_processed, "metadata_clean.rds"))
dist_bray <- readRDS(file.path(path_processed, "dist_bray.rds"))
dist_aitchison <- readRDS(file.path(path_processed, "dist_aitchison.rds"))

# --- environmental matrix + collinearity screen ---------------------------
# elevation_m is depth_m's *exact* complement here (depth_m + elevation_m =
# 1535 constant, cor = -1, confirmed numerically) -- including both would be
# perfectly collinear, so elevation_m is dropped, not just deprioritized.
env_candidates <- c("depth_m", "temperature_c", "conductivity_ms")

# Controls already dropped in 02b_controls.R -- metadata is sediment+water only.
# log_library_size travels alongside the screened env_candidates (not one of
# them -- it's a technical covariate, not an environmental driver being
# tested for collinearity with the others) for the library-size-adjusted
# db-RDA/Mantel checks below.
env_raw <- metadata |>
  dplyr::select(sample_id, dplyr::all_of(env_candidates), log_library_size)

cor_mat <- stats::cor(env_raw[, env_candidates], use = "pairwise.complete.obs")
write_result(tibble::as_tibble(cor_mat, rownames = "variable"), "07_env_correlation")

p_cor <- tibble::as_tibble(cor_mat, rownames = "var1") |>
  tidyr::pivot_longer(-var1, names_to = "var2", values_to = "r") |>
  ggplot2::ggplot(ggplot2::aes(var1, var2, fill = r)) +
  ggplot2::geom_tile() +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", r)), size = 3) +
  ggplot2::scale_fill_gradient2(low = "#2E7DA6", high = "#B0794A", mid = "white", midpoint = 0, limits = c(-1, 1)) +
  ggplot2::labs(x = NULL, y = NULL, title = "Environmental variable correlation (elevation_m excluded: cor = -1 with depth_m)")
save_plot(p_cor, "07_env_correlation", w = 5, h = 4.5)

# VIF needs complete cases and >1 predictor with variance; conductivity_ms is
# missing for most sediment (AGENTS.md), so this is necessarily a reduced-N
# check, reported as such rather than silently using casewise deletion.
env_complete <- env_raw[stats::complete.cases(env_raw), ]
vif_vals <- car::vif(stats::lm(depth_m ~ temperature_c + conductivity_ms, data = env_complete))
write_result(tibble::tibble(variable = names(vif_vals), vif = vif_vals), "07_vif")
message(sprintf("[07_environmental_drivers] VIF computed on %d/%d complete-case samples", nrow(env_complete), nrow(env_raw)))

# --- db-RDA (Bray) ---------------------------------------------------------
# conductivity_ms's 19.6% completeness (AGENTS.md) would gut casewise-complete
# N to 10/46 if it were a primary predictor -- checked, and confirmed that
# collapses db-RDA to a degenerate null model (ordiR2step correctly finds no
# significant term to add at N=10, RsquareAdj() then returns numeric(0)).
# depth_m + temperature_c alone keep N=34/46, a far better-powered primary
# model; conductivity_ms is still screened above (VIF, correlation) and gets
# a clearly-labeled secondary, reduced-N sensitivity check below it.
run_dbrda <- function(env_vars, label, condition_var = NULL) {
  all_vars <- c(env_vars, condition_var)
  env_df <- env_raw[stats::complete.cases(env_raw[, all_vars]), ]
  samples_env <- env_df$sample_id
  bray_env <- stats::as.dist(as.matrix(dist_bray)[samples_env, samples_env])
  env_scaled <- as.data.frame(scale(as.data.frame(env_df[, all_vars])))
  rownames(env_scaled) <- samples_env
  names(env_scaled) <- all_vars

  # condition_var (e.g. log_library_size), when given, is partialled out via
  # vegan's Condition() -- the standard way to add a covariate to a
  # constrained ordination -- so the reported R2/terms below are for
  # env_vars *after* adjusting for it, not instead of it.
  full_formula <- if (is.null(condition_var)) {
    stats::reformulate(env_vars)
  } else {
    stats::reformulate(c(env_vars, sprintf("Condition(%s)", condition_var)))
  }
  dbrda_full <- vegan::capscale(stats::update(full_formula, bray_env ~ .), data = env_scaled)
  dbrda_null_formula <- if (is.null(condition_var)) bray_env ~ 1 else stats::reformulate(sprintf("Condition(%s)", condition_var), response = "bray_env")
  dbrda_null <- vegan::capscale(dbrda_null_formula, data = env_scaled)
  dbrda_step <- tryCatch(
    vegan::ordiR2step(dbrda_null, scope = stats::formula(dbrda_full), direction = "forward", permutations = 999, trace = FALSE),
    error = function(e) {
      message(sprintf("[07_environmental_drivers] (%s) ordiR2step failed: %s", label, conditionMessage(e)))
      dbrda_full
    }
  )
  r2 <- vegan::RsquareAdj(dbrda_step)
  r2_val <- if (length(r2$r.squared) == 0) NA_real_ else r2$r.squared
  adj_r2_val <- if (length(r2$adj.r.squared) == 0) NA_real_ else r2$adj.r.squared

  terms_tbl <- if (length(dbrda_step$CCA$eig) > 0) {
    tibble::as_tibble(as.data.frame(vegan::anova.cca(dbrda_step, by = "margin", permutations = 999)), rownames = "term")
  } else {
    tibble::tibble(term = "(no term retained by forward selection)")
  }

  message(sprintf("[07_environmental_drivers] db-RDA (%s, N=%d): R2 = %s, adj.R2 = %s",
                   label, length(samples_env),
                   ifelse(is.na(r2_val), "NA (no constrained axes retained)", sprintf("%.3f", r2_val)),
                   ifelse(is.na(adj_r2_val), "NA", sprintf("%.3f", adj_r2_val))))

  list(step = dbrda_step, terms = terms_tbl, samples = samples_env, env_scaled = env_scaled, bray_env = bray_env,
       r2 = r2_val, adj_r2 = adj_r2_val, n = length(samples_env))
}

dbrda_primary <- run_dbrda(c("depth_m", "temperature_c"), "primary: depth_m + temperature_c")
dbrda_sensitivity <- run_dbrda(c("depth_m", "temperature_c", "conductivity_ms"),
                                "sensitivity: + conductivity_ms, reduced N")
# Library size adjusted: same depth_m + temperature_c terms as primary, but
# with log_library_size partialled out via Condition() first -- does the
# depth_m/temperature_c signal survive once library-size differences between
# samples are no longer free to explain any of the variance?
dbrda_libsize <- run_dbrda(c("depth_m", "temperature_c"), "library-size-adjusted: + Condition(log_library_size)",
                            condition_var = "log_library_size")

write_result(
  dplyr::bind_rows(
    dbrda_primary$terms |> dplyr::mutate(model = "primary"),
    dbrda_libsize$terms |> dplyr::mutate(model = "libsize_adjusted")
  ),
  "07_dbrda_terms"
)
write_result(
  tibble::tibble(model = c("primary", "sensitivity", "libsize_adjusted"),
                  n = c(dbrda_primary$n, dbrda_sensitivity$n, dbrda_libsize$n),
                  r_squared = c(dbrda_primary$r2, dbrda_sensitivity$r2, dbrda_libsize$r2),
                  adj_r_squared = c(dbrda_primary$adj_r2, dbrda_sensitivity$adj_r2, dbrda_libsize$adj_r2)),
  "07_dbrda_summary"
)

# Everything downstream (envfit, varpart, Mantel) uses the better-powered
# primary model's sample set and distance subset.
samples_env <- dbrda_primary$samples
env_df_scaled <- dbrda_primary$env_scaled
bray_env <- dbrda_primary$bray_env
env_complete <- env_raw[match(samples_env, env_raw$sample_id), ]

# --- envfit -----------------------------------------------------------
ord_coords <- readr::read_tsv("results/06_ordination_coords.tsv", show_col_types = FALSE) |>
  dplyr::filter(distance == "bray")
pcoa_mat <- as.matrix(ord_coords[match(samples_env, ord_coords$sample_id), c("PCoA1", "PCoA2")])
rownames(pcoa_mat) <- samples_env

envfit_fit <- vegan::envfit(pcoa_mat, env_df_scaled, permutations = 999, na.rm = TRUE)
write_result(tidy_envfit(envfit_fit), "07_envfit")

p_triplot <- ggplot2::ggplot(ord_coords |> dplyr::filter(sample_id %in% samples_env),
                              ggplot2::aes(PCoA1, PCoA2, color = sample_type)) +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_color_manual(values = palette_sample_type()) +
  ggplot2::geom_segment(
    data = tidy_envfit(envfit_fit),
    ggplot2::aes(x = 0, y = 0, xend = axis1 * sqrt(r2), yend = axis2 * sqrt(r2)),
    inherit.aes = FALSE, arrow = grid::arrow(length = grid::unit(0.2, "cm")), color = "black"
  ) +
  ggrepel::geom_text_repel(
    data = tidy_envfit(envfit_fit),
    ggplot2::aes(x = axis1 * sqrt(r2), y = axis2 * sqrt(r2), label = term),
    inherit.aes = FALSE, size = 3
  ) +
  ggplot2::labs(title = "PCoA (Bray) with envfit vectors", subtitle = sprintf("N = %d complete-case samples", length(samples_env)))
save_plot(p_triplot, "07_envfit_triplot", w = 6, h = 5)

# --- variance partitioning --------------------------------------------------
# Reframed per AGENTS.md departure #2: no spatial/climate blocks exist here.
# Partitioned instead between depth_m (the physical transect gradient) and
# temperature_c (the other chemistry variable with enough completeness to
# support the N=34 primary model -- conductivity_ms stays sensitivity-only,
# see 07_dbrda_summary.tsv).
varpart_fit <- vegan::varpart(bray_env, ~depth_m, ~temperature_c, data = env_df_scaled)
saveRDS(varpart_fit, file.path(path_processed, "varpart_fit.rds"))
write_result(
  tibble::as_tibble(varpart_fit$part$indfract, rownames = "fraction"),
  "07_varpart"
)

# Short Xnames: vegan's plot.varpart() sizes its venn circles/margins for
# short labels and clips anything longer ("depth_m"/"temperature_c" both ran
# off the device edge at the default). Full names are still in
# results/07_varpart.tsv; the caption goes in a message, not graphics::title
# (that clipped too, for the same reason).
draw_varpart <- function() {
  plot(varpart_fit, Xnames = c("depth", "temp"),
       bg = c(palette_sample_type()["sediment"], palette_sample_type()["water"]))
}
grDevices::pdf(file.path("plots", "07_varpart.pdf"), width = 5, height = 5)
draw_varpart()
grDevices::dev.off()
grDevices::png(file.path("plots", "07_varpart.png"), width = 5, height = 5, units = "in", res = 300)
draw_varpart()
grDevices::dev.off()

message("[07_environmental_drivers] plots/07_varpart.{pdf,png}: variance partitioning, depth vs. temperature (no spatial/climate blocks -- one cave, AGENTS.md departure #2)")

# --- Mantel / partial Mantel -------------------------------------------
# Community vs. depth distance only -- there is no geographic distance to
# partial out (every sample shares one lat/lon, AGENTS.md departure #2), so
# the classic "community ~ environment | geography" partial Mantel doesn't
# apply. Partialling temperature out of the depth relationship instead.
depth_dist <- stats::dist(env_complete$depth_m)
temp_dist <- stats::dist(env_complete$temperature_c)
libsize_dist <- stats::dist(env_complete$log_library_size)

mantel_depth <- vegan::mantel(bray_env, depth_dist, permutations = 999)
mantel_depth_partial <- tryCatch(
  vegan::mantel.partial(bray_env, depth_dist, temp_dist, permutations = 999),
  error = function(e) NULL
)
# Same idea as PERMANOVA's/db-RDA's library-size covariate above: does the
# community~depth relationship survive partialling out how far apart two
# samples are in library size?
mantel_depth_libsize_partial <- tryCatch(
  vegan::mantel.partial(bray_env, depth_dist, libsize_dist, permutations = 999),
  error = function(e) NULL
)

mantel_results <- tibble::tibble(
  test = c("community ~ depth", "community ~ depth | temperature", "community ~ depth | log_library_size"),
  statistic_r = c(mantel_depth$statistic,
                   if (!is.null(mantel_depth_partial)) mantel_depth_partial$statistic else NA,
                   if (!is.null(mantel_depth_libsize_partial)) mantel_depth_libsize_partial$statistic else NA),
  p_value = c(mantel_depth$signif,
              if (!is.null(mantel_depth_partial)) mantel_depth_partial$signif else NA,
              if (!is.null(mantel_depth_libsize_partial)) mantel_depth_libsize_partial$signif else NA)
)
write_result(mantel_results, "07_mantel")

message(sprintf(
  "[07_environmental_drivers] Mantel community~depth: r = %.3f, p = %.3f",
  mantel_depth$statistic, mantel_depth$signif
))
