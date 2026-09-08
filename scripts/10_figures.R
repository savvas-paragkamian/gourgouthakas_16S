# 10_figures.R — multi-panel publication figures from patchwork, built from
# already-computed results/*.tsv (not re-running any stats). No Fig
# "map + alpha": there's no map here (AGENTS.md departure #2) -- Fig 1 pairs
# the depth-transect view with alpha diversity instead.

source("scripts/00_setup.R")
library(patchwork)

metadata <- readRDS(file.path(path_processed, "metadata_clean.rds"))
alpha <- readr::read_tsv("results/alpha_diversity.tsv", show_col_types = FALSE)
ord_coords <- readr::read_tsv("results/06_ordination_coords.tsv", show_col_types = FALSE)
permanova <- readr::read_tsv("results/06_permanova.tsv", show_col_types = FALSE)
envfit_tbl <- readr::read_tsv("results/07_envfit.tsv", show_col_types = FALSE)
phylum_long <- readr::read_tsv("results/05_phylum_relabund.tsv", show_col_types = FALSE)
aldex_result <- readr::read_tsv("results/da_aldex2_sample_type.tsv", show_col_types = FALSE)
decay_df <- readr::read_tsv("results/09_distance_decay.tsv", show_col_types = FALSE)
varpart_fit <- readRDS(file.path(path_processed, "varpart_fit.rds"))

# --- Fig 1: depth transect + alpha diversity --------------------------------
p1a <- ggplot2::ggplot(decay_df, ggplot2::aes(depth_distance, bray_dissimilarity, color = pair_type)) +
  ggplot2::geom_point(alpha = 0.4, size = 1) +
  ggplot2::geom_smooth(method = "loess", formula = y ~ x, color = "black", se = TRUE) +
  ggplot2::labs(x = "Depth distance (m)", y = "Bray-Curtis dissimilarity", color = NULL, title = "A. Distance-decay in depth space")

p1b <- alpha |> # controls already dropped upstream in 02b_controls.R
  ggplot2::ggplot(ggplot2::aes(depth_m, shannon, color = sample_type)) +
  ggplot2::geom_point(size = 1.5, alpha = 0.7) +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  ggplot2::scale_color_manual(values = palette_sample_type()) +
  ggplot2::scale_x_continuous(breaks = scales::breaks_width(200)) +
  ggplot2::labs(x = "Depth (m)", y = "Shannon diversity", color = NULL, title = "B. Alpha diversity vs. depth")

fig1 <- p1a + p1b + patchwork::plot_layout(guides = "collect")
save_plot(fig1, "fig1_depth_transect_alpha", w = 11, h = 4.5)

# --- Fig 2: ordination + PERMANOVA ------------------------------------------
p2a <- ord_coords |>
  dplyr::filter(distance == "bray") |>
  ggplot2::ggplot(ggplot2::aes(PCoA1, PCoA2, color = sample_type)) +
  ggplot2::geom_point(size = 2, alpha = 0.8) +
  ggplot2::scale_color_manual(values = palette_sample_type()) +
  ggplot2::labs(title = "A. PCoA (Bray-Curtis)", color = NULL)

# additive_marginal only: each main effect's own marginal R2/p, the cleaner
# primary result (see 06_beta_diversity.R for why the interaction model
# can't report main effects marginally, and results/06_permanova.tsv for
# both models in full).
permanova_tbl <- permanova |>
  dplyr::filter(distance == "bray", model == "additive_marginal", !term %in% c("Residual", "Total")) |>
  dplyr::transmute(Term = term, R2 = sprintf("%.3f", R2), F = sprintf("%.2f", F), `p` = sprintf("%.3f", Pr_F_))
p2b <- gridExtra_table(permanova_tbl, "B. PERMANOVA (Bray, sediment + water)")

fig2 <- p2a + p2b
save_plot(fig2, "fig2_ordination_permanova", w = 10, h = 4.5)

# --- Fig 3: db-RDA triplot + variance partitioning --------------------------
p3a <- ord_coords |>
  dplyr::filter(distance == "bray") |>
  ggplot2::ggplot(ggplot2::aes(PCoA1, PCoA2, color = sample_type)) +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_color_manual(values = palette_sample_type()) +
  ggplot2::geom_segment(data = envfit_tbl, ggplot2::aes(x = 0, y = 0, xend = axis1 * sqrt(r2), yend = axis2 * sqrt(r2)),
                         inherit.aes = FALSE, arrow = grid::arrow(length = grid::unit(0.15, "cm"))) +
  ggrepel::geom_text_repel(data = envfit_tbl, ggplot2::aes(axis1 * sqrt(r2), axis2 * sqrt(r2), label = term),
                            inherit.aes = FALSE, size = 3) +
  ggplot2::labs(title = "A. envfit triplot", color = NULL)

# plot.varpart() is base graphics, and it doesn't compose with a ggplot in
# the same device (base layout()/print(<ggplot>) leaves the venn diagram
# blank -- ggplot ignores the base layout; grid.grabExpr() *also* came back
# blank on this specific base plot, tried and confirmed empty). Robust fix:
# reuse the already-rendered, already-correct plots/07_varpart.png (its own
# dedicated canvas, so "depth_m"/"temperature_c" aren't clipped) as a raster
# image and embed that -- a pure grid rasterGrob composes fine in patchwork,
# no base-graphics capture involved.
varpart_png <- png::readPNG(file.path("plots", "07_varpart.png"))
varpart_grob <- grid::grobTree(
  grid::rasterGrob(varpart_png, y = 0.45, height = 0.9),
  grid::textGrob("B. Variance partitioning", y = 0.97, gp = grid::gpar(fontface = "bold", fontsize = 11))
)
p3b <- patchwork::wrap_elements(varpart_grob)
fig3 <- p3a + p3b
save_plot(fig3, "fig3_drivers_varpart", w = 10, h = 5)

# --- Fig 4: composition + differential abundance ----------------------------
top_phyla <- phylum_long |> dplyr::group_by(phylum) |> dplyr::summarise(m = mean(rel_abund)) |> dplyr::slice_max(m, n = 8) |> dplyr::pull(phylum)
# Explicit sample order by sample_type then depth_m -- fct_reorder()'s
# automatic median-based ordering chokes on duplicate sample_id rows (one
# per phylum) once there are ties; see 05's fix.
sample_order <- metadata |> dplyr::arrange(sample_type, depth_m) |> dplyr::pull(sample_id)
p4a <- phylum_long |>
  dplyr::mutate(phylum2 = ifelse(phylum %in% top_phyla, phylum, "Other")) |>
  dplyr::left_join(metadata |> dplyr::select(sample_id, sample_type, depth_m), by = "sample_id") |>
  dplyr::group_by(sample_id, phylum2, sample_type, depth_m) |>
  dplyr::summarise(rel_abund = sum(rel_abund), .groups = "drop") |>
  dplyr::mutate(sample_id = factor(sample_id, levels = intersect(sample_order, sample_id))) |>
  ggplot2::ggplot(ggplot2::aes(sample_id, rel_abund, fill = phylum2)) +
  ggplot2::geom_col() +
  ggplot2::facet_grid(~sample_type, scales = "free_x", space = "free_x") +
  ggplot2::labs(x = NULL, y = "Relative abundance", fill = "Phylum", title = "A. Phylum composition") +
  ggplot2::theme(axis.text.x = ggplot2::element_blank())

p4b <- aldex_result |>
  dplyr::mutate(significant = wi.eBH < 0.05) |>
  ggplot2::ggplot(ggplot2::aes(effect, -log10(wi.eBH), color = significant)) +
  ggplot2::geom_point(alpha = 0.5, size = 1) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ggplot2::scale_color_manual(values = c(`TRUE` = "#E15759", `FALSE` = "grey70")) +
  ggplot2::labs(x = "ALDEx2 effect size", y = "-log10(q)", title = "B. Differential abundance (sediment vs. water)", color = "q<0.05")

fig4 <- p4a / p4b
save_plot(fig4, "fig4_composition_da", w = 9, h = 8)

# --- figure index -----------------------------------------------------------
figure_index <- tibble::tibble(
  figure = c("fig1_depth_transect_alpha", "fig2_ordination_permanova", "fig3_drivers_varpart", "fig4_composition_da"),
  panels = c("A: distance-decay in depth space; B: alpha diversity vs. depth",
             "A: PCoA (Bray-Curtis) by sample type; B: PERMANOVA table",
             "A: envfit triplot; B: variance partitioning (depth vs. temperature)",
             "A: phylum composition; B: ALDEx2 volcano (sediment vs. water)"),
  files_pdf = sprintf("plots/%s.pdf", c("fig1_depth_transect_alpha", "fig2_ordination_permanova", "fig3_drivers_varpart", "fig4_composition_da")),
  files_png = sprintf("plots/%s.png", c("fig1_depth_transect_alpha", "fig2_ordination_permanova", "fig3_drivers_varpart", "fig4_composition_da"))
)
write_result(figure_index, "figure_index")

message("[10_figures] wrote fig1-fig4 + figure_index.tsv")
