# 06_beta_diversity.R — Bray-Curtis + Aitchison distances (no UniFrac: no
# phylogeny, AGENTS.md departure #1), PCoA/NMDS, PERMANOVA, betadisper.

source("scripts/00_setup.R")

counts_relabund <- readRDS(file.path(path_processed, "counts_relabund.rds"))
counts_clr <- readRDS(file.path(path_processed, "counts_clr.rds"))
metadata <- readRDS(file.path(path_processed, "metadata_clean.rds"))

# Uses whichever matrix 02b_controls.R's replicate-concordance step decided on
# (results/replicate_concordance_decision.tsv). Both replicate levels are NOT
# pooled by default there, so this runs on the full 46-library (control-free,
# per AGENTS.md), unpooled matrix -- ordinations below are colored by
# tech_rep/bio_rep as a visual re-check of that call before interpreting any
# other grouping.
pooling_decision <- readr::read_tsv("results/replicate_concordance_decision.tsv", show_col_types = FALSE)
if (any(pooling_decision$pooled_by_default)) {
  warning("[06_beta_diversity] pooling_decision has pooled_by_default = TRUE somewhere -- this script still assumes the unpooled matrix; update it if that decision changes.")
}

# --- distances ------------------------------------------------------------
dist_bray <- vegan::vegdist(counts_relabund, method = "bray")
dist_aitchison <- vegan::vegdist(counts_clr, method = "euclidean")

saveRDS(dist_bray, file.path(path_processed, "dist_bray.rds"))
saveRDS(dist_aitchison, file.path(path_processed, "dist_aitchison.rds"))

# --- ordinations ------------------------------------------------------------
run_ordination <- function(dist_obj, metadata, label) {
  pcoa <- ape::pcoa(dist_obj)
  pcoa_axes <- pcoa$vectors[, 1:2]
  colnames(pcoa_axes) <- c("PCoA1", "PCoA2")
  pcoa_var <- pcoa$values$Relative_eig[1:2] * 100

  nmds <- vegan::metaMDS(dist_obj, k = 2, trymax = 100, trace = FALSE)

  coords <- tibble::as_tibble(pcoa_axes, rownames = "sample_id") |>
    dplyr::left_join(tibble::as_tibble(vegan::scores(nmds), rownames = "sample_id"), by = "sample_id") |>
    dplyr::left_join(metadata, by = "sample_id") |>
    dplyr::mutate(distance = label, nmds_stress = nmds$stress)

  list(coords = coords, pcoa_var = pcoa_var, nmds_stress = nmds$stress)
}

ord_bray <- run_ordination(dist_bray, metadata, "bray")
ord_aitchison <- run_ordination(dist_aitchison, metadata, "aitchison")

all_coords <- dplyr::bind_rows(ord_bray$coords, ord_aitchison$coords)
write_result(all_coords, "06_ordination_coords")

message(sprintf(
  "[06_beta_diversity] NMDS stress: bray = %.3f, aitchison = %.3f (>0.2 = interpret cautiously)",
  ord_bray$nmds_stress, ord_aitchison$nmds_stress
))

plot_ordination <- function(coords, x, y, color_var, palette, title) {
  p <- ggplot2::ggplot(coords, ggplot2::aes(x = .data[[x]], y = .data[[y]], color = .data[[color_var]])) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::labs(title = title, color = color_var)
  if (!is.null(palette)) p <- p + ggplot2::scale_color_manual(values = palette)
  p
}

p_pcoa_type <- plot_ordination(ord_bray$coords, "PCoA1", "PCoA2", "sample_type", palette_sample_type(),
                                sprintf("PCoA (Bray-Curtis), PC1/PC2 = %.1f%%/%.1f%% var", ord_bray$pcoa_var[1], ord_bray$pcoa_var[2]))
save_plot(p_pcoa_type, "06_pcoa_bray_sample_type", w = 6, h = 5)

p_pcoa_aitch_type <- plot_ordination(ord_aitchison$coords, "PCoA1", "PCoA2", "sample_type", palette_sample_type(),
                                      sprintf("PCoA (Aitchison), PC1/PC2 = %.1f%%/%.1f%% var", ord_aitchison$pcoa_var[1], ord_aitchison$pcoa_var[2]))
save_plot(p_pcoa_aitch_type, "06_pcoa_aitchison_sample_type", w = 6, h = 5)

# Replicate re-check: color the same Bray PCoA by tech_rep and bio_rep. If
# replicates were truly interchangeable, these should NOT separate samples
# the way sample_type/depth does -- if they visibly do, that's independent
# support for the concordance flags already written in 02b_controls.R.
p_pcoa_techrep <- plot_ordination(ord_bray$coords, "PCoA1", "PCoA2", "tech_rep", NULL,
                                   "PCoA (Bray-Curtis), colored by technical replicate")
save_plot(p_pcoa_techrep, "06_pcoa_bray_tech_rep", w = 6, h = 5)

p_pcoa_biorep <- plot_ordination(ord_bray$coords |> dplyr::filter(sample_type == "sediment"),
                                  "PCoA1", "PCoA2", "bio_rep", palette_sample_set(),
                                  "PCoA (Bray-Curtis), sediment only, colored by biological replicate arm")
save_plot(p_pcoa_biorep, "06_pcoa_bray_bio_rep", w = 6, h = 5)

p_pcoa_depth <- plot_ordination(ord_bray$coords, "PCoA1", "PCoA2", "depth_m", NULL,
                                 "PCoA (Bray-Curtis), colored by depth") + scale_color_depth()
save_plot(p_pcoa_depth, "06_pcoa_bray_depth", w = 6, h = 5)

# --- PERMANOVA + betadisper -------------------------------------------------
# Controls were already dropped in 02b_controls.R (QC artifacts, not part of
# the environmental gradient) -- `metadata` here is already sediment + water
# only, no filter needed.
eco_samples <- metadata$sample_id
meta_eco <- metadata

run_permanova <- function(dist_obj, samples, meta_eco, label) {
  d <- stats::as.dist(as.matrix(dist_obj)[samples, samples])

  # `by = "margin"` on a formula that carries an interaction only reports
  # the highest-order term -- main effects that participate in an
  # interaction aren't marginally estimable in that framework, so
  # `sample_type * depth_m` with by="margin" silently returned *just* the
  # interaction row (R2=0.027, a middling p) and dropped both main effects
  # entirely. Checked what those main effects actually look like: fit
  # separately (additive model, by="margin", both terms individually
  # meaningful) they're each highly significant (R2~0.05-0.06, p=0.001) --
  # a materially stronger and more informative result than the interaction
  # alone. Both models are reported now: the additive one as primary (each
  # main effect's own marginal contribution), the full interactive one
  # sequentially (by="terms") to also show whether sample_type's effect
  # depends on depth, which the additive model can't address.
  fit_additive <- vegan::adonis2(d ~ sample_type + depth_m, data = meta_eco, by = "margin", permutations = 999)
  fit_interaction <- vegan::adonis2(d ~ sample_type * depth_m, data = meta_eco, by = "terms", permutations = 999)

  disp <- vegan::betadisper(d, meta_eco$sample_type)
  disp_test <- vegan::permutest(disp, permutations = 999)
  list(
    adonis = dplyr::bind_rows(
      tidy_adonis2(fit_additive) |> dplyr::mutate(model = "additive_marginal"),
      tidy_adonis2(fit_interaction) |> dplyr::mutate(model = "interaction_sequential")
    ) |> dplyr::mutate(distance = label),
    betadisper_p = disp_test$tab$`Pr(>F)`[1]
  )
}

perm_bray <- run_permanova(dist_bray, eco_samples, meta_eco, "bray")
perm_aitchison <- run_permanova(dist_aitchison, eco_samples, meta_eco, "aitchison")

write_result(dplyr::bind_rows(perm_bray$adonis, perm_aitchison$adonis), "06_permanova")
write_result(
  tibble::tibble(distance = c("bray", "aitchison"),
                  betadisper_p_value = c(perm_bray$betadisper_p, perm_aitchison$betadisper_p)),
  "06_betadisper"
)

bray_additive <- perm_bray$adonis |> dplyr::filter(model == "additive_marginal", !term %in% c("Residual", "Total"))
message(sprintf(
  "[06_beta_diversity] PERMANOVA (Bray, sediment+water, additive): %s",
  paste(sprintf("%s R2=%.3f p=%s", bray_additive$term, bray_additive$R2, bray_additive$Pr_F_), collapse = "; ")
))
message(sprintf(
  "[06_beta_diversity] betadisper p = %.3f (bray), %.3f (aitchison) -- a small p means dispersion differs between groups and the PERMANOVA result may partly reflect that rather than location",
  perm_bray$betadisper_p, perm_aitchison$betadisper_p
))

# --- dissimilarity heatmap: sediment | water, clustered, library size shown -
# Same design language as 02b_controls.R's replicate-concordance figure and
# 05_taxonomic_composition.R's prevalence-abundance figure: two separate,
# bordered panels (sediment | water) rather than one shared plot, so each
# gets its own sample count and doesn't visually compress the other. Built
# with pheatmap (already available, the same tool isd_archive_scripts used
# for this exact plot type) rather than hand-rolled ggplot tiles + a
# separately-rendered dendrogram -- getting a base-graphics dendrogram to
# line up pixel-for-pixel with an independently-built ggplot heatmap is
# fragile, and pheatmap does clustering + dendrogram + annotation strip
# together, correctly, in one call.
counts_clean_hm <- readRDS(file.path(path_processed, "counts_clean.rds"))
libsize_hm <- rowSums(counts_clean_hm)
bray_mat_full <- as.matrix(dist_bray)

# library_size as 5 discrete steps (not a continuous 2-color gradient) --
# log-spaced bins across each panel's own range, since library size spans
# orders of magnitude and a linear gradient would make most cells look
# identical. Grayscale, deliberately not blue, so it can't be confused with
# the dissimilarity fill.
n_library_steps <- 5L
library_step_colors <- grDevices::colorRampPalette(c("white", "grey15"))(n_library_steps)

make_dissimilarity_heatmap <- function(type_label, out_path) {
  ids <- metadata$sample_id[metadata$sample_type == type_label]
  m <- bray_mat_full[ids, ids]
  lib <- libsize_hm[ids]
  lib_bins <- cut(log10(lib), breaks = n_library_steps, dig.lab = 10)
  levels(lib_bins) <- sprintf("%s reads", scales::label_comma()(round(10^c(
    min(log10(lib)), min(log10(lib)) + seq_len(n_library_steps) * diff(range(log10(lib))) / n_library_steps
  )[-1])))
  ann_col <- data.frame(library_size = lib_bins, row.names = ids)

  # clustering_distance_* must be set explicitly to the actual Bray-Curtis
  # values (as.dist(m)) -- left at pheatmap's default, it would instead
  # compute its own (Euclidean) distance BETWEEN ROWS OF m, clustering
  # samples by how similar their *dissimilarity profiles* are to every other
  # sample rather than by their direct pairwise Bray-Curtis distance. Average
  # linkage (UPGMA): the conventional choice for microbial community
  # dissimilarity dendrograms.
  #
  # cellwidth = cellheight (not left to auto-stretch to a fixed width/height)
  # is what actually makes the matrix cells square -- pheatmap's default
  # sizing fills whatever width/height you give it, which for a non-square
  # sample count or aspect doesn't come out square on its own. filename with
  # no explicit width/height lets pheatmap compute the canvas size FROM the
  # cell count at that fixed cell size, instead of the other way around.
  pheatmap::pheatmap(
    m,
    clustering_distance_rows = stats::as.dist(m),
    clustering_distance_cols = stats::as.dist(m),
    clustering_method = "average",
    color = grDevices::colorRampPalette(c("white", "#08306B"))(100), # white -> dark navy
    annotation_col = ann_col,
    annotation_colors = list(library_size = stats::setNames(library_step_colors, levels(lib_bins))),
    main = sprintf("%s (n=%d)", type_label, length(ids)),
    fontsize = 11, fontsize_row = 7, fontsize_col = 7,
    cellwidth = 13, cellheight = 13,
    filename = out_path
  )
}

path_hm_sediment <- file.path("plots", "_tmp_dissimilarity_sediment.png")
path_hm_water <- file.path("plots", "_tmp_dissimilarity_water.png")
make_dissimilarity_heatmap("sediment", path_hm_sediment)
make_dissimilarity_heatmap("water", path_hm_water)

# Embed both pheatmap PNGs as rasters and compose side by side (same
# png::readPNG + grid::rasterGrob technique 10_figures.R uses to embed
# 07_varpart.png -- a pure grid rasterGrob composes with ggplot/patchwork,
# a captured base-graphics device does not). cellwidth=cellheight above
# makes each PNG's OWN cells square, but patchwork's `|` still gives both
# panels equal width by default regardless of their native pixel aspect
# ratio, which would re-stretch them right back out -- plot_layout(widths=)
# set proportional to each PNG's actual pixel aspect ratio (at their shared
# row height) is what keeps the composed figure from undoing that.
embed_png <- function(path) patchwork::wrap_elements(grid::rasterGrob(png::readPNG(path), interpolate = TRUE))
img_dim <- function(path) dim(png::readPNG(path))[1:2] # c(height_px, width_px)
dim_sed <- img_dim(path_hm_sediment)
dim_wat <- img_dim(path_hm_water)
aspect_sed <- dim_sed[2] / dim_sed[1]
aspect_wat <- dim_wat[2] / dim_wat[1]

p_dissimilarity_heatmap <- (embed_png(path_hm_sediment) | embed_png(path_hm_water)) +
  patchwork::plot_layout(widths = c(aspect_sed, aspect_wat)) +
  patchwork::plot_annotation(
    title = "Bray-Curtis dissimilarity, clustered (average linkage)",
    subtitle = sprintf("Column strip = library size (%d steps, white = shallow, dark = deep)", n_library_steps),
    theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 20, face = "bold"),
                            plot.subtitle = ggplot2::element_text(size = 12))
  )
save_plot(p_dissimilarity_heatmap, "06_dissimilarity_heatmap", w = 14, h = 7.9, dpi = 400)
# The two per-panel PNGs were only ever intermediate build artifacts for the
# composed figure above -- moved into plots/_tmp/ (not deleted; this project
# doesn't use rm) so they don't clutter the main plots/ listing as if they
# were independent, publishable outputs.
dir.create(file.path("plots", "_tmp"), showWarnings = FALSE)
file.rename(path_hm_sediment, file.path("plots", "_tmp", basename(path_hm_sediment)))
file.rename(path_hm_water, file.path("plots", "_tmp", basename(path_hm_water)))

message("[06_beta_diversity] wrote plots/06_dissimilarity_heatmap.{pdf,png}")
