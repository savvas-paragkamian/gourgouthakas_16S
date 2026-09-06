# functions.R — shared helpers sourced by 00_setup.R.
#
# Kept dependency-light and phyloseq/qiime2R-free throughout, per PLAN.md.
# No tree-based helpers here: this dataset has no phylogeny (AGENTS.md
# departure #1), so there is no Faith's PD / UniFrac helper to write.

# --- I/O ---------------------------------------------------------------

#' Read the combined ASV table (counts + taxonomy in one wide file).
#'
#' The HiFi-16S-workflow output actually used here
#' (results/hifi/final/best_tax_merged_freq_tax.tsv) is not the generic
#' feature-table.tsv / taxonomy.tsv pair PLAN.md §2 sketches — it's ASV rows
#' with id/Sequence/Taxon/Confidence followed by one count column per sample.
#' This splits it into the two pieces the rest of the pipeline expects.
#' A plain TSV feature table (samples already split, no Taxon column) or a
#' .biom file both fall back to a straight matrix read with no taxonomy.
read_feature_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "biom") {
    bm <- biomformat::read_biom(path)
    mat <- as.matrix(biomformat::biom_data(bm))
    return(list(counts = t(mat), taxonomy = NULL))
  }

  tbl <- readr::read_tsv(path, show_col_types = FALSE)
  meta_cols <- intersect(c("id", "Sequence", "Taxon", "Confidence"), names(tbl))
  sample_cols <- setdiff(names(tbl), meta_cols)

  mat <- as.matrix(tbl[, sample_cols])
  storage.mode(mat) <- "integer"
  rownames(mat) <- tbl$id
  counts <- t(mat) # samples as rows, ASVs as columns

  taxonomy <- NULL
  if ("Taxon" %in% names(tbl)) {
    taxonomy <- dplyr::tibble(
      asv_id = tbl$id,
      lineage = tbl$Taxon,
      confidence = if ("Confidence" %in% names(tbl)) as.numeric(tbl$Confidence) else NA_real_,
      sequence = if ("Sequence" %in% names(tbl)) tbl$Sequence else NA_character_
    )
  }
  list(counts = counts, taxonomy = taxonomy)
}

#' Split a `d__X;p__Y;...` lineage string into rank columns, prefixes stripped.
#' Missing trailing ranks (e.g. no species call) become NA, not "".
split_taxonomy <- function(tbl, lineage_col = "lineage") {
  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  prefixes <- c("d__", "p__", "c__", "o__", "f__", "g__", "s__")

  lineage <- tbl[[lineage_col]]
  parts <- strsplit(lineage, ";", fixed = TRUE)

  rank_mat <- matrix(NA_character_, nrow = length(parts), ncol = length(ranks),
                      dimnames = list(NULL, ranks))
  for (i in seq_along(parts)) {
    p <- trimws(parts[[i]])
    for (j in seq_along(prefixes)) {
      hit <- p[startsWith(p, prefixes[j])]
      if (length(hit) == 1L) {
        val <- sub(prefixes[j], "", hit, fixed = TRUE)
        rank_mat[i, j] <- if (nzchar(val)) val else NA_character_
      }
    }
  }

  out <- dplyr::bind_cols(
    tbl[setdiff(names(tbl), ranks)],
    dplyr::as_tibble(rank_mat)
  )
  out
}

#' Stop with an informative message unless counts/taxonomy/metadata (and,
#' if present, tree) all reference exactly the same ASV/sample ids.
#' `tree = NULL` is the expected case here (no phylogeny) and skips that check.
assert_aligned <- function(counts, taxonomy = NULL, metadata = NULL, tree = NULL) {
  problems <- character(0)

  if (!is.null(metadata)) {
    samples_counts <- rownames(counts)
    samples_meta <- metadata$sample_id
    if (!setequal(samples_counts, samples_meta)) {
      problems <- c(problems, sprintf(
        "counts vs metadata sample mismatch: %d only in counts, %d only in metadata",
        length(setdiff(samples_counts, samples_meta)),
        length(setdiff(samples_meta, samples_counts))
      ))
    }
  }

  if (!is.null(taxonomy)) {
    asvs_counts <- colnames(counts)
    asvs_tax <- taxonomy$asv_id
    if (!setequal(asvs_counts, asvs_tax)) {
      problems <- c(problems, sprintf(
        "counts vs taxonomy ASV mismatch: %d only in counts, %d only in taxonomy",
        length(setdiff(asvs_counts, asvs_tax)),
        length(setdiff(asvs_tax, asvs_counts))
      ))
    }
  }

  if (!is.null(tree)) {
    asvs_counts <- colnames(counts)
    if (!setequal(asvs_counts, tree$tip.label)) {
      problems <- c(problems, "counts vs tree tip-label mismatch")
    }
  }

  if (length(problems) > 0) {
    stop("assert_aligned() failed:\n  - ", paste(problems, collapse = "\n  - "), call. = FALSE)
  }
  invisible(TRUE)
}

#' Save a ggplot as matched .pdf (vector) + .png (raster preview) in plots/.
save_plot <- function(plot, name, w = 7, h = 5, dpi = 300) {
  ggplot2::ggsave(file.path("plots", paste0(name, ".pdf")), plot, width = w, height = h, device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path("plots", paste0(name, ".png")), plot, width = w, height = h, dpi = dpi)
  invisible(plot)
}

#' Write a tidy result to results/<name>.tsv (name without extension).
write_result <- function(x, name) {
  path <- file.path("results", paste0(name, ".tsv"))
  readr::write_tsv(x, path)
  invisible(path)
}

# --- theme & palettes ----------------------------------------------------

theme_soil <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

# sample_type: the primary contrast in this dataset (AGENTS.md departure #5).
palette_sample_type <- function() {
  c(sediment = "#8C6142", water = "#2E7DA6", control = "#999999")
}

# sample_set: the two sediment biological-replicate arms + controls.
palette_sample_set <- function() {
  c(transect = "#B0794A", isolate_source = "#5A3A22", control = "#999999")
}

# Generic categorical palette for taxa (phylum/genus barplots), colorblind-safe,
# with a fixed grey for the "Other" bucket so it reads as residual, not a taxon.
palette_taxa <- function(n) {
  base <- c(
    "#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#B07AA1",
    "#76B7B2", "#EDC948", "#FF9DA7", "#9C755F", "#BAB0AC"
  )
  if (n <= length(base)) return(stats::setNames(base[seq_len(n)], NULL))
  grDevices::colorRampPalette(base)(n)
}

palette_other <- "grey70"

# Continuous depth gradient: shallow (light) -> deep (dark), single hue.
scale_color_depth <- function(...) {
  ggplot2::scale_color_viridis_c(option = "mako", direction = -1, ...)
}
scale_fill_depth <- function(...) {
  ggplot2::scale_fill_viridis_c(option = "mako", direction = -1, ...)
}

# --- tidy wrappers ---------------------------------------------------------

#' vegan::adonis2 result -> one-row-per-term tidy tibble.
tidy_adonis2 <- function(fit) {
  df <- as.data.frame(fit)
  df$term <- rownames(df)
  dplyr::as_tibble(df) |>
    dplyr::relocate(term) |>
    dplyr::rename_with(~ gsub("[^A-Za-z0-9]+", "_", .x))
}

#' Render a small data.frame/tibble as a patchwork-composable "plot" (a table
#' grob wrapped as a ggplot via patchwork::wrap_elements), for figures that
#' pair a plot panel with a compact stats table (e.g. 10_figures.R's
#' PERMANOVA panel).
gridExtra_table <- function(df, title = NULL) {
  tbl_grob <- gridExtra::tableGrob(df, rows = NULL, theme = gridExtra::ttheme_minimal(base_size = 9))
  if (!is.null(title)) {
    title_grob <- grid::textGrob(title, gp = grid::gpar(fontface = "bold", fontsize = 10))
    tbl_grob <- gridExtra::arrangeGrob(title_grob, tbl_grob, ncol = 1, heights = c(0.15, 0.85))
  }
  patchwork::wrap_elements(tbl_grob)
}

#' vegan::envfit result -> one-row-per-variable tidy tibble (vectors only;
#' extend with $factors if/when a factor is fit).
tidy_envfit <- function(fit) {
  v <- fit$vectors
  dplyr::tibble(
    term = rownames(v$arrows),
    v$arrows |> dplyr::as_tibble() |> dplyr::rename(axis1 = 1, axis2 = 2),
    r2 = v$r,
    p_value = v$pvals
  )
}
