# 05_taxonomic_composition.R — Phylum/Genus aggregation, stacked barplots,
# core taxa membership.

source("scripts/00_setup.R")

counts <- readRDS(file.path(path_processed, "counts_clean.rds"))
taxonomy <- readRDS(file.path(path_processed, "taxonomy_clean.rds"))
metadata <- readRDS(file.path(path_processed, "metadata_clean.rds"))

lib_sizes <- rowSums(counts)

counts_long <- counts |>
  as.data.frame() |>
  tibble::rownames_to_column("sample_id") |>
  tidyr::pivot_longer(-sample_id, names_to = "asv_id", values_to = "count") |>
  dplyr::filter(count > 0) |>
  dplyr::left_join(taxonomy |> dplyr::select(asv_id, phylum, genus), by = "asv_id") |>
  dplyr::left_join(metadata |> dplyr::select(sample_id, sample_type, sample_set, site, depth_m), by = "sample_id")

# --- aggregate + relative abundance (of the FULL sample, not of the
# aggregated-taxon subtotal -- every count already has a non-NA phylum,
# since unassigned-phylum ASVs were dropped in 02_qc_filter) ---------------
aggregate_taxon <- function(long, rank) {
  long |>
    dplyr::group_by(sample_id, .data[[rank]]) |>
    dplyr::summarise(count = sum(count), .groups = "drop") |>
    dplyr::left_join(tibble::tibble(sample_id = names(lib_sizes), library_size = lib_sizes), by = "sample_id") |>
    dplyr::mutate(rel_abund = count / library_size)
}

phylum_long <- aggregate_taxon(counts_long, "phylum")
genus_long <- aggregate_taxon(counts_long, "genus") |> dplyr::filter(!is.na(genus))

write_result(phylum_long, "05_phylum_relabund")
write_result(genus_long, "05_genus_relabund")

# --- top-N barplots (remainder -> "Other") --------------------------------
top_n_barplot <- function(long, rank, n_top, metadata, title) {
  top_taxa <- long |>
    dplyr::group_by(.data[[rank]]) |>
    dplyr::summarise(mean_rel = mean(rel_abund), .groups = "drop") |>
    dplyr::slice_max(mean_rel, n = n_top) |>
    dplyr::pull(.data[[rank]])

  # Order samples by sample_type then depth_m explicitly -- fct_reorder()'s
  # automatic (median-based) ordering chokes on NA depth_m (the controls,
  # which have no depth) once there are duplicate sample_id rows (one per
  # taxon), so the level order is built by hand instead.
  sample_order <- metadata |>
    dplyr::arrange(sample_type, depth_m) |>
    dplyr::pull(sample_id)

  plot_df <- long |>
    dplyr::mutate(taxon = ifelse(.data[[rank]] %in% top_taxa, .data[[rank]], "Other")) |>
    dplyr::group_by(sample_id, taxon) |>
    dplyr::summarise(rel_abund = sum(rel_abund), .groups = "drop") |>
    dplyr::left_join(metadata |> dplyr::select(sample_id, sample_type, site, depth_m), by = "sample_id") |>
    dplyr::mutate(
      taxon = forcats::fct_relevel(taxon, "Other", after = 0),
      sample_id = factor(sample_id, levels = intersect(sample_order, sample_id))
    )

  pal <- stats::setNames(c(palette_other, palette_taxa(length(top_taxa))), levels(plot_df$taxon))

  ggplot2::ggplot(plot_df, ggplot2::aes(x = sample_id, y = rel_abund, fill = taxon)) +
    ggplot2::geom_col() +
    ggplot2::facet_grid(~sample_type, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::labs(x = NULL, y = "Relative abundance", fill = rank, title = title) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
}

p_phylum <- top_n_barplot(phylum_long, "phylum", 10, metadata, "Phylum composition (top 10, ordered by depth within sample type)")
save_plot(p_phylum, "05_phylum_barplot", w = 10, h = 6)

p_genus <- top_n_barplot(genus_long, "genus", 15, metadata, "Genus composition (top 15, ordered by depth within sample type)")
save_plot(p_genus, "05_genus_barplot", w = 10, h = 6)

# --- core taxa --------------------------------------------------------
# Core within each sample_type: present (rel_abund > 0) in >= 80% of that
# type's samples, AND mean relative abundance >= 0.1% -- prevalence alone
# would let a lot of things that are barely detectable in, mean abundance
# alone would let a thing present in only 2 samples in.
core_taxa <- function(long, rank, metadata, prevalence_thresh = 0.8, mean_abund_thresh = 0.001) {
  n_samples_per_type <- metadata |> dplyr::count(sample_type, name = "n_total")

  # `long` (genus_long/phylum_long) only has rows where the taxon is present
  # (built from counts_long, which is pre-filtered to count > 0) -- absent
  # samples have no row at all, not a rel_abund = 0 row. So `mean(rel_abund)`
  # over these rows is the mean *given presence*, not the true mean across
  # every sample of that type; it silently omits the zeros from absent
  # samples and inflates the reported abundance by up to 1/prevalence_thresh
  # (e.g. 25% high right at the 80% prevalence cutoff -- confirmed
  # numerically: Parabacteroides in control read 0.218 conditional vs. 0.174
  # true). Fixed by summing (not averaging) here, then dividing by the true
  # n_total (including absent samples) after the join, below.
  long |>
    dplyr::left_join(metadata |> dplyr::select(sample_id, sample_type), by = "sample_id") |>
    dplyr::group_by(sample_type, .data[[rank]]) |>
    dplyr::summarise(
      n_present = sum(rel_abund > 0),
      sum_rel_abund = sum(rel_abund),
      .groups = "drop"
    ) |>
    dplyr::left_join(n_samples_per_type, by = "sample_type") |>
    dplyr::mutate(prevalence = n_present / n_total,
                   mean_rel_abund = sum_rel_abund / n_total, # zero-padded, true mean
                   is_core = prevalence >= prevalence_thresh & mean_rel_abund >= mean_abund_thresh) |>
    dplyr::select(-sum_rel_abund) |>
    dplyr::filter(is_core) |>
    dplyr::arrange(sample_type, dplyr::desc(mean_rel_abund))
}

core_genus <- core_taxa(genus_long, "genus", metadata)
write_result(core_genus, "05_core_genera")

p_core <- ggplot2::ggplot(core_genus, ggplot2::aes(
  x = forcats::fct_reorder(genus, mean_rel_abund), y = mean_rel_abund, fill = sample_type
)) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values = palette_sample_type()) +
  ggplot2::labs(x = NULL, y = "Mean relative abundance",
                title = sprintf("Core genera (>=80%% prevalence, >=0.1%% mean abundance) by sample type"))
save_plot(p_core, "05_core_genera", w = 7, h = max(4, 0.25 * dplyr::n_distinct(core_genus$genus)))

# Separate per-sample_type plots too: sediment (36 samples, huge depth range,
# only 1 genus clears 80% prevalence) and water (10 samples, far more
# homogeneous, 17 genera clear it) are so lopsided that one shared-scale
# plot badly serves whichever group has fewer bars -- each gets its own
# plot, scaled to its own genera and abundance range.
for (st in unique(core_genus$sample_type)) {
  core_st <- core_genus |> dplyr::filter(sample_type == st)
  write_result(core_st, sprintf("05_core_genera_%s", st))
  if (nrow(core_st) == 0) {
    message(sprintf("[05_taxonomic_composition] no core genera for sample_type = %s at these thresholds -- no plot written", st))
    next
  }
  p_core_st <- ggplot2::ggplot(core_st, ggplot2::aes(
    x = forcats::fct_reorder(genus, mean_rel_abund), y = mean_rel_abund
  )) +
    ggplot2::geom_col(fill = palette_sample_type()[[st]]) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Mean relative abundance",
                  title = sprintf("Core genera: %s (>=80%% prevalence, >=0.1%% mean abundance, N=%d)",
                                   st, sum(metadata$sample_type == st)))
  save_plot(p_core_st, sprintf("05_core_genera_%s", st), w = 7, h = max(3, 0.3 * nrow(core_st)))
}

message(sprintf(
  "[05_taxonomic_composition] %d phyla, %d genera aggregated; %d core genera identified (%s)",
  dplyr::n_distinct(phylum_long$phylum), dplyr::n_distinct(genus_long$genus), nrow(core_genus),
  paste(sprintf("%s=%d", names(table(core_genus$sample_type)), table(core_genus$sample_type)), collapse = ", ")
))
