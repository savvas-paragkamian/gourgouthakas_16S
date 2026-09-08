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

# --- ASV-level prevalence-abundance -----------------------------------------
# One point per ASV per sample_type: prevalence = how many of that type's
# samples it's present in; mean relative abundance = the TRUE mean across
# every sample of that type, absent samples counted as zero -- not the mean
# conditional on presence, which silently drops the zeros and inflates the
# reported abundance (the exact bug core_taxa() below was fixed for
# earlier this session; same fix applied here from the start).
n_samples_per_type <- metadata |> dplyr::count(sample_type, name = "n_total")

asv_prevalence_abundance <- counts_long |>
  dplyr::mutate(rel_abund = count / lib_sizes[sample_id]) |>
  dplyr::group_by(asv_id, sample_type) |>
  dplyr::summarise(prevalence = dplyr::n(), sum_rel_abund = sum(rel_abund), .groups = "drop") |>
  dplyr::left_join(n_samples_per_type, by = "sample_type") |>
  dplyr::mutate(mean_rel_abund = sum_rel_abund / n_total) |>
  dplyr::left_join(
    taxonomy |> dplyr::transmute(
      asv_id,
      # Species (GTDB-suffix-stripped, per T1) when resolved -- 80.5% of
      # ASVs overall, 46/48 of what ends up labeled here, so this is a real
      # improvement, not usually falling through. "Genus_species" ->
      # "Genus species": the underscore is a GTDB storage format, not part
      # of the name, and the italic font below is the standard binomial
      # convention -- genus alone doesn't need it as much as a real binomial
      # does. Falls back to genus, then "unclassified", when species-level
      # assignment failed.
      taxon_label = gsub("_", " ", dplyr::coalesce(species_norm, species, genus_norm, genus, "unclassified")),
      phylum
    ),
    by = "asv_id"
  ) |>
  dplyr::select(asv_id, sample_type, prevalence, n_total, mean_rel_abund, taxon_label, phylum)

write_result(asv_prevalence_abundance, "05_asv_prevalence_abundance")

# Label ASVs meeting a per-panel prevalence floor (sediment >= 10/36
# samples, water >= 6/10) -- labeled with genus (GTDB-suffix-stripped, per
# 02b_controls.R's T1; "unclassified" if genus-level assignment failed),
# not the opaque asv_id.
label_min_prevalence <- c(sediment = 10L, water = 6L)

# x axis needs the raw sample count AND that count as a % of the panel's
# own total -- sediment (n=36) and water (n=10) have different totals, so
# one shared facet_wrap scale can't label both correctly (the % denominator
# would be wrong for whichever panel didn't set it). Two separate ggplot
# panels, each with its own correct denominator, composed side by side with
# patchwork, instead of one facet_wrap -- also gives a cleaner "two separate
# plots" look than a shared panel ever would.
y_range <- range(asv_prevalence_abundance$mean_rel_abund, na.rm = TRUE)

# All labels across BOTH panels share one phylum palette/legend (built once,
# outside the per-panel function) -- a phylum present in only one panel still
# needs a stable color, and patchwork can only merge two panels' legends into
# one if their color scales are identical, breaks included.
all_labels <- asv_prevalence_abundance |>
  dplyr::filter(prevalence >= label_min_prevalence[sample_type])
label_phyla <- sort(unique(all_labels$phylum))
# A separate, more saturated palette from palette_taxa()'s muted/pastel
# default -- fine for a large filled bar area (05_phylum_barplot), too pale
# for small italic text on white that needs to read from across a room.
# Bold, print/projector-safe hues, dark enough for on-white legibility.
vibrant_taxa_colors <- c(
  "#D62728", "#1F77B4", "#2CA02C", "#E6550D", "#9467BD",
  "#8C2D04", "#C51B8A", "#00A0A0", "#0050A0", "#B8860B", "#7B3F00"
)
phylum_colors <- stats::setNames(
  if (length(label_phyla) <= length(vibrant_taxa_colors)) vibrant_taxa_colors[seq_along(label_phyla)]
  else grDevices::colorRampPalette(vibrant_taxa_colors)(length(label_phyla)),
  label_phyla
)
# The identical scale OBJECT, not just an equivalent scale_color_manual()
# call, is what patchwork's guides="collect" actually needs to recognize
# two panels' legends as mergeable -- two separately-constructed-but-
# identical scales weren't enough (confirmed: still showed as two legend
# blocks after drop=FALSE alone), reusing one object is.
shared_phylum_scale <- ggplot2::scale_color_manual(values = phylum_colors, name = "Phylum", drop = FALSE)

make_prevalence_panel <- function(type_label, show_legend = TRUE) {
  df <- asv_prevalence_abundance |> dplyr::filter(sample_type == type_label)
  n_total <- df$n_total[1]
  # phylum forced to one explicit, identically-leveled factor (not left as
  # character) via label_phyla/drop=FALSE, so both panels' legends would
  # show the identical full phylum set anyway -- patchwork's guides="collect"
  # never actually merged them into one regardless (tried both drop=FALSE
  # and a literal shared scale object; ggrepel's guide key seems to defeat
  # its merge-detection), so show_legend just keeps one copy (sediment's)
  # and drops the redundant second one, rather than fighting that further.
  labels <- df |>
    dplyr::filter(prevalence >= label_min_prevalence[[type_label]]) |>
    dplyr::mutate(phylum = factor(phylum, levels = label_phyla))
  brks <- scales::breaks_pretty()(c(0, max(df$prevalence)))
  brks <- brks[brks >= 0 & brks <= max(df$prevalence)]

  ggplot2::ggplot(df, ggplot2::aes(x = prevalence, y = mean_rel_abund)) +
    ggplot2::geom_point(alpha = 0.25, size = 1.1, color = palette_sample_type()[[type_label]]) +
    # Labeled points get a small black outline of their own -- with up to 30
    # labels crowded at the high-prevalence end, ggrepel pushes most of the
    # label *text* well away from its point to avoid overlap, and a thin
    # leader line alone is easy to lose track of; a visibly darker/thicker
    # segment plus a marked anchor point makes "this text belongs to that
    # point" traceable even in the crowd.
    ggplot2::geom_point(data = labels, shape = 21, color = "black", fill = NA, size = 2, stroke = 0.8) +
    ggrepel::geom_text_repel(
      data = labels, ggplot2::aes(label = taxon_label, color = phylum),
      size = 2.9, fontface = "bold.italic",
      max.overlaps = Inf, min.segment.length = 0,
      segment.size = 0.5, segment.color = "grey30", segment.alpha = 0.8,
      box.padding = 0.5, point.padding = 0.3,
      # Species binomials run noticeably longer than genus-only labels did
      # (e.g. "40CM-3-62-11 sp001914955") -- more repel force/iterations and
      # a wider seed spread give ggrepel more room to untangle them in the
      # already-crowded high-prevalence corner instead of stacking text.
      # water's panel is narrower in y-range (fewer, tighter-packed labels)
      # so needed the push increased further, not just left at sediment's
      # already-tuned settings, to stop stacking there specifically.
      force = 4, force_pull = 0.5, max.iter = 40000, max.time = 5, seed = 42
    ) +
    shared_phylum_scale +
    ggplot2::scale_x_continuous(
      breaks = brks,
      labels = function(x) paste0(x, "\n(", round(100 * x / n_total), "%)")
    ) +
    ggplot2::scale_y_log10(labels = scales::label_comma(), limits = y_range) +
    ggplot2::labs(x = "Prevalence (number of samples)", y = "Mean relative abundance", title = type_label) +
    (if (!show_legend) ggplot2::guides(color = "none") else NULL) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 18),
      plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
      legend.text = ggplot2::element_text(size = 15),
      legend.title = ggplot2::element_text(size = 17, face = "bold"),
      legend.key.size = ggplot2::unit(1.1, "lines")
    )
}

# Slide-ready (16:9, high dpi, larger text): ~15,457 ASVs per panel, so
# small, transparent points against overplotting; y axis log-scaled and
# shared across both panels (fixed limits) -- relative abundance spans
# several orders of magnitude, and the two panels need to stay visually
# comparable now that they're separate plots, not shared-scale facets.
# guides="collect" + legend.position="bottom" (applied to every subplot via
# `&`) puts one shared phylum legend under the whole figure instead of
# repeating it per panel.
# The larger legend.title (bumped up for readability) runs right up against
# -- and without this margin, past -- the canvas's left edge, since the
# bottom-of-figure guide box is left-aligned to the full composed width with
# no padding of its own by default.
p_prevalence_abundance <- (make_prevalence_panel("sediment") | make_prevalence_panel("water", show_legend = FALSE)) &
  ggplot2::theme(legend.position = "bottom")
p_prevalence_abundance <- p_prevalence_abundance +
  patchwork::plot_annotation(
    title = "ASV prevalence vs. mean relative abundance",
    subtitle = "Labeled: sediment ASVs in >=10/36 samples, water ASVs in >=6/10 samples (label color = phylum)",
    theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 22, face = "bold"),
                            plot.subtitle = ggplot2::element_text(size = 13),
                            plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 70))
  )
save_plot(p_prevalence_abundance, "05_asv_prevalence_abundance", w = 16, h = 9, dpi = 400)

# --- top-N barplots (remainder -> "Other") --------------------------------
# slide_quality = TRUE switches on the larger-text, 16:9-ready styling and
# the p_libsize-style plain-read-count y axis on the depth strip below;
# FALSE (default, used by the genus plot) keeps the original compact sizing.
top_n_barplot <- function(long, rank, n_top, metadata, title, slide_quality = FALSE) {
  top_taxa <- long |>
    dplyr::group_by(.data[[rank]]) |>
    dplyr::summarise(mean_rel = mean(rel_abund), .groups = "drop") |>
    dplyr::slice_max(mean_rel, n = n_top) |>
    dplyr::pull(.data[[rank]])

  # Order samples by sample_type then depth_m explicitly -- fct_reorder()'s
  # automatic (median-based) ordering chokes on duplicate sample_id rows
  # (one per taxon) once there are ties, so the level order is built by
  # hand instead.
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

  ts <- if (slide_quality) {
    list(axis_title = 18, axis_text = 14, plot_title = 20, strip = 16, legend_text = 14, legend_title = 16, depth_text = 10)
  } else {
    list(axis_title = 11, axis_text = 8, plot_title = 11, strip = 11, legend_text = 9, legend_title = 10, depth_text = 6)
  }

  p_composition <- ggplot2::ggplot(plot_df, ggplot2::aes(x = sample_id, y = rel_abund, fill = taxon)) +
    ggplot2::geom_col() +
    ggplot2::facet_grid(~sample_type, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::labs(x = NULL, y = "Relative abundance", fill = rank, title = title) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = ts$axis_text),
      axis.title = ggplot2::element_text(size = ts$axis_title),
      plot.title = ggplot2::element_text(size = ts$plot_title, face = "bold"),
      strip.text = ggplot2::element_text(size = ts$strip, face = "bold"),
      legend.text = ggplot2::element_text(size = ts$legend_text),
      legend.title = ggplot2::element_text(size = ts$legend_title)
    )

  # Read-depth strip, same x order/facets as the composition panel above --
  # a taxon at 44% of a library means something different at 300 reads than
  # at 300,000. Without this, the composition plot alone can't tell you
  # which. Log-scale bars, no fill legend (would just repeat sample_type,
  # already visible via the facet strip above). y axis matches
  # 02_qc_filter.R's p_libsize when slide_quality: plain comma-formatted
  # read counts at every power of ten, not "1e+05".
  p_depth <- tibble::tibble(sample_id = names(lib_sizes), library_size = lib_sizes) |>
    dplyr::filter(sample_id %in% levels(plot_df$sample_id)) |>
    dplyr::mutate(sample_id = factor(sample_id, levels = levels(plot_df$sample_id))) |>
    dplyr::left_join(metadata |> dplyr::select(sample_id, sample_type), by = "sample_id") |>
    ggplot2::ggplot(ggplot2::aes(x = sample_id, y = library_size)) +
    ggplot2::geom_col(fill = "grey40") +
    ggplot2::facet_grid(~sample_type, scales = "free_x", space = "free_x") +
    (if (slide_quality) {
      ggplot2::scale_y_log10(breaks = 10^(1:6), labels = scales::label_comma())
    } else {
      ggplot2::scale_y_log10()
    }) +
    ggplot2::labs(x = NULL, y = if (slide_quality) "Library size (reads)" else "Reads\n(log10)") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = ts$depth_text),
      axis.text.y = ggplot2::element_text(size = ts$depth_text),
      axis.title = ggplot2::element_text(size = ts$axis_title * 0.7),
      strip.text = ggplot2::element_blank()
    )

  p_composition / p_depth + patchwork::plot_layout(heights = c(4, 1))
}

p_phylum <- top_n_barplot(phylum_long, "phylum", 10, metadata,
                           "Phylum composition (top 10, ordered by depth within sample type)",
                           slide_quality = TRUE)
save_plot(p_phylum, "05_phylum_barplot", w = 13.333, h = 7.5, dpi = 400)

p_genus <- top_n_barplot(genus_long, "genus", 15, metadata, "Genus composition (top 15, ordered by depth within sample type)")
save_plot(p_genus, "05_genus_barplot", w = 10, h = 7)

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
