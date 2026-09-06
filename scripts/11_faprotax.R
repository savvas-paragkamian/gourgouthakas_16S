# 11_faprotax.R — functional-guild prediction from ASV taxonomy (FAPROTAX,
# Louca et al. 2016), the same tool isd_archive_scripts/isd_crete_workflow.sh
# used (lines 98-111) against the same idea (a taxonomy-string community
# matrix -> collapse_table.py -> a functional-group x sample table). New to
# this project's own 00-10 pipeline; see AGENTS.md's Downstream analysis
# notes for the naming-mismatch caveat this script measures rather than
# assumes.
#
# FAPROTAX's database is written against pre-GTDB (SILVA/NCBI-style) names
# (Proteobacteria, Firmicutes, ...); this project's taxonomy is GTDB-style
# (Pseudomonadota, Bacillota, ...). Higher-rank (phylum/class) group
# definitions will systematically under-match; genus/species-level
# definitions (the majority of specific functional guilds) should still
# resolve, since standard genus epithets are shared across naming systems.
# The unassigned-read fraction reported below is the direct, measured answer
# to how much this actually costs -- not assumed away.

source("scripts/00_setup.R")

counts <- readRDS(file.path(path_processed, "counts_clean.rds"))
taxonomy <- readRDS(file.path(path_processed, "taxonomy_clean.rds"))
metadata <- readRDS(file.path(path_processed, "metadata_clean.rds"))
lib_sizes <- rowSums(counts)

path_faprotax_root <- "tools/faprotax"
path_faprotax_dist <- file.path(path_faprotax_root, "FAPROTAX_1.2.12")
path_faprotax_script <- file.path(path_faprotax_dist, "collapse_table.py")
path_faprotax_db <- file.path(path_faprotax_dist, "FAPROTAX.txt")
path_faprotax_venv_py <- file.path(path_faprotax_root, "venv", "bin", "python3")

# --- one-time, idempotent setup: download FAPROTAX + build a local venv ----
# with numpy (the container's system python3 has neither numpy nor pip; a
# venv + ensurepip needs no root/apt). Skipped entirely if already present
# (tools/faprotax/ is .gitignore'd -- vendored tool, not source/data).
ensure_faprotax_tooling <- function() {
  if (!file.exists(path_faprotax_script) || !file.exists(path_faprotax_db)) {
    dir.create(path_faprotax_root, showWarnings = FALSE, recursive = TRUE)
    zip_path <- file.path(path_faprotax_root, "faprotax.zip")
    url <- paste0(
      "http://pages.uoregon.edu/slouca/LoucaLab/archive/FAPROTAX/",
      "SECTION_Download/MODULE_Downloads/CLASS_Latest%20release/",
      "UNIT_FAPROTAX_1.2.12/FAPROTAX_1.2.12.zip"
    )
    message("[11_faprotax] downloading FAPROTAX 1.2.12 from the Louca Lab archive...")
    status <- utils::download.file(url, zip_path, quiet = TRUE, mode = "wb")
    if (status != 0 || !file.exists(zip_path) || file.size(zip_path) < 1e5) {
      stop("[11_faprotax] FAPROTAX download failed or returned an unexpectedly small file -- ",
           "check network access / whether the Louca Lab archive URL has moved.", call. = FALSE)
    }
    utils::unzip(zip_path, exdir = path_faprotax_root)
    if (!file.exists(path_faprotax_script)) {
      stop("[11_faprotax] FAPROTAX zip did not contain the expected collapse_table.py.", call. = FALSE)
    }
  }
  if (!file.exists(path_faprotax_venv_py)) {
    message("[11_faprotax] building a local Python venv with numpy for collapse_table.py...")
    venv_status <- system2("python3", c("-m", "venv", file.path(path_faprotax_root, "venv")))
    if (venv_status != 0) stop("[11_faprotax] python3 -m venv failed.", call. = FALSE)
    pip_status <- system2(file.path(path_faprotax_root, "venv", "bin", "pip"),
                           c("install", "--quiet", "numpy"))
    if (pip_status != 0) stop("[11_faprotax] pip install numpy (into the local venv) failed -- check network access.", call. = FALSE)
  }
  invisible(TRUE)
}
ensure_faprotax_tooling()

# --- build the FAPROTAX input table -----------------------------------------
# ASVs as rows, one column per sample (raw counts) + a trailing taxonomy
# column -- the orientation collapse_table.py expects, same transpose
# 08_differential_abundance.R already does for ALDEx2. genus_norm/species_norm
# (not the raw GTDB-suffixed genus/species) are used for the last two ranks:
# they're already stripped of polyphyly-split suffixes (02b_controls.R's T1
# step), so a match that would otherwise be blocked by e.g. "_A" isn't.
faprotax_taxonomy_string <- function(tax) {
  sprintf(
    "d__%s;p__%s;c__%s;o__%s;f__%s;g__%s;s__%s",
    tax$domain, tax$phylum, tax$class, tax$order, tax$family,
    ifelse(is.na(tax$genus_norm), tax$genus, tax$genus_norm),
    ifelse(is.na(tax$species_norm), tax$species, tax$species_norm)
  )
}

faprotax_input <- taxonomy |>
  dplyr::transmute(asv_id, taxonomy_string = faprotax_taxonomy_string(taxonomy)) |>
  dplyr::left_join(
    t(counts) |> as.data.frame() |> tibble::rownames_to_column("asv_id"),
    by = "asv_id"
  ) |>
  dplyr::relocate(taxonomy_string, .after = dplyr::last_col()) |>
  dplyr::rename(taxonomy = taxonomy_string)

path_input <- "results/11_faprotax_input.tsv"
readr::write_tsv(faprotax_input, path_input)

# --- run the reference implementation ---------------------------------------
# -n none: raw summed counts per group, not FAPROTAX's own normalization --
# this project computes its own relative abundance from library_size
# everywhere else (05_taxonomic_composition.R's aggregate_taxon()), and
# should stay consistent rather than trust two different normalization
# schemes for the same read.
path_functional_table <- "results/11_faprotax_functional_table.tsv"
path_report <- "results/11_faprotax_report.txt"
path_subtables <- "results/11_faprotax_subtables"
dir.create(path_subtables, showWarnings = FALSE, recursive = TRUE)

collapse_args <- c(
  path_faprotax_script,
  "-i", path_input,
  "-o", path_functional_table,
  "-g", path_faprotax_db,
  "-d", "taxonomy",
  "--omit_columns", "0",
  "-n", "none",
  "-r", path_report,
  "--collapse_by_metadata", "taxonomy",
  "--out_sub_tables_dir", path_subtables,
  "--force" # this script must be safely re-runnable, like every other script here
)
collapse_status <- system2(path_faprotax_venv_py, collapse_args, stdout = TRUE, stderr = TRUE)
collapse_exit <- attr(collapse_status, "status")
if (!is.null(collapse_exit) && collapse_exit != 0) {
  stop("[11_faprotax] collapse_table.py failed (exit ", collapse_exit, "):\n",
       paste(utils::tail(collapse_status, 20), collapse = "\n"), call. = FALSE)
}
if (!file.exists(path_functional_table)) {
  stop("[11_faprotax] collapse_table.py ran but did not produce ", path_functional_table, call. = FALSE)
}

# --- post-process: this project's own relative-abundance convention --------
functional_wide <- readr::read_tsv(path_functional_table, show_col_types = FALSE)

functional_long <- functional_wide |>
  tidyr::pivot_longer(-group, names_to = "sample_id", values_to = "reads") |>
  dplyr::left_join(tibble::tibble(sample_id = names(lib_sizes), library_size = lib_sizes), by = "sample_id") |>
  dplyr::mutate(rel_abund = reads / library_size) |>
  dplyr::left_join(metadata |> dplyr::select(sample_id, sample_type, depth_m), by = "sample_id")

write_result(functional_long, "11_faprotax_relabund")

# Unassigned fraction: collapse_table.py's report.txt states an aggregate
# ASV/record-level figure ("N out of M records (.. %) could not be assigned
# to any group") but doesn't list which ASVs those are, and writes no
# leftover pseudo-row into the functional table itself -- confirmed by
# inspecting the actual report.txt from this run, not assumed from the
# --help text. The read-weighted, per-sample number this project actually
# wants is computed independently here: every group's sub-table
# (--out_sub_tables_dir) lists that group's member taxonomy strings as row
# names, so the union across all of them is the "assigned" set; anything in
# faprotax_input not in that union is unassigned. This is the direct,
# measured answer to the GTDB-naming-mismatch caveat above, not an assumed
# one.
assigned_taxonomy_strings <- list.files(path_subtables, pattern = "\\.tsv$", full.names = TRUE) |>
  purrr::map(~ readr::read_tsv(.x, skip = 1, show_col_types = FALSE)$record) |>
  unlist() |>
  unique()

is_assigned <- faprotax_input$taxonomy %in% assigned_taxonomy_strings
unassigned_asv_ids <- faprotax_input$asv_id[!is_assigned]

unassigned_by_sample <- tibble::tibble(
  sample_id = names(lib_sizes),
  library_size = lib_sizes,
  unassigned_reads = rowSums(counts[, unassigned_asv_ids, drop = FALSE])
) |>
  dplyr::mutate(unassigned_frac = unassigned_reads / library_size)
write_result(unassigned_by_sample, "11_faprotax_unassigned")

# Corroborating ASV/record-level figure straight from collapse_table.py's
# own report (line "N out of M records (.. %) could not be assigned to any
# group (leftovers)") -- a different denominator (ASVs, not reads) than the
# read-weighted number above, both reported rather than picking one.
report_lines <- readLines(path_report)
report_leftover_line <- grep("could not be assigned to any group", report_lines, value = TRUE)

mean_unassigned_frac <- mean(unassigned_by_sample$unassigned_frac, na.rm = TRUE)
message(sprintf(
  "[11_faprotax] %d functional groups x %d samples; mean read-weighted unassigned fraction = %.1f%%; %s",
  dplyr::n_distinct(functional_long$group),
  dplyr::n_distinct(functional_long$sample_id),
  100 * mean_unassigned_frac,
  if (length(report_leftover_line) == 1) trimws(sub("^#\\s*", "", report_leftover_line)) else "ASV-level leftover line not found in report.txt"
))

# --- sanity check: no group's raw read count should exceed the sample's ----
# total library size (a group can double-count an ASV belonging to several
# groups across DIFFERENT groups, but never within one group x sample cell)
over_capacity <- functional_long |>
  dplyr::filter(reads > library_size)
if (nrow(over_capacity) > 0) {
  stop(sprintf("[11_faprotax] %d group x sample cell(s) exceed that sample's library size -- collapse_table.py output or the join above is broken.",
               nrow(over_capacity)), call. = FALSE)
}

# --- figures -----------------------------------------------------------
# Top-N functional groups by mean relative abundance, same style as
# 05_taxonomic_composition.R's top_n_barplot(): samples ordered by
# sample_type then depth_m, faceted by sample_type.
top_groups <- functional_long |>
  dplyr::group_by(group) |>
  dplyr::summarise(mean_rel = mean(rel_abund), .groups = "drop") |>
  dplyr::slice_max(mean_rel, n = 15) |>
  dplyr::pull(group)

sample_order <- metadata |> dplyr::arrange(sample_type, depth_m) |> dplyr::pull(sample_id)

heatmap_df <- functional_long |>
  dplyr::filter(group %in% top_groups) |>
  dplyr::mutate(
    group = factor(group, levels = rev(top_groups)),
    sample_id = factor(sample_id, levels = intersect(sample_order, sample_id))
  )

p_heatmap <- ggplot2::ggplot(heatmap_df, ggplot2::aes(x = sample_id, y = group, fill = rel_abund)) +
  ggplot2::geom_tile() +
  ggplot2::facet_grid(~sample_type, scales = "free_x", space = "free_x") +
  scale_fill_depth(name = "Relative\nabundance", trans = "sqrt") +
  ggplot2::labs(x = NULL, y = NULL, title = "FAPROTAX functional groups (top 15 by mean relative abundance)") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
save_plot(p_heatmap, "11_faprotax_heatmap", w = 10, h = 6)

# Depth-transect panel for a short list of ecologically salient groups given
# this is a cave sediment/water system -- filtered to whichever of these
# candidates actually have nonzero support in this dataset, checked rather
# than assumed present.
candidate_groups <- c(
  "nitrification", "aerobic_ammonia_oxidation", "aerobic_nitrite_oxidation",
  "sulfate_respiration", "sulfur_respiration", "dark_sulfite_oxidation",
  "methanotrophy", "methanogenesis", "fermentation", "aerobic_chemoheterotrophy",
  "denitrification", "nitrogen_fixation"
)
present_groups <- functional_long |>
  dplyr::filter(group %in% candidate_groups) |>
  dplyr::group_by(group) |>
  dplyr::filter(sum(reads) > 0) |>
  dplyr::pull(group) |>
  unique()

if (length(present_groups) > 0) {
  p_depth <- functional_long |>
    dplyr::filter(group %in% present_groups, !is.na(depth_m)) |>
    ggplot2::ggplot(ggplot2::aes(depth_m, rel_abund, color = sample_type)) +
    ggplot2::geom_point(size = 1.5, alpha = 0.7) +
    ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = TRUE) +
    ggplot2::scale_color_manual(values = palette_sample_type()) +
    ggplot2::facet_wrap(~group, scales = "free_y") +
    ggplot2::labs(x = "Depth (m)", y = "Relative abundance", color = NULL,
                  title = "Ecologically salient FAPROTAX groups vs. depth")
  save_plot(p_depth, "11_faprotax_depth_trends", w = 10, h = 7)
} else {
  message("[11_faprotax] none of the candidate ecologically-salient groups (nitrification/sulfur/methane/fermentation cycling) had any nonzero support in this dataset -- no depth-trend plot written.")
}

message(sprintf(
  "[11_faprotax] wrote 11_faprotax_relabund.tsv, 11_faprotax_unassigned.tsv, 11_faprotax_heatmap.{pdf,png}%s",
  if (length(present_groups) > 0) ", 11_faprotax_depth_trends.{pdf,png}" else ""
))
