# 03_normalize.R — rarefaction curves + three normalized matrices
# (rarefied / relative abundance / CLR), each used downstream for a
# specific purpose per PLAN.md §7.

source("scripts/00_setup.R")

counts <- readRDS(file.path(path_processed, "counts_clean.rds"))
metadata <- readRDS(file.path(path_processed, "metadata_clean.rds"))

# --- rarefaction curves ----------------------------------------------------
grDevices::pdf(file.path("plots", "03_rarefaction_curves.pdf"), width = 8, height = 6)
vegan::rarecurve(counts, step = 50, label = FALSE,
                  col = ifelse(metadata$sample_type[match(rownames(counts), metadata$sample_id)] == "sediment",
                                palette_sample_type()["sediment"], palette_sample_type()["water"]))
graphics::title(main = "Rarefaction curves (post-QC ASV table)")
grDevices::dev.off()
grDevices::png(file.path("plots", "03_rarefaction_curves.png"), width = 8, height = 6, units = "in", res = 300)
vegan::rarecurve(counts, step = 50, label = FALSE,
                  col = ifelse(metadata$sample_type[match(rownames(counts), metadata$sample_id)] == "sediment",
                                palette_sample_type()["sediment"], palette_sample_type()["water"]))
graphics::title(main = "Rarefaction curves (post-QC ASV table)")
grDevices::dev.off()

# --- rarefaction depth: use the upstream pipeline's own recommendation -----
# HiFi-16S-workflow already computed this (results/hifi/final/
# rarefaction_depth_suggested.txt) -- reuse it rather than re-deriving an
# arbitrary quantile, since it comes with its own methodology from the
# workflow. (alpha_depth_suggested.txt, 116877, is a much higher *saturation*
# depth -- how deep is "deep enough" for the curve to plateau -- not a
# depth to subsample every sample to; using it here would drop all but a
# handful of the highest-biomass libraries. rarefaction_depth_suggested.txt
# is the one meant for standardizing samples.)
rarefaction_depth <- as.integer(readLines(file.path(path_hifi_final, "rarefaction_depth_suggested.txt")))
message(sprintf("[03_normalize] rarefaction depth = %d (from HiFi-16S-workflow's rarefaction_depth_suggested.txt)", rarefaction_depth))

lib_sizes <- rowSums(counts)
below_depth <- names(lib_sizes)[lib_sizes < rarefaction_depth]

rarefaction_rationale <- tibble::tibble(
  sample_id = names(lib_sizes),
  library_size = lib_sizes,
  rarefaction_depth = rarefaction_depth,
  excluded_from_rarefied = lib_sizes < rarefaction_depth
) |>
  dplyr::left_join(metadata |> dplyr::select(sample_id, sample_type, control_type), by = "sample_id") |>
  dplyr::arrange(library_size)
write_result(rarefaction_rationale, "03_rarefaction_depth_rationale")

message(sprintf(
  "[03_normalize] %d/%d samples below the rarefaction depth, excluded from the rarefied matrix (alpha diversity only -- still present in relative-abundance and CLR matrices): %s",
  length(below_depth), nrow(counts), paste(below_depth, collapse = ", ")
))

counts_for_rarefaction <- counts[!(rownames(counts) %in% below_depth), , drop = FALSE]
counts_rarefied <- vegan::rrarefy(counts_for_rarefaction, sample = rarefaction_depth)

# --- relative abundance (all samples; used for composition + unconstrained
# ordination, where a shallow sample is noisier but not invalid) -----------
counts_relabund <- vegan::decostand(counts, method = "total")

# --- CLR (all samples; zero-treatment via zCompositions, then compositions::clr) ---
#
# The full post-QC ASV table (15,460 ASVs, prevalence >=2 samples only) is
# far too sparse for compositional zero-replacement at full resolution:
# cmultRepl's default z.warning/z.delete=0.8/TRUE *silently dropped 33/51
# samples* the first time this was run here, because most ASVs are zero in
# >80% of samples at that resolution. Two independent fixes, both kept:
#   1. A CLR-specific prevalence filter, coarser than 02_qc_filter's general
#      >=2-sample floor -- >=10% of samples (>=6/51) -- brings sparsity down
#      to something compositional replacement can actually work with (717
#      ASVs) while staying far more permissive than the ~79 ASVs a 20%
#      threshold would give.
#   2. z.delete = FALSE explicitly, so if a handful of genuinely near-empty
#      samples (the mocks by design, plus the shallowest sediment libraries)
#      still exceed the zero-warning threshold even after that, cmultRepl
#      *warns* rather than *silently removing them* -- every sample present
#      in counts_clean stays present in counts_clr, full stop.
clr_prevalence_min_samples <- ceiling(0.10 * nrow(counts))
clr_asv_keep <- colSums(counts > 0) >= clr_prevalence_min_samples
counts_for_clr <- counts[, clr_asv_keep, drop = FALSE]

message(sprintf(
  "[03_normalize] CLR-specific prevalence filter (>=%d/%d samples): %d/%d ASVs kept for CLR/Aitchison",
  clr_prevalence_min_samples, nrow(counts), sum(clr_asv_keep), ncol(counts)
))

counts_zero_replaced <- zCompositions::cmultRepl(counts_for_clr, method = "CZM", output = "p-counts",
                                                  z.delete = FALSE, suppress.print = TRUE)
stopifnot(setequal(rownames(counts_zero_replaced), rownames(counts))) # no sample silently dropped

counts_clr <- compositions::clr(counts_zero_replaced)
counts_clr <- matrix(as.numeric(counts_clr), nrow = nrow(counts_clr), ncol = ncol(counts_clr),
                      dimnames = dimnames(counts_clr))

# --- save ------------------------------------------------------------------
saveRDS(counts_rarefied, file.path(path_processed, "counts_rarefied.rds"))
saveRDS(counts_relabund, file.path(path_processed, "counts_relabund.rds"))
saveRDS(counts_clr, file.path(path_processed, "counts_clr.rds"))
saveRDS(below_depth, file.path(path_processed, "rarefaction_excluded_samples.rds"))

message(sprintf(
  "[03_normalize] rarefied: %d x %d | relabund: %d x %d | clr: %d x %d",
  nrow(counts_rarefied), ncol(counts_rarefied),
  nrow(counts_relabund), ncol(counts_relabund),
  nrow(counts_clr), ncol(counts_clr)
))
