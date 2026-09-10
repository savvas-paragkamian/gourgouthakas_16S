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
  dplyr::left_join(metadata |> dplyr::select(sample_id, sample_type), by = "sample_id") |>
  dplyr::arrange(library_size)
write_result(rarefaction_rationale, "03_rarefaction_depth_rationale")

# --- T8: audit the upstream depth, report only, don't act on it ------------
# HiFi-16S-workflow's final_stats.R computes rarefaction_depth_suggested.txt
# by sorting EVERY sample's final read count -- no sample_type awareness --
# and picking the value at rank floor(0.8*N). Controls (including ctr_EB=41,
# ctr_MM=70 reads) were part of that sort when 499 was generated, so the
# chosen rank's value is pulled down by the near-empty blanks. This is a
# stop-and-ask finding, not a unilateral fix: AGENTS.md records reusing the
# pipeline's own value as a deliberate choice, and re-deriving it changes
# which samples get excluded from the rarefied matrix -- report the
# recomputed alternative for a decision, keep using 499 either way.
recomputed_depth <- {
  sorted <- sort(lib_sizes, decreasing = TRUE)
  n <- length(sorted)
  rank <- max(1L, as.integer(n * 0.8))
  as.integer(floor(sorted[rank]))
}
write_result(
  tibble::tibble(
    value = c("upstream_rarefaction_depth_suggested", "recomputed_on_46_real_samples"),
    depth = c(rarefaction_depth, recomputed_depth),
    note = c("computed by HiFi-16S-workflow/bin/final_stats.R over all 51 libraries including controls",
             "same rank-floor(0.8*N) algorithm, recomputed here over the 46 control-free libraries -- NOT applied, reported for a decision")
  ),
  "T8_rarefaction_depth_audit"
)
message(sprintf(
  "[03_normalize] T8 audit: upstream depth = %d (computed over 51 libraries incl. controls); recomputed over 46 real libraries = %d. Still using %d (stop-and-ask, not changed here).",
  rarefaction_depth, recomputed_depth, rarefaction_depth
))

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
# The full post-QC ASV table is far too sparse for compositional
# zero-replacement at full resolution: cmultRepl's default
# z.warning/z.delete=0.8/TRUE *silently dropped samples* the first time this
# was run here (33/51, back when controls were still in counts_clean), since
# most ASVs are zero in >80% of samples at that resolution. Two independent
# fixes, both kept:
#   1. A CLR-specific prevalence filter, coarser than 02b_controls.R's
#      general >=2-sample floor -- >=10% of samples. Re-checked after
#      controls moved out of counts_clean (PLAN.md's filter-order fix):
#      the *same* 10% threshold now recovers ~1,200 ASVs on 46 samples,
#      up from 717 on the old 51-sample (control-inflated-sparsity) table --
#      controls, especially the log-distributed mock's ~578 low-abundance
#      ASVs, were dragging down apparent prevalence for real-sample ASVs
#      too. A more aggressive (lower) threshold was tried (5% -> ~5,400
#      ASVs) but not adopted: cmultRepl's runtime scales with ASV count,
#      and there's no clear analytical benefit at that resolution to justify
#      the extra compute here -- 10% stays the default, re-derived rather
#      than re-affirmed blindly, per PLAN.md.
#   2. z.delete = FALSE explicitly, so if a handful of genuinely near-empty
#      samples (the shallowest sediment/water libraries) still exceed the
#      zero-warning threshold even after that, cmultRepl *warns* rather than
#      *silently removing them* -- every sample present in counts_clean
#      stays present in counts_clr, full stop.
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

# ============================================================================
# Genus-level equivalents (Part 3c of golden-napping-breeze.md) -- same
# three transforms, same rarefaction_depth (499, from the upstream
# workflow), same re-derived CLR prevalence LOGIC (>=10% of samples) but
# re-applied to the genus matrix's own column count and sparsity, not
# reusing the ASV-level threshold's absolute number. Genus-level becomes
# primary for 04/06/07; ASV-level (above) stays supplementary throughout.
# ============================================================================
counts_genus <- readRDS(file.path(path_processed, "counts_genus_clean.rds"))

# below_depth is recomputed from counts_genus's OWN row sums, not reused
# from the ASV-level `below_depth` above -- a sample's genus-level total is
# lower than its ASV-level total (genus-NA reads excluded, ~4.6% overall,
# but not uniformly per sample), so the same 499-read depth can exclude a
# different set of samples here.
lib_sizes_genus <- rowSums(counts_genus)
below_depth_genus <- names(lib_sizes_genus)[lib_sizes_genus < rarefaction_depth]
message(sprintf(
  "[03_normalize] genus-level: %d/%d samples below the rarefaction depth (%d), excluded from counts_genus_rarefied: %s",
  length(below_depth_genus), nrow(counts_genus), rarefaction_depth, paste(below_depth_genus, collapse = ", ")
))

counts_genus_for_rarefaction <- counts_genus[!(rownames(counts_genus) %in% below_depth_genus), , drop = FALSE]
counts_genus_rarefied <- vegan::rrarefy(counts_genus_for_rarefaction, sample = rarefaction_depth)

counts_genus_relabund <- vegan::decostand(counts_genus, method = "total")

clr_prevalence_min_samples_genus <- ceiling(0.10 * nrow(counts_genus))
clr_genus_keep <- colSums(counts_genus > 0) >= clr_prevalence_min_samples_genus
counts_genus_for_clr <- counts_genus[, clr_genus_keep, drop = FALSE]
message(sprintf(
  "[03_normalize] genus-level CLR-specific prevalence filter (>=%d/%d samples): %d/%d genera kept for CLR/Aitchison",
  clr_prevalence_min_samples_genus, nrow(counts_genus), sum(clr_genus_keep), ncol(counts_genus)
))

counts_genus_zero_replaced <- zCompositions::cmultRepl(counts_genus_for_clr, method = "CZM", output = "p-counts",
                                                        z.delete = FALSE, suppress.print = TRUE)
stopifnot(setequal(rownames(counts_genus_zero_replaced), rownames(counts_genus)))

counts_genus_clr <- compositions::clr(counts_genus_zero_replaced)
counts_genus_clr <- matrix(as.numeric(counts_genus_clr), nrow = nrow(counts_genus_clr), ncol = ncol(counts_genus_clr),
                            dimnames = dimnames(counts_genus_clr))

saveRDS(counts_genus_rarefied, file.path(path_processed, "counts_genus_rarefied.rds"))
saveRDS(counts_genus_relabund, file.path(path_processed, "counts_genus_relabund.rds"))
saveRDS(counts_genus_clr, file.path(path_processed, "counts_genus_clr.rds"))
saveRDS(below_depth_genus, file.path(path_processed, "rarefaction_excluded_samples_genus.rds"))

message(sprintf(
  "[03_normalize] genus-level rarefied: %d x %d | relabund: %d x %d | clr: %d x %d",
  nrow(counts_genus_rarefied), ncol(counts_genus_rarefied),
  nrow(counts_genus_relabund), ncol(counts_genus_relabund),
  nrow(counts_genus_clr), ncol(counts_genus_clr)
))
