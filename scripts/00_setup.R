# 00_setup.R — libraries, seed, paths, theme, palettes.
#
# Every script in this pipeline starts with:
#   source("scripts/00_setup.R")
#
# renv activates from the project root (.Rprofile), so this assumes R was
# started at the repo root (/work), not from inside scripts/.

# --- fail loudly on any missing package -------------------------------------
required_pkgs <- c(
  # tidyverse (loaded as a unit below, listed here so the check is uniform)
  "tidyverse",
  # core ecology / stats
  "vegan", "ape",
  # compositional normalization
  "compositions", "zCompositions",
  # contamination removal
  "decontam",
  # differential abundance
  "ALDEx2", "Maaslin2",
  # flat-file I/O (qiime-free)
  "biomformat", "Biostrings",
  # figures
  "patchwork", "ggpubr", "ggrepel", "scales", "here", "gridExtra", "png", "pheatmap",
  # collinearity screening (§07)
  "car"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required package(s): ", paste(missing_pkgs, collapse = ", "),
    ". Run renv::restore() (inside the container, so renv activates from ",
    "the project root) before continuing.",
    call. = FALSE
  )
}

# picante (Faith's PD) and GUniFrac (UniFrac) are deliberately NOT required:
# this dataset has no phylogeny (AGENTS.md departure #1). sf/terra are
# likewise not required: there is no spatial gradient to map, only depth
# (AGENTS.md departure #2). Do not add them back without re-reading AGENTS.md.

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ape)
})

source("scripts/functions.R")

set.seed(42) # rarefaction, NMDS, permutations, ALDEx2 MC sampling

# --- paths -------------------------------------------------------------

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("plots", showWarnings = FALSE, recursive = TRUE)

path_hifi_final <- "results/hifi/final"
path_asv_table <- file.path(path_hifi_final, "best_tax_merged_freq_tax.tsv")
path_tracking <- file.path(path_hifi_final, "dada2_tracking.tsv")
path_metadata <- "data/metadata.tsv"

path_processed <- "data/processed"

# --- theme + palettes ---------------------------------------------------

ggplot2::theme_set(theme_soil())

# --- session marker so 09/etc. can tell they're running standalone ---------
if (!exists(".soil_setup_done")) {
  .soil_setup_done <- TRUE
  message(sprintf("[00_setup] R %s, seed 42, wd = %s", getRversion(), getwd()))
}
