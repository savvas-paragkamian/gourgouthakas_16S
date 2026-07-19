#!/usr/bin/env Rscript
#
# Build the two inputs the HiFi-16S-workflow requires:
#   data/samplesheet.tsv  — columns 'sample-id', 'filepath'
#   data/metadata.tsv     — column 'sample_name' + the study design
#
# The pipeline's inspect_metadata process diffs the sorted 'sample-id' and
# 'sample_name' columns and dies on any mismatch, so both files are generated
# from the same barcode table to keep them in lockstep.
#
# One row per *library* (51), not per site. The two technical replicates of
# each site stay separate here and are only pooled downstream in R, so that
# replicate concordance can be measured before it is averaged away.
#
# Base R only: this runs on the host before the analysis container exists.

suppressWarnings(options(stringsAsFactors = FALSE))

root      <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
fastq_dir <- file.path(root, "data", "PB482_SP")
out_dir   <- file.path(root, "data")

read_mixs <- function(path) {
  # MIxS/ENA checklist sheets carry two junk header rows above the data:
  # row 1 is the checklist banner, row 3 is a units line.
  hdr <- scan(path, what = "", sep = "\t", nlines = 1, skip = 1, quiet = TRUE)
  df  <- read.delim(path, sep = "\t", skip = 3, header = FALSE,
                    quote = "", check.names = FALSE)
  names(df) <- hdr[seq_len(ncol(df))]
  df[nzchar(trimws(df$sample_alias)), , drop = FALSE]
}

num <- function(x) suppressWarnings(as.numeric(trimws(x)))

# --- barcode -> biosample ---------------------------------------------------
bc <- read.csv(file.path(fastq_dir, "user_biosamples.csv"),
               header = TRUE, fileEncoding = "UTF-8-BOM")
names(bc) <- c("barcode", "biosample")
bc$biosample <- trimws(bc$biosample)

# --- resolve each barcode to its FASTQ --------------------------------------
# This is the reason a samplesheet exists at all: user_biosamples.csv knows the
# barcode but not the file, and the FASTQ names bury the barcode mid-string.
#
# Paths are written RELATIVE to the repo root, not absolute. The pipeline runs
# inside the container with the repo bind-mounted at /work, so a host path like
# /mnt/data/... does not resolve there. Nextflow resolves relative paths against
# the launch directory, which is /work.
fastqs <- list.files(fastq_dir, pattern = "\\.fastq\\.gz$")
bc$filepath <- vapply(bc$barcode, function(b) {
  hit <- fastqs[grepl(b, fastqs, fixed = TRUE)]
  if (length(hit) != 1L) {
    stop(sprintf("Expected exactly 1 FASTQ for barcode %s, found %d", b, length(hit)))
  }
  file.path("data", "PB482_SP", hit)
}, character(1))

# --- decompose the biosample name -------------------------------------------
# Biological libraries are "<site><rep>_A", e.g. C1_1_A, C9I_2_A, W3_1_A.
# Controls are "ctr_*" and have no replicate structure.
is_ctrl <- grepl("^ctr_", bc$biosample)

bc$replicate <- NA_integer_
bc$site      <- NA_character_

bio <- !is_ctrl
parts <- regmatches(bc$biosample[bio],
                    regexec("^([A-Z]+[0-9]+I?)_([0-9]+)_A$", bc$biosample[bio]))
stopifnot(all(lengths(parts) == 3L))

# "C9I" in the barcode sheet is written "C9_I" in the MIxS sheet — normalise.
bc$site[bio]      <- sub("^([A-Z]+[0-9]+)I$", "\\1_I", vapply(parts, `[`, "", 2))
bc$replicate[bio] <- as.integer(vapply(parts, `[`, "", 3))
bc$site[is_ctrl]  <- bc$biosample[is_ctrl]

# --- sample metadata ---------------------------------------------------------
sed <- read_mixs(file.path(root, "data", "gourgouthakas_metagenomic_samples-sediment.tsv"))
wat <- read_mixs(file.path(root, "data", "gourgouthakas_metagenomic_samples-water.tsv"))

sed_meta <- data.frame(
  site            = trimws(sed$sample_alias),
  sample_type     = "sediment",
  depth_m         = num(sed$depth),
  elevation_m     = num(sed$elevation),
  temperature_c   = num(sed$temperature),
  conductivity_ms = NA_real_,
  lat             = num(sed$`geographic location (latitude)`),
  lon             = num(sed$`geographic location (longitude)`),
  collection_date = trimws(sed$`collection date`),
  description     = trimws(sed$sample_description)
)

# NOTE: W5's conductivity (8.4) and temperature (195.9) are transposed in the
# source sheet — 195.9 degC is impossible in a cave stream and 8.4 mS/cm is far
# off the ~188 of every other water sample. Carried through verbatim here and
# flagged, rather than silently "fixed"; see CLAUDE.md.
wat_meta <- data.frame(
  site            = trimws(wat$sample_alias),
  sample_type     = "water",
  depth_m         = num(wat$depth),
  elevation_m     = num(wat$altitude),
  temperature_c   = num(wat$temperature),
  conductivity_ms = num(wat$conductivity),
  lat             = num(wat$`geographic location (latitude)`),
  lon             = num(wat$`geographic location (longitude)`),
  collection_date = trimws(wat$`collection date`),
  description     = trimws(wat$sample_description)
)

site_meta <- rbind(sed_meta, wat_meta)

meta <- merge(bc, site_meta, by = "site", all.x = TRUE, sort = FALSE)

# Controls carry no environmental data; label what kind each one is so that
# 02_qc_filter.R can route blanks to decontam and mocks to a accuracy check.
meta$sample_type[is.na(meta$sample_type)] <- "control"
meta$control_type <- NA_character_
meta$control_type[meta$biosample %in% c("ctr_zymo_com", "ctr_zymo_log", "ctr_msa_3001")] <- "mock"
meta$control_type[meta$biosample == "ctr_EB"] <- "extraction_blank"
meta$control_type[meta$biosample == "ctr_MM"] <- "mastermix_blank"

# sample_set separates the three arms of the design: the main depth transect,
# the parallel sediments that also seeded culture isolates, and the controls.
meta$sample_set <- ifelse(meta$sample_type == "control", "control",
                   ifelse(grepl("_I$", meta$site), "isolate_source", "transect"))

meta$sample_name <- meta$biosample

meta <- meta[order(meta$sample_type, meta$depth_m, meta$site, meta$replicate), ]

meta_cols <- c("sample_name", "site", "replicate", "barcode", "sample_type",
               "sample_set", "control_type", "depth_m", "elevation_m",
               "temperature_c", "conductivity_ms", "lat", "lon",
               "collection_date", "description")

# --- write -------------------------------------------------------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

samplesheet <- data.frame(`sample-id` = meta$sample_name,
                          filepath = meta$filepath,
                          check.names = FALSE)

write.table(samplesheet, file.path(out_dir, "samplesheet.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(meta[, meta_cols], file.path(out_dir, "metadata.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# --- the same check inspect_metadata will run --------------------------------
stopifnot(identical(sort(samplesheet$`sample-id`), sort(meta$sample_name)))
stopifnot(!anyDuplicated(meta$sample_name))
stopifnot(all(file.exists(file.path(root, samplesheet$filepath))))

cat(sprintf("samplesheet.tsv: %d libraries\n", nrow(samplesheet)))
cat(sprintf("metadata.tsv:    %d rows, %d columns\n", nrow(meta), length(meta_cols)))
print(table(meta$sample_type, meta$sample_set))
