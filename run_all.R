# run_all.R — sources scripts/00 -> 10 in order, with logging and a runtime
# summary, then writes results/session_info.txt. Each script is independently
# runnable (sources its own scripts/00_setup.R), so this is a thin wrapper,
# not a special orchestration mechanism.

scripts <- c(
  "00_setup.R",
  "01_import.R",
  "02_qc_filter.R",
  "02b_controls.R",
  "03_normalize.R",
  "04_alpha_diversity.R",
  "05_taxonomic_composition.R",
  "06_beta_diversity.R",
  "07_environmental_drivers.R",
  "08_differential_abundance.R",
  "09_spatial_analysis.R",
  "10_figures.R",
  "11_faprotax.R"
)

timing <- data.frame(script = character(0), seconds = numeric(0), status = character(0))

for (s in scripts) {
  path <- file.path("scripts", s)
  cat(sprintf("\n=== %s ===\n", s))
  t0 <- Sys.time()
  status <- tryCatch({
    source(path, echo = FALSE)
    "OK"
  }, error = function(e) {
    message(sprintf("[run_all] %s FAILED: %s", s, conditionMessage(e)))
    "FAILED"
  })
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  timing <- rbind(timing, data.frame(script = s, seconds = elapsed, status = status))
  cat(sprintf("--- %s: %s (%.1fs) ---\n", s, status, elapsed))
  if (status == "FAILED") {
    warning(sprintf("run_all.R stopped at %s -- fix and rerun (downstream scripts read from data/processed/ and results/, so completed steps do not need to be redone).", s))
    break
  }
}

dir.create("results", showWarnings = FALSE)
readr::write_tsv(timing, "results/run_all_timing.tsv")

writeLines(utils::capture.output(utils::sessionInfo()), "results/session_info.txt")

cat("\n=== run_all.R summary ===\n")
print(timing)
cat(sprintf("\nTotal: %.1fs\n", sum(timing$seconds)))
