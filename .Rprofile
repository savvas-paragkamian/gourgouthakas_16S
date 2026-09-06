# renv's autoloader is disabled here — see AGENTS.md "environment traps" (#4)
# for why: its Bioconductor bootstrap step fails outright in this image (it
# looks for BiocManager in a per-project/user cache that was never actually
# populated for this container's HOME), which would otherwise turn every
# plain `Rscript` invocation into a crash before any script code runs.
# renv.lock is a version *manifest* here, not a live package source: all 357
# pinned versions already match what's installed in the system library
# (/usr/lib64/R/library), which is what every script actually loads from.
Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")
source("renv/activate.R")
