#!/usr/bin/env Rscript

# 99_session_info.R
# Run this script in the environment used for the final analysis and commit
# the resulting text file with the tagged release.

dir.create("results", showWarnings = FALSE, recursive = TRUE)

info <- capture.output({
  cat("R version:\n")
  print(R.version.string)
  cat("\nSession information:\n")
  print(sessionInfo())
})

writeLines(info, "results/sessionInfo.txt")
cat(paste(info, collapse = "\n"))
