# Post-process the overnight log file into a one-page summary.
# Run as: Rscript inst/scripts/overnight_summarize.R [log_file]

args <- commandArgs(trailingOnly = TRUE)
log_files <- if (length(args) > 0L) args[1L] else {
  reports_dir <- "C:/Users/ammonsj/DID GPU/didgpu/inst/overnight_reports"
  fs <- list.files(reports_dir, pattern = "^overnight_.*\\.log$",
                    full.names = TRUE)
  if (length(fs) == 0L) {
    cat("No overnight logs found in ", reports_dir, "\n", sep = "")
    quit(status = 1L)
  }
  fs <- fs[order(file.mtime(fs), decreasing = TRUE)]
  fs[1L]
}

cat("=========================================================\n")
cat(" Overnight summary: ", basename(log_files), "\n", sep = "")
cat("=========================================================\n")

lines <- readLines(log_files, warn = FALSE)

# Extract test-result lines (typically "test-name: ..............")
test_lines <- grep("^[a-z][a-zA-Z0-9_-]+:\\s+[.WSF1-9]+$", lines, value = TRUE)
fail_lines <- grep("Failed|FAILED|^── ", lines, value = TRUE)
skip_lines <- grep("Skipped", lines, value = TRUE)

cat("\n## Test-file results\n\n")
if (length(test_lines) > 0L) {
  for (l in test_lines) cat("  ", l, "\n", sep = "")
} else {
  cat("  (no test-file lines detected)\n")
}

cat("\n## Failures and warnings\n\n")
if (length(fail_lines) > 0L) {
  for (l in head(fail_lines, 30L)) cat("  ", l, "\n", sep = "")
  if (length(fail_lines) > 30L) cat("  ... ", length(fail_lines) - 30L,
                                       " more lines suppressed\n", sep = "")
} else {
  cat("  ALL CLEAN. No failures or test-warning headers in the log.\n")
}

cat("\n## Skips\n\n")
if (length(skip_lines) > 0L) {
  for (l in head(skip_lines, 10L)) cat("  ", l, "\n", sep = "")
} else {
  cat("  (none)\n")
}

# Pull elapsed timing if present.
cat("\n## Timing\n\n")
elapsed_lines <- grep("elapsed:|finished:|started:|TOTAL", lines, value = TRUE)
for (l in elapsed_lines) cat("  ", l, "\n", sep = "")

cat("\n=========================================================\n")
cat(" Full log at: ", log_files, "\n", sep = "")
cat("=========================================================\n")
