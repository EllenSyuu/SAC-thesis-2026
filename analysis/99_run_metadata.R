# Step 99 - record what produced the results in results/
source(file.path(rprojroot::find_root("DESCRIPTION"), "analysis", "_setup.R"))

git <- function(...) {
  out <- tryCatch(system2("git", c("-C", ROOT, ...), stdout = TRUE, stderr = FALSE),
                  error = function(e) NA_character_, warning = function(w) NA_character_)
  if (length(out) == 0) NA_character_ else out[1]
}
sha <- git("rev-parse", "HEAD")
dirty <- !is.na(sha) && length(suppressWarnings(
  system2("git", c("-C", ROOT, "status", "--porcelain", "--", "R", "analysis", "data/raw"),
          stdout = TRUE))) > 0

pkgs <- c("covarchimax", "extRemes", "copula", "dplyr", "tidyr", "readr",
          "ggplot2", "patchwork", "ragg", "sf", "rnaturalearth")
ver <- vapply(pkgs, function(p) {
  if (p == "covarchimax") return(read.dcf(file.path(ROOT, "DESCRIPTION"), "Version")[1])
  tryCatch(as.character(utils::packageVersion(p)), error = function(e) "not installed")
}, character(1))

raw <- list.files(RAW, full.names = TRUE)
steps <- if (file.exists(RUN_LOG)) utils::read.csv(RUN_LOG, colClasses = c(step = "character")) else data.frame()

lines <- c(
  "# Run metadata (written by analysis/99_run_metadata.R)",
  sprintf("timestamp:        %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  sprintf("git_sha:          %s%s", sha, if (dirty) "  (uncommitted changes in R/, analysis/ or data/raw/)" else ""),
  sprintf("r_version:        %s", R.version.string),
  sprintf("platform:         %s", R.version$platform),
  sprintf("cores:            %d", N_CORES),
  "",
  "## parameters",
  sprintf("seed:             %d", SEED),
  sprintf("quick:            %s", QUICK),
  sprintf("years:            %d-%d", YEARS[1], YEARS[2]),
  sprintf("period_cut:       %d", PERIOD_CUT),
  sprintf("min_days:         %d", MIN_DAYS),
  sprintf("bootstrap_B:      %d", BOOT_B),
  sprintf("omega_grid:       %g..%g by %g (%d points)", min(W_GRID), max(W_GRID),
          diff(W_GRID)[1], length(W_GRID)),
  "",
  "## input checksums (md5)",
  sprintf("%-40s %s", basename(raw), unname(tools::md5sum(raw))),
  "",
  "## package versions",
  sprintf("%-16s %s", names(ver), ver),
  "",
  "## step timings (seconds, this run)",
  if (nrow(steps)) sprintf("%-6s %-8s %8.1f  quick=%s pobs=%s min_days=%s", steps$step,
                           substr(steps$finished, 12, 19), steps$seconds,
                           steps$quick, steps$pobs, steps$min_days) else "(no step log)",
  if (nrow(steps)) sprintf("total  %17.1f", sum(steps$seconds)) else NULL
)
writeLines(lines, file.path(RES, "run_metadata.txt"))
log_step("99 run metadata written")
