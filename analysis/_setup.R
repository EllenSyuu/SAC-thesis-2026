# Shared configuration for the analysis pipeline. Sourced by every step.
# A UTF-8 locale is needed for accented station names in figures.
if (!l10n_info()$`UTF-8`) {
  for (loc in c("C.UTF-8", "en_US.UTF-8")) {
    if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break
  }
}
suppressPackageStartupMessages({
  pkgload::load_all(rprojroot::find_root("DESCRIPTION"), quiet = TRUE,
                    export_all = FALSE)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# QUICK=1 shrinks the Monte-Carlo study, bootstrap and omega-grid for a fast
# smoke run and writes to results/quick/ so the real results are untouched.
QUICK   <- identical(Sys.getenv("QUICK"), "1")

ROOT    <- rprojroot::find_root("DESCRIPTION")
RAW     <- file.path(ROOT, "data", "raw")
PROC    <- file.path(ROOT, "data", "processed", if (QUICK) "quick" else "")
RES     <- file.path(ROOT, "results", if (QUICK) "quick" else "")
FIG     <- file.path(RES, "figures")
TAB     <- file.path(RES, "tables")
for (d in c(PROC, FIG, TAB)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
SEED    <- 20260408
N_CORES <- max(1L, parallel::detectCores() - 0L)

MIN_DAYS     <- as.integer(Sys.getenv("MIN_DAYS", "0"))  # 0 = thesis (no threshold)
BOOT_B       <- as.integer(Sys.getenv("BOOT_B", if (QUICK) "20" else "200"))
YEARS        <- c(1930, 2017)   # 2018-19 dropped: ECCC manual->automatic gap
PERIOD_CUT   <- 1975
W_GRID       <- if (QUICK) seq(0.05, 0.95, by = 0.05) else seq(0.01, 0.99, by = 0.01)

# Station order defines u1..u6 throughout.
STATIONS <- data.frame(
  key   = c("DRUMMONDVILLE", "L'ASSOMPTION", "LES CEDRES",
            "ST CHARLES DE MANDEVILLE", "LENNOXVILLE", "NICOLET"),
  label = c("Drummondville", "L'Assomption", "Les Cèdres",
            "St-Charles-de-Mandeville", "Lennoxville", "Nicolet"),
  short = c("Drum.", "L'Ass.", "Les Cèd.", "St-Ch.", "Lennox.", "Nicolet"),
  abbr  = c("D", "A", "C", "S", "L", "N"),
  id    = c(5348, 5237, 5239, 5262, 5397, 5426),
  lat   = c(45.88, 45.80, 45.30, 46.40, 45.37, 46.22),
  lon   = c(-72.48, -73.42, -74.05, -73.35, -71.82, -72.62),
  elev  = c(82.3, 21.0, 47.2, 168.0, 181.0, 30.4),
  stringsAsFactors = FALSE
)
U_COLS <- paste0("u", seq_len(nrow(STATIONS)))

SEASONS <- c("AM", "JJ", "AS", "ON")
PAL_SEASON <- c(AM = "#0072B2", JJ = "#E69F00", AS = "#009E73", ON = "#CC79A7")
PAL_PERIOD <- c(pre = "#2166AC", post = "#B2182B")
PAL_METHOD <- c(Pickands = "#999999", CFG = "#E69F00", MLE = "#0072B2",
                `MLE (link)` = "#D55E00", True = "black")

theme_set(
  theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "grey80", fill = NA),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold"),
          plot.title.position = "plot",
          legend.position = "bottom")
)

save_fig <- function(p, name, width = 8, height = 5) {
  dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = width,
         height = height, dpi = 200, bg = "white", device = dev)
  invisible(p)
}

save_tab <- function(df, name, digits = 4) {
  # round only genuinely fractional columns; counts stay exact
  frac <- vapply(df, function(x) is.numeric(x) && any(abs(x - round(x)) > 1e-9, na.rm = TRUE),
                 logical(1))
  df[frac] <- lapply(df[frac], signif, digits = digits)
  readr::write_csv(df, file.path(TAB, paste0(name, ".csv")))
  invisible(df)
}

pair_index <- function(d = nrow(STATIONS)) {
  idx <- utils::combn(d, 2)
  data.frame(i = idx[1, ], j = idx[2, ],
             pair = paste(STATIONS$short[idx[1, ]], "&", STATIONS$short[idx[2, ]]),
             stringsAsFactors = FALSE)
}

load_pobs <- function(source = c("raw", "gumbel")) {
  f <- c(raw = "pseudo_obs.csv", gumbel = "pseudo_obs_gumbel.csv")[match.arg(source)]
  readr::read_csv(file.path(PROC, f), show_col_types = FALSE) |>
    mutate(season = factor(season, levels = SEASONS),
           period = factor(period, levels = c("pre", "post")))
}

log_step <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

# ---- run bookkeeping ---------------------------------------------------------
STEP_T0 <- Sys.time()
RUN_LOG <- file.path(RES, "run_steps.csv")
finish_step <- function(id) {
  secs <- as.numeric(difftime(Sys.time(), STEP_T0, units = "secs"))
  row <- data.frame(step = id, finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
                    seconds = round(secs, 1), quick = QUICK,
                    pobs = Sys.getenv("POBS", "raw"), min_days = MIN_DAYS)
  utils::write.table(row, RUN_LOG, sep = ",", row.names = FALSE,
                     col.names = !file.exists(RUN_LOG), append = file.exists(RUN_LOG))
  log_step(id, " done (", round(secs), " s)")
}

# ---- result cache ------------------------------------------------------------
# Expensive results (bootstrap, Monte-Carlo, EV test) are cached under
# results/cache/. The key hashes the inputs, the parameters and the package
# source, so any change invalidates the cache; FORCE=1 recomputes regardless.
CACHE <- file.path(RES, "cache")          # QUICK runs keep a separate cache
dir.create(CACHE, showWarnings = FALSE)
code_hash <- function() {
  f <- sort(list.files(file.path(ROOT, "R"), full.names = TRUE))
  unname(tools::md5sum(f))
}
cached <- function(name, key, expr) {
  path <- file.path(CACHE, paste0(name, ".rds"))
  key_hash <- digest_key(list(key, code_hash()))
  if (!identical(Sys.getenv("FORCE"), "1") && file.exists(path)) {
    obj <- readRDS(path)
    if (identical(attr(obj, "cache_key"), key_hash)) {
      log_step("  cache hit: ", name)
      return(obj)
    }
  }
  t0 <- Sys.time()
  obj <- force(expr)
  attr(obj, "cache_key") <- key_hash
  attr(obj, "compute_seconds") <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  saveRDS(obj, path)
  log_step("  computed and cached: ", name, " (",
           round(attr(obj, "compute_seconds")), " s)")
  obj
}
digest_key <- function(x) {
  f <- tempfile()
  saveRDS(x, f, version = 2)
  on.exit(unlink(f))
  unname(tools::md5sum(f))
}
