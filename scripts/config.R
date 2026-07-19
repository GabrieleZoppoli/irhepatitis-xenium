# ============================================================================
# config.R — centralised paths for the irHepatitis publication pipeline
# ============================================================================
# Every R script in this folder sources this file. Paths resolve relative
# to the bundle root (the directory containing this scripts/ folder).
#
# Override at runtime by setting environment variables:
#   IRHEP_BUNDLE_ROOT   — bundle root (default: detected via here::here())
#   IRHEP_RAW_DIR       — raw Xenium data dir (default: bundle/data/raw)
#   IRHEP_PROCESSED_DIR — processed data dir (default: bundle/data/processed)
# ============================================================================

# Detect bundle root: prefer IRHEP_BUNDLE_ROOT env, else here::here()
.detect_bundle_root <- function() {
  env_root <- Sys.getenv("IRHEP_BUNDLE_ROOT", unset = NA)
  if (!is.na(env_root) && nzchar(env_root)) return(normalizePath(env_root))
  if (requireNamespace("here", quietly = TRUE)) return(here::here())
  # Fallback: assume cwd is at scripts/ or bundle root
  cwd <- normalizePath(getwd())
  if (basename(cwd) == "scripts") return(dirname(cwd))
  cwd
}

BUNDLE_ROOT <- .detect_bundle_root()

# Resolve all dependent paths
RAW_DIR       <- Sys.getenv("IRHEP_RAW_DIR",
                            unset = file.path(BUNDLE_ROOT, "data", "raw"))
PROCESSED_DIR <- Sys.getenv("IRHEP_PROCESSED_DIR",
                            unset = file.path(BUNDLE_ROOT, "data", "processed"))
FIGURES_DIR   <- file.path(BUNDLE_ROOT, "figures")
TABLES_DIR    <- file.path(BUNDLE_ROOT, "tables")
SCRIPTS_DIR   <- file.path(BUNDLE_ROOT, "scripts")

# Standard analysis-internal paths used by Fig 2 scripts
SESSION_SNAPSHOT  <- file.path(PROCESSED_DIR, "session_snapshot.RData")
CLUSTER_CONFIG    <- file.path(PROCESSED_DIR, "cluster_config.rds")
CLUSTER_MARKERS   <- file.path(PROCESSED_DIR, "clustering", "cluster_markers.csv")
CELLCHAT_COHERENT <- file.path(PROCESSED_DIR, "cellchat",
                               "cellchat_coherent_interactions.csv")

# Sample identifiers
SAMPLES <- list(
  irHep = c("18321_a", "24774_d"),
  AIH   = c("54023_b", "56784_c")
)

# Per-sample x/y shifts used in multi-sample composite plotting
SAMPLE_SHIFTS <- data.frame(
  sample_id = c("18321_a", "24774_d", "54023_b", "56784_c"),
  x_shift   = c(0, 12000, 0, 12000),
  y_shift   = c(0, 0, -8000, -8000),
  condition = c("irHep", "irHep", "AIH", "AIH"),
  stringsAsFactors = FALSE
)

# Helper: locate per-sample raw directory
sample_raw_dir <- function(sample_id) {
  file.path(RAW_DIR, paste0("0027420_", sample_id))
}

# Confirm bundle root exists; warn if expected dirs missing
if (!dir.exists(BUNDLE_ROOT))
  stop(sprintf("BUNDLE_ROOT not found: %s — set IRHEP_BUNDLE_ROOT", BUNDLE_ROOT))
if (!dir.exists(RAW_DIR))
  warning(sprintf("RAW_DIR not found: %s — populate or set IRHEP_RAW_DIR", RAW_DIR))
if (!dir.exists(PROCESSED_DIR))
  warning(sprintf("PROCESSED_DIR not found: %s", PROCESSED_DIR))

invisible(NULL)
