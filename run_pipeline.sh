#!/usr/bin/env bash
# =============================================================================
# Script:      run_pipeline.sh
# Pipeline:    16S Microbiota Pipeline
# Author:      Roberto C. Torres, PhD. <torres.roberto.c@gmail.com>
# Institution: Bioinformatics Lab, Infectious Diseases Research Unit,
#              CMN SXXI, IMSS, Mexico City
# Description: Run the full pipeline or from a specific step
# =============================================================================

set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SCRIPTS_DIR="${SCRIPT_DIR}/scripts"
export SCRIPTS_DIR
DEFAULT_CONFIG="${SCRIPT_DIR}/config/biopsia_saliva_V4V5_fullrun_20260505.R"
CONDA_RSCRIPT="${SCRIPT_DIR}/../../../../../../miniconda3/bin/Rscript"
ALT_CONDA_RSCRIPT="/Users/rtorres/miniconda3/bin/Rscript"

print_help() {
  cat <<EOF
Usage: run_pipeline.sh [OPTIONS]

Options:
  --config <path>   Path to config R file (required)
  --outdir <path>   Existing results_* run directory (required when --from > 1;
                    ignored when --from = 1, where the directory is created here)
  --from <N>        Start from step N (1–12, default: 1)
  --help            Show this help message

Steps:
  1  01_qc.R
  2  02_demux.R
  3  03_dada2.R
  4  04_collapse.R
  5  05_tax.R
  6  06_filter.R
  7  07_alpha.R
  8  08_beta.R
  9  09_composition.R
  10 10_differential.R
  11 11_informative.R
  12 12_report.R
EOF
}

detect_rscript() {
  local detected=""

  if detected="$(which Rscript 2>/dev/null)"; then
    if [ "$detected" = "$CONDA_RSCRIPT" ] || [ "$detected" = "$ALT_CONDA_RSCRIPT" ]; then
      printf '%s\n' "$detected"
      return 0
    fi
    if [ -x "$ALT_CONDA_RSCRIPT" ]; then
      printf '%s\n' "$ALT_CONDA_RSCRIPT"
      return 0
    fi
    printf '%s\n' "$detected"
    return 0
  fi

  if [ -x "$ALT_CONDA_RSCRIPT" ]; then
    printf '%s\n' "$ALT_CONDA_RSCRIPT"
    return 0
  fi

  if [ -x "$CONDA_RSCRIPT" ]; then
    printf '%s\n' "$CONDA_RSCRIPT"
    return 0
  fi

  return 1
}

CONFIG="$DEFAULT_CONFIG"
OUTDIR=""
FROM_STEP=1

while [ $# -gt 0 ]; do
  case "$1" in
    --config)
      [ $# -ge 2 ] || { echo "ERROR: --config requires a path."; exit 1; }
      CONFIG="$2"
      shift 2
      ;;
    --outdir)
      [ $# -ge 2 ] || { echo "ERROR: --outdir requires a path."; exit 1; }
      OUTDIR="$2"
      shift 2
      ;;
    --from)
      [ $# -ge 2 ] || { echo "ERROR: --from requires an integer."; exit 1; }
      FROM_STEP="$2"
      shift 2
      ;;
    --help)
      print_help
      exit 0
      ;;
    *)
      echo "ERROR: unknown option '$1'."
      print_help
      exit 1
      ;;
  esac
done

case "$FROM_STEP" in
  ''|*[!0-9]*)
    echo "ERROR: --from must be an integer between 1 and 12."
    exit 1
    ;;
esac

if [ "$FROM_STEP" -lt 1 ] || [ "$FROM_STEP" -gt 12 ]; then
  echo "ERROR: --from must be an integer between 1 and 12."
  exit 1
fi

if [ ! -f "$CONFIG" ]; then
  echo "ERROR: config file not found: $CONFIG"
  exit 1
fi

RSCRIPT_BIN="$(detect_rscript)" || {
  echo "ERROR: Rscript not found."
  exit 1
}

STEPS=(
  "01_qc.R"
  "02_demux.R"
  "03_dada2.R"
  "04_collapse.R"
  "05_tax.R"
  "06_filter.R"
  "07_alpha.R"
  "08_beta.R"
  "09_composition.R"
  "10_differential.R"
  "11_informative.R"
  "12_report.R"
)

TOTAL_STEPS="${#STEPS[@]}"

# --- Determine the run directory ---
if [ "$FROM_STEP" -eq 1 ]; then
  # Extract RESULTS_BASE from the config and create a fresh timestamped run dir.
  RESULTS_BASE=$("$RSCRIPT_BIN" --vanilla -e "source('${CONFIG}'); cat(RESULTS_BASE)" 2>/dev/null)
  if [ -z "$RESULTS_BASE" ]; then
    echo "ERROR: RESULTS_BASE is not defined in config: ${CONFIG}"
    exit 1
  fi
  RUN_DIR="${RESULTS_BASE}/results_$(date +%Y%m%d_%H%M%S)"
  mkdir -p "$RUN_DIR" || { echo "ERROR: cannot create run directory: ${RUN_DIR}"; exit 1; }
  echo "Run directory: ${RUN_DIR}"
else
  # Resuming from a later step — the run directory must be supplied.
  if [ -z "$OUTDIR" ]; then
    echo "ERROR: --outdir is required when using --from > 1"
    exit 1
  fi
  RUN_DIR="$OUTDIR"
  if [ ! -d "$RUN_DIR" ]; then
    echo "ERROR: outdir does not exist: ${RUN_DIR}"
    exit 1
  fi
fi

for ((idx = FROM_STEP - 1; idx < TOTAL_STEPS; idx++)); do
  step_num=$((idx + 1))
  script_name="${STEPS[$idx]}"

  echo "=== [Step ${step_num}/12] ${script_name} ==="

  "$RSCRIPT_BIN" "${SCRIPTS_DIR}/${script_name}" "$CONFIG" --outdir "$RUN_DIR"

  if [ $? -ne 0 ]; then
    echo "ERROR: step ${step_num} (${script_name}) failed. Aborting."
    exit 1
  fi
done

echo "Pipeline completed successfully. Steps ${FROM_STEP}–12 done."
