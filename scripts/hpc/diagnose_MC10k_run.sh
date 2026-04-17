#!/bin/bash
set -euo pipefail

# Diagnose status of a chunked MC run folder produced by scripts/hpc/submit_MC_10k.sh and scripts/hpc/run_MC_10k_array.sh
# Usage:
#   bash scripts/hpc/diagnose_MC10k_run.sh MC_results/MC10k_chunks
#   ROOT=MC_results/MC_dry_run bash scripts/hpc/diagnose_MC10k_run.sh

ROOT="${1:-${ROOT:-}}"
if [[ -z "${ROOT}" ]]; then
  echo "Usage: bash diagnose_MC10k_run.sh <run_root>"
  exit 1
fi

if [[ ! -d "${ROOT}" ]]; then
  echo "ERROR: run root not found: ${ROOT}"
  exit 1
fi

count_files() {
  local pattern="$1"
  find "${ROOT}" -type f -name "${pattern}" | wc -l | tr -d ' '
}

TASK_DIRS=$(find "${ROOT}" -maxdepth 1 -type d -name "task_*" | sort)
TOTAL_TASKS=$(echo "${TASK_DIRS}" | sed '/^$/d' | wc -l | tr -d ' ')
SUCCESS_TASKS=$(find "${ROOT}" -type f -name "TASK_SUCCESS" | wc -l | tr -d ' ')
FAILED_TASKS=$(find "${ROOT}" -type f -name "TASK_FAILED" | wc -l | tr -d ' ')
CSV_COUNT=$(count_files "MC_chunk_outputs_*.csv")
MAT_COUNT=$(count_files "MC_5FU_PK_chunk_*.mat")
STATUS_JSON_COUNT=$(count_files "task_status.json")
LOG_COUNT=$(count_files "task_runtime.log")

echo "MC10k run diagnostics"
echo "  ROOT                : ${ROOT}"
echo "  TASK DIRS           : ${TOTAL_TASKS}"
echo "  TASK_SUCCESS files  : ${SUCCESS_TASKS}"
echo "  TASK_FAILED files   : ${FAILED_TASKS}"
echo "  chunk CSV files     : ${CSV_COUNT}"
echo "  chunk MAT files     : ${MAT_COUNT}"
echo "  task_status.json    : ${STATUS_JSON_COUNT}"
echo "  task_runtime.log    : ${LOG_COUNT}"

echo
if [[ "${FAILED_TASKS}" -gt 0 ]]; then
  echo "Failed tasks (up to first 20):"
  find "${ROOT}" -type f -name "TASK_FAILED" | sort | head -n 20
  echo
  IDS=$(find "${ROOT}" -type f -name "TASK_FAILED" | sed -n 's#.*task_\([0-9][0-9]*\)/TASK_FAILED#\1#p' | sed 's/^0*//' | sed '/^$/d' | sort -n | paste -sd, -)
  if [[ -n "${IDS}" ]]; then
    echo "Suggested selective rerun command:"
    echo "  sbatch --array=${IDS} --export=\"ALL,CHUNK_SIZE=<chunk>,SEED_BASE=<seed_base>,OUT_ROOT=<same_root>,RETRY_ATTEMPTS=<n>\" scripts/hpc/run_MC_10k_array.sh"
  fi
else
  echo "No failed task markers detected."
fi

echo
if [[ "${CSV_COUNT}" -eq 0 ]]; then
  echo "WARNING: no chunk CSV outputs found. Check task_runtime.log and matlab_attempt_*.log files."
fi

if [[ "${CSV_COUNT}" -gt 0 ]]; then
  echo "Merge command:"
  echo "  python scripts/python/merge_MC_chunk_results.py --root ${ROOT}"
fi
