#!/bin/bash
# Helper submitter: auto-detects input count and submits matching array range.
# Usage:
#   bash submit_auc_tune_wp2.sh
#   INPUT_LIST=scripts/hpc/slurm_auc_inputs_wp2.txt TARGET_AUC=30 bash submit_auc_tune_wp2.sh

set -euo pipefail

INPUT_LIST="${INPUT_LIST:-scripts/hpc/slurm_auc_inputs_wp2.txt}"
JOB_SCRIPT="${JOB_SCRIPT:-scripts/hpc/run_auc_tune_array.sh}"

if [[ ! -f "$INPUT_LIST" ]]; then
  echo "ERROR: Input list not found: $INPUT_LIST"
  exit 1
fi

if [[ ! -f "$JOB_SCRIPT" ]]; then
  echo "ERROR: Job script not found: $JOB_SCRIPT"
  exit 1
fi

N_FILES=$(grep -v '^\s*#' "$INPUT_LIST" | sed '/^\s*$/d' | wc -l | tr -d '[:space:]')
if [[ -z "$N_FILES" || "$N_FILES" -lt 1 ]]; then
  echo "ERROR: No CSV inputs found in $INPUT_LIST"
  exit 1
fi

echo "Submitting $JOB_SCRIPT with array range 1-$N_FILES using $INPUT_LIST"
sbatch --array="1-${N_FILES}" --export=ALL,INPUT_LIST="$INPUT_LIST" "$JOB_SCRIPT"
