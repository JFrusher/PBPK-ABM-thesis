#!/bin/bash
#SBATCH --job-name=auc_tune_wp2
#SBATCH --output=output_%A_%a.txt
#SBATCH --error=error_%A_%a.txt
#SBATCH --array=1-8
#SBATCH --account=student
#SBATCH --partition=amd_student
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=08:00:00

set -euo pipefail

# For dynamic sizing based on the input list, submit via submit_auc_tune_wp2.sh

module purge
module load matlab/2025a

cd "$SLURM_SUBMIT_DIR"

# Tunable runtime settings (override with: sbatch --export=ALL,TARGET_AUC=30,...)
TARGET_AUC="${TARGET_AUC:-25.0}"
TOLERANCE="${TOLERANCE:-0.1}"
MAX_ITERATIONS="${MAX_ITERATIONS:-20}"
DISABLE_PLOTS="${DISABLE_PLOTS:-1}"
CLEANUP_ON_CONVERGED="${CLEANUP_ON_CONVERGED:-1}"
INPUT_LIST="${INPUT_LIST:-scripts/hpc/slurm_auc_inputs_wp2.txt}"

if [[ ! -f "$INPUT_LIST" ]]; then
  echo "ERROR: Input list not found: $INPUT_LIST"
  exit 1
fi

mapfile -t CSV_LIST < <(grep -v '^\s*#' "$INPUT_LIST" | sed '/^\s*$/d')
N_FILES="${#CSV_LIST[@]}"

if [[ "$N_FILES" -eq 0 ]]; then
  echo "ERROR: No CSV inputs found in $INPUT_LIST"
  exit 1
fi

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
if (( TASK_ID < 1 || TASK_ID > N_FILES )); then
  echo "ERROR: SLURM_ARRAY_TASK_ID=${TASK_ID} is out of bounds for ${N_FILES} inputs"
  exit 1
fi

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "WARNING: SLURM_ARRAY_TASK_ID is not set; running first input only."
fi

INPUT_CSV="${CSV_LIST[$((TASK_ID-1))]}"
if [[ ! -f "$INPUT_CSV" ]]; then
  echo "ERROR: CSV file not found: $INPUT_CSV"
  exit 1
fi

CSV_BASE="$(basename "$INPUT_CSV" .csv)"
RUN_ROOT="AUC_tuning_runs/${CSV_BASE}"
mkdir -p "$RUN_ROOT"

OUTPUT_PREFIX="$RUN_ROOT/${CSV_BASE}_AUCtuned"
SUMMARY_FILE="$RUN_ROOT/${CSV_BASE}_summary.txt"

echo "=== AUC tuning task start ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array task: ${TASK_ID}/${N_FILES}"
echo "Input CSV: ${INPUT_CSV}"
echo "Target AUC: ${TARGET_AUC}"
echo "Tolerance: ${TOLERANCE}"
echo "Max iterations: ${MAX_ITERATIONS}"
echo "Output prefix: ${OUTPUT_PREFIX}"

MATLAB_CMD="try, addpath(genpath(pwd)); [optCsv, s] = tune_csv_to_auc_target('${INPUT_CSV}', 'TargetAUC', ${TARGET_AUC}, 'Tolerance', ${TOLERANCE}, 'MaxIterations', ${MAX_ITERATIONS}, 'OutputPrefix', '${OUTPUT_PREFIX}', 'DisablePlots', logical(${DISABLE_PLOTS}), 'CleanupOnConverged', logical(${CLEANUP_ON_CONVERGED})); fid=fopen('${SUMMARY_FILE}','w'); fprintf(fid,'inputCsv=%s\\n', s.inputCsv); fprintf(fid,'optimizedCsv=%s\\n', s.optimizedCsv); fprintf(fid,'converged=%d\\n', s.converged); fprintf(fid,'iterations=%d\\n', s.iterations); fprintf(fid,'finalAUC=%.8f\\n', s.finalAUC); fprintf(fid,'finalError=%.8f\\n', s.finalError); fprintf(fid,'cumulativeScale=%.8f\\n', s.cumulativeScale); fclose(fid); disp(optCsv); catch ME, disp(getReport(ME,'extended')); exit(1); end; exit(0);"

matlab -batch "$MATLAB_CMD"

echo "=== AUC tuning task complete ==="
