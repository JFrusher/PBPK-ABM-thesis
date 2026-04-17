#!/bin/bash
set -euo pipefail

# Usage:
#   bash submit_MC_10k.sh
#   TOTAL_SAMPLES=10000 CHUNK_SIZE=100 MAX_CONCURRENT=25 bash submit_MC_10k.sh
# Optional runtime knobs:
#   DRY_RUN=1                 # print/record submission plan without calling sbatch
#   RETRY_ATTEMPTS=1          # per-task MATLAB retries handled in run_MC_10k_array.sh
#   LOG_ROOT=<path>           # where submission manifests are saved (default: <OUT_ROOT>/_submission)

TOTAL_SAMPLES="${TOTAL_SAMPLES:-10000}"
CHUNK_SIZE="${CHUNK_SIZE:-100}"
MAX_CONCURRENT="${MAX_CONCURRENT:-20}"
SEED_BASE="${SEED_BASE:-100000}"
OUT_ROOT="${OUT_ROOT:-$PWD/MC_results/MC10k_chunks}"
DRY_RUN="${DRY_RUN:-0}"
RETRY_ATTEMPTS="${RETRY_ATTEMPTS:-1}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_SCRIPT="${SCRIPT_DIR}/run_MC_10k_array.sh"

die() {
  echo "ERROR: $*" >&2
  exit 1
}

is_pos_int() {
  [[ "$1" =~ ^[0-9]+$ ]] && [[ "$1" -ge 1 ]]
}

is_nonneg_int() {
  [[ "$1" =~ ^[0-9]+$ ]]
}

is_pos_int "${TOTAL_SAMPLES}" || die "TOTAL_SAMPLES must be a positive integer (got: ${TOTAL_SAMPLES})"
is_pos_int "${CHUNK_SIZE}" || die "CHUNK_SIZE must be a positive integer (got: ${CHUNK_SIZE})"
is_pos_int "${MAX_CONCURRENT}" || die "MAX_CONCURRENT must be a positive integer (got: ${MAX_CONCURRENT})"
is_nonneg_int "${SEED_BASE}" || die "SEED_BASE must be a non-negative integer (got: ${SEED_BASE})"
is_nonneg_int "${RETRY_ATTEMPTS}" || die "RETRY_ATTEMPTS must be a non-negative integer (got: ${RETRY_ATTEMPTS})"

[[ -f "${RUN_SCRIPT}" ]] || die "run script not found: ${RUN_SCRIPT}"

if [[ "${DRY_RUN}" != "1" ]]; then
  command -v sbatch >/dev/null 2>&1 || die "sbatch not found in PATH"
fi

ARRAY_TASKS=$(( (TOTAL_SAMPLES + CHUNK_SIZE - 1) / CHUNK_SIZE ))
ARRAY_SPEC="1-${ARRAY_TASKS}%${MAX_CONCURRENT}"

mkdir -p "${OUT_ROOT}"

LOG_ROOT="${LOG_ROOT:-${OUT_ROOT}/_submission}"
mkdir -p "${LOG_ROOT}"
STAMP="$(date +%Y%m%d_%H%M%S)"
MANIFEST="${LOG_ROOT}/submission_${STAMP}.txt"

{
  echo "timestamp=${STAMP}"
  echo "pwd=${PWD}"
  echo "script_dir=${SCRIPT_DIR}"
  echo "run_script=${RUN_SCRIPT}"
  echo "TOTAL_SAMPLES=${TOTAL_SAMPLES}"
  echo "CHUNK_SIZE=${CHUNK_SIZE}"
  echo "MAX_CONCURRENT=${MAX_CONCURRENT}"
  echo "ARRAY_TASKS=${ARRAY_TASKS}"
  echo "ARRAY_SPEC=${ARRAY_SPEC}"
  echo "SEED_BASE=${SEED_BASE}"
  echo "RETRY_ATTEMPTS=${RETRY_ATTEMPTS}"
  echo "OUT_ROOT=${OUT_ROOT}"
  echo "DRY_RUN=${DRY_RUN}"
} > "${MANIFEST}"

echo "Submitting MC chunk array"
echo "  TOTAL_SAMPLES : ${TOTAL_SAMPLES}"
echo "  CHUNK_SIZE    : ${CHUNK_SIZE}"
echo "  ARRAY_TASKS   : ${ARRAY_TASKS}"
echo "  MAX_CONCURRENT: ${MAX_CONCURRENT}"
echo "  ARRAY_SPEC    : ${ARRAY_SPEC}"
echo "  OUT_ROOT      : ${OUT_ROOT}"
echo "  RETRY_ATTEMPTS: ${RETRY_ATTEMPTS}"
echo "  MANIFEST      : ${MANIFEST}"

SBATCH_CMD=(
  sbatch
  "--array=${ARRAY_SPEC}"
  "--export=ALL,CHUNK_SIZE=${CHUNK_SIZE},SEED_BASE=${SEED_BASE},OUT_ROOT=${OUT_ROOT},RETRY_ATTEMPTS=${RETRY_ATTEMPTS}"
  "${RUN_SCRIPT}"
)

printf "sbatch_command=" >> "${MANIFEST}"
printf "%q " "${SBATCH_CMD[@]}" >> "${MANIFEST}"
printf "\n" >> "${MANIFEST}"

if [[ "${DRY_RUN}" == "1" ]]; then
  echo "DRY_RUN=1 -> not submitting job"
  printf "Would run: "
  printf "%q " "${SBATCH_CMD[@]}"
  printf "\n"
  exit 0
fi

set +e
SBATCH_OUT="$(${SBATCH_CMD[@]} 2>&1)"
SBATCH_RC=$?
set -e

echo "sbatch_output=${SBATCH_OUT}" >> "${MANIFEST}"
echo "sbatch_rc=${SBATCH_RC}" >> "${MANIFEST}"

if [[ ${SBATCH_RC} -ne 0 ]]; then
  echo "ERROR: sbatch submission failed"
  echo "${SBATCH_OUT}"
  exit ${SBATCH_RC}
fi

echo "${SBATCH_OUT}"

JOB_ID="$(echo "${SBATCH_OUT}" | sed -n 's/.*Submitted batch job \([0-9][0-9]*\).*/\1/p' | head -n 1)"
if [[ -n "${JOB_ID}" ]]; then
  echo "job_id=${JOB_ID}" >> "${MANIFEST}"
  echo "Submission succeeded: job_id=${JOB_ID}"
else
  echo "WARNING: could not parse job id from sbatch output"
fi

echo "Submission metadata saved to: ${MANIFEST}"
