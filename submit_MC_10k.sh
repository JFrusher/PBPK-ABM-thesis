#!/bin/bash
set -euo pipefail

# Usage:
#   bash submit_MC_10k.sh
#   TOTAL_SAMPLES=10000 CHUNK_SIZE=100 MAX_CONCURRENT=25 bash submit_MC_10k.sh

TOTAL_SAMPLES="${TOTAL_SAMPLES:-10000}"
CHUNK_SIZE="${CHUNK_SIZE:-100}"
MAX_CONCURRENT="${MAX_CONCURRENT:-20}"
SEED_BASE="${SEED_BASE:-100000}"
OUT_ROOT="${OUT_ROOT:-$PWD/MC_results/MC10k_chunks}"

if [[ "${CHUNK_SIZE}" -lt 1 ]]; then
  echo "ERROR: CHUNK_SIZE must be >= 1"
  exit 1
fi

ARRAY_TASKS=$(( (TOTAL_SAMPLES + CHUNK_SIZE - 1) / CHUNK_SIZE ))
ARRAY_SPEC="1-${ARRAY_TASKS}%${MAX_CONCURRENT}"

mkdir -p "${OUT_ROOT}"

echo "Submitting MC chunk array"
echo "  TOTAL_SAMPLES : ${TOTAL_SAMPLES}"
echo "  CHUNK_SIZE    : ${CHUNK_SIZE}"
echo "  ARRAY_TASKS   : ${ARRAY_TASKS}"
echo "  MAX_CONCURRENT: ${MAX_CONCURRENT}"
echo "  ARRAY_SPEC    : ${ARRAY_SPEC}"
echo "  OUT_ROOT      : ${OUT_ROOT}"

sbatch \
  --array="${ARRAY_SPEC}" \
  --export="ALL,CHUNK_SIZE=${CHUNK_SIZE},SEED_BASE=${SEED_BASE},OUT_ROOT=${OUT_ROOT}" \
  run_MC_10k_array.sh
