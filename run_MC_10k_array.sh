#!/bin/bash
#SBATCH --job-name=mc5fu_chunk
#SBATCH --output=slurm_mc_chunk_%A_%a.out
#SBATCH --error=slurm_mc_chunk_%A_%a.err
#SBATCH --account=student
#SBATCH --partition=amd_student
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --array=1-100%20

set -euo pipefail

module purge
module load matlab/2025a

cd "$SLURM_SUBMIT_DIR"

# Runtime knobs (override with sbatch --export=ALL,VAR=...)
CHUNK_SIZE="${CHUNK_SIZE:-100}"
SEED_BASE="${SEED_BASE:-100000}"
OUT_ROOT="${OUT_ROOT:-$SLURM_SUBMIT_DIR/MC_results/MC10k_chunks}"

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
TASK_SEED=$((SEED_BASE + TASK_ID - 1))
TASK_OUT="${OUT_ROOT}/task_$(printf "%05d" "${TASK_ID}")"
mkdir -p "${TASK_OUT}"

echo "[MC-CHUNK] Job=${SLURM_JOB_ID} Task=${TASK_ID} ChunkSize=${CHUNK_SIZE} Seed=${TASK_SEED}"
echo "[MC-CHUNK] Output=${TASK_OUT}"

matlab -batch "try, MC_5FU_PK_sensitivity(${CHUNK_SIZE}, '${TASK_OUT}', true, ${TASK_SEED}, true); catch ME, disp(getReport(ME,'extended')); exit(1); end; exit(0);"

echo "[MC-CHUNK] Done task ${TASK_ID}"
