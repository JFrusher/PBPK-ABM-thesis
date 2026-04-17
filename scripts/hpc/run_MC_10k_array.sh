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

timestamp() {
	date +"%Y-%m-%d %H:%M:%S"
}

log() {
	local level="$1"
	shift
	local msg="$*"
	echo "[$(timestamp)] [${level}] ${msg}"
	if [[ -n "${TASK_LOG:-}" ]]; then
		echo "[$(timestamp)] [${level}] ${msg}" >> "${TASK_LOG}"
	fi
}

die() {
	log "ERROR" "$*"
	exit 1
}

is_pos_int() {
	[[ "$1" =~ ^[0-9]+$ ]] && [[ "$1" -ge 1 ]]
}

is_nonneg_int() {
	[[ "$1" =~ ^[0-9]+$ ]]
}

mark_failed() {
	local reason="$1"
	if [[ -n "${TASK_OUT:-}" ]]; then
		echo "${reason}" > "${TASK_OUT}/TASK_FAILED"
	fi
	if [[ -n "${TASK_STATUS_FILE:-}" ]]; then
		printf '{"status":"failed","reason":"%s","finished_at":"%s"}\n' "${reason//\"/\\\"}" "$(timestamp)" > "${TASK_STATUS_FILE}"
	fi
}

on_err() {
	local rc=$?
	mark_failed "script_error_rc_${rc}"
	log "ERROR" "Unhandled failure (exit code ${rc})"
	exit ${rc}
}

trap on_err ERR

module purge
module load matlab/2025a

cd "$SLURM_SUBMIT_DIR"

# Runtime knobs (override with sbatch --export=ALL,VAR=...)
CHUNK_SIZE="${CHUNK_SIZE:-100}"
SEED_BASE="${SEED_BASE:-100000}"
OUT_ROOT="${OUT_ROOT:-$SLURM_SUBMIT_DIR/MC_results/MC10k_chunks}"
RETRY_ATTEMPTS="${RETRY_ATTEMPTS:-1}"

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
TASK_SEED=$((SEED_BASE + TASK_ID - 1))
TASK_OUT="${OUT_ROOT}/task_$(printf "%05d" "${TASK_ID}")"
mkdir -p "${TASK_OUT}"

TASK_LOG="${TASK_OUT}/task_runtime.log"
TASK_STATUS_FILE="${TASK_OUT}/task_status.json"
TASK_ENV_FILE="${TASK_OUT}/task_env.txt"

is_pos_int "${CHUNK_SIZE}" || die "CHUNK_SIZE must be a positive integer (got ${CHUNK_SIZE})"
is_nonneg_int "${SEED_BASE}" || die "SEED_BASE must be a non-negative integer (got ${SEED_BASE})"
is_nonneg_int "${RETRY_ATTEMPTS}" || die "RETRY_ATTEMPTS must be a non-negative integer (got ${RETRY_ATTEMPTS})"

[[ -f "matlab/analysis/MC_5FU_PK_sensitivity.m" ]] || die "matlab/analysis/MC_5FU_PK_sensitivity.m not found in ${SLURM_SUBMIT_DIR}"

{
	echo "job_id=${SLURM_JOB_ID:-na}"
	echo "job_name=${SLURM_JOB_NAME:-na}"
	echo "array_job_id=${SLURM_ARRAY_JOB_ID:-na}"
	echo "array_task_id=${TASK_ID}"
	echo "node=${SLURMD_NODENAME:-na}"
	echo "cpus_per_task=${SLURM_CPUS_PER_TASK:-na}"
	echo "chunk_size=${CHUNK_SIZE}"
	echo "seed_base=${SEED_BASE}"
	echo "task_seed=${TASK_SEED}"
	echo "retry_attempts=${RETRY_ATTEMPTS}"
	echo "submit_dir=${SLURM_SUBMIT_DIR:-na}"
	echo "task_out=${TASK_OUT}"
	echo "started_at=$(timestamp)"
} > "${TASK_ENV_FILE}"

printf '{"status":"running","started_at":"%s","task_id":%s,"seed":%s}\n' "$(timestamp)" "${TASK_ID}" "${TASK_SEED}" > "${TASK_STATUS_FILE}"

log "INFO" "MC-CHUNK start: job=${SLURM_JOB_ID:-na} task=${TASK_ID} chunk=${CHUNK_SIZE} seed=${TASK_SEED}"
log "INFO" "Output folder: ${TASK_OUT}"

ATTEMPT=0
MAX_ATTEMPTS=$((RETRY_ATTEMPTS + 1))
MATLAB_RC=1
while [[ ${ATTEMPT} -lt ${MAX_ATTEMPTS} ]]; do
	ATTEMPT=$((ATTEMPT + 1))
	ATTEMPT_LOG="${TASK_OUT}/matlab_attempt_${ATTEMPT}.log"
	log "INFO" "Starting MATLAB attempt ${ATTEMPT}/${MAX_ATTEMPTS}"

	set +e
	matlab -batch "try, addpath(genpath(fullfile(pwd,'matlab'))); MC_5FU_PK_sensitivity(${CHUNK_SIZE}, '${TASK_OUT}', true, ${TASK_SEED}, true); catch ME, disp(getReport(ME,'extended')); exit(1); end; exit(0);" > "${ATTEMPT_LOG}" 2>&1
	MATLAB_RC=$?
	set -e

	if [[ ${MATLAB_RC} -eq 0 ]]; then
		log "INFO" "MATLAB attempt ${ATTEMPT} succeeded"
		break
	fi

	log "WARN" "MATLAB attempt ${ATTEMPT} failed (rc=${MATLAB_RC}). Log: ${ATTEMPT_LOG}"
done

if [[ ${MATLAB_RC} -ne 0 ]]; then
	mark_failed "matlab_failed_after_${MAX_ATTEMPTS}_attempts"
	die "MATLAB failed after ${MAX_ATTEMPTS} attempt(s)"
fi

CSV_COUNT=$(find "${TASK_OUT}" -type f -name "MC_chunk_outputs_*.csv" | wc -l | tr -d ' ')
MAT_COUNT=$(find "${TASK_OUT}" -type f -name "MC_5FU_PK_chunk_*.mat" | wc -l | tr -d ' ')

if [[ "${CSV_COUNT}" -lt 1 || "${MAT_COUNT}" -lt 1 ]]; then
	mark_failed "missing_chunk_outputs_csv_${CSV_COUNT}_mat_${MAT_COUNT}"
	die "Expected chunk outputs not found under ${TASK_OUT} (csv=${CSV_COUNT}, mat=${MAT_COUNT})"
fi

rm -f "${TASK_OUT}/TASK_FAILED"
touch "${TASK_OUT}/TASK_SUCCESS"
printf '{"status":"success","finished_at":"%s","csv_count":%s,"mat_count":%s,"attempts_used":%s}\n' \
	"$(timestamp)" "${CSV_COUNT}" "${MAT_COUNT}" "${ATTEMPT}" > "${TASK_STATUS_FILE}"

log "INFO" "Chunk outputs detected: csv=${CSV_COUNT}, mat=${MAT_COUNT}"
log "INFO" "MC-CHUNK done: task=${TASK_ID}"
