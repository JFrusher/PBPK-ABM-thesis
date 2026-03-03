#!/bin/bash
#SBATCH --job-name=mc5fu_single
#SBATCH --output=slurm_MC_total.out
#SBATCH --error=slurm_MC_total.err
#SBATCH --account=student
#SBATCH --partition=amd_student
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G              # Optimized from 64G based on 30.75G usage
#SBATCH --time=16:00:00        # Optimized from 12h based on 16m usage

set -euo pipefail

module purge
module load matlab/2025a 

cd "$SLURM_SUBMIT_DIR"

# Set samples to 1000
N_SAMPLES=1000
RUN_OUT="$SLURM_SUBMIT_DIR/MC_results/Total_Batch"
mkdir -p "$RUN_OUT"

# Call the function ONCE with the 1000 argument
matlab -batch "try, MC_5FU_PK_sensitivity(${N_SAMPLES}, '${RUN_OUT}', true); catch ME, disp(getReport(ME,'extended')); exit(1); end; exit(0);"
