#!/bin/bash
#SBATCH --job-name=PhysiCell_Post
#SBATCH --output=post_%j.txt
#SBATCH --error=post_err_%j.txt
#SBATCH --account=student
#SBATCH --partition=amd_student
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4      # Reduced: 4 is plenty for image/excel work
#SBATCH --mem=8G               # Increased: Image processing can be RAM-heavy
#SBATCH --time=02:00:00        # Shortened: 30 mins is usually enough

# 1. Load Environment
module load python/3.12.4
source $(conda info --base)/etc/profile.d/conda.sh
conda activate my_imaging_env

# 2. Simple Threading (for 'make jpeg' if it uses Magick++)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# 3. Execution
echo "Converting snapshots to JPEGs..."
make jpeg

echo "Exporting to Excel..."
# Using ${SLURM_JOB_ID} is more reliable in bash than %j
python3 export_outputs_to_excel.py --folder output --outdir Sim_results_${SLURM_JOB_ID}_csv

conda deactivate
make movie
echo "Post-processing complete."
