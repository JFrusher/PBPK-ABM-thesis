#!/bin/bash
#SBATCH --job-name=PhysiCell_Array
#SBATCH --output=output_%A_%a.txt
#SBATCH --error=error_%A_%a.txt
#SBATCH --array=1-10               # This creates 10 sub-jobs
#SBATCH --account=student
#SBATCH --partition=amd_student
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --cpus-per-task=16
#SBATCH --time=03:00:00            # Shorter time needed since they run at once

# 1. Load Environment
module load gcc/15.2.0

# 2. Thread Binding
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

# 3. Create a unique output folder for THIS array instance
# $SLURM_ARRAY_TASK_ID will be 1, 2, 3... up to 10
OUT_DIR="output_run_$SLURM_ARRAY_TASK_ID"
mkdir -p $OUT_DIR

# 4. Run the simulation
# Note: We point the project to a specific folder so they don't overwrite each other
./project ./config/PhysiCell_settings.xml

# 5. Move results and Zip
mv output/* $OUT_DIR/
zip -r "${OUT_DIR}.zip" $OUT_DIR
rm -rf $OUT_DIR

echo "Array Task $SLURM_ARRAY_TASK_ID finished."