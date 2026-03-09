#!/bin/bash
#SBATCH --job-name=PhysiCell_Opt
#SBATCH --output=output_%j.txt
#SBATCH --error=error_%j.txt
#SBATCH --account=student
#SBATCH --partition=amd_student
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --cpus-per-task=16
#SBATCH --time=20:00:00

# 1. Load Compiler Environment
module load gcc/15.2.0

# 2. OPTIMIZATION: Re-compile for AMD EPYC Architecture
# Using -march=native ensures the binary uses the specific AVX instructions available on the node
echo "Compiling optimized binary..."
make clean
CFLAGS="-O3 -march=native -fopenmp -fno-trapping-math -funroll-loops" \
CXXFLAGS="-O3 -march=native -fopenmp -fno-trapping-math -funroll-loops" \
make -j$SLURM_CPUS_PER_TASK

# 3. Thread Binding (Critical for AMD Chiplet Performance)
# Spreading threads across cores maximizes memory bandwidth for the solver
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

# 4. Run Main Simulation
echo "Starting Simulation..."
./project ./config/PhysiCell_settings.xml

echo "Job Finished."


