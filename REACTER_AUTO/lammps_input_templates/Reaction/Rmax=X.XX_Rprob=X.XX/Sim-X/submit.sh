#!/bin/bash
#SBATCH --job-name=<<job-name>>
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=6:00:00
#SBATCH --export=ALL
#SBATCH --mem=60G

module purge 
module load openmpi/4.0.6-intel-2021.4.0

LMP=/data/sahmed73/LAMMPS/LATEST_KOKKOS_PYTHON/lammps-29Sep2021/build/lmp_kk

NSLOTS=$(($SLURM_NNODES * $SLURM_NTASKS_PER_NODE))
OMP_NUM_THREADS=1

# === Log Start Time ===
start_time=$(date +%s)
echo "Job started at: $(date)"

# === Run Simulation ===
mpirun -np $NSLOTS $LMP -in input.in > output.out

# === Log End Time and Runtime ===
end_time=$(date +%s)
echo "Job ended at: $(date)"

# === Calculate Elapsed Time ===
elapsed=$((end_time - start_time))
hours=$((elapsed / 3600))
minutes=$(((elapsed % 3600) / 60))
seconds=$((elapsed % 60))
printf "Total runtime: %02d:%02d:%02d (hh:mm:ss)\n" $hours $minutes $seconds