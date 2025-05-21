#!/usr/bin/env bash
#SBATCH --job-name=REACTER_AUTOMATION_ML
#SBATCH --partition=pi.amartini
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --time=10-00:00:00
#SBATCH --export=ALL
#SBATCH --mem=10G
#SBATCH --output=automation_log.out
#SBATCH --error=automation_error.out

# Activate environment
source ~/.bashrc
conda activate saenv

# Source required functions
source /mnt/borgstore/amartini/sahmed73/data/automation_scripts/auto_job_submitter/job_submission_functions.sh
source FUNCTION_SCRIPTS/reacter_input_generation_functions.sh

# User Inputs
ROOT_DIR="/mnt/borgstore/amartini/sahmed73/data/automation_scripts/REACTER_AUTO_MACHINE_LEARNING"
WORKING_DIR="${ROOT_DIR}/test_delete_afterwards"
CSV_FILE="/mnt/borgstore/amartini/sahmed73/data/SMILES/test.csv"

mkdir -p "$WORKING_DIR" || { echo "Directory already exists. Exiting." >&2; exit 1; }

# Ensure CSV has a newline at end
sed -i -e '$a\' "$CSV_FILE"
cp -n "$CSV_FILE" "$WORKING_DIR/"

echo "===== Parsed CSV lines (excluding header) ====="
awk -F',' 'NR > 1 && NF == 2 { printf "%s|%s\n", $1, $2 }' "$CSV_FILE"
echo "==============================================="

# Export required environment variables and functions for parallel
export WORKING_DIR
source FUNCTION_SCRIPTS/export_reacter_input_generation_functions_env.sh # exporting

PARALLEL_JOBS=${SLURM_NTASKS_PER_NODE:-4}  # fallback to 4 if not set
echo "Executing $PARALLEL_JOBS parallel jobs..."
awk -F',' 'NR>1 && NF==2' "$CSV_FILE" | parallel -j "$PARALLEL_JOBS" --colsep ',' \
    run_reacter_in_gen \
        --mole1 "PAO-OO" \
        --smiles1 "CCCCCCCCCCC" \
        --mole2 {1} \
        --smiles2 {2} \
        --n_mole1 "100" \
        --n_mole2 "50" \
        --density "0.15" \
        --forcefield "PCFF-IFF" \
        --ff_path "frc_files/pcff_interface_v1_6mBN.frc" \
        --output_dir "$WORKING_DIR" \
        --runtime "10" \
        --Rmax "2.0" \
        --Rprob "0.1" \
        --nframe "1000" \
        --nsim "2"

echo "ALL REACTER INPUT FILES CREATED!!"