#!/usr/bin/env bash
#SBATCH --job-name=REACTER_AUTOMATION_PEROXY
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=3-00:00:00 #10-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --export=ALL
#SBATCH --mem=60G
#SBATCH --output=automation_log.out
#SBATCH --error=automation_error.out

# activate conda environment
source ~/.bashrc
conda activate saenv

# source my_functions.sh
source /mnt/borgstore/amartini/sahmed73/data/automation_scripts/auto_job_submitter/job_submission_functions.sh
source reacter_input_generation_functions.sh

# User Inputs
ROOT_DIR="/mnt/borgstore/amartini/sahmed73/data/automation_scripts/REACTER_AUTO_PEROXY_RADICAL"

#=============================================================================
#                          Create Simulations Inputs
#=============================================================================

WORKING_DIR="${ROOT_DIR}/test_delete_afterwards"
mkdir "${WORKING_DIR}" || { echo "Directory already exists. Exiting."; exit 1; }


CSV_FILE="/mnt/borgstore/amartini/sahmed73/data/SMILES/12_AO_SMILES.csv"
cp -n "$CSV_FILE" "$WORKING_DIR/"

create_reacter_datafiles_from_smiles \
    --mole1 "PAO-OO" \
    --smiles1 "CCCCCCCCCCC(O[O])(CCCCCCCC)CC(C)CCCCCCCC" \
    --mole2 "BHT" \
    --smiles2 "Cc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1" \
    --n_mole1 100 \
    --n_mole2 50 \
    --density 0.20 \
    --forcefield "PCFF-IFF" \
    --ff_path "frc_files/pcff_interface_v1_6mBN.frc" \
    --output_dir "$WORKING_DIR"

echo "ALL REACTER INPUT FILES CREATED!!"

conda deactivate