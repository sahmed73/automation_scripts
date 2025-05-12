#!/usr/bin/env bash
#SBATCH --job-name=REACTER_AUTOMATION
#SBATCH --partition=pi.amartini
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=10-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --export=ALL
#SBATCH --mem=60G
#SBATCH --output=automation_log.out
#SBATCH --error=automation_error.out

# activate conda environment
source ~/.bashrc
conda activate saenv

# source my_functions.sh
source job_submission_functions.sh
source reacter_input_generation_functions.sh

# User Inputs
ROOT_DIR="/mnt/borgstore/amartini/sahmed73/data/REACTER"

#=============================================================================
#                          Create Simulations Inputs
#=============================================================================

# WORKING_DIR=$(create_next_folder "$ROOT_DIR" "Batch")

WORKING_DIR="${ROOT_DIR}/12_AO_20ns"
mkdir "${WORKING_DIR}" || { echo "Directory already exists. Exiting."; exit 1; }


CSV_FILE="/mnt/borgstore/amartini/sahmed73/data/SMILES/12_AO_SMILES.csv"
cp -n "$CSV_FILE" "$WORKING_DIR/"

# Export working dir and functions so parallel can access them
export WORKING_DIR
export -f create_reacter_inputs_from_smiles
export -f create_lammps_input_from_template
export -f create_next_folder

# Export all variables
export SCRIPT_DIR
export PYTHON
export LUNAR_DIR
export TEMPLATE_DIR
export PRE_PAO_RADICAL_MOL
export POST_PAO_RADICAL_MOL
export PAO_RADICAL_TYPED_PCFF
export PAO_RADICAL_MOL
export FORCE_FIELD
export PY_MERGE_MOL
export PY_DELETE_PHENOLIC_HYDROGEN
export PY_BOXSIZE_FROM_DENSITY
export PY_UPDATE_MAPFILE
export PY_ALL2LMP_CUSTOMS
export TARGET_DENSITY
export N_MOL1
export N_MOL2
export RXN_N_SIM

# Export associative array RXN_PARAM
export RXN_PARAM_EXPORT="$(declare -p RXN_PARAM)"

# Determine number of cores (fallback to 4 if not set)
CORES=${SLURM_NTASKS_PER_NODE:-4}

# Count total number of jobs
TOTAL=$(awk -F',' 'NR>1 && NF==2' "$CSV_FILE" | wc -l)

echo "Running parallel with ${CORES} cores... Total jobs: ${TOTAL}"

# Run parallel processing (skip header row)
awk -F',' 'NR>1 && NF==2' "$CSV_FILE" | parallel -j "$CORES" --colsep ',' "
    source <(printf '%s' \"\$RXN_PARAM_EXPORT\")
    create_reacter_inputs_from_smiles {1} {2} \"$WORKING_DIR\"
"

echo "ALL REACTER INPUT FILES CREATED!!"


#=============================================================================
#                              Submitting Jobs
#=============================================================================

job_submission_manager "${WORKING_DIR}"


# deactivate conda env
conda deactivate