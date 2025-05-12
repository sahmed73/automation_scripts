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

# source my_functions.sh
source job_submission_functions.sh
source reacter_input_generation_functions.sh

# User Inputs
ROOT_DIR="/mnt/borgstore/amartini/sahmed73/data/Solubility/Solubility_Automation/Solubility005"

#=============================================================================
#                              Submitting Jobs
#=============================================================================

job_submission_manager "${ROOT_DIR}"