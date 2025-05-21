#!/usr/bin/env bash
#SBATCH --job-name=auto_job_submitter
#SBATCH --partition=pi.amartini
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=10-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --export=ALL
#SBATCH --mem=10G
#SBATCH --output=automation_log.out
#SBATCH --error=automation_error.out

# source my_functions.sh
source job_submission_functions.sh

# User Inputs
ROOT_DIR="/mnt/borgstore/amartini/sahmed73/data/REACTER/12_AO_20ns"

#=============================================================================
#                              Submitting Jobs
#=============================================================================

job_submission_manager "${ROOT_DIR}"
