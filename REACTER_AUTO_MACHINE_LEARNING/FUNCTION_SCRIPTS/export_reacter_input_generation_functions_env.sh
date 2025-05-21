#!/bin/bash

# Load variables and functions
source FUNCTION_SCRIPTS/reacter_input_generation_functions.sh

# Export all variables
export PYTHON
export LUNAR_DIR
export PY_SINGLE_SMILES_TO_MOL
export PY_MERGE_MOL
export PY_DELETE_PHENOLIC_HYDROGEN
export PY_BOXSIZE_FROM_DENSITY
export PY_UPDATE_MAPFILE
export TEMPLATE_DIR
export POST_MOL1_SMILES
export ERROR

# Export necessary functions
export -f run_reacter_in_gen
export -f create_reacter_datafiles_from_smiles
export -f create_lammps_input_from_template
