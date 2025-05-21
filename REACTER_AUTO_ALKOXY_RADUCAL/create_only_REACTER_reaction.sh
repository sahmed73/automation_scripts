#!/bin/bash

#====== Inputs ======
PARENT_DIR="/mnt/borgstore/amartini/sahmed73/data/REACTER/12_AO_20ns"
SCRIPT_DIR="/mnt/borgstore/amartini/sahmed73/data/automation_scripts/REACTER_AUTO"
PYTHON="/home/sahmed73/anaconda3/envs/saenv/bin/python3"
LUNAR_DIR="/mnt/borgstore/amartini/sahmed73/data/LUNAR"
TEMPLATE_DIR="${SCRIPT_DIR}/lammps_input_templates"

# Static parameters
declare -A RXN_PARAM
RXN_PARAM["sim_run_time"]="30"
RXN_PARAM["Rmax"]="1.90"
RXN_PARAM["Rprob"]="0.1"
RXN_PARAM["ff"]="PCFF-IFF"
RXN_PARAM["nframe"]="4000"

RXN_N_SIM=5  # Number of reaction simulations

#====== Function ======
create_lammps_input_from_template() {
    local parent_dir="$1"
    local template_dir="$2"
    local number_of_simulations="$3"

    echo "Creating LAMMPS inputs for: $parent_dir"
    echo "Template directory: $template_dir"
    echo "Simulations to generate: $number_of_simulations"
    
    local antioxidant_name
    antioxidant_name=$(basename "$parent_dir")
    echo "Antioxidant name: $antioxidant_name"

    for ((i = 1; i <= number_of_simulations; i++)); do
        local sim_dir="$parent_dir/Reaction/R=${RXN_PARAM[Rmax]}_Rp=${RXN_PARAM[Rprob]}_${RXN_PARAM[sim_run_time]}ns/Sim-${i}"

        if [[ -d "$sim_dir" ]]; then
            echo "Skipping Sim-${i}: already exists"
            continue
        fi

        mkdir -p "$sim_dir"

        # Input.in
        local input_template="$template_dir/Reaction/Rmax=X.XX_Rprob=X.XX/Sim-X/input.in"
        if [[ ! -f "$input_template" ]]; then
            echo "ERROR: Missing input template: $input_template"
            continue
        fi

        cp "$input_template" "$sim_dir/input.in"

        for parameter in "${!RXN_PARAM[@]}"; do
            value="${RXN_PARAM[$parameter]}"
            sed -i "s|<<$parameter>>|$value|g" "$sim_dir/input.in"
        done

        sed -i "s|<<sim_number>>|$i|g" "$sim_dir/input.in"

        # submit.sh
        local submit_template="$template_dir/Reaction/Rmax=X.XX_Rprob=X.XX/Sim-X/submit.sh"
        if [[ ! -f "$submit_template" ]]; then
            echo "ERROR: Missing submit template: $submit_template"
            continue
        fi

        cp "$submit_template" "$sim_dir/submit.sh"
        sed -i "s|<<job-name>>|RxS${i}-${antioxidant_name}_R=${RXN_PARAM[Rmax]}_Rp=${RXN_PARAM[Rprob]}_${RXN_PARAM[sim_run_time]}ns|g" "$sim_dir/submit.sh"
    done

    echo "All requested simulations processed."
}


for AO_DIR in "$PARENT_DIR"/*; do
    if [[ -d "$AO_DIR/Reaction" ]]; then
        create_lammps_input_from_template "$AO_DIR" "$TEMPLATE_DIR" "$RXN_N_SIM"
    fi
done

