#!/bin/bash

# Here I am using actual LUNAR code. I am not modifying anything

#====== Inputs ====== 
SCRIPT_DIR="/mnt/borgstore/amartini/sahmed73/data/automation_scripts/REACTER_AUTO_PEROXY_RADICAL"
PYTHON="/home/sahmed73/anaconda3/envs/saenv/bin/python3"
LUNAR_DIR="/mnt/borgstore/amartini/sahmed73/data/LUNAR"

TEMPLATE_DIR="${SCRIPT_DIR}/lammps_input_templates"
FORCE_FIELD="PCFF-IFF"

# Python codes
PY_MERGE_MOL="${SCRIPT_DIR}/CL_merge_into_single_mol_file.py"
PY_DELETE_PHENOLIC_HYDROGEN="${SCRIPT_DIR}/CL_delete_single_phenolic_hydrogen_from_molfile.py"
PY_BOXSIZE_FROM_DENSITY="${SCRIPT_DIR}/CL_boxSize_from_density.py"
PY_UPDATE_MAPFILE="${SCRIPT_DIR}/CL_get_initiator_ids_and_update_mapfile.py"


# define products (post reaction mole1, which is fixed for all simulations)
POST_MOL1_SMILES="CCCCCCCCCCC(OO)(CCCCCCCC)CC(C)CCCCCCCC" # PAO Hydroperoxide


run_reacter_in_gen() {
    # Default values
    local mole1="" mole2=""
    local smiles1="" smiles2=""
    local n_mole1=0 n_mole2=0
    local density=0
    local runtime=0
    local Rmax=0 Rprob=0
    local ff=""
    local nframe=0 nsim=0
    local output_dir=""

    # Argument parser
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -m1|--mole1) mole1="$2"; shift 2 ;;
            -m2|--mole2) mole2="$2"; shift 2 ;;
            -s1|--smiles1) smiles1="$2"; shift 2 ;;
            -s2|--smiles2) smiles2="$2"; shift 2 ;;
            -n1|--n_mole1) n_mole1="$2"; shift 2 ;;
            -n2|--n_mole2) n_mole2="$2"; shift 2 ;;
            -d|--density) density="$2"; shift 2 ;;
            -r|--runtime) runtime="$2"; shift 2 ;;
            -R|--Rmax) Rmax="$2"; shift 2 ;;
            -rp|--Rprob) Rprob="$2"; shift 2 ;;
            -ff|--forcefield) ff="$2"; shift 2 ;;
            -fp|--ff_path) ff_path="$2"; shift 2 ;;
            -nf|--nframe) nframe="$2"; shift 2 ;;
            -ns|--nsim) nsim="$2"; shift 2 ;;
            -o|--output_dir) output_dir="$2"; shift 2 ;;
            -h|--help)
                echo "Usage: create_reacter_datafiles_from_smiles [options]"
                echo "  -m1|--mole1        Molecule 1 name"
                echo "  -m2|--mole2        Molecule 2 name"
                echo "  -s1|--smiles1      SMILES string for molecule 1"
                echo "  -s2|--smiles2      SMILES string for molecule 2"
                echo "  -n1|--n_mole1      Number of molecules of type 1"
                echo "  -n2|--n_mole2      Number of molecules of type 2"
                echo "  -d|--density       Density in g/cm³"
                echo "  -r|--runtime       Runtime in ns"
                echo "  -R|--Rmax          Reaction cutoff distance"
                echo "  -rp|--Rprob        Reaction probability"
                echo "  -ff|--forcefield   Force field name"
                echo "  -fp|--ff_path      Force filed path"
                echo "  -nf|--nframe       Number of frames"
                echo "  -ns|--nsim         Number of simulations"
                return 0
                ;;
            *)
                echo "Unknown option: $1"
                return 1
                ;;
        esac
    done


    create_reacter_datafiles_from_smiles \
        --mole1 "$mole1" \
        --smiles1 "$smiles1" \
        --mole2 "$mole2" \
        --smiles2 "$smiles2" \
        --n_mole1 "$n_mole1" \
        --n_mole2 "$n_mole2" \
        --density "$density" \
        --forcefield "$ff" \
        --ff_path "$ff_path" \
        --output_dir "$output_dir"

}





create_reacter_datafiles_from_smiles() {
    # Default values
    local mole1="" mole2=""
    local smiles1="" smiles2=""
    local n_mole1=0 n_mole2=0
    local density=0
    local ff=""
    local fp=""
    local output_dir=""

    # Argument parser
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -m1|--mole1) mole1="$2"; shift 2 ;;
            -m2|--mole2) mole2="$2"; shift 2 ;;
            -s1|--smiles1) smiles1="$2"; shift 2 ;;
            -s2|--smiles2) smiles2="$2"; shift 2 ;;
            -n1|--n_mole1) n_mole1="$2"; shift 2 ;;
            -n2|--n_mole2) n_mole2="$2"; shift 2 ;;
            -d|--density) density="$2"; shift 2 ;;
            -ff|--forcefield) ff="$2"; shift 2 ;;
            -fp|--ff_path) ff_path="$2"; shift 2 ;;
            -o|--output_dir) output_dir="$2"; shift 2 ;;
            -h|--help)
                echo "Usage: create_reacter_datafiles_from_smiles [options]"
                echo "  -m1|--mole1        Molecule 1 name"
                echo "  -m2|--mole2        Molecule 2 name"
                echo "  -s1|--smiles1      SMILES string for molecule 1"
                echo "  -s2|--smiles2      SMILES string for molecule 2"
                echo "  -n1|--n_mole1      Number of molecules of type 1"
                echo "  -n2|--n_mole2      Number of molecules of type 2"
                echo "  -d|--density       Density in g/cm³"
                echo "  -ff|--forcefield   Force field name e.g. PCFF, PCFF-IFF"
                echo "  -fp|--ff_path      Force filed path"
                echo "  -o|--output_dir    Output directory"
                return 0
                ;;
            *)
                echo "Unknown option: $1"
                return 1
                ;;
        esac
    done

    # Define all directories relative to output_dir and molecule name (mole2 is the antioxidant)
    mole_dir="${output_dir}/${mole2}"                             # Directory for this specific molecule
    data_dir="${output_dir}/${mole2}/DataFile"                    # Subdir for data files
    single_molecules_dir="${output_dir}/${mole2}/DataFile/Single_Molecules"  # Where individual .mol files will go

    # Create necessary directories (all parents included with -p)
    mkdir -p "${single_molecules_dir}"

    # Create individual .mol files from SMILES for each molecule
    "$PYTHON" CL_single_smiles_to_mol.py "$mole1" "$smiles1" "${single_molecules_dir}"  # Reactant 1 (e.g., PAO-OO)
    "$PYTHON" CL_single_smiles_to_mol.py "$mole2" "$smiles2" "${single_molecules_dir}"  # Reactant 2 (e.g., A1 antioxidant)

    # Create post-reaction mole1 mol file from modified SMILES (e.g., peroxy radical after reaction)
    # This will be used to create the post reaction structure, delete afterwards
    "$PYTHON" CL_single_smiles_to_mol.py "post_${mole2}" "${POST_MOL1_SMILES}" "${single_molecules_dir}"

    # Save the full file paths of the individual molecules
    mole1_mol="${single_molecules_dir}/${mole1}.mol"
    mole2_mol="${single_molecules_dir}/${mole2}.mol"
    post_mole1_mol="${single_molecules_dir}/post_${mole2}.mol"

    # Create pre-reaction combined .mol file
    # Assumes there's only one reactive site (e.g., one phenolic OH group)
    pre_reaction_mol="${single_molecules_dir}/pre_reaction_1.mol"
    "$PYTHON" "${PY_MERGE_MOL}"  "${mole1_mol}" "${mole2_mol}" "${pre_reaction_mol}"

    # Create post-reaction combined .mol file using the same antioxidant structure (needs to be changed)
    post_reaction_mol="${single_molecules_dir}/post_reaction_1.mol"
    "$PYTHON" "${PY_MERGE_MOL}"  "${post_mole1_mol}" "${mole2_mol}" "${post_reaction_mol}"

    # Remove the phenolic hydrogen from the antioxidant in the *post-reaction* structure
    # Script must be able to identify phenolic OH and remove *only one* hydrogen
    "$PYTHON" "${PY_DELETE_PHENOLIC_HYDROGEN}" "${post_reaction_mol}"

    # Cleanup temporary file (the intermediate post-antioxidant molecule used only for merging)
    rm "$post_mole1_mol"

    
    #=============================================================================
    #                                 LUNAR
    #=============================================================================
    
    # save in an array
    mole_files=("$mole1_mol" "$mole2_mol" "$pre_reaction_mol" "$post_reaction_mol")

    ############################
    #    Run atom_typing.py    #
    ############################
    echo "Running atom_typing.py for $mole2" # for antioxdiant in this case

    pushd "$LUNAR_DIR" > /dev/null || exit 1

    atom_typing_dir="${data_dir}/from_LUNAR/from_atom_typing"

    mkdir -p "${ATOM_TYPING_DIR}" # do i need this?

    # Run atom_typing.py for each .mol file
    for topo in "${mole_files[@]}"; do
        echo "Processing: $topo"
        
        $PYTHON atom_typing.py \
        -dir "${atom_typing_dir}" \
        -topo "$topo" \
        -bond "n.u." \
        -charge-file "frc_files/Gasteiger_parameters.txt" # no-need \
        -newfile ":_typed" \
        -ff "$ff" \
        -nta-comments F \
        -reset-charges F \
        -vdw-scale 1.1 \
        -pdb skip \
        -boundary p-p-p \
        -bond-reset F
    done

        
    ########################
    #    Run all2lmp.py    #
    ########################
    echo "Running all2lmp.py for each _typed.data file"

    # Define input force field file
    all2lmp_dir="${data_dir}/from_LUNAR/from_all2lmp"

    mkdir -p "${all2lmp_dir}"

    # Loop through all _typed.data files and process them individually
    for topo in "${atom_typing_dir}"/*_typed.data; do

        # Find the corresponding _typed.nta file
        nta_file="${topo%.data}.nta"

        echo "Processing: $topo with all2lmp.py"

        # Run all2lmp.py
        $PYTHON all2lmp.py -dir "${all2lmp_dir}" \
          -class 2 \
          -frc "$ff_path" \
          -newfile ":_${ff}" \
          -topo "$topo" \
          -nta "$nta_file" \
          -ignore T \
          -type-labels F \
          -reset-charges T \
          -atomstyle full \
          -auto-equivs T \
          -asm "frc_files/general_assumed_equivs.coeffs" \
          -assumed-equivs F \
          -reset-molids T \
          -write-comments T \
          -write-bond-react F \
          -morse-bond F \
          -rx 0 -ry 0 -rz 0 \
          -sx 0.0 -sy 0.0 -sz 0.0 \
          -add2box 0.0
        sleep 2
    done

    ##################################
    #     Run bond_react_merge.py    #
    ##################################
    
    bond_react_merge_dir="${data_dir}/from_LUNAR/from_bond_react_merge"

    mkdir -p "${bond_react_merge_dir}"

    data1="${all2lmp_dir}/${mole1}_typed_${ff}.data"
    data2="${all2lmp_dir}/${mole2}_typed_${ff}.data"
    pre1="${all2lmp_dir}/pre_reaction_1_typed_${ff}.data"
    post1="${all2lmp_dir}/post_reaction_1_typed_${ff}.data"

    $PYTHON bond_react_merge.py \
        -files data1:"${data1}",data2:"${data2}",pre1:"${pre1}",post1:"${post1}" \
        -dir "${bond_react_merge_dir}" \
        -newfile ':_merged' \
        -atomstyle full \
        -map T \
        -write-rxn-mol2files T \
        -write-moleculefiles T \
        -write-rxn-datafiles F \
        -type-labels F \
        -edge F
    sleep 2



    ########################
    #  Run cell_bulder.py  #
    ########################
    echo "Running cell_bulder.py"
    
    bulk_dir="${data_dir}/Bulk"

    mkdir -p "${bulk_dir}"

    mole1_typed_ff_merged="${bond_react_merge_dir}/${mole1}_typed_${ff}_merged.data"
    mole2_typed_ff_merged="${bond_react_merge_dir}/${mole2}_typed_${ff}_merged.data"

    file_string="${mole1_typed_ff_merged}:${n_mole1},$mole2_typed_ff_merged:${n_mole2}"

    box_length=$($PYTHON "$PY_BOXSIZE_FROM_DENSITY" "$file_string" "$density")
    box_size="${box_length}Ax${box_length}Ax${box_length}A"
    echo "Calculated Box Size: $box_size"
    output_datafile="${n_mole1}_${mole1}_${n_mole2}_${mole2}_density=${density}"

    $PYTHON cell_builder.py \
        -files "$file_string" \
        -dir "$bulk_dir" \
        -newfile "${output_datafile}" \
        -duplicate 1 \
        -dist-scale 1.2 \
        -atomstyle full \
        -type-labels F \
        -reset-molids clusters \
        -unwrap T \
        -grp-mono F \
        -seed 12345 \
        -rall 360.0 \
        -domain "$box_size" \
        -boundary p-p-p \
        -maxtry 100 \
        -tolerance 2.0 \
        -mixing tolerance \
        -ff-join none
    sleep 2

    popd > /dev/null
}


create_lammps_input_from_template() {
	local parent_dir="$1"       # Expects subdirectories: DataFile, Eq, Reaction
	local template_dir="$2"
	local number_of_simulations="$3"
	local antioxidant_name
	local sim_dir
	antioxidant_name=$(basename "$parent_dir")

	local datafile
	datafile=$(find "$parent_dir/DataFile/Bulk" -maxdepth 1 -name "*.data" | head -n 1)

	for ((i = 1; i <= number_of_simulations; i++)); do

		######## Equilibration ######

		sim_dir="$parent_dir/Eq/Sim-${i}"

        if [[ ! -d "$sim_dir" ]]; then
    		mkdir -p "$sim_dir"

    		# --- Eq INPUT.IN ---
    		cp "$template_dir/Eq/Sim-X/input.in" "$sim_dir/input.in"
    		sed -i "s|<<datafile>>|$datafile|g" "$sim_dir/input.in"
    		sed -i "s|<<sim_number>>|$i|g" "$sim_dir/input.in"

    		# --- Eq SUBMIT.SH ---
    		cp "$template_dir/Eq/Sim-X/submit.sh" "$sim_dir/submit.sh"
    		sed -i "s|<<job-name>>|Eq-S${i}-${antioxidant_name}|g" "$sim_dir/submit.sh"
        fi

		######## Reaction ######

		sim_dir="$parent_dir/Reaction/R=${RXN_PARAM[Rmax]}_Rp=${RXN_PARAM[Rprob]}_${RXN_PARAM[sim_run_time]}ns/Sim-${i}"
        if [[ ! -d "$sim_dir" ]]; then
    		mkdir -p "$sim_dir"

    		# --- Reaction INPUT.IN ---
    		cp "$template_dir/Reaction/Rmax=X.XX_Rprob=X.XX/Sim-X/input.in" "$sim_dir/input.in"

            # replacing static parameters
            for parameter in "${!RXN_PARAM[@]}"; do
                value="${RXN_PARAM[$parameter]}"
                sed -i "s|<<$parameter>>|$value|g" "$sim_dir/input.in"
            done

            # replacing dynamic parameters
    		sed -i "s|<<sim_number>>|$i|g" "$sim_dir/input.in"

    		# --- Reaction SUBMIT.SH ---
    		cp "$template_dir/Reaction/Rmax=X.XX_Rprob=X.XX/Sim-X/submit.sh" "$sim_dir/submit.sh"
    		sed -i "s|<<job-name>>|RxS${i}-${antioxidant_name}_R=${RXN_PARAM[Rmax]}_Rp=${RXN_PARAM[Rprob]}_${RXN_PARAM[sim_run_time]}ns|g" "$sim_dir/submit.sh"
        fi
	done
}