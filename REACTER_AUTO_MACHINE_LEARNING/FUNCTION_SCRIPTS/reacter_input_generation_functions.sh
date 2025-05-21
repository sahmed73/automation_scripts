#!/bin/bash

# Here I am using actual LUNAR code. I am not modifying anything

#====== Inputs ====== 
PYTHON="/home/sahmed73/anaconda3/envs/saenv/bin/python3"
LUNAR_DIR="/mnt/borgstore/amartini/sahmed73/data/LUNAR"

# Python codes
PY_SINGLE_SMILES_TO_MOL="PYTHON_SCRIPTS/CL_single_smiles_to_mol.py"
PY_MERGE_MOL="PYTHON_SCRIPTS/CL_merge_into_single_mol_file.py"
PY_DELETE_PHENOLIC_HYDROGEN="PYTHON_SCRIPTS/CL_delete_single_phenolic_hydrogen_from_molfile.py"
PY_BOXSIZE_FROM_DENSITY="PYTHON_SCRIPTS/CL_boxSize_from_density.py"
PY_UPDATE_MAPFILE="PYTHON_SCRIPTS/CL_get_initiator_ids_and_update_mapfile.py"

# Lammps input template files
TEMPLATE_DIR="LAMMPS_INPUT_TEMPLATE"


# define products (post reaction mole1, which is fixed for all simulations)
POST_MOL1_SMILES="CCCCCCCCCCC(OO)(CCCCCCCC)CC(C)CCCCCCCC" # PAO Hydroperoxide

ERROR="false"


run_reacter_in_gen() {
    local this_func_name="${FUNCNAME[0]}"
    echo "Running $this_func_name ..."
    # all arguments are required
    # Arguments for datafile generation
    local mole1
    local mole2               # for both
    local smiles1 smiles2
    local n_mole1 n_mole2
    local density 
    local ff ff_path 
    local output_dir          # for both, # relative to LUNAR direcotry

    # Arguments for LAMMPS input generation
    local runtime Rmax Rprob nframe nsim


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
                echo "[${this_func_name} @ line ${LINENO}] Error: Unknown option '$1'" >&2
                return 1
                ;;
        esac
    done

    # Validate required arguments
    if [[ -z "$mole1" || -z "$mole2" || -z "$smiles1" || -z "$smiles2" ||
          -z "$n_mole1" || -z "$n_mole2" || -z "$density" || -z "$runtime" ||
          -z "$Rmax" || -z "$Rprob" || -z "$ff" || -z "$ff_path" ||
          -z "$nframe" || -z "$nsim" || -z "$output_dir" ]]; then
        echo "[${this_func_name} @ line ${LINENO}] ERROR: Missing one or more required arguments." >&2
        echo "Use -h or --help to see required options." >&2
        return 1
    fi

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

    create_lammps_input_from_template \
    --runtime "$runtime" \
    --Rmax "$Rmax" \
    --Rprob "$Rprob" \
    --nframe "$nframe" \
    --nsim "$nsim" \
    --output_dir "$output_dir" \
    --mole2 "$mole2" \
    --forcefield "$ff"


    ## error handlings
    if [[ "$ERROR" == "true" ]]; then
        mv "${output_dir}/${mole2}" "${output_dir}/ERROR_${mole2}"
    fi
}


create_reacter_datafiles_from_smiles() {
    local this_func_name="${FUNCNAME[0]}"

    # All arguments are required
    local mole1 mole2
    local smiles1 smiles2
    local n_mole1 n_mole2
    local density ff ff_path
    local output_dir

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
                echo "  -fp|--ff_path      Force field path"
                echo "  -o|--output_dir    Output directory"
                return 0
                ;;
            *)
                echo "[${this_func_name} @ line ${LINENO}] ERROR: Unknown option: $1" >&2
                return 1
                ;;
        esac
    done

    # Validation block
    if [[ -z "$mole1" || -z "$mole2" || -z "$smiles1" || -z "$smiles2" ||
          -z "$n_mole1" || -z "$n_mole2" || -z "$density" ||
          -z "$ff" || -z "$ff_path" || -z "$output_dir" ]]; then
        echo "[${this_func_name} @ line ${LINENO}] ERROR: Missing one or more required arguments." >&2
        echo "Use -h or --help to see usage." >&2
        return 1
    fi

    echo ""
    echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%% Processing $mole2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    echo ""

    # Define all directories relative to output_dir and molecule name (mole2 is the antioxidant)
    data_dir="${output_dir}/${mole2}/DataFile"                    # Subdir for data files
    single_molecules_dir="${output_dir}/${mole2}/DataFile/Single_Molecules"  # Where individual .mol files will go

    # Create necessary directories (all parents included with -p)
    mkdir -p "${single_molecules_dir}"

    # Create individual .mol files from SMILES for each molecule
    "$PYTHON" "$PY_SINGLE_SMILES_TO_MOL" "$mole1" "$smiles1" "${single_molecules_dir}"  # Reactant 1 (e.g., PAO-OO)
    "$PYTHON" "$PY_SINGLE_SMILES_TO_MOL" "$mole2" "$smiles2" "${single_molecules_dir}"  # Reactant 2 (e.g., A1 antioxidant)

    # Create post-reaction mole1 mol file from modified SMILES (e.g., peroxy radical after reaction)
    # This will be used to create the post reaction structure, delete afterwards
    "$PYTHON" "$PY_SINGLE_SMILES_TO_MOL" "post_${mole2}" "${POST_MOL1_SMILES}" "${single_molecules_dir}"

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

    pushd "$LUNAR_DIR" > /dev/null || exit 1
    
    # save in an array
    mole_files=("$mole1_mol" "$mole2_mol" "$pre_reaction_mol" "$post_reaction_mol")

    ############################
    #    Run atom_typing.py    #
    ############################
    echo "============================================"
    echo "Running atom_typing.py for $mole2" # for antioxdiant in this case
    echo "============================================"

    atom_typing_dir="${data_dir}/from_LUNAR/from_atom_typing"

    mkdir -p "${atom_typing_dir}" 

    # Run atom_typing.py for each .mol file
    for topo in "${mole_files[@]}"; do
        echo "  -Processing: $topo"
        
        $PYTHON "atom_typing.py" \
        -dir "${atom_typing_dir}" \
        -topo "$topo" \
        -bond "n.u." \
        -charge-file "frc_files/Gasteiger_parameters.txt" \
        -newfile ":_typed" \
        -ff "$ff" \
        -nta-comments F \
        -reset-charges F \
        -vdw-scale 1.1 \
        -pdb skip \
        -boundary p-p-p \
        -bond-reset F > /dev/null
    done

        
    ########################
    #    Run all2lmp.py    #
    ########################
    echo "============================================"
    echo "Running all2lmp.py for each _typed.data file"
    echo "============================================"

    # Define input force field file
    all2lmp_dir="${data_dir}/from_LUNAR/from_all2lmp"

    mkdir -p "${all2lmp_dir}"

    # Loop through all _typed.data files and process them individually
    for topo in "${atom_typing_dir}"/*_typed.data; do
        echo "  -Processing: $topo"
        # Find the corresponding _typed.nta file
        nta_file="${topo%.data}.nta"
        # Run all2lmp.py
        $PYTHON "all2lmp.py" \
            -dir "${all2lmp_dir}" \
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
            -add2box 0.0 > /dev/null
    done

    ##################################
    #     Run bond_react_merge.py    #
    ##################################
    echo "============================================"
    echo "Running bond_react_merge.py"
    echo "============================================"
    
    bond_react_merge_dir="${data_dir}/from_LUNAR/from_bond_react_merge"

    mkdir -p "${bond_react_merge_dir}"

    data1="${all2lmp_dir}/${mole1}_typed_${ff}.data"
    data2="${all2lmp_dir}/${mole2}_typed_${ff}.data"
    pre1="${all2lmp_dir}/pre_reaction_1_typed_${ff}.data"
    post1="${all2lmp_dir}/post_reaction_1_typed_${ff}.data"

    echo "  -data1=${all2lmp_dir}/${mole1}_typed_${ff}.data"
    echo "  -data2=${all2lmp_dir}/${mole2}_typed_${ff}.data"
    echo "  -pre1=${all2lmp_dir}/pre_reaction_1_typed_${ff}.data"
    echo "  -post1=${all2lmp_dir}/post_reaction_1_typed_${ff}.data"

    $PYTHON "bond_react_merge.py" \
        -files data1:"${data1}",data2:"${data2}",pre1:"${pre1}",post1:"${post1}" \
        -dir "${bond_react_merge_dir}" \
        -newfile ':_merged' \
        -atomstyle full \
        -map T \
        -write-rxn-mol2files T \
        -write-moleculefiles T \
        -write-rxn-datafiles F \
        -type-labels F \
        -edge F > /dev/null
    sleep 2


    popd > /dev/null
    ########################
    #  Run cell_bulder.py  #
    ########################
    echo "============================================"
    echo "Running cell_bulder.py"
    echo "============================================"
    
    bulk_dir="${data_dir}/Bulk"

    mkdir -p "${bulk_dir}"

    mole1_typed_ff_merged="${bond_react_merge_dir}/${mole1}_typed_${ff}_merged.data"
    mole2_typed_ff_merged="${bond_react_merge_dir}/${mole2}_typed_${ff}_merged.data"

    file_string="${mole1_typed_ff_merged}:${n_mole1},$mole2_typed_ff_merged:${n_mole2}"

    box_length=$($PYTHON "$PY_BOXSIZE_FROM_DENSITY" "$file_string" "$density")
    box_size="${box_length}Ax${box_length}Ax${box_length}A"
    echo "  -Targey density: $density"
    echo "  -Calculated Box Size: $box_size"
    output_datafile="${n_mole1}_${mole1}_${n_mole2}_${mole2}_density=${density}"

    pushd "$LUNAR_DIR" > /dev/null || exit 1
    $PYTHON "cell_builder.py" \
        -files "$file_string" \
        -dir "$bulk_dir" \
        -newfile "${output_datafile}" \
        -duplicate 1 \
        -dist-scale 1.2 \
        -atomstyle full \
        -type-labels F \
        -reset-molids offset \
        -unwrap T \
        -grp-mono F \
        -seed 12345 \
        -rall 360.0 \
        -domain "$box_size" \
        -boundary p-p-p \
        -maxtry 100 \
        -tolerance 2.0 \
        -mixing tolerance \
        -ff-join none > /dev/null
    
    # check if the cell_builder works correctly
    lunar_logfile=$(find "$bulk_dir" -maxdepth 1 -name "*.log.lunar" | head -n 1)
    if ! [[ -f "$lunar_logfile" ]] || grep -q "WARNING inserted" "$lunar_logfile"; then
        echo "[${this_func_name} @ line ${LINENO}] Error: .log.lunar not found, or perfect packing may unsuccessful... renaming ${mole2} to ERROR_${mole2}" >&2
        ERROR="true"
    fi


    popd > /dev/null
}


create_lammps_input_from_template() {
    local this_func_name="${FUNCNAME[0]}"

    local runtime Rmax Rprob nframe nsim
    local output_dir mole2
    local ff

    # Argument parser
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -r|--runtime) runtime="$2"; shift 2 ;;
            -R|--Rmax) Rmax="$2"; shift 2 ;;
            -rp|--Rprob) Rprob="$2"; shift 2 ;;
            -nf|--nframe) nframe="$2"; shift 2 ;;
            -ns|--nsim) nsim="$2"; shift 2 ;;
            -o|--output_dir) output_dir="$2"; shift 2 ;;
            -m2|--mole2) mole2="$2"; shift 2 ;;
            -ff|--forcefield) ff="$2"; shift 2 ;;
            -h|--help)
                echo "Usage: create_lammps_input_from_template [options]"
                echo "  -r|--runtime       Runtime in nanoseconds"
                echo "  -R|--Rmax          Reaction cutoff distance"
                echo "  -rp|--Rprob        Reaction probability"
                echo "  -nf|--nframe       Number of frames to save"
                echo "  -ns|--nsim         Number of simulations"
                echo "  -o|--output_dir    Output directory"
                echo "  -m2|--mole2        Name of molecule 2 (e.g., antioxidant name)"
                echo "  -ff|--forcefield   Force field name"
                return 0
                ;;
            *)
                echo "[${this_func_name} @ line ${LINENO}] Error: Unknown option: $1" >&2
                return 1
                ;;
        esac
    done

    # Validation block
    if [[ -z "$runtime" || -z "$Rmax" || -z "$Rprob" || -z "$ff" ||
          -z "$nframe" || -z "$nsim" || -z "$output_dir" || -z "$mole2" ]]; then
        echo "[${this_func_name} @ line ${LINENO}] Error: Missing one or more required arguments." >&2
        echo "Use -h or --help to see required options." >&2
        return 1
    fi


    mole_dir="${output_dir}/${mole2}"
	local datafile
	datafile=$(find "${mole_dir}/DataFile/Bulk" -maxdepth 1 -name "*.data" | head -n 1)

	for ((i = 1; i <= nsim; i++)); do

		######## Equilibration ######

		sim_dir="${mole_dir}/Eq/Sim-${i}"

        if [[ ! -d "$sim_dir" ]]; then
    		mkdir -p "$sim_dir"

    		# --- Eq INPUT.IN ---
    		cp "$TEMPLATE_DIR/Eq/Sim-X/input.in" "$sim_dir/input.in"
    		sed -i "s|<<datafile>>|$datafile|g" "$sim_dir/input.in"
    		sed -i "s|<<sim_number>>|$i|g" "$sim_dir/input.in"

    		# --- Eq SUBMIT.SH ---
    		cp "$TEMPLATE_DIR/Eq/Sim-X/submit.sh" "$sim_dir/submit.sh"
    		sed -i "s|<<job-name>>|Eq-S${i}-${mole2}|g" "$sim_dir/submit.sh"
        fi

		######## Reaction ######

		sim_dir="${mole_dir}/Reaction/R=${Rmax}_Rp=${Rprob}_${runtime}ns/Sim-${i}"
        if [[ ! -d "$sim_dir" ]]; then
    		mkdir -p "$sim_dir"

    		# --- Reaction INPUT.IN ---
    		cp "$TEMPLATE_DIR/Reaction/Rmax=X.XX_Rprob=X.XX/Sim-X/input.in" "$sim_dir/input.in"

            # replacing static parameters
            sed -i "s|<<Rmax>>|$Rmax|g" "$sim_dir/input.in"
            sed -i "s|<<Rprob>>|$Rprob|g" "$sim_dir/input.in"
            sed -i "s|<<sim_run_time>>|$runtime|g" "$sim_dir/input.in"
            sed -i "s|<<ff>>|$ff|g" "$sim_dir/input.in"
            sed -i "s|<<nframe>>|$nframe|g" "$sim_dir/input.in"

            # replacing dynamic parameters
    		sed -i "s|<<sim_number>>|$i|g" "$sim_dir/input.in"

    		# --- Reaction SUBMIT.SH ---
    		cp "$TEMPLATE_DIR/Reaction/Rmax=X.XX_Rprob=X.XX/Sim-X/submit.sh" "$sim_dir/submit.sh"
    		sed -i "s|<<job-name>>|RxS${i}-${mole2}_R=${Rmax}_Rp=${Rprob}_${runtime}ns|g" "$sim_dir/submit.sh"
        fi
	done
}