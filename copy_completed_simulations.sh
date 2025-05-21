#!/bin/bash

ROOT_DIR="$1"
DEST_DIR="$2"

if [[ -z "$ROOT_DIR" || -z "$DEST_DIR" ]]; then
    echo "Usage: $0 <root_directory> <destination_directory>"
    exit 1
fi

copy_if_completed() {
    local sim_dir="$1"
    local log_file="$sim_dir/log.lammps"

    if [[ -f "$log_file" ]] && grep -q "Total wall time" "$log_file"; then
        # Get relative path from ROOT_DIR to sim_dir
        local rel_path="${sim_dir#$ROOT_DIR/}"
        local target_path="$DEST_DIR/$rel_path"

        mkdir -p "$target_path"
        cp "$log_file" "$target_path/"
        [ -f "$sim_dir/msd_AO_all.txt" ] && cp "$sim_dir/msd_AO_all.txt" "$target_path/"
    fi
}

# Reaction simulations (may include subfolders before Sim-*)
find "$ROOT_DIR" -type d -path "*/Reaction/*/Sim-*" | while read -r sim; do
    copy_if_completed "$sim"
done

# Eq simulations (direct Sim-* under Eq)
find "$ROOT_DIR" -type d -path "*/Eq/Sim-*" | while read -r sim; do
    copy_if_completed "$sim"
done
