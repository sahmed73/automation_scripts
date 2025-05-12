#!/bin/bash

parent_dir="/mnt/borgstore/amartini/sahmed73/data/REACTER/12_AO_20ns"
log_file="$parent_dir/directorywise_completion.log"

# Initialize log file with timestamp
echo "$(date '+%Y-%m-%d %H:%M:%S')" > "$log_file"
echo "" >> "$log_file"

# Initialize directory-level and simulation-level counters
folder_complete_count=0
folder_incomplete_count=0

count_complete=0
count_incomplete=0
count_missing=0

for antioxidant_dir in "$parent_dir"/*; do
    [[ -d "$antioxidant_dir" ]] || continue  # Skip non-directory items

    antioxidant_name=$(basename "$antioxidant_dir")
    status_overall="complete"
    sub_log=""

    while IFS= read -r sim_dir; do
        outfile="$sim_dir/output.out"
        rel_path="${sim_dir#$parent_dir/}"  # Relative path for readability

        if [[ -f "$outfile" ]]; then
            if grep -q "Total wall time" "$outfile"; then
                sim_status="complete"
                ((count_complete++))
            else
                sim_status="incomplete"
                ((count_incomplete++))
                status_overall="incomplete"
            fi
        else
            sim_status="missing"
            ((count_missing++))
            status_overall="incomplete"
        fi

        sub_log+="    $rel_path: $sim_status"$'\n'

    done < <(find "$antioxidant_dir" -type f -name "input.in" -exec dirname {} \;)

    # Count the folder-level status
    if [[ "$status_overall" == "complete" ]]; then
        ((folder_complete_count++))
    else
        ((folder_incomplete_count++))
    fi

    {
        echo "$antioxidant_name: $status_overall"
        echo "$sub_log"
        echo "----------------------------------------"
    } >> "$log_file"
done

# Append summary to log
{
    echo ""
    echo "========== Summary =========="
    echo "Simulation status:"
    echo "  Complete simulations:   $count_complete"
    echo "  Incomplete simulations: $count_incomplete"
    echo "  Missing simulations:    $count_missing"
    echo ""
    echo "Directory status:"
    echo "  Folders marked complete:   $folder_complete_count"
    echo "  Folders marked incomplete: $folder_incomplete_count"
    echo "============================="
} >> "$log_file"

echo "Full summary written to: $log_file"
