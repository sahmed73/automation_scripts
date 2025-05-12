#!/usr/bin/env bash
#SBATCH --job-name=auto_job_submitter
#SBATCH --partition=pi.amartini
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=10-00:00:00
#SBATCH --export=ALL
#SBATCH --mem=60G
#SBATCH --output=automation_log.out
#SBATCH --error=automation_error.out


# User Inputs
root_directory="/mnt/borgstore/amartini/sahmed73/data/REACTER/Batch001"
sleep_time=300

source "/mnt/borgstore/amartini/sahmed73/data/ARLOI/Global_Scripts/functions.sh"
log_file="$root_directory/simulation_log.txt"
mypy="/home/sahmed73/anaconda3/envs/saenv/bin/python3"


echo "Simulation Log - $(date)" > "$log_file"
echo "Root directory: $root_directory" >> "$log_file"
echo "=============================" >> "$log_file"

log_message() {
    local message="$1"
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $message" >> "$log_file"
}

# Priority order: highest to lowest
partition_priority=("short" "compute" "medium" "long")

# Job limits per partition
short_limit=12
long_limit=3
medium_limit=6
compute_limit=5

# Time and core settings
declare -A partition_time=(
    ["short"]="6:00:00"
    ["long"]="72:00:00"
    ["medium"]="24:00:00"
    ["compute"]="120:00:00"
)

declare -A partition_cores=(
    ["short"]=48
    ["long"]=48
    ["medium"]=48
    ["compute"]=32
)


submit_job() {
    local sim_dir="$1"
    local submit_file="$sim_dir/submit.sh"
    local job_submitted=false

    dos2unix "$submit_file"

    for partition in "${partition_priority[@]}"; do
        # Get limit and queue info
        local limit_var="${partition}_limit"
        local limit=${!limit_var}
        local core_count="${partition_cores[$partition]}"
        local wall_time="${partition_time[$partition]}"

        # Get queue count
        if [[ "$partition" == "compute" ]]; then
            count=$(squeue -M merced -u sahmed73 | grep -cw "$partition")
        else
            count=$(squeue -u sahmed73 | grep -cw "$partition")
        fi

        if (( count < limit )); then
            # Update SLURM script
            sed -i "/^#SBATCH --partition=/c\#SBATCH --partition=$partition" "$submit_file"
            sed -i "/^#SBATCH --ntasks-per-node=/c\#SBATCH --ntasks-per-node=$core_count" "$submit_file"
            sed -i "/^#SBATCH --time=/c\#SBATCH --time=$wall_time" "$submit_file"

            # Submit job
            if [[ "$partition" == "compute" ]]; then
                sbatch -M merced "$submit_file" > /dev/null && job_submitted=true
            else
                sbatch "$submit_file" > /dev/null && job_submitted=true
            fi

            sleep 2
            break
        fi
    done

    echo "$job_submitted"
}


submit_job_OLD() {
    local sim_dir="$1"
    local submit_file="$sim_dir/submit.sh"
    local partition=$(grep -oP '(?<=^#SBATCH --partition=).*' "$submit_file" | xargs)
    local job_submitted=false

    if [[ "$partition" == "short" && $(squeue -u sahmed73 | grep -cw 'short') -lt $short_limit ]]; then
	
		# change the partition
		sed -i '/^#SBATCH --partition=/c\#SBATCH --partition=short' "$submit_file"
		sed -i '/^#SBATCH --ntasks-per-node=/c\#SBATCH --ntasks-per-node=48' "$submit_file"
		sed -i '/^#SBATCH --time=/c\#SBATCH --time=6:00:00' "$submit_file"
		
        sbatch "$submit_file" > /dev/null && job_submitted=true
        sleep ${break_time}
		
    elif [[ "$partition" == "short" && $(squeue -u sahmed73 | grep -cw 'long') -lt $long_limit ]]; then
        
		# change the partition
		sed -i '/^#SBATCH --partition=/c\#SBATCH --partition=long' "$submit_file"
		sed -i '/^#SBATCH --ntasks-per-node=/c\#SBATCH --ntasks-per-node=48' "$submit_file"
		sed -i '/^#SBATCH --time=/c\#SBATCH --time=72:00:00' "$submit_file"
		
		sbatch "$submit_file" > /dev/null && job_submitted=true
        sleep ${break_time}
		
    elif [[ "$partition" == "long" && $(squeue -u sahmed73 | grep -cw 'long') -lt $long_limit ]]; then
        sbatch "$submit_file" > /dev/null && job_submitted=true
        sleep ${break_time}
		
    elif [[ $(squeue -M merced -u sahmed73 | grep -cw 'compute') -lt $compute_limit ]]; then
		dos2unix "$sim_dir/submit.merced.sh"
        sbatch -M merced "$sim_dir/submit.merced.sh" > /dev/null && job_submitted=true
        sleep ${break_time}
    fi

    echo "$job_submitted"
}


check_dependencies() {
    local sim_dir="$1"
    local output=$(timeout 10s $mypy -c "import dependency_checker; print(dependency_checker.check_dependencies('$sim_dir'))")
    echo "$output"
}


process_simulation_dirs() {
    # Initialize sim_dirs with simulations that have not yet completed
    mapfile -t sim_dirs < <(find "$1" -type f -name "input.in" -exec dirname {} \; | while read -r sim_dir; do
        check_simulation_status "$sim_dir"
		status=$?
        [[ $status -ne 1 ]] && echo "$sim_dir"
    done)

    loop=1
    echo "Number of simulations: ${#sim_dirs[@]}"
    while [ ${#sim_dirs[@]} -gt 0 ]; do
		
        log_message "Submission loop: ${loop}"
        echo "-----------------------------------------------------------" >> "$log_file"
        ((loop++))

        echo "Current simulations in the queue:" >> "$log_file"
        count=1 
        for dir in "${sim_dirs[@]}"; do
            echo "$count. ${dir#"$root_directory/"}" >> "$log_file"
            ((count++))
        done

        local indices_to_remove=()

        for i in "${!sim_dirs[@]}"; do
            echo ""
            echo ">>> Processing sim_dirs[$i] = ${sim_dirs[i]}"
            
            local sim_dir="${sim_dirs[i]}"
            
            dependency=$(check_dependencies "$sim_dir")
            echo ">>> Dependency: $dependency"
            
            if [[ -f "$dependency" ]]; then
                echo ">>> Valid dependency file found."
                
                pushd "$sim_dir" > /dev/null
                
                local job_submitted=$(submit_job "$sim_dir")
                job_submitted=$(echo "$job_submitted" | tr -d '\r')  # clean CR if any
                
                echo ">>> Job submitted flag: '$job_submitted'"
                
                if [[ "$job_submitted" == "true" ]]; then
                    log_message "Job submitted: '${sim_dir#"$root_directory/"}'"
                    indices_to_remove+=("$i")
                fi
                
                popd > /dev/null
            else
                echo ">>> Skipping: Dependency file does not exist."
            fi

            echo ">>> Reached end of inner loop for index $i"
            echo ""
        done


        echo ""
        echo "============================"
        echo "2. I am here! Bro!"
        echo "============================"
        echo ""

        # Remove submitted jobs from sim_dirs by creating a filtered list
        for index in "${indices_to_remove[@]}"; do
            echo ""
            echo "============================"
            echo "3. I am here! Bro!"
            echo "============================"
            echo ""

            unset 'sim_dirs[index]'
        done
		
        sim_dirs=("${sim_dirs[@]}")

        
        echo "Waiting ${sleep_time} secs..." >> "$log_file"
		echo "" >> "$log_file"
        sleep $((sleep_time * 1))
    done
}

process_simulation_dirs "${root_directory}"

log_message "Script completed."
