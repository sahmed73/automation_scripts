#!/bin/bash


# ======== functions ========
# 1. job_submission_manager
# 2. submit_job
# 3. check_dependencies
# 4. update_job_dir_array
# 5. check_simulation_status
# 6. check_slurm_job_status
# 7. check_lammps_output_file_status
# 8. job_id_from_working_dir
# 9. generate_job_csv_log
# ===========================


# Helper modules
PYTHON="/home/sahmed73/anaconda3/envs/saenv/bin/python3"

# Static Global Variables
CLUSTERS=("pinnacles" "merced")
USERNAME="sahmed73"

# Dynamic Global Variables
JOB_DIR_ARRAY=()




# ------------------------------------------------------------------------------
# job_submission_manager
#
# Description:
#   Automates the submission of simulation jobs from a given parent directory.
#   Continuously scans for eligible simulations, checks for dependency readiness,
#   submits jobs to SLURM clusters based on defined partition priorities,
#   and tracks submission status until all jobs are processed or declared stuck.
#
# Inputs:
#   $1 - parent_dir: Root directory containing simulation setup folders.
#
# Main Workflow:
#   - Discover all simulation folders containing 'input.in' files.
#   - For each pending simulation:
#       - Check if its dependency file is ready (exists).
#       - Count existing submission attempts (slurm-*.out files).
#       - If dependency is ready and submission attempts are within limits:
#           - Submit the job to SLURM using available partition priorities.
#           - On successful submission or reaching max attempts:
#               - Remove the job from the tracking list.
#   - Monitor SLURM queues across clusters to prevent queue flooding.
#   - If no jobs can be submitted due to unmet dependencies, detect stuck conditions.
#   - Generate and update 'simulation_job_report.csv' to log statuses.
#   - Continue submission cycles until:
#       - All jobs are submitted and completed, OR
#       - Jobs remain stuck for consecutive loops and are declared unprocessable.
#
# Outputs:
#   - Submits jobs using sbatch to SLURM clusters.
#   - Generates 'simulation_job_report.csv' in the parent directory.
#   - Writes detailed submission log with timestamps.
#
# Dependencies (helper functions required):
#   - update_job_dir_array()   # Collects simulation jobs to process
#   - check_dependencies()     # Checks if job dependencies are ready
#   - submit_job()             # Submits job to SLURM partitions
#   - generate_job_csv_log()   # Logs all jobs and statuses in CSV format
#
# Notes:
#   - COMPLETED, RUNNING, and PENDING jobs are skipped. (can be changed)
#   - A maximum of 5 submission attempts are allowed per job directory. (can be changed)
#   - A small sleep is inserted after job submission to avoid scheduler flooding.
#   - Stuck jobs are detected if no progress is made for 3 consecutive loops.
# ------------------------------------------------------------------------------
job_submission_manager() {

    # Parent directory to process (passed as argument)
    local parent_dir="$1"

    # Other local variables
    local job_dir
    local dependency
    local job_submitted
    local max_submission_attempts=5  # Max submission attempts per job
    local log_file="${parent_dir}/job_submission.log"
    local submission_loop_count=0
    local jobs_in_slurm_queue

    # Cluster partition and priority definition
    local partitions_priority="pinnacles.short:0 \
                               pinnacles.long:3 \
                               pinnacles.medium:0 \
                               pinnacles.bigmem:2 \
                               pinnacles.pi.amartini:4 \
                               merced.compute:3 \
                               merced.bigmem:0 \
                               pinnacles.cenvalarc.compute:3 \
                               pinnacles.cenvalarc.bigmem:2"

    # Extract unique cluster names (before first dot) and store in array
    local clusters
    readarray -t clusters < <(printf "%s\n" $partitions_priority | awk -F. '{print $1}' | sort -u)

    # Start logging
    echo "=== Job Submission Log Started at $(date '+%Y-%m-%d %H:%M:%S') ===" > "$log_file"

    # Count how many consecutive loops had stuck jobs (no progress due to dependency issues)
    local consecutive_stuck_loops=0

    # Infinite loop: will break based on exit conditions below
    while :; do
        submission_loop_count=$((submission_loop_count + 1))

        # Update the job list (filter out COMPLETED, RUNNING, PENDING jobs)
        update_job_dir_array "$parent_dir" "Exclude:COMPLETED RUNNING PENDING"
        local simulations_submitted=0

        # Log loop start
        echo "===========================================================================================" >> "$log_file"
        echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting submission loop #$submission_loop_count with $initial_simulation_count simulations pending." >> "$log_file"
        echo "===========================================================================================" >> "$log_file"

        # Count jobs in this loop that have unmet dependencies
        local dependency_not_met_count=0

        # Iterate over all job directories that need submission
        for i in "${!JOB_DIR_ARRAY[@]}"; do
            job_dir="${JOB_DIR_ARRAY[${i}]}"
            dependency=$(check_dependencies "$job_dir")
            local short_path="${job_dir#$parent_dir/}"

            # Check if dependency file exists (i.e., dependencies satisfied)
            if [[ -f "$dependency" ]]; then
                echo "Dependency satisfied for: $short_path" >> "$log_file"

                pushd "$job_dir" > /dev/null

                # Count number of previous submission attempts by checking slurm output files
                local submission_attempt_count
                submission_attempt_count=$(find "$job_dir" -maxdepth 1 -name "slurm-*.out" | wc -l)

                # If max attempts reached, skip this job
                if (( submission_attempt_count >= max_submission_attempts )); then
                    echo "Max submission attempts reached for: $short_path (Attempts: $submission_attempt_count). Skipping ..." >> "$log_file"
                    unset JOB_DIR_ARRAY[i]
                    popd > /dev/null
                    continue
                fi

                # Submit the job using available partitions
                job_submitted=$(submit_job "$job_dir" "$partitions_priority")

                # If submission successful
                if [[ "$job_submitted" == "true" ]]; then
                    simulations_submitted=$((simulations_submitted + 1))
                    echo "${simulations_submitted}. Job submitted successfully: $short_path" >> "$log_file"
                    unset JOB_DIR_ARRAY[i]
                    sleep 1  # Sleep to avoid overloading submission
                fi

                popd > /dev/null
            else
                # Dependency not ready → count it
                dependency_not_met_count=$((dependency_not_met_count + 1))
                echo "Dependency: $dependency" >> "$log_file"
                echo "Dependency not ready for: $short_path" >> "$log_file"
            fi
        done

        # Generate and update CSV status report after this loop
        generate_job_csv_log "$parent_dir"

        # Calculate total jobs running in SLURM queues (sum of all clusters)
        total_jobs_in_slurm_queue=0
        for cluster in "${clusters[@]}"; do
            job_count=$(squeue -M "$cluster" -u "$USERNAME" --noheader 2>/dev/null | grep -c .)

            if [[ $? -ne 0 ]]; then
                echo "Cluster $cluster: Error accessing squeue"
            else
                echo "Cluster $cluster: $job_count jobs"
                total_jobs_in_slurm_queue=$((total_jobs_in_slurm_queue + job_count))
            fi
        done

        # Calculate how many jobs are still left to submit
        local simulations_left=${#JOB_DIR_ARRAY[@]}

        # If no jobs could be submitted in this loop, queue is almost empty → count this as stuck loop
        if [[ "$simulations_left" -eq "$dependency_not_met_count" && "$total_jobs_in_slurm_queue" -le 1 ]]; then
            consecutive_stuck_loops=$((consecutive_stuck_loops + 1))
        fi

        # Log end of this loop
        echo "$(date '+%Y-%m-%d %H:%M:%S') - Loop #$submission_loop_count finished: $simulations_submitted submitted, $simulations_left still pending." >> "$log_file"
        echo "" >> "$log_file"

        # Exit if no simulations left
        if [[ "$simulations_left" -eq 0 ]]; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') - All jobs processed after $submission_loop_count loops. Exiting." >> "$log_file"
            break
        fi

        # Exit if stuck for 3 consecutive loops (nothing progressed)
        if [[ "$consecutive_stuck_loops" -eq 3 ]]; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') - Remaining $simulations_left jobs failed to meet dependencies after three attempts. Exiting." >> "$log_file"
            break
        fi

        # Sleep before the next loop
        echo "Sleeping for 30 sec....." >> "$log_file"
        sleep 30
    done
}





#------------------------------------------------------------------------------
# submit_job
#
# Description:
#   Submits a simulation job by modifying its SLURM submission script to match
#   the available cluster partition (short, medium, long, compute, or pi.amartini),
#   according to real-time queue limits. If partition availability is found,
#   updates the SLURM parameters (partition, core count, wall time) and submits
#   the job using sbatch. Priority is given to partitions in a predefined order.
#
# Inputs:
#   $1 - job_dir: Path to the simulation directory containing 'submit.sh'.
#   $2 - partitions (Optional): Space-separated "partition:limit" pairs
#        (default: "short:12 long:3 medium:6 compute:0 pi.amartini:0").
#
# Key Operations:
#   - Define partition priorities (e.g., short > compute > medium > long > pi.amartini).
#   - Define maximum core counts and wall times for each partition.
#   - Check the current number of running or pending jobs in each partition.
#   - Find an available partition (under the allowed limit).
#   - Update the 'submit.sh' file to:
#       - Set correct partition name.
#       - Set correct number of cores.
#       - Set correct maximum wall time.
#   - Submit the job using sbatch (normal for Pinnacle, `-M merced` for compute partition).
#   - If submission is successful, return 'true'; otherwise 'false'.
#
# Outputs:
#   - Echoes "true" if the job was successfully submitted, "false" otherwise.
#
# Dependencies (Helper Functions Required):
#   - None

# Notes:
#   - Automatically uses dos2unix to fix Windows line endings in 'submit.sh' before editing.
#   - Only submits to a partition if the current job count is below its limit.
#   - Priority order ensures faster queues (like 'short') are attempted first.
#   - Designed for SLURM job management on mixed-cluster environments (Pinnacle and Merced).
#------------------------------------------------------------------------------

submit_job() {
    local job_dir="$1"

    # Default cluster.partition:limit map
    local default_partitions_priority="pinnacles.short:12 \
        pinnacles.long:3 \
        pinnacles.medium:6 \
        pinnacles.bigmem:2 \
        pinnacles.cenvalarc.bigmem:2 \
        pinnacles.cenvalarc.compute:3 \
        merced.compute:6 \
        merced.bigmem:6 \
        pinnacles.pi.amartini:0"

    local partitions="${2:-$default_partitions_priority}"

    # Time and core settings
    declare -A partition_time=(
        ["pinnacles.short"]="6:00:00"
        ["pinnacles.long"]="3-00:00:00"
        ["pinnacles.medium"]="1-00:00:00"
        ["pinnacles.bigmem"]="3-00:00:00"
        ["pinnacles.cenvalarc.bigmem"]="3-00:00:00"
        ["pinnacles.cenvalarc.compute"]="3-00:00:00"
        ["merced.compute"]="5-00:00:00"
        ["merced.bigmem"]="5-00:00:00"
        ["pinnacles.pi.amartini"]="3-00:00:00"
    )

    declare -A partition_cores=(
        ["pinnacles.short"]=48
        ["pinnacles.long"]=48
        ["pinnacles.medium"]=48
        ["pinnacles.bigmem"]=48
        ["pinnacles.cenvalarc.bigmem"]=48
        ["pinnacles.cenvalarc.compute"]=48
        ["merced.compute"]=40
        ["merced.bigmem"]=24
        ["pinnacles.pi.amartini"]=48
    )

    local submit_file="$job_dir/submit.sh"
    local job_submitted=false

    local this_function="${FUNCNAME[0]}"

    if [[ ! -f "$submit_file" ]]; then
        echo "[$this_function] ERROR: Missing submit file: $submit_file" >&2
        echo "false"
        return
    fi

    dos2unix "$submit_file" > /dev/null 2>&1

    local partition_priority=()
    declare -A partition_limit

    for item in $partitions; do
        item=$(echo "$item" | xargs)
        local name=${item%%:*}
        local limit=${item##*:}
        partition_priority+=("$name")
        partition_limit["$name"]="$limit"
    done

    for cluster_partition in "${partition_priority[@]}"; do
        local cluster="${cluster_partition%%.*}"
        local partition="${cluster_partition#*.}"
        local core_count="${partition_cores[$cluster_partition]}"
        local wall_time="${partition_time[$cluster_partition]}"

        if [[ -z "$core_count" || -z "$wall_time" ]]; then
            echo "[$this_function] WARNING: No core/time config for $cluster_partition, skipping." >&2
            continue
        fi

        # Get running+pending job count on this partition
        local count
        count=$(squeue -M "$cluster" --user "$USERNAME" --noheader -o "%i|%P" \
                | awk -F'|' -v p="$partition" '$2 == p' | wc -l)

        if (( count < partition_limit["$cluster_partition"] )); then
            sed -i "/^#SBATCH --partition=/c\\#SBATCH --partition=$partition" "$submit_file"
            sed -i "/^#SBATCH --ntasks-per-node=/c\\#SBATCH --ntasks-per-node=$core_count" "$submit_file"
            sed -i "/^#SBATCH --time=/c\\#SBATCH --time=$wall_time" "$submit_file"

            local output
            output=$(sbatch -M "$cluster" "$submit_file" 2>&1)
            local status=$?

            if [[ $status -eq 0 ]]; then
                job_submitted=true
            else
                echo "[$this_function] ERROR submitting job: $job_dir" >&2
                echo "[$this_function] sbatch error: $output" >&2
                job_submitted=false
            fi
            break
        else
            echo "[$this_function] Partition full: $cluster_partition (current: $count / limit: ${partition_limit[$cluster_partition]})" >&2
        fi
    done

    echo "$job_submitted"
}



#------------------------------------------------------------------------------
# check_dependencies
#
# Description:
#   Checks whether the necessary dependency files for a given simulation job
#   are ready. Runs a Python script (dependency_checker.py) that performs the 
#   actual dependency check logic, with a timeout to avoid hanging indefinitely.
#
# Inputs:
#   $1 - job_dir: Path to the simulation directory to check for dependencies.
#
# Key Operations:
#   - Calls an external Python script ('dependency_checker.py') with job_dir as input.
#   - Limits the Python call to 10 seconds using 'timeout' command.
#   - Captures and returns the output from the Python script.
#
# Outputs:
#   - Echoes the dependency status (file path of satisfied dependency, or empty string if missing).
#
# Dependencies (Helper Functions Required):
#   - None inside Bash.
#   - Requires 'dependency_checker.py' Python script to be available and properly working.
#
# Notes:
#   - If the dependency check Python script exceeds 10 seconds, the command will fail.
#   - It is assumed that dependency_checker.py prints a filepath if successful, or nothing otherwise.
#   - Requires the environment variable $PYTHON to point to the correct Python executable.
#------------------------------------------------------------------------------
check_dependencies() {
    local job_dir="$1"
    local output
    output=$(timeout 10s "$PYTHON" -c "import dependency_checker; print(dependency_checker.check_dependencies('$job_dir'))")
    echo "$output"
}

update_job_dir_array() {
    local parent_dir="$1"
    local except_string="${2:-Exclude:COMPLETED RUNNING PENDING}"
    local status

    echo "Updating the job-dir array..." >> "$log_file"

    JOB_DIR_ARRAY=()  # Clear previous entries

    while IFS= read -r sim_dir; do
        status=$(check_simulation_status "$sim_dir")
        if ! grep -qw "$status" <<< "${except_string#Exclude:}"; then
            JOB_DIR_ARRAY+=("$sim_dir")
        fi
    done < <(find "$parent_dir" -type f -name "input.in" -exec dirname {} \;)

    echo "Update Completed! ${#JOB_DIR_ARRAY[@]} jobs added." >> "$log_file"
}



#------------------------------------------------------------------------------
# check_simulation_status
#
# Description:
#   Determines the current status of a simulation directory based on:
#   (1) the LAMMPS output file and 
#   (2) the SLURM job status. 
#   It intelligently distinguishes between successfully completed, incomplete, failed, or not submitted jobs.
#
# Inputs:
#   $1 - job_dir: The directory path containing simulation files, including 'output.out' and 'slurm-*.out'.
#
# Key Operations:
#   - Checks if 'output.out' exists and contains the string "Total wall time".
#     - If found → marks as 'complete'.
#     - If not found → marks as 'incomplete'.
#     - If missing → marks as 'missing'.
#     - If 'output.out' is not missing:
#     - Finds the latest 'slurm-*.out' file and extracts the job ID.
#     - Uses 'sacct' to query SLURM for the job status.
#     - First tries the default cluster; if no info, tries the 'merced' cluster.
#   - Combines the LAMMPS output check and SLURM status to decide the final status:
#     - 'HARD_FAIL' : SLURM job failed.
#     - 'COMPLETED' : SLURM completed + output file complete.
#     - 'SOFT_FAIL' : SLURM completed + output file incomplete.
#     - Other SLURM statuses returned directly (e.g., PENDING, RUNNING).
#
# Outputs:
#   - Prints the determined simulation status to stdout.
#
# Dependencies (Helper functions Used):
#   - None
#
# Notes:
#   - Assumes job output files follow the naming 'slurm-*.out'.
#   - Assumes 'sacct' is available and configured for both Pinnacle and Merced clusters.
#------------------------------------------------------------------------------
check_simulation_status() {
    local job_dir="$1"

    
    # Check LAMMPS output file
    local completion=$(check_lammps_output_file_status "$job_dir")
    local slurm_status=$(check_slurm_job_status "$job_dir")

    # Decide final simulation status
    if [[ "$slurm_status" == "COMPLETED" && "$completion" == "complete" ]]; then
        echo "COMPLETED"
    elif [[ "$slurm_status" == "COMPLETED" && "$completion" == "incomplete" ]]; then
        echo "SOFT_FAIL"
    else
        echo "$slurm_status"
    fi

    # logging 
    if [[ "$slurm_status" != "COMPLETED" ]]; then
        echo "    completion: $completion, slurm_status: $slurm_status" >> "$log_file"
    fi
}


# check sacct command here
check_slurm_job_status() {  
    local job_dir="$1"
    local job_id
    local slurm_status=""

    job_id=$(job_id_from_working_dir "$job_dir")

    if [[ -n "$job_id" ]]; then
        for cluster in "${CLUSTERS[@]}"; do
            slurm_status=$(sacct -j "$job_id" -M "$cluster" --format=State --noheader --parsable2 \
                            | awk -F'|' 'NR==1{split($1,a," "); gsub(/\+/, "", a[1]); print a[1]}')

            [[ -n "$slurm_status" ]] && break
        done
    fi

    if [[ "$slurm_status" != "COMPLETED" ]]; then
        echo "$job_dir" >> "$log_file"
        echo "    JobID: $job_id" >> "$log_file"
    fi

    echo "$slurm_status"
}



check_lammps_output_file_status() {
    local job_dir="$1"
    local lammps_output_file="$job_dir/output.out"
    local output_file_status

    if [ -f "$lammps_output_file" ]; then
        if grep -q "Total wall time" "$lammps_output_file"; then
            output_file_status="complete"
        else
            output_file_status="incomplete"
        fi
    else
        output_file_status="missing"
    fi

    echo "$output_file_status"
}

job_id_from_working_dir() {
    local job_dir="$1"
    local all_entries=""
    local job_id=""

    for cluster in "${CLUSTERS[@]}"; do
        entries=$(sacct -M "$cluster" --user="$USER" \
                    --start="$(date -d '1 month ago' +%F)" \
                    --format=JobID,Submit,WorkDir \
                    --noheader --parsable2 2>/dev/null \
                | grep -vE '\.batch|\.extern' \
                | awk -F'|' -v wd="$job_dir" '$3 == wd')
        all_entries+=$'\n'"$entries"
    done

    job_id=$(echo "$all_entries" \
             | grep -v '^$' \
             | sort -t '|' -k2,2r \
             | head -n 1 \
             | cut -d '|' -f1)

    echo "$job_id"
}


#------------------------------------------------------------------------------
# generate_job_csv_log
#
# Description:
#   Generates a detailed CSV report summarizing the status and resource usage
#   of all simulation jobs within a given parent directory. It collects
#   information from LAMMPS output files and SLURM job records (using sacct).
#
# Inputs:
#   $1 - parent_dir: The root directory containing simulation setup folders.
#
# Key Operations:
#   - Searches for all simulation directories containing 'input.in' files.
#   - For each simulation:
#     - Checks the simulation status using 'check_simulation_status'.
#     - Extracts job metadata (JobID, Partition, Start, End, etc.) using 'sacct'.
#     - Captures any LAMMPS error message from 'output.out'.
#     - Falls back to default values ('N/A') if fields are unavailable.
#   - Writes the collected information into a CSV file named 'simulation_job_report.csv'.
#
# Outputs:
#   - A CSV file ('simulation_job_report.csv') created inside the parent directory
#     with the following columns:
#     JobID, JobName, Status, Partition, Submit, Start, End, Elapsed,
#     AllocCPUs, NodeList, MaxRSS, ReqMem, SimDir, Error
#
# Dependencies (Helper Functions/Commands Used):
#   - check_simulation_status()   # To determine job state
#
# Notes:
#   - If the job info is missing on the default cluster (Pinnacle), it retries using the Merced cluster.
#   - Empty or missing fields are recorded as 'N/A' to maintain consistency.
#   - Automatically captures and logs any detected LAMMPS runtime errors.
#------------------------------------------------------------------------------
generate_job_csv_log() {
    local parent_dir="$1"
    local csv_file="${parent_dir}/simulation_job_report.csv"
    local temp_csv_file="${parent_dir}/.simulation_job_report_tmp.csv"
    local sim_dir status job_id sacct_output

    echo "JobID,JobName,Status,Partition,Submit,Start,End,Elapsed,AllocCPUs,NodeList,MaxRSS,ReqMem,SimDir,Error" > "$temp_csv_file"

    while IFS= read -r sim_dir; do
        status=$(check_simulation_status "$sim_dir")

        # Default values
        job_id="N/A"; jobname="N/A"; partition="N/A"; submit="N/A"
        start="N/A"; end="N/A"; elapsed="N/A"; alloccpus="N/A"
        nodelist="N/A"; maxrss="N/A"; reqmem="N/A"; error="None"

        # Capture LAMMPS error output, joined in a single line
        if [[ -f "${sim_dir}/output.out" ]]; then
            error=$(grep -i "error" "${sim_dir}/output.out" | tr '\n' ';' || echo "None")
        fi

        if [[ -n "$status" ]]; then
            job_id=$(job_id_from_working_dir "$sim_dir")

            if [[ -n "$job_id" ]]; then
                for cluster in "${CLUSTERS[@]}"; do
                    sacct_output=$(sacct -M "$cluster" -j "$job_id" \
                        --format=JobID,JobName,Start,End,Elapsed,Submit,AllocCPUs,NodeList,MaxRSS,ReqMem,Partition \
                        --noheader --parsable2 | grep -E "^${job_id}\|")
                    if [[ -n "$sacct_output" ]]; then
                        IFS="|" read -r job_id jobname start end elapsed submit alloccpus nodelist maxrss reqmem partition <<< "$sacct_output"
                        break
                    fi
                done
            fi
        fi

        printf '"%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s"\n' \
            "$job_id" "$jobname" "$status" "$partition" "$submit" "$start" "$end" "$elapsed" \
            "$alloccpus" "$nodelist" "$maxrss" "$reqmem" "$sim_dir" "$error" >> "$temp_csv_file"

    done < <(find "$parent_dir" -type f -name "input.in" -exec dirname {} \;)

    mv "$temp_csv_file" "$csv_file"
}



