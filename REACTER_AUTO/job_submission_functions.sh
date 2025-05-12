#!/bin/bash


# ======== functions ========
# 1. job_submission_manager
# 2. submit_job
# 3. check_dependencies
# 4. update_job_dir_array
# 5. check_simulation_status
# 6. generate_job_csv_log
# ===========================


# Helper modules
PYTHON="/home/sahmed73/anaconda3/envs/saenv/bin/python3"

# Global variables
JOB_DIR_ARRAY=()
USERNAME="sahmed73"


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
        update_job_dir_array "$parent_dir" "Except:COMPLETED RUNNING PENDING"
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
                    sleep 2  # Sleep to avoid overloading submission
                fi

                popd > /dev/null
            else
                # Dependency not ready → count it
                dependency_not_met_count=$((dependency_not_met_count + 1))
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
        echo "Sleeping for 1 min....." >> "$log_file"
        sleep 60
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
    # cluster_name.partition_name:max_job_submit_limit
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

    
    local submit_file="$job_dir/submit.sh"
    local job_submitted=false

    # Priority order: highest to lowest
    local partition_priority=()
    declare -A partition_limit

    # Fill arrays
    for item in $partitions; do   # cannot quotes the partitions since looping
        item=$(echo "$item" | xargs)  # Trim leading/trailing whitespace
        local name=${item%%:*}
        local limit=${item##*:}

        partition_priority+=("$name")
        partition_limit["$name"]="$limit"
    done


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
        ["pinnacles.pi.amartini"]="3-00:00:00" # pi.amartini has no time limit
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

    dos2unix "$submit_file" > /dev/null 2>&1

    for cluster_partition in "${partition_priority[@]}"; do

        cluster=$(cut -d'.' -f1 <<< "$cluster_partition") # Cluster is everything before the first dot
        partition=${cluster_partition#*.} # Partition is everything after the first dot

        # Get limit and queue info
        local core_count="${partition_cores[$cluster_partition]}"
        local wall_time="${partition_time[$cluster_partition]}"

        # Get queue count
        local count
        count=$(squeue -M "$cluster" --user "$USERNAME" -o "%.18i %.30P" | awk -v p="$partition" '$2 == p' | wc -l)

        if (( count < "${partition_limit[$cluster_partition]}" )); then
            # Update SLURM script
            sed -i "/^#SBATCH --partition=/c\\#SBATCH --partition=$partition" "$submit_file"
            sed -i "/^#SBATCH --ntasks-per-node=/c\\#SBATCH --ntasks-per-node=$core_count" "$submit_file"
            sed -i "/^#SBATCH --time=/c\\#SBATCH --time=$wall_time" "$submit_file"

            # Submit job
            sbatch -M "$cluster" "$submit_file" > /dev/null && job_submitted=true

            break
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


#------------------------------------------------------------------------------
# update_job_dir_array
#
# Description:
#   Updates the global JOB_DIR_ARRAY with simulation directories that are 
#   eligible for submission. It filters out directories whose simulation status 
#   matches any of the excluded statuses provided.
#
# Inputs:
#   $1 - parent_dir: The root directory where all simulation subfolders are located.
#   $2 - except_string (Optional): Space-separated statuses to exclude.
#        Default: "Except:COMPLETED RUNNING PENDING"
#
# Key Operations:
#   - Finds all simulation directories that contain an 'input.in' file.
#   - Checks the current status of each simulation by calling check_simulation_status().
#   - Excludes directories whose status matches any provided in except_string.
#   - Populates the JOB_DIR_ARRAY with only eligible simulation directories.
#
# Outputs:
#   - Updates the global array variable JOB_DIR_ARRAY.
#
# Dependencies (Helper Functions Required):
#   - check_simulation_status()   # Determines the status of a simulation directory
#
# Notes:
#   - The exclusion statuses are space-separated after the prefix "Except:".
#   - The JOB_DIR_ARRAY is cleared and repopulated fresh each time this function runs.
#   - Only simulation directories containing 'input.in' files are considered.
#------------------------------------------------------------------------------
# update_job_dir_array() {
# 	# except_string="Except:COMPLETED RUNNING PENDING" # space separted
#     local parent_dir="$1"
#     local except_string="${2:-Except:COMPLETED RUNNING PENDING}"

#     JOB_DIR_ARRAY=()  # Clear previous entries

#     while IFS= read -r sim_dir; do
#         local status
#         status=$(check_simulation_status "$sim_dir")

#         if [[ ! " ${except_string#Except:} " =~ " $status " ]]; then
#             JOB_DIR_ARRAY+=("$sim_dir")
#         fi
#     done < <(find "$parent_dir" -type f -name "input.in" -exec dirname {} \;)
# }

update_job_dir_array() {
    local parent_dir="$1"
    local except_string="${2:-Except:COMPLETED RUNNING PENDING}"

    JOB_DIR_ARRAY=()  # Clear previous entries

    echo ""
    echo "========================================================================"
    echo ">>> Starting update_job_dir_array"
    echo ">>> Parent directory: $parent_dir"
    echo ">>> Except statuses: ${except_string#Except:}"

    while IFS= read -r sim_dir; do
        local status
        status=$(check_simulation_status "$sim_dir")
        echo "--- Checking simulation dir: $sim_dir"
        echo "    Status: $status"

        if ! grep -qw "$status" <<< "${except_string#Except:}"; then
            echo "    -> Status not in except list. Adding to JOB_DIR_ARRAY."
            JOB_DIR_ARRAY+=("$sim_dir")
        else
            echo "    -> Status in except list. Skipping."
        fi
    done < <(find "$parent_dir" -type f -name "input.in" -exec dirname {} \;)

    
    echo ">>> Total directories added: ${#JOB_DIR_ARRAY[@]}"
    echo ">>> Finished update_job_dir_array"
    echo "========================================================================"
    echo ""
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
    local completion
    local lammps_output_file="$job_dir/output.out"
    if [ -f "$lammps_output_file" ]; then
        if grep -q "Total wall time" "$lammps_output_file"; then
            completion="complete"
        else
            completion="incomplete"
        fi
    else
        completion="missing"
    fi

    # Check SLURM job status
    local slurm_status
    if [[ "$completion" != "missing" ]]; then
        if ls "$job_dir"/slurm-*.out &> /dev/null; then
            local slurm_file
            slurm_file=$(ls -t "$job_dir"/slurm-*.out | head -n 1)
            local job_id
            job_id=$(basename "$slurm_file" | sed -E 's/slurm-([0-9]+)\.out/\1/')

            slurm_status=$(sacct -j "$job_id" --format=State --noheader | awk 'NR==1{gsub(/\+/, "", $1); print $1}') # try pinnacle
            if [[ -z "$slurm_status" ]]; then
                slurm_status=$(sacct -M merced -j "$job_id" --format=State --noheader | awk 'NR==1{gsub(/\+/, "", $1); print $1}') # then try merced
            fi
        else
            slurm_status="NOT_SUBMITTED"
        fi
    else
        slurm_status="NOT_SUBMITTED"
    fi

    # Decide final simulation status
    if [[ "$slurm_status" == "FAILED" ]]; then
        echo "HARD_FAIL"
    elif [[ "$slurm_status" == "COMPLETED" && "$completion" == "complete" ]]; then
        echo "COMPLETED"
    elif [[ "$slurm_status" == "COMPLETED" && "$completion" == "incomplete" ]]; then
        echo "SOFT_FAIL"
    else
        echo "$slurm_status"
    fi
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
    local sim_dir status job_id
    local sacct_output

    # Write header line to temporary file
    echo "JobID,JobName,Status,Partition,Submit,Start,End,Elapsed,AllocCPUs,NodeList,MaxRSS,ReqMem,SimDir,Error" > "$temp_csv_file"

    # Loop over all simulation directories
    while IFS= read -r sim_dir; do
        status=$(check_simulation_status "$sim_dir")

        # Default values
        local jobid="N/A" jobname="N/A" partition="N/A" submit="N/A"
        local start="N/A" end="N/A" elapsed="N/A" alloccpus="N/A"
        local nodelist="N/A" maxrss="N/A" reqmem="N/A" error="None"

        # Getting lammps error (if any)
        if [[ -f "${sim_dir}/output.out" ]]; then
            error=$(grep -i "error" "${sim_dir}/output.out" || echo "None")
        fi

        if [[ "$status" != "NOT_SUBMITTED" ]]; then
            if ls "$sim_dir"/slurm-*.out &> /dev/null; then
                local slurm_file
                slurm_file=$(ls -t "$sim_dir"/slurm-*.out | head -n 1)
                job_id=$(basename "$slurm_file" | sed -E 's/slurm-([0-9]+)\.out/\1/')

                # Try Pinnacle cluster first
                sacct_output=$(sacct -j "$job_id" --format=JobID,JobName,Start,End,Elapsed,Submit,AllocCPUs,NodeList,MaxRSS,ReqMem,Partition --noheader --parsable2 | grep -E "^${job_id}\|")

                # If not found, try Merced cluster
                if [[ -z "$sacct_output" ]]; then
                    sacct_output=$(sacct -M merced -j "$job_id" --format=JobID,JobName,Start,End,Elapsed,Submit,AllocCPUs,NodeList,MaxRSS,ReqMem,Partition --noheader --parsable2 | grep -E "^${job_id}\|")
                fi

                # Parse sacct output if available
                if [[ -n "$sacct_output" ]]; then
                    IFS="|" read -r jobid jobname start end elapsed submit alloccpus nodelist maxrss reqmem partition <<< "$sacct_output"
                fi
            fi
        fi

        # Write the final line to temporary CSV
        echo "\"$jobid\",\"$jobname\",\"$status\",\"$partition\",\"$submit\",\"$start\",\"$end\",\"$elapsed\",\"$alloccpus\",\"$nodelist\",\"$maxrss\",\"$reqmem\",\"$sim_dir\",\"$error\"" >> "$temp_csv_file"

    done < <(find "$parent_dir" -type f -name "input.in" -exec dirname {} \;)

    # Once everything is done, replace the old CSV file
    mv "$temp_csv_file" "$csv_file"
}

