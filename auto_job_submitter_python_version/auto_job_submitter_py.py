# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sun Jun  1 15:47:42 2025

Full python version of auto slurm job submitter
"""

import os
import re
import subprocess
import time
from datetime import datetime
import concurrent.futures

def run_slurm_command(command, config):
    """
    Execute a SLURM command and handle errors gracefully.
    
    Args:
        command (list): The SLURM command to run as a list of strings.
        config (dict): Configuration dictionary containing logger.
        
    Returns:
        subprocess.CompletedProcess or None: The result of the command or None if it failed.
    """
    log = config.get("logger", print)
    label = command[0]
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        return result
    except subprocess.CalledProcessError as e:
        prefix = f"\t\t-- [SLURM Error] '{label}' access error: "
        indent = " " * len(prefix)
        formatted_error = prefix + e.stderr.strip().replace("\n", f"\n{indent}")
        log(formatted_error, timestamp=True)
        return None



# 1
def extract_lammps_variables(input_file_path):   
    """
    Extract variable definitions from a LAMMPS input file.
    
    Args:
        input_file_path (str): Path to the LAMMPS input file.
        
    Returns:
        dict: Dictionary of variable names and their values.
    """
    variables = {}
    # only can handle literal-string
    styles = "(string|equal)"  # Only handling a subset of variable styles
    variable_pattern = re.compile(
        rf"^\s*variable\s+(\w+)\s+{styles}(?:\s+(.*))?", re.IGNORECASE
    )

    with open(input_file_path, 'r') as file:
        for line in file:
            line = line.split('#', 1)[0].strip()  # Remove comments
            if not line:
                continue
            match = variable_pattern.match(line)
            if match:
                var_name = match.group(1)
                var_value = match.group(3).strip() if match.group(3) else ""
                variables[var_name] = var_value
    
    return variables

# 2
def get_dependency_file_path(dirr, config):
    """
    Extracts dependency file paths from LAMMPS input files.
    Looks for read_restart or read_data commands and resolves variable references.
    
    Args:
        dirr (str): Directory containing the LAMMPS input file.
        config (dict): Configuration dictionary with job settings.
        
    Returns:
        str or None: Path to the dependency file if found, None otherwise.
    """
    job_dir_identifying_filename = config["job_dir_identifying_filename"]
    
    input_file_path = os.path.join(dirr, job_dir_identifying_filename)
    variables = extract_lammps_variables(input_file_path)
        
    max_try = 20
    with open(input_file_path, 'r') as file:
        for line in file: 
            if line.strip().startswith('read_restart') or line.strip().startswith('read_data'):
                parts = line.split()
                if len(parts) > 1:
                    file_location = parts[1].strip()
                    pattern = r"\$\{?(\w+)\}?"
                    try_count = 0
                    while try_count < max_try:
                        try_count+=1
                        match = re.search(pattern, file_location)
                        if match:
                            var_name = match.group(1)
                            var_value = variables.get(var_name, f"${{{var_name}}}")
                            file_location = re.sub(pattern, var_value, file_location, count=1)
                        else:
                            break

                    if not os.path.isabs(file_location) and file_location:
                        file_location = os.path.abspath(os.path.join(dirr, file_location))

                    return file_location
    return None

# 3
def read_job_id_from_text_file(job_dir):
    """
    Read the job ID from the job_id.txt file in the job directory.
    
    Args:
        job_dir (str): Path to the job directory.
        
    Returns:
        str or None: The job ID if found, None otherwise.
    """
    job_id_path = os.path.join(job_dir, "job_id.txt")

    if os.path.isfile(job_id_path):
        with open(job_id_path, 'r') as f:
            return f.readline().strip()
    else:
        return None
    
    
# 4
def check_simulation_output_file_status(job_dir, config):
    """
    Check the simulation output file to determine if the job completed successfully.
    
    Args:
        job_dir (str): Path to the job directory.
        config (dict): Configuration containing output filename and completion keyword.
        
    Returns:
        str: Status of the output file - "complete", "incomplete", or "missing".
    """
    job_output_filename = config["job_output_filename"]
    job_completion_keyword = config["job_completion_keyword"]

    output_filepath = os.path.join(job_dir, job_output_filename)

    if os.path.isfile(output_filepath):
        with open(output_filepath, 'r') as file:
            if any(job_completion_keyword in line for line in file):
                return "complete"
            
            # will delete later
            elif "eq.nvt.data" in os.listdir(job_dir):
                short_path = job_dir[job_dir.find("/", 65)+1:]
                log = config.get("logger", print)
                log(f"\t\t-- No wall-time but has eq.nvt.data: {short_path}")
                return "complete"
            
            else:
                return "incomplete"
    else:
        return "missing"


# 5
def check_slurm_job_status(job_dir, config):
    """
    Check the status of a SLURM job using its job ID.
    
    Args:
        job_dir (str): Path to the job directory.
        config (dict): Configuration containing clusters information.
        
    Returns:
        str or None: The SLURM job status (e.g., "COMPLETED", "RUNNING") if found, None otherwise.
    """
    clusters = config["clusters"]
    
    job_id = read_job_id_from_text_file(job_dir)
    slurm_status = None

    if job_id:
        for cluster in clusters:
            command = ["sacct", "-j", job_id, "-M", cluster, "--format=State", "--noheader", "--parsable2"]
            result = run_slurm_command(command, config)
            if result:
                output = result.stdout.strip()
                if output:
                    first_line = output.splitlines()[0]
                    state = first_line.split('|')[0].split()[0].replace('+', '')
                    slurm_status = state
                    break

    return slurm_status


def check_simulation_status(job_dir, config):
    """
    Determine overall simulation status by combining output file status and SLURM job status.
    
    Args:
        job_dir (str): Path to the job directory.
        config (dict): Configuration dictionary.
        
    Returns:
        str: Combined status of the simulation ("COMPLETED", "SOFT_FAIL", "NOT_YET_SUBMITTED", or SLURM status).
    """
    
    sim_output_file_status = check_simulation_output_file_status(job_dir, config)
    slurm_status = check_slurm_job_status(job_dir, config)
    
    if slurm_status == "COMPLETED" and sim_output_file_status == "complete":
        return "COMPLETED"
    elif slurm_status == "COMPLETED" and sim_output_file_status == "incomplete":
        return "SOFT_FAIL"
    elif slurm_status == None and sim_output_file_status == "missing":
        return "NOT_YET_SUBMITTED"
    else:
        return slurm_status

def update_job_dir_array(root_dir, config):
    """
    Find all job directories that are not in excluded states.
    
    Args:
        root_dir (str): Root directory to search for job directories.
        config (dict): Configuration dictionary containing job directory identification criteria.
        
    Returns:
        list: List of job directories that are not in excluded states.
    """
    job_dir_identifying_filename = config["job_dir_identifying_filename"]
    exclude_status = config["exclude_status"]

    # Collect candidate directories first
    candidate_dirs = []
    for dirpath, _, files in os.walk(root_dir):
        if job_dir_identifying_filename in files:
            candidate_dirs.append(dirpath)
    
    def check_and_filter(dirpath):
        sim_status = check_simulation_status(dirpath, config)
        if sim_status not in exclude_status:
            return dirpath
        return None

    job_dir_array = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=16) as executor:
        for result in executor.map(check_and_filter, candidate_dirs):
            if result:
                job_dir_array.append(result)

    return job_dir_array


def get_jobs_ready_for_submission(job_dir_array, config):    
    """
    Filters job directories that are ready for submission based on dependency availability.
    
    Args:
        job_dir_array (list): List of job directory paths to check.
        config (dict): Configuration dictionary with job settings.
        
    Returns:
        list: List of job directories that are ready for submission.
    """
    jobs_ready_for_submission_array = []
    for job_dir in job_dir_array:
        dependency_file_path = get_dependency_file_path(job_dir, config)
        if dependency_file_path:
            # Add only if the dependency file exists
            if os.path.isfile(dependency_file_path):
                jobs_ready_for_submission_array.append(job_dir)
        else:
            # No dependency specified; consider ready for submission
            jobs_ready_for_submission_array.append(job_dir)
            
    return jobs_ready_for_submission_array

def update_slurm_resources_in_submit_file(job_dir, partition, ntasks_per_node, time_limit, config):
    """
    Update SLURM submit script in job_dir with the desired partition, cores, and time,
    """
    submit_filepath = os.path.join(job_dir, config["job_submit_filename"])

    # Read existing submit file
    with open(submit_filepath, 'r') as f:
        lines = f.readlines()

    # Replace or insert SLURM directives
    def replace_or_insert(lines, key, value):
        pattern = re.compile(rf'^#SBATCH\s+{re.escape(key)}=.*')
        replaced = False
        for i, line in enumerate(lines):
            if pattern.match(line):
                lines[i] = f"#SBATCH {key}={value}\n"
                replaced = True
                break
        if not replaced:
            lines.insert(2, f"#SBATCH {key}={value}\n")
        return lines

    lines = replace_or_insert(lines, "--partition", partition)
    lines = replace_or_insert(lines, "--ntasks-per-node", str(ntasks_per_node))
    lines = replace_or_insert(lines, "--time", time_limit)

    # Write back to file
    with open(submit_filepath, 'w') as f:
        f.writelines(lines)

    
def submit_job(job_dir, config):   
    '''
        Max resources:
            
        "pinnacles|short            |12|6:00:00   |56",
        "pinnacles|long             |3 |3-00:00:00|56",
        "pinnacles|medium           |6 |1-00:00:00|56",
        "pinnacles|bigmem           |2 |3-00:00:00|56",
        "pinnacles|pi.amartini      |2 |3-00:00:00|56",
        "merced   |compute          |3 |5-00:00:00|40",
        "merced   |bigmem           |3 |5-00:00:00|24",
        "pinnacles|cenvalarc.compute|3 |3-00:00:00|56", 
        "pinnacles|cenvalarc.bigmem |2 |3-00:00:00|56"
    '''
    username = config["username"]
    job_submit_filename = config["job_submit_filename"]
    cluster_partition_config = config["cluster_partition_config"]
    logger = config.get("logger", print)
    
    short_path = job_dir[job_dir.find("/", 65)+1:]  # Hard Code

    # logger(f"[DEBUG] Starting job submission attempt in directory: {job_dir}")

    for cf in cluster_partition_config:
        cluster, partition, max_jobs, time_limit, num_cores = [s.strip() for s in cf.split('|')]
        max_jobs = int(max_jobs)
        num_cores = int(num_cores)

        # logger(f"[DEBUG] Checking queue on cluster='{cluster}', partition='{partition}' (max_jobs={max_jobs})")

        # Check current job count in queue
        command = ["squeue", "-M", cluster, "--user", username, "--partition", partition, "--noheader"]
        result = run_slurm_command(command, config)
        njob_queue = len(result.stdout.strip().splitlines()) if result else max_jobs  # prevent submission to a failed partition

        # logger(f"[DEBUG] Current job count in queue: {njob_queue}/{max_jobs}")

        if njob_queue < max_jobs:
            # logger(f"[DEBUG] Preparing to submit job to partition '{partition}' on cluster '{cluster}'")

            # Prepare and submit job
            update_slurm_resources_in_submit_file(job_dir, partition, num_cores, time_limit, config)
            submit_file = os.path.join(job_dir, job_submit_filename)

            current_dir = os.getcwd()
            os.chdir(job_dir)

            try:
                command = ["sbatch", "-M", cluster, submit_file]
                # logger(f"[DEBUG] Submitting job with command: {' '.join(command)}")
                result = run_slurm_command(command, config)
                if result and result.returncode == 0 and "Submitted batch job" in result.stdout:
                    job_id = result.stdout.strip().split()[3]
                    job_id_file = os.path.join(job_dir, "job_id.txt")

                    try: old = open(job_id_file).readlines()
                    except: old = []

                    with open(job_id_file, 'w') as f:
                        f.write(f"{job_id}\n")
                        f.writelines(old)

                    logger(f"\t\t-- JOB SUBMITTED: {short_path} on {cluster}.{partition}")
                    return True
            finally:
                os.chdir(current_dir)
    
    logger(f"\t\t-- JOB NOT SUBMITTED: {short_path} (Queue full or error)")
    return False


def setup_logger(log_file_base):
    timestamp = time.strftime('%Y%m%d-%H%M%S')
    base, ext = os.path.splitext(log_file_base)
    log_file = f"{base}_{timestamp}{ext}"  # Makes: auto_submit_20250604_142100.log

    def log(message, timestamp=False, print_to_console=True):
        parts = []
        parts.append(str(message))
        if timestamp:
            parts.append(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}]")

        log_line = " ".join(parts) + "\n"
        with open(log_file, 'a') as f:  # Append mode
            f.write(log_line)
        if print_to_console:
            print(log_line.strip())

    return log     


def initialize_jobid_files_for_legacy_runs(root_dir, config):
    """
    Create job_id.txt files for simulations from legacy runs (those without job_id.txt).
    It looks for slurm-<jobid>.out files, checks submit times via sacct,
    and writes the most recent job ID into job_id.txt.
    """
    log = config.get("logger", print)
    clusters = config["clusters"]
    job_id_pattern = re.compile(r"slurm-(\d+)\.out")
    job_dir_identifying_filename = config["job_dir_identifying_filename"]
    
    start_time = time.time()
    
    log("initialization of job_id.txt files for legacy runs...", timestamp=True)


    total_dirs = 0
    created_count = 0
    skipped_count = 0

    for dirpath, _, files in os.walk(root_dir):
        if job_dir_identifying_filename not in files: continue
    
        total_dirs += 1
        if "job_id.txt" in files:
            skipped_count += 1
            continue

        slurm_job_ids = [match.group(1) for f in files if (match := job_id_pattern.match(f))]

        job_submit_times = []
        for job_id in slurm_job_ids:
            for cluster in clusters:
                command = ["sacct", "-j", job_id, "-M", cluster, "--format=Submit", "--noheader", "--parsable2"]
                result = run_slurm_command(command, config)
                if result and result.stdout.strip():
                    try:
                        submit_time = datetime.strptime(result.stdout.strip().splitlines()[0], "%Y-%m-%dT%H:%M:%S")
                        job_submit_times.append((job_id, submit_time))
                    except ValueError:
                        continue

        if job_submit_times:
            latest_job_id = sorted(job_submit_times, key=lambda x: x[1], reverse=True)[0][0]
            job_id_path = os.path.join(dirpath, "job_id.txt")
            with open(job_id_path, 'w') as f:
                f.write(f"{latest_job_id}\n")
                
            relative_path = os.path.relpath(dirpath, root_dir)
            log(f"\t\t- Created job_id.txt in '{relative_path}' with SLURM job ID: {latest_job_id}")
            created_count += 1
        else:
            skipped_count += 1

    log(f"\t- Directories scanned: {total_dirs} | job_id.txt created: {created_count} | Skipped: {skipped_count}", timestamp=True)
    
    end_time = time.time()
    extime = end_time - start_time
    log(f"\t- Runtime: {extime}")

       

def main():
    # inputs
    root_dir = "/mnt/borgstore/amartini/sahmed73/data/REACTER/TEST_REVERSE_REAC_V2"
    log_filepath = os.path.join(root_dir, "auto_submit.log")
    log = setup_logger(log_filepath)
    
    
    config = {
        "username": "sahmed73",
        "job_dir_identifying_filename": "input.in",
        "clusters": ["pinnacles"],
        "exclude_status": ["COMPLETED", "RUNNING", "PENDING"],
        "job_submit_filename": "submit.sh",
        "job_output_filename": "output.out",
        "job_completion_keyword": "Total wall time", # "Normal termination of Gaussian"
        "logger": log,
        "cluster_partition_config": [
            "pinnacles|test              |1 |1:00:00   |52",
            "pinnacles|short             |12|6:00:00   |52",
            "pinnacles|long              |3 |6:00:00   |52",
            "pinnacles|medium            |6 |6:00:00   |52",
            "pinnacles|bigmem            |2 |6:00:00   |52",
            "pinnacles|pi.amartini       |6 |6:00:00   |52",
            "pinnacles|cenvalarc.compute |3 |6:00:00   |52",
            "pinnacles|cenvalarc.bigmem  |2 |6:00:00   |52"
        ]

    }


    # initialize_jobid_files_for_legacy_runs(root_dir, config)
    
    log("Auto job submission script started.", timestamp=True)
    job_dir_array = update_job_dir_array(root_dir, config)
    log(f"Detected {len(job_dir_array)} job directories in the root directory.")
    log("="*70)
    log("")

    
    submission_loop_count = 0
    while job_dir_array:
        start_time = time.time()
        submission_loop_count += 1
        log(f"\nSubmission loop: {submission_loop_count}", timestamp=True)
        
        jobs_ready_for_submission_array = get_jobs_ready_for_submission(job_dir_array, config)
        log(f"\t- {len(jobs_ready_for_submission_array)}/{len(job_dir_array)} jobs ready for submission.")
        
        job_submitted_count = 0
    
        for job_dir in jobs_ready_for_submission_array:
            submitted = submit_job(job_dir, config)
            if submitted:
                job_submitted_count += 1
                if job_dir in job_dir_array:
                    job_dir_array.remove(job_dir)
                time.sleep(1)  # wait to avoid submission overlap
            else:
                break
    
        runtime = time.time() - start_time
        log(f"\t- {job_submitted_count}/{len(jobs_ready_for_submission_array)} jobs submitted successfully")
        log(f"\t- Runtime: {runtime:.2f} seconds")
        log("\t- Waiting 60 seconds before next loop ...")
        time.sleep(60)


    
if __name__ == "__main__":
    main()












