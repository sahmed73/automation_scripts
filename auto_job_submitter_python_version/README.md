# SLURM Auto Job Submitter

A Python script for automatically monitoring and submitting LAMMPS simulation jobs to SLURM clusters.

## Features

- Automatically finds simulation directories with available dependencies
- Monitors job status and completion 
- Submits jobs to appropriate SLURM partitions when resources are available
- Handles job dependencies
- Multi-threaded scanning for improved performance
- Detailed logging with timestamps

## Usage

1. Configure the script by editing the parameters in the `main()` function
2. Run the script:
   ```bash
   python auto_job_submitter_py.py
   ```

## Configuration

The main configuration parameters include:

- `root_dir`: Root directory to search for job directories
- `clusters`: List of SLURM clusters to use
- `exclude_status`: Job statuses to exclude from submission attempts
- `job_dir_identifying_filename`: File that identifies a directory as a job directory
- `job_completion_keyword`: String in output file indicating job completion
- `cluster_partition_config`: List of cluster/partition configurations with resource limits

## Requirements

- Python 3.8+
- SLURM environment with `squeue`, `sacct`, and `sbatch` commands available
- LAMMPS input files with proper structure
