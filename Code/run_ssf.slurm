#!/bin/bash
#SBATCH --job-name=fit_steps
#SBATCH --array=1-500          # Array job for 500 rasters
#SBATCH --output=logs/ssf_%A_%a.out
#SBATCH --error=logs/ssf_%A_%a.err
#SBATCH --time=00:10:00         # Adjust time as needed
#SBATCH --cpus-per-task=4       # Number of CPUs per task
#SBATCH --mem=4G                # Memory per node
#SBATCH --partition=scavenge    # Partition name
#SBATCH --requeue               # rerun if stopped

# Load R module (if on a cluster that requires it)
module load R/4.2

# Define the input script, sims file, and output directory
SCRIPT="base_fxns.R"
SIMS_FILE="sims.RDS"
STEP_FILE="output_steps/used_${SLURM_ARRAY_TASK_ID}.rds"
OUTDIR="output_ssfs"
mkdir -p $OUTDIR

# Create the output filename based on the SLURM array task ID
FILENAME="ssfs_${SLURM_ARRAY_TASK_ID}.rds"

# Run the R script
Rscript -e "source('$SCRIPT'); \
            set.seed(100); \
             result <- generate_ssf(surfaces = readRDS('$SIMS_FILE'), movement = readRDS('$STEP_FILE'), index = ${SLURM_ARRAY_TASK_ID}); \
             saveRDS(result, file=file.path('$OUTDIR', '$FILENAME'))"