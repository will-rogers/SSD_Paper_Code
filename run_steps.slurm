#!/bin/bash
#SBATCH --job-name=simulate_steps
#SBATCH --array=1-500          # Array job for 500 rasters
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --error=logs/sim_%A_%a.err
#SBATCH --time=00:05:00         # Adjust time as needed
#SBATCH --cpus-per-task=4       # Number of CPUs per task
#SBATCH --mem=4G                # Memory per node
#SBATCH --partition=scavenge    # Partition name
#SBATCH --requeue               # rerun if stopped

# Load R module (if on a cluster that requires it)
module load R/4.2

# Define the input script, sims file, and output directory
SCRIPT="base_fxns.R"
SIMS_FILE="sur.val.packed.RDS"
OUTDIR="output_steps"
mkdir -p $OUTDIR

# Create the output filename based on the SLURM array task ID
FILENAME="used_${SLURM_ARRAY_TASK_ID}.rds"

# Run the R script
Rscript -e "source('$SCRIPT'); \
             raster <- rast(readRDS('$SIMS_FILE')[[${SLURM_ARRAY_TASK_ID}]]); \
             result <- simulate_movements(raster, step_no = 1E4, teleport = NULL); \
             saveRDS(result, file=file.path('$OUTDIR', '$FILENAME'))"