#!/bin/bash
#SBATCH --job-name=run_rsf_maps
#SBATCH --array=1-500          # Array job for 500 rasters
#SBATCH --output=logs/rsfm_%A_%a.out
#SBATCH --error=logs/rsfm_%A_%a.err
#SBATCH --time=00:15:00         # Adjust time as needed
#SBATCH --cpus-per-task=4       # Number of CPUs per task
#SBATCH --mem=16G                # Memory per node
#SBATCH --partition=scavenge    # Partition name
#SBATCH --requeue               # rerun if stopped

# Load R module (if on a cluster that requires it)
module load R/4.2

# Define the input script, sims file, and output directory
SCRIPT="base_fxns.R"
SIMS_FILE="sims.RDS"
SSF_FILE="output_ssfs/ssfs_${SLURM_ARRAY_TASK_ID}.rds"
RSF_FILE="output_rsfs/rsfs_${SLURM_ARRAY_TASK_ID}.rds"
OUTDIR="output_rsf_maps"
mkdir -p $OUTDIR

# Create the output filename based on the SLURM array task ID
FILENAME="rsfm_${SLURM_ARRAY_TASK_ID}.rds"

# Run the R script
Rscript -e "source('$SCRIPT'); \
            set.seed(100); \
            sim <- select_right_surface(readRDS('$SIMS_FILE'), index = ${SLURM_ARRAY_TASK_ID}); \
            issf.fit <- readRDS('$SSF_FILE')[[1]]; \
            rsf.fit <- readRDS('$RSF_FILE'); \
            mock.surface <- create_mock_surface(sim, F, list(x = 1, y = 1)); \
            pred.data <- get_cells(issf.fit, mock.surface, sim); \
            rsf.prob.raster <- generate_rsf_map(pred.data = pred.data, model = rsf.fit, mock.surface); \
            saveRDS(rsf.prob.raster, file=file.path('$OUTDIR', '$FILENAME'))"