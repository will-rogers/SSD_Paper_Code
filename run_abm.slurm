#!/bin/bash
#SBATCH --job-name=run_abm
#SBATCH --array=1-10          # Array job for 500 rasters
#SBATCH --output=logs/ssd_%A_%a.out
#SBATCH --error=logs/ssd_%A_%a.err
#SBATCH --time=00:10:00         # Adjust time as needed
#SBATCH --cpus-per-task=1       # Number of CPUs per task
#SBATCH --mem=4G                # Memory per node
#SBATCH --partition=scavenge    # Partition name
#SBATCH --requeue               # rerun if stopped

# Load R module (if on a cluster that requires it)
module load R/4.2

# Define the input script, sims file, and output directory
SCRIPT="base_fxns.R"
SIMS_FILE="sims.RDS"
OUTDIR="output_abms"
mkdir -p $OUTDIR

# Create the output filename based on the SLURM array task ID
FILENAME="abms_1_${SLURM_ARRAY_TASK_ID}.rds"

# Run the R script
Rscript -e "source('$SCRIPT'); \
            index <- ${SLURM_ARRAY_TASK_ID}
            index.scenario <- read.csv('index.sheet.csv')[index,]; \
            sim <- rast(readRDS('$SIMS_FILE')[[index.scenario[,'landscape']]]); \
            issf.fit <- readRDS(paste0('output_ssfs/ssfs_',index.scenario[,'scenarios'],'.rds')); \
            locations.df <- as.data.frame(readRDS(paste0('output_steps/used_',index.scenario[,'scenarios'],'.rds'))); \
            abm.prob.raster <- run_abm_2_convergence(issf.fit, sim, locations.df, iterator = index);\
            saveRDS(abm.prob.raster, file=file.path('$OUTDIR', '$FILENAME'))"