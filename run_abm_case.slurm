#!/bin/bash
#SBATCH --job-name=abm
#SBATCH --array=0-9999 #distributed over 10000 jobs
#SBATCH --output=logs/abm_%A.out
#SBATCH --error=logs/abm_%A.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=scavenge
#SBATCH --requeue

module load R/4.2

# === Configuration ===
SCRIPT="base_fxns.R"
OUTDIR="output_abm_case/fisher_focal/1" #directory
BATCH_SIZE=10                    # Number of tasks per array job
START_INDEX=$((SLURM_ARRAY_TASK_ID * BATCH_SIZE))
END_INDEX=$((START_INDEX + BATCH_SIZE - 1))

mkdir -p $OUTDIR

# === R Loop ===
Rscript - <<EOF
library(terra)

source('$SCRIPT')
setwd("~/project/movement")
total <- readRDS("fisher.nonfocal.abm.RDS")

print("setting conditions")
spatial.data <- rast(total[[3]]) # rasters
model <- total[[2]][[1]] # movement model (second list item is the model number)
ssf.dat <- total[[1]] # mock ssf data to init simulation

# Create a temporary multi-layer stack for each batch of 10
stack_list <- list()
print("starting simulations")
for (INDEX in $START_INDEX:$END_INDEX) {
   
    abm.prob.raster <- rast(run_abm_2_convergence2(model, spatial.data, ssf.dat, iterator = INDEX)) 
    
    stack_list[[length(stack_list) + 1]] <- (abm.prob.raster)
}
print("merging")
# === Final Raster Stack ===
if (length(stack_list) > 0) {
    output_stack <- rast(stack_list)

    # Write to disk with LZW compression
    output_file <- file.path('$OUTDIR', paste0('abm_stack_', $SLURM_ARRAY_TASK_ID, '.tif'))
    writeRaster(output_stack, output_file, overwrite = TRUE)
    cat("Saved:", output_file, "\n")
} else {
    cat("No valid rasters to write for this batch.\n")
}

EOF
