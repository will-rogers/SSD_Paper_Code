#!/bin/bash
#SBATCH --job-name=run_abm
#SBATCH --array=0-9999          # 10,000 Array Jobs
#SBATCH --time=12:00:00
#SBATCH --output=logs/abm_%A.out
#SBATCH --error=logs/abm_%A.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=scavenge
#SBATCH --requeue

module load R/4.2

# === Configuration ===
SCRIPT="base_fxns.R"
OUTDIR="output_abm_size"
BATCH_SIZE=10                    # Number of tasks per array job
START_INDEX=$((SLURM_ARRAY_TASK_ID * BATCH_SIZE))
END_INDEX=$((START_INDEX + BATCH_SIZE - 1))

mkdir -p $OUTDIR

# === R Loop ===
Rscript - <<EOF
library(terra)

source('$SCRIPT')
setwd("~/project/movement")
total <- readRDS("abm.example.RDS")

print("setting conditions")

spatial.data <- rast(total[[3]]) # rasters
model <- total[[2]] # movement model (second list item is the model number)
ssf.dat <- total[[1]] %>% mutate(V1 = x1_, V2 = x2_) %>% dplyr::select('V1','V2') # mock ssf data to init simulation


sizes <- round(sqrt(seq(10*10,100*100, length.out = 10)))
size.choice <- ($SLURM_ARRAY_TASK_ID %/% 1000) + 1

spatial.data <- crop(spatial.data, c(0,sizes[size.choice],0,sizes[size.choice]))

# Create a temporary multi-layer stack
stack_list <- list()

for (INDEX in $START_INDEX:$END_INDEX) {
    cat("Processing index:", INDEX, "\n")
    
    set.seed(INDEX)
    
    abm.prob.raster <- rast(run_abm_2_convergence(model, spatial.data, ssf.dat, iterator = INDEX))
    
    stack_list[[length(stack_list) + 1]] <- (abm.prob.raster)
    
}

# === Final Raster Stack ===
if (length(stack_list) > 0) {
    output_stack <- rast(stack_list)


    # Write to disk with LZW compression
    output_file <- paste0('output_abm_size/abm_size_',sizes[size.choice],'_', $SLURM_ARRAY_TASK_ID, '.tif')
    writeRaster(output_stack, output_file, overwrite = TRUE)
    cat("Saved:", output_file, "\n")
    
} else {

    cat("No valid rasters to write for this batch.\n")
    
}

EOF
