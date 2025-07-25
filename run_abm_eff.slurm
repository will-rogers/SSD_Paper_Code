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
SIMS_FILE="sims.RDS"
OUTDIR="output_abm"
BATCH_SIZE=250                    # Number of tasks per array job
START_INDEX=$((SLURM_ARRAY_TASK_ID * BATCH_SIZE))
END_INDEX=$((START_INDEX + BATCH_SIZE - 1))

mkdir -p $OUTDIR

# === R Loop ===
Rscript - <<EOF
library(terra)

source('$SCRIPT')
sims <- readRDS('$SIMS_FILE')
index_sheet <- read.csv('index.sheet.csv')

# Create a temporary multi-layer stack
stack_list <- list()

for (INDEX in $START_INDEX:$END_INDEX) {
    cat("Processing index:", INDEX, "\n")
    
    index.scenario <- index_sheet[INDEX + 1, ]
    
    # ==== Validate Inputs ====
    if (is.na(index.scenario[['landscape']]) || is.na(index.scenario[['scenarios']])) {
        cat("ERROR: Missing scenario or landscape for index", INDEX, "\n")
        next
    }
    
    # ==== Load Rasters ====
    sim <- try(rast(sims[[index.scenario[['landscape']]]]), silent = TRUE)
    if (inherits(sim, "try-error")) {
        cat("ERROR: Failed to load raster for landscape", index.scenario[['landscape']], "\n")
        next
    }
    
    issf_path <- paste0('output_ssfs/ssfs_', index.scenario[['scenarios']], '.rds')
    if (!file.exists(issf_path)) {
        cat("ERROR: Missing ISSF file for", issf_path, "\n")
        next
    }
    issf.fit <- readRDS(issf_path)[[2]]

    locations_path <- paste0('output_steps/used_', index.scenario[['scenarios']], '.rds')
    if (!file.exists(locations_path)) {
        cat("ERROR: Missing locations file for", locations_path, "\n")
        next
    }
    locations.df <- as.data.frame(readRDS(locations_path))
    
    # ==== Run Simulation ====
    abm.prob.raster <- try(rast(run_abm_2_convergence(issf.fit, sim, locations.df, iterator = INDEX)), silent = TRUE)
    
    if (inherits(abm.prob.raster, "try-error")) {
        cat("ERROR: ABM failed for index", INDEX, "\n")
        next
    }
    
    # ==== Validate Output ====
    if (!inherits(abm.prob.raster, "SpatRaster")) {
        cat("ERROR: Output is not a SpatRaster for index", INDEX, "\n")
        next
    }
    
    # ==== Add to Stack ====
    stack_list[[length(stack_list) + 1]] <- (abm.prob.raster)
}

# === Final Raster Stack ===
if (length(stack_list) > 0) {
    output_stack <- rast(stack_list)

    # Prepare output directory and save
    scenario_dir <- file.path('$OUTDIR', index.scenario[['scenarios']])
    if (!dir.exists(scenario_dir)) dir.create(scenario_dir, recursive = TRUE)

    # Write to disk with LZW compression
    output_file <- file.path(scenario_dir, paste0('abm_stack_', $SLURM_ARRAY_TASK_ID, '.tif'))
    writeRaster(output_stack, output_file, overwrite = TRUE)
    cat("Saved:", output_file, "\n")
} else {
    cat("No valid rasters to write for this batch.\n")
}

EOF
