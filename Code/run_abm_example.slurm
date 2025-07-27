#!/bin/bash
#SBATCH --job-name=run_abm
#SBATCH --array=0-999          # 10,000 Array Jobs
#SBATCH --time=01:30:00
#SBATCH --output=logs/abm_%A.out
#SBATCH --error=logs/abm_%A.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=scavenge
#SBATCH --requeue

module load R/4.2

# === Configuration ===
SCRIPT="base_fxns.R"
OUTDIR="output_abm_example"
mkdir -p $OUTDIR

# === R Loop ===
Rscript - <<EOF
library(terra)

source('$SCRIPT')
total <- readRDS("abm.example.RDS")

spatial.data <- rast(total[[3]]) # rasters
model <- total[[2]] # movement model
ssf.dat <- total[[1]] # mock ssf data to init simulation

set.seed( $SLURM_ARRAY_TASK_ID)

abm.steps <- generate.abm(model = model, spatial.data = spatial.data, ssf.dat = ssf.dat, steps = 1E5, n.control = 100)

saveRDS(abm.steps, paste0('output_abm_example/abm_steps_', $SLURM_ARRAY_TASK_ID, '.rds'))

EOF