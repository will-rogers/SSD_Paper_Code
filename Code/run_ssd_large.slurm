#!/bin/bash
#SBATCH --job-name=run_ssd_larger
#SBATCH --array=1-500          # Array job for 500 rasters
#SBATCH --output=logs/ssd_l_%A_%a.out
#SBATCH --error=logs/ssd_l_%A_%a.err
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
STEP_FILE="output_steps/used_${SLURM_ARRAY_TASK_ID}.rds"
OUTDIR="output_ssds_l"
mkdir -p $OUTDIR

# Create the output filename based on the SLURM array task ID
FILENAME="ssds_l_${SLURM_ARRAY_TASK_ID}.rds"

# Run the R script
Rscript -e "source('$SCRIPT'); \
            set.seed(100); \
            sim <- select_right_surface(readRDS('$SIMS_FILE'), index = ${SLURM_ARRAY_TASK_ID}); \
            issf.fit <- readRDS('$SSF_FILE')[[1]]; \
            mock.surface <- create_mock_surface(sim, F, list(x = 1, y = 1)); \
            pred.data <- get_cells(issf.fit, mock.surface, sim); \
            cell.data <- get_cell_data(issf.fit, pred.data); \
            neighbors.found <- neighbor_lookup(mock.surface, cell.data); \
            limit = 25; \
            sparse.neighbors <- neighbor_finder(issf.fit, cell.data, neighbors.found, distance.override = limit); \
            print('a'); \
            ssf.comparisons <- compile_ssf_comparisons(sparse.neighbors, cell.data); \
            print('b'); \
            ssf.comparisons <- lapply(ssf.comparisons, function(x) {list(.for = x$.for, .given = x$.given) }); \
            print('c'); \
            surface <- predict_ssf_comparisons(issf.fit, ssf.comparisons); \
            print('d'); \
            ssd.prob.raster <- wrap(eigen_decompose_to_raster(surface, mock.surface)); \
            saveRDS(ssd.prob.raster, file=file.path('$OUTDIR', '$FILENAME'))"