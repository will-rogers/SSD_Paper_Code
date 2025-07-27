setwd("~/project/movement")
source("base_fxns.R")
library(fs)
library(terra)
library(foreach)
library(doParallel)

# Directory where rasters are stored
raster_dir <- "/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/fisher_nonfocal/2"
raster_files <- list.files(raster_dir, pattern = "\\.tif$", full.names = TRUE)
length(raster_files)

setwd("/gpfs/gibbs/project/ezenwa/wr254/movement/abm_compiled")
# Directory with 10,000 rasters (each 10 layers)


store.convergence <- matrix(NA, nrow = 10000, ncol = c(1+5+7+7+6))

colnames(store.convergence) <- c("red1", 
                                 "roe1", "roe2", "roe3", "roe4", "roe5",
                                 "fisher.focal1","fisher.focal2","fisher.focal3","fisher.focal4","fisher.focal5","fisher.focal6","fisher.focal7",
                                 "fisher.nonfocal1","fisher.nonfocal2","fisher.nonfocal3","fisher.nonfocal4","fisher.nonfocal5","fisher.nonfocal6","fisher.nonfocal7",
                                 "buffalo1", "buffalo2", "buffalo3", "buffalo4", "buffalo5", "buffalo6")


setwd("/gpfs/gibbs/project/ezenwa/wr254/movement")
for(i in 1:26){
  if(i == 1) {
    raster_dir <- paste0("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/red/",i,"/")
    rast100k <- rast(paste0("abm_compiled/final_red_",i,".tif"))/1E5
  }
  
  if(i %in% 2:6) {
    raster_dir <- paste0("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/roe/",i-1,"/")
    rast100k <- rast(paste0("abm_compiled/final_roe_",i-1,".tif"))/1E5
  }
  
  if(i %in% 7:13) {
    raster_dir <- paste0("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/fisher_focal/",i-6,"/")
    rast100k <- rast(paste0("abm_compiled/final_ff_",i-6,".tif"))/1E5
  }
  
  if(i %in% 14:20) {
    raster_dir <- paste0("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/fisher_nonfocal/",i-13,"/")
    rast100k <- rast(paste0("abm_compiled/final_fnf_",i-13,".tif"))/1E5
  }
  
  if(i %in% 21:26) {
    raster_dir <- paste0("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/buffalo/",i-20,"/")
    rast100k <- rast(paste0("abm_compiled/final_sum_b",i-20,".tif"))/1E5
  }
  
  raster_files <- dir_ls(raster_dir, glob = "*.tif")
  
  # Initialize sum raster with first file
  r0 <- rast(raster_files[1])
  r_sum <- sum(r0)  # sum over layers in first file
  
  store.convergence[1, i] <- sum(sqrt(values((r_sum/sum(values(r_sum))))*values(rast100k)))
  
  length.through <- length(raster_files)
  # Loop through remaining files and incrementally add to output
  pb = txtProgressBar(min = 0, max = length.through, initial = 0)
  for (f in 2:length.through) {
    r <- rast(raster_files[f])
    r <- sum(r)  # sum layers within the file
    r_sum <- r_sum + r
    
    store.convergence[f, i] <- sum(sqrt(values((r_sum/sum(values(r_sum))))*values(rast100k)))
    setTxtProgressBar(pb,f)
  }
  close(pb)
  print(i)
}

saveRDS(store.convergence, "store.convergence.case.RDS")



store.convergence.sim <- matrix(NA, nrow = 20, ncol = c(500))

colnames(store.convergence.sim) <- c(1:500)


setwd("/gpfs/gibbs/project/ezenwa/wr254/movement")
for(i in 1:500){
  raster_dir <- paste0("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm/",i,"/")
  rast.used <- read.RDS(paste0("abm_rasters/used_",i,".rds"))
  
  rast.used <- rast.used/sum(values(rast.used))
  
  raster_files <- dir_ls(raster_dir, glob = "*.tif")
  
  # Initialize sum raster with first file
  r0 <- rast(raster_files[1])
  r_sum <- sum(r0)  # sum over layers in first file
  
  store.convergence.sim[1, i] <- sum(sqrt(values((r_sum/sum(values(r_sum))))*values(rast.used)))
  
  length.through <- length(raster_files)
  # Loop through remaining files and incrementally add to output
  pb = txtProgressBar(min = 0, max = length.through, initial = 0)
  for (f in 2:length.through) {
    r <- rast(raster_files[f])
    r <- sum(r)  # sum layers within the file
    r_sum <- r_sum + r
    
    store.convergence.sim[f, i] <- sum(sqrt(values((r_sum/sum(values(r_sum))))*values(rast.used)))
    setTxtProgressBar(pb,f)
  }
  close(pb)
  print(i)
}

saveRDS(store.convergence.sim, "store.convergence.simulations.RDS")


store.convergence.sim <- readRDS("store.convergence.simulations.RDS")

as.data.frame(store.convergence.sim) %>% 
  mutate(iter = (1:20)*250) %>% 
  pivot_longer(colnames(store.convergence.sim)) %>% 
  ggplot(aes(x= iter, y = value, color = name)) +
  geom_path() +
  scale_color_discrete(guide = "none") +
  scale_x_log10()





###########

library(terra)

source('base_fxns.R')
setwd("~/project/movement")

raster_dir <- "/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_example/"
raster_files <- list.files(raster_dir, full.names = TRUE)
length(raster_files)

example.steps <- lapply(raster_files, readRDS)


example.steps <- lapply(1:101, function(x) {
  example.steps[[x]]$iteration <- x
  example.steps[[x]]$step <- 1:nrow(example.steps[[1]])
  example.steps[[x]]
})

bound <- rbindlist(example.steps)

bound


readRDS(raster_files[[2]]) %>% data.frame() %>% 
  mutate(t = 1:n()) %>% 
  ggplot(aes(x = x_, y = y_, color = t)) +
  geom_path() +
  scale_color_viridis_c()






library(terra)

source('$SCRIPT')
setwd("~/project/movement")
total <- readRDS("abm.example.RDS")

print("setting conditions")

spatial.data <- rast(total[[3]]) # rasters
model <- total[[2]] # movement model (second list item is the model number)
ssf.dat <- total[[1]] %>% mutate(V1 = x1_, V2 = x2_) %>% dplyr::select('V1','V2') # mock ssf data to init simulation

sizes <- round(sqrt(seq(10*10,100*100, length.out = 10)))
size.choice <- (1 %/% 1000) + 1

spatial.data <- crop(spatial.data, c(0,sizes[size.choice],0,sizes[size.choice]))

# Create a temporary multi-layer stack
stack_list <- list()

for (INDEX in 11:20) {
  cat("Processing index:", INDEX, "\n")
  
  set.seed(INDEX)
  
  ssf.dat$V1 <- ssf.dat$x1_
  ssf.dat$V2 <- ssf.dat$y1_
  ssf.dat$t <- ssf.dat$t1_
  
  abm.prob.raster <- try(rast(run_abm_2_convergence(model, spatial.data, ssf.dat, iterator = INDEX)), silent = TRUE)
  
  stack_list[[length(stack_list) + 1]] <- (abm.prob.raster)
}

# === Final Raster Stack ===
if (length(stack_list) > 0) {
  output_stack <- rast(stack_list)
  
  # Prepare output directory and save
  scenario_dir <- file.path('$OUTDIR', index.scenario[['scenarios']])
  if (!dir.exists(scenario_dir)) dir.create(scenario_dir, recursive = TRUE)
  
  # Write to disk with LZW compression
  output_file <- file.path(scenario_dir, paste0('abm_size_',sizes[size.choice],'_', $SLURM_ARRAY_TASK_ID, '.tif'))
  writeRaster(output_stack, output_file, overwrite = TRUE)
  cat("Saved:", output_file, "\n")
} else {
  cat("No valid rasters to write for this batch.\n")
}




