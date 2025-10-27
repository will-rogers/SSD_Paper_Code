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

# Directory with 10,000 rasters (each 10 layers)

d <- read.in.stack("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/fisher_nonfocal/2/")


merge.abm.case <- function(x) {
  pb = txtProgressBar(min = 0, max = length(x)*10, initial = 0, style = 3) 
  
  r_sum <- NULL
  count <- 0
  for (f in x) {
    r <- rast(f)               # load the raster (10 layers)
    r_sum <- if (is.null(r_sum)) {
      app(r, sum)              # sum across layers for the first file
    } else {
      r_sum + app(r, sum)      # add layer-sum from each file
    }
    count <- count + nlyr(r)   # increment total layer count
    # print(count)
    # rm(r); gc()
    setTxtProgressBar(pb,count)
  }
  
  # Compute mean across all pixels and layers
  r_mean <- r_sum / count
  close(pb)
  return(r_mean)
}


setwd("/gpfs/gibbs/project/ezenwa/wr254/movement/abm_compiled")
files <- list.files("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/buffalo/1/", pattern = "\\.tif$", full.names = TRUE)


final_red_1 <- merge.abm.case(list.files("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/red/1/", pattern = "\\.tif$", full.names = TRUE))
writeRaster(final_red_1, "final_red_1.tif", overwrite=T)



final_roe_1 <- merge.abm.case(list.files("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/roe/1/", pattern = "\\.tif$", full.names = TRUE))
writeRaster(final_roe_1, "final_roe_1.tif", overwrite=T)

final_roe_2 <- merge.abm.case(list.files("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/roe/2/", pattern = "\\.tif$", full.names = TRUE))
writeRaster(final_roe_2, "final_roe_2.tif", overwrite=T)

final_roe_3 <- merge.abm.case(list.files("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/roe/3/", pattern = "\\.tif$", full.names = TRUE))
writeRaster(final_roe_3, "final_roe_3.tif", overwrite=T)

final_roe_4 <- merge.abm.case(list.files("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/roe/4/", pattern = "\\.tif$", full.names = TRUE))
writeRaster(final_roe_4, "final_roe_4.tif", overwrite=T)

final_roe_5 <- merge.abm.case(list.files("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/roe/5/", pattern = "\\.tif$", full.names = TRUE))
writeRaster(final_roe_5, "final_roe_5.tif", overwrite=T)



final_ff_1 <- merge.abm.case(list.files("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/fisher_focal/1/", pattern = "\\.tif$", full.names = TRUE))
writeRaster(final_ff_1, "final_ff_1.tif", overwrite=T)

final_ff_2 <- merge.abm.case(list.files("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/fisher_focal/2/", pattern = "\\.tif$", full.names = TRUE))
writeRaster(final_ff_2, "final_ff_2.tif", overwrite=T)

final_ff_3 <- merge.abm.case(list.files("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/fisher_focal/3/", pattern = "\\.tif$", full.names = TRUE))
writeRaster(final_ff_3, "final_ff_3.tif", overwrite=T)

final_ff_4 <- merge.abm.case(list.files("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/fisher_focal/4/", pattern = "\\.tif$", full.names = TRUE))
writeRaster(final_ff_4, "final_ff_4.tif", overwrite=T)

final_ff_5 <- merge.abm.case(list.files("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/fisher_focal/5/", pattern = "\\.tif$", full.names = TRUE))
writeRaster(final_ff_5, "final_ff_5.tif", overwrite=T)

final_ff_6 <- merge.abm.case(list.files("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/fisher_focal/6/", pattern = "\\.tif$", full.names = TRUE))
writeRaster(final_ff_6, "final_ff_6.tif", overwrite=T)

final_ff_7 <- merge.abm.case(list.files("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/fisher_focal/7/", pattern = "\\.tif$", full.names = TRUE))
writeRaster(final_ff_7, "final_ff_7.tif", overwrite=T)



final_fnf_1 <- merge.abm.case(list.files("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/fisher_nonfocal/1/", pattern = "\\.tif$", full.names = TRUE))
writeRaster(final_fnf_1, "final_fnf_1.tif", overwrite=T)

final_fnf_2 <- merge.abm.case(list.files("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/fisher_nonfocal/2/", pattern = "\\.tif$", full.names = TRUE))
writeRaster(final_fnf_2, "final_fnf_2.tif", overwrite=T)

final_fnf_3 <- merge.abm.case(list.files("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/fisher_nonfocal/3/", pattern = "\\.tif$", full.names = TRUE))
writeRaster(final_fnf_3, "final_fnf_3.tif", overwrite=T)

final_fnf_4 <- merge.abm.case(list.files("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/fisher_nonfocal/4/", pattern = "\\.tif$", full.names = TRUE))
writeRaster(final_fnf_4, "final_fnf_4.tif", overwrite=T)

final_fnf_5 <- merge.abm.case(list.files("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/fisher_nonfocal/5/", pattern = "\\.tif$", full.names = TRUE))
writeRaster(final_fnf_5, "final_fnf_5.tif", overwrite=T)

final_fnf_6 <- merge.abm.case(list.files("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/fisher_nonfocal/6/", pattern = "\\.tif$", full.names = TRUE))
writeRaster(final_fnf_6, "final_fnf_6.tif", overwrite=T)

final_fnf_7 <- merge.abm.case(list.files("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/fisher_nonfocal/7/", pattern = "\\.tif$", full.names = TRUE))
writeRaster(final_fnf_7, "final_fnf_7.tif", overwrite=T)



final_sum_b1 <- merge.abm.case(list.files("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/buffalo/1/", pattern = "\\.tif$", full.names = TRUE))
writeRaster(final_sum_b1, "final_sum_b1.tif", overwrite=T)

final_sum_b2 <- merge.abm.case(list.files("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/buffalo/2/", pattern = "\\.tif$", full.names = TRUE))
writeRaster(final_sum_b2, "final_sum_b2.tif", overwrite=T)

final_sum_b3 <- merge.abm.case(list.files("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/buffalo/3/", pattern = "\\.tif$", full.names = TRUE))
writeRaster(final_sum_b3, "final_sum_b3.tif", overwrite=T)

final_sum_b4 <- merge.abm.case(list.files("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/buffalo/4/", pattern = "\\.tif$", full.names = TRUE))
writeRaster(final_sum_b4, "final_sum_b4.tif", overwrite=T)

final_sum_b5 <- merge.abm.case(list.files("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/buffalo/5/", pattern = "\\.tif$", full.names = TRUE))
writeRaster(final_sum_b5, "final_sum_b5.tif", overwrite=T)

final_sum_b6 <- merge.abm.case(list.files("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/buffalo/6/", pattern = "\\.tif$", full.names = TRUE))
writeRaster(final_sum_b6, "final_sum_b6.tif", overwrite=T)







read.in.stack <- function(dir) {
  take <- list.files(dir)
  take <- take[which(!str_detect(take, ".tif."))]
  stack <- rast(lapply(take, function(x){
    rast(paste0(dir,x))
  }))
  
  return(stack)
}

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
    rast100k <- rast(paste0("abm_compiled/final_red_",i,".tif"))
  }
  
  if(i %in% 2:6) {
    raster_dir <- paste0("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/roe/",i-1,"/")
    rast100k <- rast(paste0("abm_compiled/final_roe_",i-1,".tif"))
  }
  
  if(i %in% 7:13) {
    raster_dir <- paste0("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/fisher_focal/",i-6,"/")
    rast100k <- rast(paste0("abm_compiled/final_ff_",i-6,".tif"))
  }
  
  if(i %in% 14:20) {
    raster_dir <- paste0("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/fisher_nonfocal/",i-13,"/")
    rast100k <- rast(paste0("abm_compiled/final_fnf_",i-13,".tif"))
  }
  
  if(i %in% 21:26) {
    raster_dir <- paste0("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm_case/buffalo/",i-20,"/")
    rast100k <- rast(paste0("abm_compiled/final_sum_b",i-20,".tif"))
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
for(i in 379:500){
  raster_dir <- paste0("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm/",i,"/")
  rast.used <- readRDS(paste0("output_rasters/used_",i,".rds"))
  
  rast.used$n <- rast.used$n/sum(rast.used$n)
  rast.used <- rast(rast.used[which(!is.na(rast.used$y)),], type = "xyz")
  values(rast.used$n) <- ifelse(is.na(values(rast.used$n)), 0, values(rast.used$n))
  
  raster_files <- dir_ls(raster_dir, glob = "*.tif")
  
  # Initialize sum raster with first file
  r0 <- rast(raster_files[1])
  r_sum <- sum(r0)  # sum over layers in first file
  
  store.convergence.sim[1, i] <- sum(sqrt(values((r_sum/sum(values(r_sum))))*values(rast.used$n)))
  
  length.through <- length(raster_files)
  # Loop through remaining files and incrementally add to output
  pb = txtProgressBar(min = 0, max = length.through, initial = 0)
  for (f in 2:length.through) {
    r <- rast(raster_files[f])
    r <- sum(r)  # sum layers within the file
    r_sum <- r_sum + r
    
    store.convergence.sim[f, i] <- sum(sqrt(values((r_sum/sum(values(r_sum))))*values(rast.used$n)))
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

store.convergence.sim <- readRDS("old/store.convergence.simulations.RDS")

as.data.frame(store.convergence.sim) %>% 
  mutate(iter = (1:20)*250) %>% 
  pivot_longer(colnames(store.convergence.sim)) %>% 
  ggplot(aes(x= iter, y = value, color = name)) +
  geom_path() +
  scale_color_discrete(guide = "none") +
  scale_x_log10()




