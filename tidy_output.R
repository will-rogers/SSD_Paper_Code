setwd("~/project/movement")
source("base_fxns.R")


read.in.stack <- function(dir) {
  take <- list.files(dir)
  take <- take[which(!str_detect(take, ".tif."))]
  stack <- rast(lapply(take, function(x){
    rast(paste0(dir,x))
  }))

  return(stack)
}

for (i in 1:500) {
  combined <- mean(read.in.stack(paste0("/gpfs/gibbs/project/ezenwa/wr254/movement/output_abm/",i,"/")))
  writeRaster(combined, paste0("output_abm_sim/sim",i,".tif"), overwrite=T)
}

require(amt)
ssf.coef <- rep(list(NA), 500)
ssf.est <- for(i in 1:500){
  ssf.fit <- readRDS(paste0("/gpfs/gibbs/project/ezenwa/wr254/movement/output_ssfs/ssfs_",i,".rds"))

  a <- summary(ssf.fit$base.ssf$model)$coefficients

  ssf.coef[[i]] <- data.frame(term = c("z.1","z.2","z.3","z.4"),
                   coef = a[1:4,1],
                   se = a[1:4,3],
                   group = i)
  print(i)
}

saveRDS(ssf.coef, "ssf.coefs.RDS")

rsf.coef <- rep(list(NA), 500)
for(i in 1:500){
  rsf.fit <- readRDS(paste0("/gpfs/gibbs/project/ezenwa/wr254/movement/output_rsfs/rsfs_",i,".rds"))
  
  a <- summary(rsf.fit)$coefficients
  
  rsf.coef[[i]] <- data.frame(term = c("z.1","z.2","z.3","z.4"),
                   coef = a[1:4,1],
                   se = a[1:4,3],
                   group = i)
  print(i)
}

saveRDS(rsf.coef, "rsf.coefs.RDS")


