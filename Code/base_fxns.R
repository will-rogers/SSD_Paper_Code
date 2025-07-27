
#### Paper Functions

.libPaths(c("/home/wr254/R/4.2", .libPaths()))

# devtools::install_github("jmsigner/amt", force = T)
# install.packages("tidyverse", "/home/wr254/R/4.2")
# install.packages("amt", "/home/wr254/R/4.2")
# install.packages("parallel", "/home/wr254/R/4.2")
# install.packages("data.table", "/home/wr254/R/4.2")
# install.packages("Matrix", "/home/wr254/R/4.2")
# install.packages("progress", "/home/wr254/R/4.2")
# install.packages("pbapply", "/home/wr254/R/4.2")
# install.packages("pbmcapply", "/home/wr254/R/4.2")
# install.packages("igraph", "/home/wr254/R/4.2")
# install.packages("RSpectra", "/home/wr254/R/4.2")
# install.packages("terra", "/home/wr254/R/4.2")

library("readr",lib.loc = "/home/wr254/R/4.2")
library("tidyverse",lib.loc = "/home/wr254/R/4.2")
library("sf",lib.loc = "/home/wr254/R/4.2")
library("amt",lib.loc = "/home/wr254/R/4.2")
library("parallel", lib.loc ="/home/wr254/R/4.2")
library("data.table",lib.loc = "/home/wr254/R/4.2")
library("Matrix")
library("progress", lib.loc ="/home/wr254/R/4.2")
library("pbapply",lib.loc = "/home/wr254/R/4.2")
library("pbmcapply", lib.loc ="/home/wr254/R/4.2")
library("igraph",lib.loc = "/home/wr254/R/4.2")
library("RSpectra",lib.loc = "/home/wr254/R/4.2")
library("terra")

create_mock_surface <- function(raster.list, multiple.extents = F, resolution = list(x = 100, y = 100)){
  
  # If rasters are all of the same extent, take the extent
  if(multiple.extents == F){
    xmin <- ext(raster.list)[1]
    xmax <- ext(raster.list)[2]
    ymin <- ext(raster.list)[3]
    ymax <- ext(raster.list)[4]
  }
  
  # If rasters are of different extent, take the overlap extent
  if(multiple.extents == T){
    xmin <- max(unlist(lapply(raster.list, function(x) {ext(x)[1]})))
    xmax <- max(unlist(lapply(raster.list, function(x) {ext(x)[2]})))
    ymin <- max(unlist(lapply(raster.list, function(x) {ext(x)[3]})))
    ymax <- max(unlist(lapply(raster.list, function(x) {ext(x)[4]})))
  }
  
  # Create new raster 
  mock.surface <- rast(
    ncol=(xmax-xmin)/resolution$x, # raster automatically rounds, total cols
    nrow=(ymax-ymin)/resolution$y, # raster automatically rounds, total rows
    xmin=xmin, # min x exent
    xmax=xmax, # max x exent
    ymin=ymin, # min y exent
    ymax=ymax, # max y exent
    crs = crs(raster.list[[1]])) # take CRS from first raster, requires that rasters match in CRS
  
  values(mock.surface) <- 1:(ncol(mock.surface)*nrow(mock.surface)) # not important, just for visualization
  
  return(mock.surface)
}

check_ssf <- function(ssf.obj) {
  inherits(ssf.obj, c("fit_clogit")) | inherits(ssf.obj, c("gam")) # not broad enough, basically just from amt at the moment
  
}

step_distance <- function(ssf.obj, quantile) {
  if(!check_ssf(ssf.obj)) stop("Check that SSF model is valid")
  
  if(ssf.obj$sl_$name == "gamma") {
    step <- qgamma(quantile, # user specified quantile
                   shape = ssf.obj$sl_$params$shape, # estimated from amt
                   scale = ssf.obj$sl_$params$scale) # estimated from amt
  }
  
  if(ssf.obj$sl_$name == "exp") {
    step <- qexp(quantile, # user specified quantile
                 rate = ssf.obj$sl_$params$rate) # estimated from amt
  }
  
  if(ssf.obj$sl_$name == "unif") {
    step <- quantile(ssf.obj$model$model$sl_[which(as.character(ssf.obj$model$model[,1])=="1")], 0.95) # estimated from amt
  }
  
  return(step)
}

get_cells <- function(ssf.obj, mock.surface, raster, accessory.x.preds = NULL){
  if(!check_ssf(ssf.obj)) stop("Check that SSF model is valid")
  
  pred.xy <- crds(mock.surface) # get coordinates from grid we created
  
  predict.data <- data.frame(cbind(pred.xy, extract(raster, pred.xy))) # makes raster values a data frame
  
  predict.data$step_id_unique = ssf.obj$model$xlevels$`strata(step_id_)`[1] # fix the strata to something reasonable
  
  if(!is.null(accessory.x.preds)) {
    predict.data <- cbind(predict.data, accessory.x.preds) # this adds extraneous x values that are not in matrix
  }
  
  cells <- nrow(pred.xy) # number of cells
  
  predict.data$cellnr <- 1:cells # assign cell numbers, redundant of ID
  
  return(predict.data)
}


get_cell_data <- function(ssf.obj, pred.data){
  if(!check_ssf(ssf.obj)) stop("Check that SSF model is valid")
  
  cells <- nrow(pred.data) # number of cells
  
  sample <- pred.data %>% 
    drop_na() %>% 
    sample_n(1)
  
  pred.data$sl_ <- sqrt((sample$x - pred.data$x)^2 + (sample$y - pred.data$y)^2) 
  sample$sl_ <- 0
  
  pred.data$x2_ <- pred.data$x
  pred.data$y2_ <- pred.data$y
  
  sample$x2_ <- sample$x
  sample$y2_ <- sample$y
  
  # relying on amt step-selection models, we can predict log-RSS
  # there are better ways to predict, but this is simple
  
  if("gam" %in% class(ssf.obj)) {
    log.rss <- log(c(predict(ssf.obj, newdata = pred.data))/ c(predict(ssf.obj, newdata = sample)))
    full.raster.data <- data.frame(pred.data, lRSS = log.rss)
  } else {
    
    log.rss <- amt::log_rss(ssf.obj, # the model
                            pred.data, # the raster data (including missing values)
                            sample,  # a row of the raster data (excluding missing values)
                            ci = NA) 
    
    if(c(Inf) %in% abs(log.rss$df$log_rss)) {
      pred.data$sl_[which(pred.data$sl_ == 0)] <- 0.01
      sample$sl_[which(sample$sl_ == 0)] <- 0.01
      log.rss <- amt::log_rss(ssf.obj, # the model
                              pred.data, # the raster data (including missing values)
                              sample,  # a row of the raster data (excluding missing values)
                              ci = NA)
    }
    
    full.raster.data <- data.frame(pred.data, lRSS = log.rss$df$log_rss)
    
  }
  # bind predictions to the original data
  
  return(full.raster.data)
  
}

neighbor_lookup <- function(mock.surface, cell.data, cell.data.list = NULL){
  cols <- ncol(mock.surface) # columns in prediction
  rows <- nrow(mock.surface) # rows in prediction
  cells <- cols*rows # number of cells
  index <- 1:cells # all index values in our prediction raster
  
  # create a matrix for each column in the first row and its comparisons to distance to all other cells 
  
  # if(is.null(cell.data.list)){
  #   print("Splitting cell.data into list")
  #   cell.data.list. <- split(cell.data, cell.data$cellnr) # split the prediction data into row-wise lists to use lapply
  # }
  # 
  # if(!is.null(cell.data.list)){
  #   print("Using inputted list of cell data")
  #   cell.data.list. <- cell.data.list # split the prediction data into row-wise lists to use lapply
  # }
  
  # this is a progress bar we can use in a for loop
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]", total = cols, complete = "=", incomplete = "-", current = ">", clear = FALSE, width = 100)
  
  for(i in 1:cols) { # step through columns
    dist <- distance(as.matrix(cell.data[i,c("x","y")]), # choose the first row cell by column (raster indices are row-wise)
                     as.matrix(cell.data[i:cells,c("x","y")]), # choose all other cells
                     lonlat = F) # we are using UTM 
    if(i == 1) mat.dist <- matrix(dist, ncol = 1)
    if(i > 1) mat.dist <- cbind(mat.dist, c(dist, rep(NA, i-1))) # for each additional column, there are i-1 comparisons that are repeated (unnecessary)
    pb$tick() # for progress bar
  }
  
  return(mat.dist)
}

neighbor_finder <- function(ssf.obj = m2, cell.data, neighbors.found, quantile = 0.95, cell.data.list = NULL, distance.override = NULL){
  
  if(is.null(distance.override)) neighborhood.distance <- step_distance(ssf.obj, quantile) # take the X% step distance as your neighborhood 
  
  if(!is.null(distance.override)) neighborhood.distance <- distance.override
  
  cols <- ncol(neighbors.found) # columns of our call-up table
  differences <- nrow(neighbors.found) # number of differences in index values 
  
  print("Creating neighbor comparisons")
  vector <- c(neighbors.found) # convert the neighbors to a vector we can index later, this is the purpose of all those NA's earlier based on i-1 unique distances
  
  print("Finding valid comparisons")
  valid <- vector < neighborhood.distance # T/F whether those neighborhood distances are less than our threshold
  
  if(is.null(cell.data.list)){
    print("Splitting cell.data into list")
    cell.data.list. <- split(cell.data, cell.data$cellnr) # split the prediction data into row-wise lists to use lapply
  }
  
  if(!is.null(cell.data.list)){
    print("Using inputted list of cell data")
    cell.data.list. <- cell.data.list # split the prediction data into row-wise lists to use lapply
  }
  
  print("Running comparisons")
  neighbor.mat <- pblapply(cell.data.list., function(x){ # step through each row (see list split above)
    focal <- as.numeric(x$cellnr) # value of row cell number 
    delta <- abs(focal - cell.data$cellnr) # difference in row cell number vs all others
    index <- ifelse(focal < cell.data$cellnr, focal, cell.data$cellnr) # report the minimum cell index (based on our call-up structure)
    
    index.col <- index%%cols # use the remainder function to get the column number (see below, we have to force zeros to the column number because the remainder of the final column is zero)
    
    df <- data.frame(difference = delta + 1, # we have to add one because differences of 0 are stored in row 1, differences of 1 in row 2, etc. 
                     col = ifelse(index.col == 0, cols, index.col), # forcing remainders of zero the number of columns
                     cell.nr = cell.data$cellnr) # just tracking the cell number we are comparing against for use later
    
    # filter the data frame based on whether our call-up values are less than the neighborhood
    df <- df %>% 
      filter(valid[difference+((col-1)*differences)]) 
    
    # filter the data set to unique rows and columns for call-up (to accommodate memory issues)
    df.distinct <- df %>% distinct(difference, col) 
    
    # find the unique distances (trims time down)
    df.distinct$distances <- vector[df.distinct$difference + ((df.distinct$col - 1)*differences)]
    
    # throw the unique distances back to the full data set 
    df <- merge(df, df.distinct, by = c("difference","col"))
    
    # package into a nice data frame for export
    data.frame(row = focal, column = df$cell.nr, distance = df$distances)
  })
  
  # bind the output list
  neighbors <- rbindlist(neighbor.mat)
  
  # create a sparse matrix based on the focal id, alternate id, and distance
  sparse.neighbors <- sparseMatrix(i = neighbors$row, j = neighbors$column, x = neighbors$distance)
  
  # return both the matrix and the unbound list of neighbor cells
  return(list(matrix = sparse.neighbors, by.cell = neighbor.mat))
}


compile_ssf_comparisons <- function(sparse.neighbors, cell.data) {
  
  # this is why the export of the neighbors as individual lists was important
  ssf.comparisons <- lapply(sparse.neighbors$by.cell, function(x){ 
    baseline <- cell.data[x$column[which(x$distance == 0)],] # baseline will have a distance of zero (focal)
    
    baseline$sl_ <- 0 # create a variable "step" that records this zero distance
    
    alternate <- cell.data[x$column,] # grab all the other cell.data for neighboring cells (including focal cell)
    
    alternate$sl_ <- x$distance # force distance to this new variable step
    
    list(.given = baseline, .for = alternate) # return a list of focal and neigboring cell data
    
  })
  
  return(ssf.comparisons)
}


predict_ssf_comparisons <- function(ssf.obj, ssf.comparisons = ssf.comparisons) {
  
  print("Estimating probability surface")
  
  ## this is straight from amt - i just dont want to recalculate the uncenter term every exposure because it saves about 90% of the time of this function
  uncenter <- sum(coef(ssf.obj$model) * ssf.obj$model$means, na.rm=TRUE)
  
  prediction.list <- pblapply(ssf.comparisons, function(x){ # step through the list of SSF objects
    x1_dummy <- x$.for
    x2_dummy <- x$.given
    
    x1_dummy$step_id_ = ssf.obj$model$model$`strata(step_id_)`[1]
    x2_dummy$step_id_ = ssf.obj$model$model$`strata(step_id_)`[1]
    
    x1_dummy$sl_[which(x1_dummy$sl_ == 0)] <- 0.001
    x2_dummy$sl_[which(x2_dummy$sl_ == 0)] <- 0.001
    
    #Calculate y_x
    pred_x1 <- predict(ssf.obj$model, newdata = x1_dummy, type = "lp", reference = "sample",
                       se.fit = F)
    pred_x2 <- predict(ssf.obj$model, newdata = x2_dummy, type = "lp", reference = "sample",
                       se.fit = F)
    
    y_x1 <- pred_x1 + uncenter
    y_x2 <- pred_x2 + uncenter
    
    log_rss <- unname(y_x1 - y_x2)
    x$.for$Prob <- exp(log_rss)*(1/sum(exp(log_rss), na.rm = T)) # exponentiate and multiply against relative risk
    
    x$.for # return the data frame with probabilities
  })
  
  print("Compiling probability surface")
  
  
  
  rows <- unlist(lapply(prediction.list, nrow))
  
  rows2 <- cumsum(rows)
  rows1 <- c(1,rows2[-length(rows2)] + 1)
  
  focal <- cellnr <- prob <- rep(NA, sum(unlist(lapply(prediction.list, nrow))))
  for(i in 1:length(prediction.list)){
    focal[rows1[i]:rows2[i]] <- i # specify the focal cell for each comparison
    cellnr[rows1[i]:rows2[i]] <- prediction.list[[i]]$cellnr # specify the alt id cell for each comparison
    prob[rows1[i]:rows2[i]] <- prediction.list[[i]]$Prob # specify the prob for each comparison
  } 
  
  print("Making sparse matrix for transitions")
  # bound <- rbindlist(prediction.list) # bind all data frames 
  
  
  
  # use indexing to make a massive sparse matrix quickly 
  # Sparse.Matrix.lrss <- sparseMatrix(bound$focal.cell, bound$cellnr, x = bound$log_rss)
  # Sparse.Matrix.rss <- sparseMatrix(bound$focal.cell, bound$cellnr, x = bound$rss)
  Sparse.Matrix <- sparseMatrix(focal, cellnr, x = prob)
  # Sparse.Matrix.l <- sparseMatrix(bound$focal.cell, bound$cellnr, x = bound$Prob.l) 
  # Sparse.Matrix.h <- sparseMatrix(bound$focal.cell, bound$cellnr, x = bound$Prob.h) 
  
  # return the prediction list and sparse matrix
  return(list(#prob.surface = prediction.list, 
    # lrss.matrix = Sparse.Matrix.lrss, rss.matrix = Sparse.Matrix.rss, 
    prob.matrix = Sparse.Matrix
    # , sparse.matrix = Sparse.Matrix, sparse.matrix.l = Sparse.Matrix.l, sparse.matrix.h = Sparse.Matrix.h
  ))  
}


eigen_decompose_to_raster <- function(surface, mock.surface){
    A <- surface$prob.matrix
    d <- eigs(t(A), 1)
    d.1 <- Re(d$vectors[,1])
    prob.d <- d.1/sum(d.1)
    ssd.prob.raster <- setValues(mock.surface, prob.d)
    return(ssd.prob.raster)
}



simulate_movements <- function(raster, lambda = 1, kappa = 1, step_no = 1E4, teleport = NULL){
  
  set.seed(100)
  
  locations <- matrix(NA, nrow = step_no, ncol = 2) # storage matrix
  
  spts <- cbind(terra::xyFromCell(raster, 1:(100*100)),terra::extract(raster, 1:(100*100))) # get points from matrix for simulation
  
  # Simulate path
  alpha_x = 0 # initial turn angle
  
  grid <- expand.grid(x = 1:100, y = 1:100) # 100 x 100 landscape
  grid$z <- 1
  sim <- terra::rast(grid)
  
  if(is.null(teleport)) {
    for(step in 1:step_no) # go through steps
    {
      if(step == 1) {
        newxy<-sample(1:nrow(spts), 1, prob=exp(spts$z.1)) # if its the first step, derive random start based on suitability
        locations[step,] <- xyFromCell(sim, newxy) # put this cell in storage and jump to the next step
        next
      }
      
      alpha_z <- atan2(spts$y-locations[step-1,2],spts$x-locations[step-1,1]) # calculate difference in turn angle for all cells based on possible steps
      
      unnorm_mk<-exp(-lambda*sqrt((spts$x-locations[step-1,1])^2+(spts$y-locations[step-1,2])^2) + # negatively weighted distance
                       spts$z.1 + # directly weighted suitability surface
                       kappa*cos(alpha_x-alpha_z)) # Positively weighted small differences in turn angle
      
      mk<-unnorm_mk/sum(unnorm_mk) # make the surface a kernel
      
      newxy<-sample(1:nrow(spts), 1, prob=mk) # draw sample
      
      locations[step,] <- xyFromCell(sim, newxy) # store sampled location
      
      alpha_x<-atan2(spts$y-locations[step-1,2],spts$x-locations[step-1,1]) # calculate difference in turn angle based on the step taken
    }
  } else {
    for(step in 1:step_no) # go through steps
    {
      if(step %in% seq(1,step_no, by = teleport)) { # randomly restart the track every 1000 steps
        newxy<-sample(1:nrow(spts), 1, prob=exp(spts$z.1)) # if its the first or restart steps, derive random starts based on suitability
        locations[step,] <- xyFromCell(sim, newxy) # put this cell in storage and jump to the next step
        alpha_x <- 0 # reset turn angle to 0
        next
      }
      
      alpha_z <- atan2(spts$y-locations[step-1,2],spts$x-locations[step-1,1]) # calculate difference in turn angle for all cells based on possible steps
      
      unnorm_mk<-exp(-lambda*sqrt((spts$x-locations[step-1,1])^2+(spts$y-locations[step-1,2])^2) + # negatively weighted distance
                       spts$z.1 + # directly weighted suitability surface
                       kappa*cos(alpha_x-alpha_z)) # Positively weighted small differences in turn angle
      
      mk<-unnorm_mk/sum(unnorm_mk) # make the surface a kernel
      
      newxy<-sample(1:nrow(spts), 1, prob=mk) # draw sample
      
      locations[step,] <- xyFromCell(sim, newxy) # store sampled location
      
      alpha_x<-atan2(spts$y-locations[step-1,2],spts$x-locations[step-1,1]) # calculate difference in turn angle based on the step taken
    }
  }
  
  return(locations)
}


simulate_movements_tallied <- function(raster, lambda = 1, kappa = 1, step_no = 1E4, teleport = NULL){
  
  set.seed(100)
  
  locations <- matrix(NA, nrow = step_no, ncol = 2) # storage matrix
  
  spts <- cbind(terra::xyFromCell(raster, 1:(100*100)),terra::extract(raster, 1:(100*100))) # get points from matrix for simulation
  
  # Simulate path
  alpha_x = 0 # initial turn angle
  
  grid <- expand.grid(x = 1:100, y = 1:100) # 100 x 100 landscape
  grid$z <- 1
  sim <- terra::rast(grid)
  
  if(is.null(teleport)) {
    for(step in 1:step_no) # go through steps
    {
      if(step == 1) {
        newxy<-sample(1:nrow(spts), 1, prob=exp(spts$z.1)) # if its the first step, derive random start based on suitability
        locations[step,] <- xyFromCell(sim, newxy) # put this cell in storage and jump to the next step
        next
      }
      
      alpha_z <- atan2(spts$y-locations[step-1,2],spts$x-locations[step-1,1]) # calculate difference in turn angle for all cells based on possible steps
      
      unnorm_mk<-exp(-lambda*sqrt((spts$x-locations[step-1,1])^2+(spts$y-locations[step-1,2])^2) + # negatively weighted distance
                       spts$z.1 + # directly weighted suitability surface
                       kappa*cos(alpha_x-alpha_z)) # Positively weighted small differences in turn angle
      
      mk<-unnorm_mk/sum(unnorm_mk) # make the surface a kernel
      
      newxy<-sample(1:nrow(spts), 1, prob=mk) # draw sample
      
      locations[step,] <- xyFromCell(sim, newxy) # store sampled location
      
      alpha_x<-atan2(spts$y-locations[step-1,2],spts$x-locations[step-1,1]) # calculate difference in turn angle based on the step taken
    }
  } else {
    for(step in 1:step_no) # go through steps
    {
      if(step %in% seq(1,step_no, by = teleport)) { # randomly restart the track every 1000 steps
        newxy<-sample(1:nrow(spts), 1, prob=exp(spts$z.1)) # if its the first or restart steps, derive random starts based on suitability
        locations[step,] <- xyFromCell(sim, newxy) # put this cell in storage and jump to the next step
        alpha_x <- 0 # reset turn angle to 0
        next
      }
      
      alpha_z <- atan2(spts$y-locations[step-1,2],spts$x-locations[step-1,1]) # calculate difference in turn angle for all cells based on possible steps
      
      unnorm_mk<-exp(-lambda*sqrt((spts$x-locations[step-1,1])^2+(spts$y-locations[step-1,2])^2) + # negatively weighted distance
                       spts$z.1 + # directly weighted suitability surface
                       kappa*cos(alpha_x-alpha_z)) # Positively weighted small differences in turn angle
      
      mk<-unnorm_mk/sum(unnorm_mk) # make the surface a kernel
      
      newxy<-sample(1:nrow(spts), 1, prob=mk) # draw sample
      
      locations[step,] <- xyFromCell(sim, newxy) # store sampled location
      
      alpha_x<-atan2(spts$y-locations[step-1,2],spts$x-locations[step-1,1]) # calculate difference in turn angle based on the step taken
    }
  }
  
  indexed <- setValues(sim, 1:ncell(sim)) # store raster index as an alternate raster
  names(indexed) <- "index"
  
  indexed <- as.data.frame(indexed, xy = TRUE)
  colnames(indexed) <- c("x","y", "index")
  
  locations <- as.data.frame(locations)
  colnames(locations) <- c("x","y")
  
  merged <- merge(indexed, locations, by = c("x", "y")) # merge map with use
  
  tallied <- merged %>%
    group_by(x, y, index) %>%
    summarize(n = n()) # observed selected steps per cell
  
  not.used <- indexed$index[which(!indexed$index %in% tallied$index)] # cells are not used, meaning zero counts
  not.used.df <- tibble(index = not.used, n = 0) # we can make a mock df of these unused points
  
  tallied <- rbind(tallied, not.used.df)
  
  return(tallied)
}


select_right_surface <- function(sims, index){
   landscape.to.use <- index %% 50
   landscape.to.use <- ifelse(landscape.to.use == 0, 50, landscape.to.use)
   sim <- rast(sims[[landscape.to.use]])
   return(sim)
}

select_right_surface2 <- function(sims, index){
   landscape.to.use <- index %% 50
   landscape.to.use <- ifelse(landscape.to.use == 0, 50, landscape.to.use)
   sim <- rast(sims[[landscape.to.use]])
   return(sim)
}

generate_ssf <- function(surfaces, movement, index) {

    locations.df <- as.data.frame(movement)
    
   locations.df$t <- as.Date(1:nrow(locations.df), origin = "1970-01-01")
   trk <- amt::make_track(locations.df, .x = V1, .y = V2, .t = t)

   ssf.dat <- trk %>%
     mutate(x_ = x_ + runif(nrow(trk), -0.001, 0.001),
            y_ = y_ + runif(nrow(trk), -0.001, 0.001)) %>%
     steps()

   fatter.unif <- original.unif <- fit_distr(ssf.dat$sl_, "unif") # get the original unif distribution
   fatter.unif$params$min <- 0 # take the min to 0
   fatter.unif$params$max <-  fatter.unif$params$max*1.25 # bump up the max
   attributes(fatter.unif)$class[2] <- "sl_distr"


   landscape.to.use <- index %% 50
   landscape.to.use <- ifelse(landscape.to.use == 0, 50, landscape.to.use)

   sim <- rast(surfaces[[landscape.to.use]])

   ssf.dat. <- ssf.dat %>%
     random_steps(100,
                  sl_distr = fatter.unif) %>% # implement the fatter uniform distribution for random steps
     extract_covariates(sim)

   ssf.fit <- fit_issf(ssf.dat. %>%
                        filter(!is.na(z.1),
                               !is.na(z.1),
                               !is.na(z.3),
                               !is.na(z.4)) , case_ ~ (z.1 + z.2 + z.3 + z.4) + sl_ + log(sl_) + strata(step_id_), model = T)
                               
   ssf.fit.ta <- fit_issf(ssf.dat. %>%
                        filter(!is.na(z.1),
                               !is.na(z.1),
                               !is.na(z.3),
                               !is.na(z.4)) , case_ ~ (z.1 + z.2 + z.3 + z.4) + sl_ + log(sl_) + cos(ta_) + strata(step_id_), model = T)
                               
    return(list(base.ssf = ssf.fit, ta.ssf = ssf.fit.ta))
}


generate_rsf <- function(surfaces, movement, index) {

    locations.df <- as.data.frame(movement)
    
   locations.df$t <- as.Date(1:nrow(locations.df), origin = "1970-01-01")
   trk <- amt::make_track(locations.df, .x = V1, .y = V2, .t = t)

   landscape.to.use <- index %% 50
   landscape.to.use <- ifelse(landscape.to.use == 0, 50, landscape.to.use)

   sim <- rast(surfaces[[landscape.to.use]])

   rsf <- trk %>% 
      random_points() 
      
    rsf <- cbind(rsf, extract(sim, rsf[,2:3]))

   rsf.fit <- glm(case_ ~ z.1 + z.2 + z.3 + z.4, data = rsf, family = "binomial")
                               
    return(rsf.fit)
}




generate_ssf_map <- function(pred.data, movement, model, mock.surface) {

   locations.df <- as.data.frame(movement)
    
   locations.df$t <- as.Date(1:nrow(locations.df), origin = "1970-01-01")
   trk <- amt::make_track(locations.df, .x = V1, .y = V2, .t = t)

   ssf.dat <- trk %>%
     mutate(x_ = x_ + runif(nrow(trk), -0.001, 0.001),
            y_ = y_ + runif(nrow(trk), -0.001, 0.001)) %>%
     steps()

   a.data <- pred.data #grab the prediction data we made
    a.data$sl_ <- mean(ssf.dat$sl_,na.rm = T) #use the mean step-length for our model
    b.data <- a.data[1,]
    log.rss <- amt::log_rss(model, # the model
                              a.data, # the raster data (including missing values)
                              b.data,  # a row of the raster data (excluding missing values)
                              ci = NA)
    ssf. <-  exp(log.rss$df$log_rss)/sum(exp(log.rss$df$log_rss)) # take out of link-scale and turn into kernel
    ssf.prob.raster <- setValues(mock.surface, ssf.)
                               
    return(wrap(ssf.prob.raster))
}

generate_rsf_map <- function(pred.data, model, mock.surface) {


   probabilities <- exp(predict(model, newdata = pred.data))/(1+exp(predict(model, newdata = pred.data))) # take the rsf output into probabilities

    probabilities.kern <- probabilities/sum(probabilities) # turn into kernel
    rsf.prob.raster <- setValues(mock.surface, probabilities.kern)
                               
    return(wrap(rsf.prob.raster))
}



generate_ssf_quantile <- function(ssf.movement, level = 0.95) {
    locations.df <- as.data.frame(ssf.movement)
    
   locations.df$t <- as.Date(1:nrow(locations.df), origin = "1970-01-01")
   trk <- amt::make_track(locations.df, .x = V1, .y = V2, .t = t)

   ssf.dat <- trk %>%
     mutate(x_ = x_ + runif(nrow(trk), -0.001, 0.001),
            y_ = y_ + runif(nrow(trk), -0.001, 0.001)) %>%
     steps()
     
     limit <- quantile(ssf.dat$sl_, level)
    return(limit)
}





generate.abm <- function(model = issf.fit2, spatial.data = sim, ssf.dat = ssf.dat., steps = 1E3, n.control = 50) {
  
  #find valid cells
  valid <- !is.na(spatial.data[[1]])
  
  #convert to shapefiles
  outline <- as.polygons(valid)
  outline <- outline[which(outline[[1]]==1),]
  
  #grab a random location within
  startingxy <- st_coordinates(st_sample(st_as_sf(outline), 1))
  
  # fix to starting value
  start <- make_start(ssf.dat[1,])
  start$x_ <- startingxy[1]
  start$y_ <- startingxy[2]
  
  n.control = n.control
  
  # set redistribution kernel 
  k1 <- redistribution_kernel(model, map = spatial.data, start = start,
                              landscape = "continuous", tolerance.outside = 0.99, 
                              n.control = n.control,fun = function(xy, map) {
                                extract_covariates(xy, map, where = "end")
                              })
  
  s <- simulate_path(k1, n.steps = steps)
  
  return(tracks = s[-c(1:100),1:2])
}

aggregate.abm <- function(abm.list, spatial.data) {
  
  bound <- rbindlist(abm.list)
  
  bound.agg <- bound %>% 
    as_track() %>% 
    steps() %>% 
    extract_covariates(spatial.data) %>% 
    group_by(index) %>% 
    summarise(n = n())
  
  simulation <- setValues(spatial.data$index, 0)
  values(simulation)[bound.agg$index] <- bound.agg$n/sum(bound.agg$n)
  return(simulation)
}


generate.abm.r <- function(model = issf.fit2, spatial.data = sim, ssf.dat = ssf.dat., steps = 5E3, n.control = 100) {
  
  #find valid cells
  valid <- !is.na(spatial.data[[1]])
  
  #convert to shape files
  outline <- as.polygons(valid)
  outline <- outline[which(outline[[1]]==1),]
  
  #grab a random location within
  startingxy <- st_coordinates(st_sample(st_as_sf(outline), 1))
  
  # fix to starting value
  start <- make_start(ssf.dat[2,])
  start$x_ <- startingxy[1]
  start$y_ <- startingxy[2]
  
  n.control = n.control
  
  # set redistribution kernel 
  k1 <- redistribution_kernel(model, map = spatial.data, start = start,
                              landscape = "continuous", tolerance.outside = 0.99, 
                              n.control = n.control, fun = function(xy, map) {
                                extract_covariates(xy, map, where = "end")
                              })
  
  s <- simulate_path(k1, n.steps = steps)
  
  s.agg <- s %>% 
    slice(-c(1:100)) %>% 
    as_track() %>% 
    steps() %>% 
    extract_covariates(spatial.data) %>% 
    group_by(index) %>% 
    summarise(n = n())
  simulation <- setValues(spatial.data$index, 0)
  values(simulation)[s.agg$index] <- s.agg$n/sum(s.agg$n)
  return(raster = wrap(simulation))
}



run_abm_2_convergence <- function(model = issf.fit, spatial.data = sim, ssf.dat = ssf.movement, iterator){

    set.seed(iterator)
    locations.df <- as.data.frame(ssf.dat)

   locations.df$t <- as.Date(1:nrow(locations.df), origin = "1970-01-01")
   trk <- amt::make_track(locations.df, .x = V1, .y = V2, .t = t)

   ssf.dat <- trk %>%
     mutate(x_ = x_ + runif(nrow(trk), -0.001, 0.001),
            y_ = y_ + runif(nrow(trk), -0.001, 0.001)) %>%
     steps()

    raster <- spatial.data
    raster$index <- setValues(raster[[1]], 1:ncell(raster))
    #raster.wrap <- wrap(raster)
    
    model <- model
    
    
    #print("making cluster")
    #cl <- makeCluster(24)
    #clusterEvalQ(cl, source("base_fxns.R"))
    #clusterExport(cl, varlist = c("model","raster.wrap","ssf.dat","steps","n.control"), environment())
    #print("cluster made")
    #raster.list <- parLapply(cl = cl, 1:tracks, function(x){
    #    set.seed(x)
    #    raster <- rast(raster.wrap)
    #    return(tryCatch(generate.abm.r(model, raster, ssf.dat, steps = steps, n.control = n.control), error=function(e) NULL))
    #    print(x)
    #})
    #stopCluster(cl)


    #return(abm.raster = raster.list)
    
    return(generate.abm.r(model, raster, ssf.dat, steps = 5000, n.control = 100))
    
}




redistribution_kernel2 <- function (x = make_issf_model(), start = make_start(), map, outline, fun = function(xy, 
    map) {
    extract_covariates(xy, map, where = "both")
}, covars = NULL, max.dist = get_max_dist(x), n.control = 1e+06, 
    n.sample = 1, landscape = "continuous", compensate.movement = landscape == 
        "discrete", normalize = TRUE, interpolate = FALSE, as.rast = FALSE, 
    tolerance.outside = 0) 
{
    arguments <- as.list(environment())
    checkmate::assert_class(start, "sim_start")
    if (!landscape %in% c("continuous", "discrete")) {
        stop("Argument `landscape` is invalid. Valid values are 'continuous' or 'discrete'.")
    }
    if (landscape == "continuous") {
        xy <- random_steps_simple(start, sl_model = x$sl_, ta_model = x$ta_, 
            n.control = n.control)
    }
    else {
        xy <- kernel_setup(map, max.dist, start, covars)
    }

    possible <- vect(cbind(xy$x2_, xy$y2_), crs = crs(outline))

    fraction.outside <- mean(is.na(terra::extract(outline, possible)[,2]))
    
    if (fraction.outside > tolerance.outside) {
        warning(paste0(round(fraction.outside * 100, 3), "% of steps are ending outside the study area but only ", 
            round(tolerance.outside * 100, 3), "% is allowed. ", 
            "Terminating simulations here."))
        return(NULL)
    }
    xy$t1_ <- start$t_
    xy$t2_ <- start$t_ + start$dt
    xy <- fun(xy, map)
    w <- amt:::ssf_weights(xy, x, compensate.movement = compensate.movement)
    r <- if (!as.rast) {
        dplyr::select(xy[sample.int(nrow(xy), size = n.sample, 
            prob = w), ], x_ = x2_, y_ = y2_, t2_)
    }
    else {
        if (landscape == "continuous") {
            stop("`as.rast` not implemented for `landscape = 'continuous'`")
        }
        else {
            terra::rast(data.frame(xy[, c("x2_", "y2_")], w))
        }
    }
    if (as.rast & normalize) {
        r <- normalize(r)
    }
    res <- list(args = arguments, redistribution.kernel = r)
    class(res) <- c("redistribution_kernel", "list")
    res
}


generate.abm2 <- function(model = issf.fit2, spatial.data = sim, ssf.dat = ssf.dat., steps = 1E3, n.control = 50) {
  
  #find valid cells
  valid <- !is.na(spatial.data[[1]])
  
  #convert to shapefiles
  outline <- as.polygons(valid)
  outline <- outline[which(outline[[1]]==1),]
  
  #grab a random location within
  startingxy <- st_coordinates(st_sample(st_as_sf(outline), 1))
  
  # fix to starting value
  start <- make_start(ssf.dat[2,])
  start$x_ <- startingxy[1]
  start$y_ <- startingxy[2]
  
  n.control = n.control
  
  # set redistribution kernel 
  k1 <- redistribution_kernel2(model, map = spatial.data, start = start, outline = outline, 
                              landscape = "continuous", tolerance.outside = 0.99, 
                              n.control = n.control,fun = function(xy, map) {
                                extract_covariates(xy, map, where = "end")
                              })
  
  s <- simulate_path(k1, n.steps = steps)
  
  return(tracks = s[-c(1:100),1:2])
}


generate.abm.r2 <- function(model = issf.fit2, spatial.data = sim, ssf.dat = ssf.dat., steps = 5E3, n.control = 100) {
  
  #find valid cells
  valid <- !is.na(spatial.data[[1]])
  
  #convert to shape files
  outline <- as.polygons(valid)
  outline <- outline[which(outline[[1]]==1),]
  
  #grab a random location within
  startingxy <- st_coordinates(st_sample(st_as_sf(outline), 1))
  
  # fix to starting value
  start <- make_start(ssf.dat[2,])
  start$x_ <- startingxy[1]
  start$y_ <- startingxy[2]
  
  n.control = n.control
  
  # set redistribution kernel 
  k1 <- redistribution_kernel2(model, map = spatial.data, start = start, outline = outline,
                              landscape = "continuous", tolerance.outside = 0.99, 
                              n.control = n.control, fun = function(xy, map) {
                                extract_covariates(xy, map, where = "end")
                              })
  
  s <- simulate_path(k1, n.steps = steps)
  
  s.agg <- s %>% 
    slice(-c(1:100)) %>% 
    as_track() %>% 
    steps() %>% 
    extract_covariates(spatial.data) %>% 
    group_by(index) %>% 
    summarise(n = n())
  simulation <- setValues(spatial.data$index, 0)
  values(simulation)[s.agg$index] <- s.agg$n/sum(s.agg$n)
  return(raster = wrap(simulation))
}


run_abm_2_convergence2 <- function(model = issf.fit, spatial.data = sim, ssf.dat = ssf.movement, iterator){

    set.seed(iterator)
    
    raster <- (spatial.data)
    raster$index <- setValues(raster[[1]], 1:ncell(raster))
    
    return(generate.abm.r2(model, raster, ssf.dat, steps = 5000, n.control = 100))
    
}



compare_spearman <- function(data = ssf1.train, rstack = deer1.rstack, bootstrap = 100, bins = 10){
  rasters <- names(rstack) # layers of predicted surface
  intersected <- data %>% 
    extract_covariates(rstack) %>% 
    drop_na() # pull out all the predictions per point
  
  predictions <- values(rstack) # grab the matrix of predicted surfaces
  
  interval <- apply(predictions, 2, function(x) {
    quantile(x, seq(0, 1, length.out = bins + 1), na.rm = T)
    # seq(min(x), max(x), length.out = bins + 1)
  }) # generate centiles of the predicted surfaces
  
  
  unique.bins <- apply(interval, 2, function(x) length(unique(x)))<(bins + 1)
  
  if(T %in% unique.bins){
    print("unique bins less than requested bins, decreasing bins to minimum")
    
    index <- which(unique.bins)
    
    for(i in index){
      interval[which(duplicated(interval[,i])),i] <- NA
    }
    
  }
  
  # df.area <- expand.grid(bin = 1:bins,
  #                        layer = 1:length(rasters)) # make a storage data frame
  # df.area$area <- NA # add the outcome 
  # total.size <- nrow(predictions)
  # for(i in 1:length(rasters)){
  #   for(j in 2:(bins + 1)){
  #     df.area[j-1,i] <- sum(predictions[,i] <= interval[j,i] & predictions[,i] > interval[j-1,i])/total.size
  #   }
  # }
  
  df <- expand.grid(bin = 1:bins,
                    name = 1:length(rasters),
                    bootstrap = 1:bootstrap) # make a storage data frame
  df$intensity <- NA # add the outcome 
  df$bin.numeric <- NA # add the outcome 
  
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]", total = bootstrap, complete = "=", incomplete = "-", current = ">", clear = FALSE, width = 100)
  
  for (i in 1:bootstrap) { # for each iteration
    sampled <- intersected %>% 
      sample_n(nrow(intersected), replace = T) # resample with replacement to the full dataset size
    for (j in 2:(bins + 1)){ # step through centiles
      for (k in 1:length(rasters)){ # step through predicted surfaces
        true <- sampled[,rasters[k]] <= interval[j,k] & sampled[,rasters[k]] > interval[j-1,k] # check whether each point is in the focal centile
        df[which(
          df$bin == (j-1) &
            df$name == k &
            df$bootstrap == i
        ),"intensity"] <- sum(true, na.rm = T)/(nrow(intersected)) # store the percentage of bootstrapped points in interval, and the proportion of the landscape in that interval
        
        df[which(
          df$bin == (j-1) &
            df$name == k &
            df$bootstrap == i
        ),"bin.numeric"]  <- interval[j,k]
      }
    }
    pb$tick()
  }
  
  df <- df %>% 
    filter(!is.na(bin.numeric)) %>% 
    group_by(name, bootstrap) %>%  
    summarize(measure = cor(bin, intensity, method = "spearman")) %>% 
    mutate(name = factor(name))
  return (df) # return dataframe
}


compare_mean <- function(data = ssf1.train, rstack = deer1.rstack, bootstrap = 100){
  rasters <- names(rstack) # layers of predicted surface
  
  intersected <- data %>% 
    extract_covariates(rstack) %>% 
    drop_na() 
  
  store <- rep(list(NA), length = bootstrap)
  for(i in 1:bootstrap){
    store[[i]] <- intersected %>% 
      sample_n(nrow(intersected), replace = T) %>%
      pivot_longer(all_of(rasters)) %>% 
      mutate(bootstrap = i) %>% 
      group_by(name, bootstrap) %>% 
      summarize(measure = exp(mean(log(value))), .groups = "keep")%>% 
      mutate(name = factor(name))
  }
  
  return (rbindlist(store)) # return dataframe
}

sqrt.both.sides <- function(x) {
  new <- ifelse(x == 0, 0, sqrt(abs(x)) * x/abs(x))
}
rev.sqrt.both.sides <- function(x) {
  new <- ifelse(x == 0, 0, x^2 * x/abs(x))
}
transform_both.sqrt.trans_trans <- function() trans_new("both.sqrt.trans", sqrt.both.sides, rev.sqrt.both.sides)


large.power <- function(x) {
  new <- x^5
  }
rev.large.power <- function(x) {
  new <- x^(1/5)
  }
transform_large.power_trans <- function() trans_new("large.power", large.power, rev.large.power)
