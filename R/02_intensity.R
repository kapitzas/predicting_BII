# Preproccesing

# This script estimates parameters of the relationship between land use type and intensity at 0.5 degree resolution with population intensity and fractional land use cover (following Newbold et al., 2015). The modelling is done using gam. The parameters are then applied to 0.5 arcmin resolution population density and fractional cover, to obtain a global map of land use type and intensity at that resoltuion. 



rm(list = ls())

#----------------------------#
#### 0. packages, folders ####
#----------------------------#


# Packages
# devtools::install_github("kapitzas/gdaltools") # package to directly link R to rgdal
# devtools::install_github("kapitzas/flutes") # integerify function

require(devtools)
require(WorldClimTiles)
require(rgdal)
require(gdaltools)
require(raster)
require(stringr)
require(rgeos)
require(sf)
require(mgcv)
library(BBmisc)
require(data.table)
library("doParallel")


# Folders
# setwd("/home/student.unimelb.edu.au/kapitzas/ch3") only on boab
raw_path <-  file.path(path.expand("~"), "OneDrive - The University of Melbourne", "PhD - Large Files", "PhD - Raw Data")
data_path <- file.path(".", "data")
temp_path <- file.path(data_path, "temp")
int_path <- file.path(data_path, "lu_intensity")

#-----------------------#
#### 1. Build models ####
#----------------------## 

# 1. a Create modelling dataframe

# Load mask, extract non-NA indices 
mask_30min <- raster(file.path(processed_path, "mask_30min.tif"))
inds_30min <- which(!is.na(mask_30min[]))


# Load 0.5 degree data for modelling into data frame
files_30min <- list.files(processed_path, pattern = "30min", full.names = TRUE)
files_30min <- grep(paste(c("predicts", "unsubregions", "base_pop"),collapse="|"), files_30min, value=TRUE)

gam_data <- as.data.frame(stack(files_30min))[inds_30min,]
landuses <- c("cropland", "pasture", "primary", "secondary", "Urban")
intensities <- c("minimal", "light", "intense")
colnames(gam_data)
# turn response columns (fractions of land use type and intensity classes) into integers for gam

for(i in 1:length(landuses)){
  lu_inds <- grep(paste(paste(landuses[i], intensities, sep = "_"), collapse = "|"), colnames(gam_data))
  gam_data[landuses[i]] <- rowSums(gam_data[,lu_inds])
  gam_data[,lu_inds] <- integerify(gam_data[,lu_inds], resolution = 1000)
}

# Combine pacific micro states into one region (otherwise there is not enough data to build models on some un subregions)
gam_data[is.na(gam_data)] <- 0
gam_data$unsubregions_30min[gam_data$unsubregions_30min%in%c(57, 54, 61,29,155)] <- 99
gam_data$unsubregions_30min <- as.factor(gam_data$unsubregions_30min) # this is random effect, need to trun into factor


# 1. b Build models

# Prepare CPU cluster
cl <- makeCluster(5)
registerDoParallel(cl)

# Log file to keep track
logfile <- file.path(data_path, paste0("gam_log.txt"))
writeLines(c(""), logfile)

results <- foreach(i = 1:5, .packages = c("mgcv")) %dopar% {
  
  cat(paste(landuses[i],"\n"), file = logfile, append = T)
  
  lus <- landuses[i] # land use classes types
  
  # Two different formulas, one for minimal or not, one for lihgt or intense.
  f_minimal_or_not <- as.formula(paste("minimal_or_not ~", paste(c("unsubregions_30min", paste0("s(", c("base_pop_30min", lus), ")"), paste0("s(", c(lus, "base_pop_30min"), ",by = unsubregions_30min)")), collapse = "+")))
  f_light_or_intense <- as.formula(paste("light_or_intense ~", paste(c("unsubregions_30min", paste0("s(", c("base_pop_30min", lus), ")"), paste0("s(", c(lus, "base_pop_30min"), ",by = unsubregions_30min)")), collapse = "+")))
  
  # Get columns matching land use type and intensity classes
  lu_inds <- grep(paste(paste(landuses[i], intensities, sep = "_"), collapse = "|"), colnames(gam_data))
  
  
  # Modelling probability for minimial or not minimal (first model)
  
  cat(paste("minimal or not: ",landuses[i],"\n"), file = logfile, append = T)
  
  # Response data matrix: When we only two intensities for a lnad use type, get minimal and not minimal (whichever it is).
  if(length(lu_inds) == 2){
    minimal_or_not <- cbind(gam_data[,lu_inds][,1], gam_data[,lu_inds][, 2])
  }
  
  # When we all three intensities for a lnad use type, get minimal and the sum of those that are not not minimal to model against
  if(length(lu_inds) == 3){
    minimal_or_not <- cbind(gam_data[,lu_inds][,1], rowSums(gam_data[,lu_inds][,-1]))
  }
  
  m_minimal_or_not <- tryCatch(mgcv::gam(f_minimal_or_not, family = binomial, data = gam_data), error = function(e) NA)
  if(is.na(m_minimal_or_not)) {cat(paste("returned NA - minimal or not: ",landuses[i],"\n"), file = logfile, append = T)}
  
  # Modelling probability for light or intense, when not minimal (second model)
  
  # only when we have all three classes in a type we need to do this
  if(length(lu_inds)==3){
    cat(paste("light or intense: ",landuses[i],"\n"), file = logfile, append = T)
    
    # Response data matrix: all intensities except minimal
    light_or_intense <- as.matrix(gam_data[,lu_inds][,-1])
    
    m_light_or_intense <- tryCatch(mgcv::gam(f_light_or_intense, family = binomial, data = gam_data), error = function(e) NA)
    if(is.na(m_light_or_intense)) {cat(paste("returned NA - light or intense: ",landuses[i],"\n"), file = logfile, append = T)}
  } else {
    m_light_or_intense <- NA 
  }
  
  # Return the intensity models for each land use type
  list(m_minimal_or_not, m_light_or_intense, landuses[i])
}

# Save
saveRDS(results, file.path(data_path, "lumodels_2.rds"))

stopCluster(cl)

#---------------------#
#### 2. Validation ####
#---------------------#

# 2. a Load intensity models
result <- readRDS(file.path(data_path, "lumodels_2.rds"))

valid_df <- list()

# 2. b Predict intensity for each land use type
for(i in 1:length(landuses)){
  
  # Load models for current land use type
  m_minimal_or_not <- result[[i]][[1]]
  m_light_or_intense <- result[[i]][[2]]
  
  # predict first model (probability of being minimal)
  minimal <- predict(m_minimal_or_not, type = "response", data = gam_data)
  not_minimal <- 1 - minimal
  
  # predict second model (probability of being light, *if not minimal*)
  if(length(m_light_or_intense) != 1){
    light_if_not_minimal <- predict(m_light_or_intense, type = "response", data = gam_data)
    intense_if_not_minimal <- 1 - light_if_not_minimal
    light <- light_if_not_minimal * not_minimal
    intense <- intense_if_not_minimal * not_minimal
    predicted <- cbind(minimal, light, intense)
  } else {
    predicted <- cbind(minimal, not_minimal)
  }
  
  # Now we have the probability of being minimal, light or intense, given our model covariates pop dens and fractional land use type (wihtout intnesity)
  
  # multiply each of the land use type columns with the according two or three predicted probability columns to get land use type and intensity map
  valid_out <- sweep(predicted, 1, as.matrix(gam_data[landuses[i]]), FUN = "*")
  colnames(valid_out) <- paste0(landuses[i], "_", colnames(valid_out))
  valid_df[[i]] <- valid_out
}

# combine to data frame
valid_df <- do.call("cbind", valid_df)

# 2. c Map errors

# Get observed data (the data set we created from gls and harmonized lu data)
obs_df <- as.data.frame(readRDS("/Users/simon/OneDrive - The University of Melbourne/PhD/chapter3/data/lu_predicts_30min.rds"))[inds_30min,]
obs_df - valid_df

# Caluclate and map RMSE
error <- as.matrix(obs_df) - as.matrix(valid_df)
error <- sqrt(rowMeans(error^2))
error_map <- mask_30min
error_map[inds_30min] <- error
spplot(test)


#-------------------------------------#
#### 3. Prediction of fine res map ####
#-------------------------------------#

# 3. a Make prediction data frame with high res data

# Load 0.5 arcmin mask and extract indices
mask_30sec <- readRDS(file.path(data_path, "mask_30sec.rds"))
inds_30sec <- which(!is.na(getValues(mask_30sec)))

# Since there are 2E10^8 data points, we have to predict the models in digestable chunks to avoid memory issues. 
# Create chunks
inds_30sec_breaks <- chunk(inds_30sec, n.chunks = 500)

# Load fine res data
files_30sec <- list.files(data_path, pattern = "30sec", full.names = TRUE)
files_30sec <- files_30sec[which(grepl(paste(c("lu", "2000pop", "unsubregions_30sec") , collapse = "|"), files_30sec))]

# Write very large files into the smaller chunks so we can read them in from disk step by step
for(i in 1:length(files_30sec)){
  r <- readRDS(files_30sec[[i]])
  nm <- names(r)
  for(j in 1:length(inds_30sec_breaks)){
    saveRDS(r[inds_30sec_breaks[[j]]], file = file.path(temp_path, paste0(nm, "_chunk", j, "_nona.rds")))
    print(paste(i, j))
  }
  rm(r)
  removeTmpFiles(h = 0)
}

landuses <- c("Cropland", "Pasture", "Primary", "Secondary", "Urban")
intensities <- c("minimal", "light", "intense")

# Load models
result <- readRDS(file.path(data_path, "lumodels_2.rds"))

# 3. b Predict, processing chunks on parallel cores

# create cluster
cl <- makeCluster(10)
registerDoParallel(cl)

logfile <- file.path(data_path, paste0("gam_log.txt"))
writeLines(c(""), logfile)

# system("ps")
# system("pkill -f R")
landuses <- c("Cropland", "Pasture", "Primary", "Secondary", "Urban")
intensities <- c("minimal", "light", "intense")

lu_int_names <- as.vector(t(outer(landuses, intensities, paste, sep="_")))[-c(4,14)]

foreach(j = 1:length(inds_30sec_breaks), .packages = c("mgcv", "foreach")) %dopar% {
  
  cat(paste("chunk", j,"\n"), file = logfile, append = T)
  
  # Load covariates of current chunk into data frame
  gam_data <- lapply(list.files(temp_path, pattern = paste0("chunk", j, "_nona.rds"), full.names = TRUE), FUN = function(x) {readRDS(x)})
  gam_data <- as.data.frame(do.call("cbind", gam_data))
  colnames(gam_data) <- c("Cropland", "Pasture", "pop_2000pop_30min", "Primary", "Secondary", "unsubregions_30min", "Urban")
  chunk_inds <- inds_30sec_breaks[[j]]
  gam_data$unsubregions_30min <- as.factor(gam_data$unsubregions_30min)
  
  # Predict the models (as in validation)
  results <- foreach(i = 1:5, .packages = c("mgcv")) %do% {
    
    cat(paste(landuses[i], "\n"), file = logfile, append = T)
    
    # Load estimated model objects
    m_minimal_or_not <- result[[i]][[1]]
    m_light_or_intense <- result[[i]][[2]]
    
    # get predictions from first model (probability of being low)
    minimal <- predict(m_minimal_or_not, type = "response", newdata = gam_data)
    not_minimal <- 1 - minimal
    
    # and second model (probability of being medium *if not low*)
    if(length(m_light_or_intense) != 1){
      light_if_not_minimal <- predict(m_light_or_intense, type = "response", newdata = gam_data)
      intense_if_not_minimal <- 1 - light_if_not_minimal
      light <- light_if_not_minimal * not_minimal
      intense <- intense_if_not_minimal * not_minimal
      predicted <- cbind(minimal, light, intense)
    } else {
      predicted <- cbind(minimal, not_minimal)
    }
    
    # Apply modelled intensity probabilities to land use type fractions
    predicted_out <- sweep(predicted, 1, as.matrix(gam_data[landuses[i]]), FUN = "*")
    colnames(predicted_out) <- paste0(landuses[i], "_", colnames(predicted_out))
    predicted_out
  }
  predicted_out <- as.data.frame(do.call('cbind', results))
  colnames(predicted_out) <- lu_int_names
  
  # Save each predicted chunk indiividually to avoid memory issues, we can load them back in and merge into map later
  saveRDS(predicted_out, file = file.path(temp_path, paste0("predicted_chunk", j, "_nona.rds")))
}
stopCluster(cl)

# Recombine predicted chunks to one file for eaach land use type and intensity
landuses <- c("Cropland", "Pasture", "Primary", "Secondary", "Urban")
intensities <- c("minimal", "light", "intense")

lu_int_names <- as.vector(t(outer(landuses, intensities, paste, sep="_")))[-c(4,14)]

cl <- makeCluster(7)
registerDoParallel(cl)

logfile <- file.path(data_path, paste0("log.txt"))
writeLines(c(""), logfile)
j <- 3
foreach(j = 1:length(lu_int_names)) %do% {
  out <- numeric()
  lu_int <- lu_int_names[j]
  
  for(i in 1:500){
    cat(paste(paste(j, i), "\n"), file = logfile, append = T)
    predicted_chunk <- readRDS(list.files(temp_path, pattern = paste0("predicted_chunk",i, ".rds"), full.names = TRUE))
    out <- c(out, predicted_chunk[,which(colnames(predicted_chunk) == lu_int)])
  }
  saveRDS(out, file.path(out_path, paste0(lu_int, "_30sec", ".rds")))
}

stopCluster(cl)