#----------------------------#
#### 0. Packages, folders ####
#----------------------------#

rm(list = ls())
require(devtools)
require(WorldClimTiles)
require(rgdal)
require(gdaltools)
require(raster)
require(stringr)
require(rgeos)
require(sf)
require(flutes)
library("doParallel")
require(mgcv)

raw_path <-  file.path(path.expand("~"), "OneDrive - The University of Melbourne", "PhD - Large Files", "PhD - Raw Data")
data_path <- file.path(".", "data")
temp_path <- file.path(data_path, "temp")

# Create temp folder
if(!dir.exists(temp_path)){
   dir.create(temp_path)
}

#---------------#
#### 1. Mask ####
#---------------#

# Global mask at ca 1km
input <- file.path(data_path, "TM_WORLD_BORDERS_SIMPL-0", "TM_WORLD_BORDERS_SIMPL-0.3.shp")
output <- file.path(data_path, "mask_30sec.tif")
rasterize_shp(input, output, res =  0.0083, c(-180, 180, -90, 90), no_data = NA)
subr <- raster(output)
saveRDS(readAll(subr), file.path(data_path, "mask_30sec.rds"))
unlink(output)

# GLobal mask at ca 0.5 degrees
output <- file.path(data_path, "mask_30min.tif")
rasterize_shp(un_subregions, output, res =  0.5, c(-180, 180, -90, 90), no_data = NA)
subr <- raster(output)
saveRDS(readAll(subr), file.path(data_path, "mask_30min.rds"))
unlink(output)

output <- file.path(data_path, "unsubregions_30sec.tif")
gdaltools::rasterize_shp(un_subregions, output, res = 0.0083, c(-180, 180, -90, 90), attribute = "SUBREGION")
subr <- raster(output)
saveRDS(readAll(subr), file.path(data_path, "unsubregions_30sec.rds"))
unlink(output)

output <- file.path(data_path, "unsubregions_30min.tif")
gdaltools::rasterize_shp(un_subregions, output, res = 0.5, c(-180, 180, -90, 90), attribute = "SUBREGION")
subr <- raster(output)
saveRDS(readAll(subr), file.path(data_path, "unsubregions_30min.rds"))
unlink(output)


#-------------------------------#
#### 2. Prepare land use data####
#-------------------------------#

mask_30sec <- readRDS(file.path(data_path, "mask_30sec.rds"))
mask_30min <- readRDS(file.path(data_path, "mask_30min.rds"))

# 2. a GLS Data

# Assign CRS and write back to disk
GLS <- raster(file.path(raw_path, "Global", "GLS_data", "land_systems.asc"))
crs(GLS) <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
writeRaster(GLS, file.path(raw_path, "Global", "GLS_data", "land_systems.tif"), driver = "GTiff", overwrite = TRUE)

input <- file.path(raw_path, "Global", "GLS_data", "land_systems.tif")
output <- file.path(temp_path, "gls_5min.tif")
reproj_ras(input, output, crs = crs(mask_30sec), res = res(mask_30sec)*10, method = "near", ext = extent(mask_30sec))

gls <- raster(output)
gls_layers <- layerize(gls)
gls_layers <- mask(gls_layers, gls)
gls_layers <- aggregate(gls_layers, fac = 6, fun = "mean")
gls_layers <- resample(gls_layers, mask_30min, method = "ngb")

gls_classes <- gsub("X", "gls", names(gls_layers))
names(gls_layers) <- gls_classes
saveRDS(readAll(gls_layers), file.path(data_path, "gls_30min.rds"))
unlink(output)

# 2. b Harmonized land use data

# load data
harmonized <- list.files(file.path(raw_path, "Global", "harmonised land use downscaled"), pattern = ".bil", recursive = TRUE, full.name = TRUE)

# Get rid of ice class (no existe en Aus)

# harmonized <- harmonized[-which(grepl("ICE", harmonized))]
lu_names <- c("Cropland", "Pasture", "Primary", "Secondary", "Urban")

# reproject harmonised to match gls and for downscaling
for (i in 1:length(harmonized)){
  
  lu_path <- harmonized[i]
  
  out <-  file.path(temp_path, paste0(lu_names[i], "_lu_30sec.tif"))
  reproj_ras(lu_path, out, crs = crs(mask_30sec), res = res(mask_30sec), method = "near", ext = extent(mask_30sec))
  subr <- raster(out)
  print("masking")
  subr <- mask(subr, mask_30sec)
  saveRDS(readAll(subr), file.path(data_path, paste0(lu_names[i], "_lu_30sec.rds")))
  unlink(out)
  
  out <- file.path(temp_path, paste0(lu_names[i], "_lu_30min.tif"))
  reproj_ras(lu_path, out, crs = crs(mask_30min), res = res(mask_30min), method = "near", ext = extent(mask_30min))
  subr <- raster(out)
  subr <- mask(subr, mask_30min)
  saveRDS(readAll(subr), file.path(data_path, paste0(lu_names[i], "_lu_30min.rds")))
  unlink(out)
  removeTmpFiles(h=0)
}

names(gls_layers) <- gsub("X", "gls", names(gls_layers))
saveRDS(readAll(lu_30min), file.path(".", "data", "harmonized_30min.rds"))
saveRDS(readAll(gls_layers), file.path(".", "data", "gls_30min.rds"))

# 2. c PREDICTS land use at 0.5 degree
# Make conversion table
conversion_table <- data.frame("GLS_class" = paste0("gls", 0:29))
conversion_table$Cropland <- c(rep("minimal", 3), rep("light", 3), rep("intense", 5), "minimal", "light", rep("intense", 2), "minimal", "light", "intense", rep(NA, 12))
conversion_table$Pasture <- c("light", rep("intense", 2), "light", rep("intense", 2), "light", rep("intense", 4), rep("light", 3), "intense", rep("light", 3), NA, "light", "intense", rep(NA,3), "light", "intense", NA, "light", rep(NA, 2))
conversion_table$Primary <- c(rep(NA, 18), c("minimal", "light", "intense", rep("minimal", 3)), rep(NA, 2), "minimal", rep(NA, 3))
conversion_table$Secondary <- conversion_table$Primary
conversion_table$Urban <- c(rep(NA, 28), "minimal", "intense")

write.csv(conversion_table, file = "./data/GLS_conversion_table.csv")

glc_har_restr <- list()
intensities <- c("minimal", "light", "intense")
landuses <- colnames(conversion_table)[-1]

# Get gls classes per land use type and intensity class
for (i in 1:5){
  har_restr <- list() 
  for (j in 1:3){
    har_restr[[j]] <- conversion_table$GLS_class[which(conversion_table[,i + 1] == intensities[j])]
  }
  names(har_restr) <- intensities
  glc_har_restr[[i]] <- har_restr
}

names(glc_har_restr) <- landuses

gls_layers <- readRDS("/Users/simon/OneDrive - The University of Melbourne/PhD/chapter3/data/gls_30min.rds")


# Sum gls classes per assigned land use and intensity level to find the fraction occupied by different intensities within each type

nl2 <- list()
i <- j <- 1

for (i in 1:5){
  nl1 <- list()
  for(j in 1:3){
    ly_ind <- glc_har_restr[[i]][[j]]
    name <- paste(landuses[i], intensities[j], sep = "_")
    if(!length(ly_ind)  == 0){
      r <- gls_layers[[which(names(gls_layers) %in% ly_ind)]]
      if(length(ly_ind) > 1){
        r <- sum(r)
      }
      names(r) <- name
      nl1[[j]] <- r
    }
  }
  nl1 <- stack(unlist(nl1))
  name2 <- names(nl1)
  names(nl1) <- name2
  nl2[[i]] <- nl1
}

nl2 <- stack(unlist(nl2))

type_int_classes <- names(nl2)

lu_30min <- list()
lu_files <- list.files(data_path, pattern = "lu_30min", full.names = TRUE)
for (i in 1:length(lu_files)){
  lu_30min[[i]] <- readRDS(lu_files[[i]])
}
lu_30min <- stack(lu_30min)

for(i in 1:length(landuses)){
  r1 <- nl2[[which(grepl(landuses[i], names(nl2)))]]
  r2 <- lu_30min[[which(grepl(landuses[i], names(lu_30min)))]]
  r <- r1 * r2
  names(r) <- names(nl2)[which(grepl(landuses[i], names(nl2)))]
  final_lu30min[[i]] <- r
}

final_lu30min <- stack(final_lu30min)
plot(sum(final_lu30min) > 0)
names <- names(final_lu30min)
spplot(sum(final_lu30min))
spplot(final_lu30min)
final_lu30min <- final_lu30min/sum(final_lu30min) # Here we are rescaling so cells sum to 1. Sometimes the harmonized classes in a cell aren't represented by any of the gls classes. In those areas, we assume that the gls class is more precise and assume no cover of the according harmoinzed class. Unallocated areas are distributed proportionately into existing harmonized classes.

names(final_lu30min) <- names

saveRDS(readAll(final_lu30min), file = file.path(temp_path, "lu_predicts_30min.rds"))

# 2. c Load pop dens data

popdens <- list.files(file.path(raw_path, "Global", "popdens"), pattern = ".tif$", recursive = TRUE, full.name = TRUE)

popdens <- c(
  popdens[grepl("total_2000", popdens)],
  popdens[grepl("total_2050", popdens)],
  popdens[grepl("total_2100", popdens)]
)

# reproject to match mask
pop_names <- c("pop_2000", "pop_2050", "pop_2100")
for (i in 2:length(popdens)){
  input <- popdens[i]
  output <- file.path(temp_path, paste0(pop_names[i], "pop_30sec.tif"))
  reproj_ras(input, output, crs = crs(mask_30sec), res = res(mask_30sec), method = "near", ext = extent(mask_30sec))
  subr <- raster(output)
  subr <- mask(subr, mask_30sec)
  saveRDS(readAll(subr), file.path(data_path, paste0(pop_names[i], "pop_30sec.rds")))
  unlink(output)
  
  output <- file.path(temp_path, paste0(pop_names[i], "pop_30min.tif"))
  reproj_ras(input, output, crs = crs(mask_30min), res = res(mask_30min), method = "near", ext = extent(mask_30min))
  subr <- raster(output)
  subr <- mask(subr, mask_30min)
  saveRDS(readAll(subr), file.path(data_path, paste0(pop_names[i], "pop_30min.rds")))
  unlink(output)
  removeTmpFiles(h=0)
}

# synch NA

files_30min <- list.files(data_path, pattern = "30min", full.names = TRUE)
for(i in 1:length(files_30min)){
  r <- readRDS(files_30min[[i]])
  mask_30min[is.na(r[])] <- NA
  print(i)
}
saveRDS(mask_30min, file.path(data_path, "mask_30min.rds"))

for(i in 1:length(files_30min)){
  r <- readRDS(files_30min[[i]])
  r[is.na(mask_30min[])] <- NA
  saveRDS(r, files_30min[[i]])
  print(i)
}


files_30sec <- list.files(data_path, pattern = "30sec", full.names = TRUE)
for(i in 1:length(files_30sec)){
  r <- readRDS(files_30sec[[i]])
  mask_30sec[is.na(r[])] <- NA
  print(i)
}
saveRDS(mask_30sec, file.path(data_path, "mask_30sec.rds"))

for(i in 1:length(files_30sec)){
  r <- readRDS(files_30sec[[i]])
  r[is.na(mask_30sec[])] <- NA
  saveRDS(r, files_30sec[[i]])
  print(i)
}

# Modelling
# make modelling dataframe
mask_30min <- readRDS(file.path(data_path, "mask_30min.rds"))
inds_30min <- which(!is.na(mask_30min[]))
files_30min <- list.files(data_path, pattern = "30min", full.names = TRUE)
files_30min <- grep(paste(c("lu_predicts", "unsubregions", "pop_2000"),collapse="|"), files_30min, value=TRUE)

df_list <- list()
for( i in 1:length(files_30min)){
  r <- readRDS(files_30min[[i]])
  df_list[[i]] <- r
}


gam_data <- as.data.frame(stack(df_list))[inds_30min,]
landuses <- c("Cropland", "Pasture", "Primary", "Secondary", "Urban")
intensities <- c("minimal", "light", "intense")

for(i in 1:length(landuses)){
  lu_inds <- grep(paste(paste(landuses[i], intensities, sep = "_"), collapse = "|"), colnames(gam_data))
  gam_data[landuses[i]] <- rowSums(gam_data[,lu_inds])
  gam_data[,lu_inds] <- integerify(gam_data[,lu_inds], resolution = 1000)
}

# Combine pacific micro states into one region
gam_data[is.na(gam_data)] <- 0

gam_data$unsubregions_30min[gam_data$unsubregions_30min%in%c(57, 54, 61,29,155)] <- 99

gam_data$unsubregions_30min <- as.factor(gam_data$unsubregions_30min)

stopCluster(cl)
cl <- makeCluster(5)
registerDoParallel(cl)

logfile <- file.path(data_path, paste0("gam_log.txt"))
writeLines(c(""), logfile)

mask_30min <- readRDS("./data/mask_30min.rds")

results <- foreach(i = 1:5, .packages = c("mgcv")) %dopar% {
  
  lus <- landuses[i]
  cat(paste(landuses[i],"\n"), file = logfile, append = T)
  f_minimal_or_not <- as.formula(paste("minimal_or_not ~", paste(c("unsubregions_30min", paste0("s(", c("pop_2000pop_30min", lus), ")"), paste0("s(", c(lus, "pop_2000pop_30min"), ",by = unsubregions_30min)")), collapse = "+")))
  f_light_or_intense <- as.formula(paste("light_or_intense ~", paste(c("unsubregions_30min", paste0("s(", c("pop_2000pop_30min", lus), ")"), paste0("s(", c(lus, "pop_2000pop_30min"), ",by = unsubregions_30min)")), collapse = "+")))
  
  # Get columns matching land use type and intensity classes
  lu_inds <- grep(paste(paste(landuses[i], intensities, sep = "_"), collapse = "|"), colnames(gam_data))
  
  # Get the data 
  if(length(lu_inds) == 2){
    minimal_or_not <- cbind(gam_data[,lu_inds][,1], gam_data[,lu_inds][, 2])
  }
  
  if(length(lu_inds) == 3){
    minimal_or_not <- cbind(gam_data[,lu_inds][,1], rowSums(gam_data[,lu_inds][,-1]))
  }
  
  cat(paste("minimal or not: ",landuses[i],"\n"), file = logfile, append = T)
  m_minimal_or_not <- tryCatch(mgcv::gam(f_minimal_or_not, family = binomial, data = gam_data), error = function(e) NA)
  if(is.na(m_minimal_or_not)) {cat(paste("returned NA - minimal or not: ",landuses[i],"\n"), file = logfile, append = T)}
  if(length(lu_inds)==3){
    light_or_intense <- as.matrix(gam_data[,lu_inds][,-1])
    cat(paste("light or intense: ",landuses[i],"\n"), file = logfile, append = T)
    m_light_or_intense <- tryCatch(mgcv::gam(f_light_or_intense, family = binomial, data = gam_data), error = function(e) NA)
    if(is.na(m_light_or_intense)) {cat(paste("returned NA - light or intense: ",landuses[i],"\n"), file = logfile, append = T)}
  } else {
    m_light_or_intense <- NA 
  }
  list(m_minimal_or_not, m_light_or_intense, landuses[i])
}
saveRDS(results, file.path(data_path, "lumodels_2.rds"))



# Predicting for validation
result <- readRDS(file.path(data_path, "lumodels_2.rds"))


valid_df <- list()
for(i in 1:length(landuses)){
  
  # Calculate error
  
  m_minimal_or_not <- result[[i]][[1]]
  m_light_or_intense <- result[[i]][[2]]
  
  # get predictions from first model (probability of being low)
  minimal <- predict(m_minimal_or_not, type = "response", data = gam_data)
  not_minimal <- 1 - minimal
  
  # and second model (probability of being medium *if not low*)
  if(length(m_light_or_intense) != 1){
    light_if_not_minimal <- predict(m_light_or_intense, type = "response", data = gam_data)
    intense_if_not_minimal <- 1 - light_if_not_minimal
    light <- light_if_not_minimal * not_minimal
    intense <- intense_if_not_minimal * not_minimal
    predicted <- cbind(minimal, light, intense)
  } else {
    predicted <- cbind(minimal, not_minimal)
  }
  valid_out <- sweep(predicted, 1, as.matrix(gam_data[landuses[i]]), FUN = "*")
  colnames(valid_out) <- paste0(landuses[i], "_", colnames(valid_out))
  valid_df[[i]] <- valid_out
}

summary(result[[1]][[1]])

valid_df <- do.call("cbind", valid_df)
obs_df <- as.data.frame(readRDS("/Users/simon/OneDrive - The University of Melbourne/PhD/chapter3/data/lu_predicts_30min.rds"))[inds_30min,]
obs_df - valid_df

error <- as.matrix(obs_df) - as.matrix(valid_df)
error <- sqrt(rowMeans(error^2))

test <- mask_30min

test[inds_30min] <- gam_data$pop_2000pop_30min
test[inds_30min] <- error
rasterVis::levelplot(test)


# Predicting high res lu intensity map
result <- readRDS(file.path(data_path, "lumodels_2.rds"))
str(gam_data)

files_30sec <- list.files(data_path, pattern = "30sec", full.names = TRUE)[c(1,2,3,4,7,8,9,10)]

data_30sec <- list()
for( i in 1:length(files_30sec)){
  r <- readRDS(files_30sec[[i]])
  data_30sec[[i]] <- r
}

gam_data <- na.omit(as.data.frame(stack(data_30sec)))

valid_df <- list()
for(i in 1:length(landuses)){
  
  # Calculate error
  
  m_minimal_or_not <- result[[i]][[1]]
  m_light_or_intense <- result[[i]][[2]]
  
  # get predictions from first model (probability of being low)
  minimal <- predict(m_minimal_or_not, type = "response", data = gam_data)
  not_minimal <- 1 - minimal
  
  # and second model (probability of being medium *if not low*)
  if(length(m_light_or_intense) != 1){
    light_if_not_minimal <- predict(m_light_or_intense, type = "response", data = gam_data)
    intense_if_not_minimal <- 1 - light_if_not_minimal
    light <- light_if_not_minimal * not_minimal
    intense <- intense_if_not_minimal * not_minimal
    predicted <- cbind(minimal, light, intense)
  } else {
    predicted <- cbind(minimal, not_minimal)
  }
  valid_out <- sweep(predicted, 1, as.matrix(gam_data[landuses[i]]), FUN = "*")
  colnames(valid_out) <- paste0(landuses[i], "_", colnames(valid_out))
  valid_df[[i]] <- valid_out
}

