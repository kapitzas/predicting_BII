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

raw_path <-  file.path("/Volumes", "external", "OneDrive - The University of Melbourne", "PhD - Large Files", "raw data")
processed_path <- file.path(raw_path, "Global", "processed rasters")
temp_path <- file.path("/Volumes", "external", "c3 processing", "temp")

# data_path <- file.path(".", "data")

# Create temp folder
if(!dir.exists(temp_path)){
  dir.create(temp_path)
}

#---------------#
#### 1. Mask ####
#---------------#

# Global mask at ca 1km
input <- file.path(raw_path, "Global", "TM_WORLD_BORDERS_SIMPL-0", "TM_WORLD_BORDERS_SIMPL-0.3.shp")
output <- file.path(temp_path, "mask_30sec.tif")
rasterize_shp(input, output, res =  0.0083, c(-180, 180, -90, 90), no_data = NA)
# subr <- raster(output)
# saveRDS(readAll(subr), file.path(data_path, "mask_30sec.rds"))
# unlink(output)

# GLobal mask at ca 0.5 degrees
output <- file.path(temp_path, "mask_30min.tif")
rasterize_shp(input, output, res =  0.5, c(-180, 180, -90, 90), no_data = NA)
# subr <- raster(output)
# saveRDS(readAll(subr), file.path(data_path, "mask_30min.rds"))
# unlink(output)

output <- file.path(temp_path, "unsubregions_30sec.tif")
gdaltools::rasterize_shp(input, output, res = 0.0083, c(-180, 180, -90, 90), attribute = "SUBREGION")
# subr <- raster(output)
# saveRDS(readAll(subr), file.path(data_path, "unsubregions_30sec.rds"))
# unlink(output)

output <- file.path(temp_path, "unsubregions_30min.tif")
gdaltools::rasterize_shp(input, output, res = 0.5, c(-180, 180, -90, 90), attribute = "SUBREGION")
# subr <- raster(output)
# saveRDS(readAll(subr), file.path(data_path, "unsubregions_30min.rds"))
# unlink(output)

#-----------------#
#### 2. Layers ####
#-----------------#

# 2. a Process pop dens data

popdens <- list.files(file.path(raw_path, "Global", "popdens"), pattern = ".tif$", recursive = TRUE, full.name = TRUE)

# Get SSP 2, SSP 7 and base year
popdens <- c(
  popdens[grepl("total_2000", popdens)],
  popdens[grepl("total_2070", popdens)]
)

# reproject to match mask
pop_names <- c("base", "ssp2", "ssp5")
for (i in 1:length(popdens)){
  input <- popdens[i]
  output <- file.path(temp_path, paste0(pop_names[i], "_pop_30sec.tif"))
  reproj_ras(input, output, crs = crs(mask_30sec), res = res(mask_30sec), method = "near", ext = extent(mask_30sec))
  
  
  output <- file.path(temp_path, paste0(pop_names[i], "_pop_30min.tif"))
  reproj_ras(input, output, crs = crs(mask_30min), res = res(mask_30min), method = "near", ext = extent(mask_30min))
  removeTmpFiles(h=0)
}

# 2.b Making distance rasters for land use model (need these to be global)

mask_30sec <- raster(file.path(temp_path, "mask_30sec.tif"))

# Roads
roads <- st_read(paste0(file.path(raw_path, "Global", "groads-v1-global-gdb", "gROADS_v1.gdb")))

rid <- matrix(c(0:7, "hwy", "pri", "sec", "tert", "loc", "trail", "priv", "unspec"), ncol = 2, nrow = 8)
rid <- data.frame(rid)
road_classes <- sort(unique(roads$FCLASS))

for (i in road_classes[1:7]){
  
  print(paste0("Writing subset ", i))
  roads_sub <- roads[which(roads$FCLASS == i),]
  outfile <- file.path(temp_path, paste0("diro.shp"))
  unlink(list.files(temp_path, pattern = "diro.", full.names = TRUE))
  st_write(roads_sub, outfile)
  
  print(paste0("Rasterizing subset ", i))
  infile <- file.path(temp_path, paste0("diro.shp"))
  outfile <- file.path(temp_path, paste0("roads_temp.tif"))
  rasterize_shp(infile, outfile, res = res(mask_30sec)[1], ext = extent(mask_30sec), no_data = NA)
  
  print(paste0("Proximity of subset ", i))
  infile <- file.path(temp_path, paste0("roads_temp.tif"))
  outfile <- file.path(temp_path, paste0("diro_", rid[i+1,2], "_30sec.tif"))
  proximity_ras(infile, outfile)
  
  unlink(infile)
  unlink(outfile)
  rm(out)
  raster::removeTmpFiles(h=0)
}


# Built-up areas
infile <- file.path(raw_path, "Global", "Global Built up areas", "bltupa.shp")
unlink(file.path(temp_path, "builtup_raster.tif"))
outfile <- file.path(temp_path, "builtup_raster.tif")
rasterize_shp(infile, outfile, res = res(mask_30sec), ext = extent(mask_30sec))

infile <- file.path(temp_path, "builtup_raster.tif")
outfile <- file.path(temp_path, "dibu_30sec.tif")
proximity_ras(infile, outfile)
unlink(infile)

# Lakes
infile <- file.path(raw_path, "Global", "GSHHS data", "GSHHS_lakes_L2-L4.shp")
unlink(file.path(temp_path, "dila_raster.tif"))
outfile <- file.path(temp_path, "dila_raster.tif")
rasterize_shp(infile, outfile, res = res(mask_30sec), ext = extent(mask_30sec))

infile <- file.path(temp_path, "dila_raster.tif")
outfile <- file.path(temp_path, "dila_30sec.tif")
proximity_ras(infile, outfile)
unlink(infile)

# Rivers
infile <- file.path(raw_path, "Global", "GSHHS data", "WDBII_rivers_global_L2-L9.shp")
unlink(file.path(temp_path, "diri_raster.tif"))
outfile <- file.path(temp_path, "diri_raster.tif")
rasterize_shp(infile, outfile, res = res(mask_30sec), ext = extent(mask_30sec))

infile <- file.path(temp_path, "diri_raster.tif")
outfile <- file.path(temp_path, "diri_30sec.tif")
proximity_ras(infile, outfile)
unlink(infile)

# 2. c Protected areas
infile <- file.path(raw_path, "Global", "Global Protected Areas", "WDPA_Mar2018-shapefile-polygons.shp")
require(rgeos)
pa <- sf::st_read(infile)
pa <- pa[which(pa$IUCN_CAT%in%c("Ia", "Ib", "II")),]
pa2 <- as(pa, 'Spatial')
pa2 <- gSimplify(pa2, tol = 0.00083)
pa2 <- st_as_sf(pa2)
outfile <- file.path(temp_path, "PA_IaIbII.shp")
sf::st_write(pa2, outfile)
pa <- st_read(outfile)

infile <- file.path(temp_path, "PA_IaIbII.shp")
outfile <- file.path(temp_path, "pa_30sec.tif")

rasterize_shp(infile, outfile, res = res(mask_30sec), ext = extent(mask_30sec), no_data = NA)
out <- raster(outfile)

out[which(!is.na(mask_30sec[]) & is.na(out[]))] <- 0
writeRaster(out, file.path(temp_path, paste0("pa_30sec.tif")), format = "GTiff", overwrite = TRUE)
rm(out, pa2, pa)
unlink(infile)

# 2. d ELevation
mask_30sec <- raster(file.path(temp_path, "mask_30sec.tif"))
infile <-file.path(raw_path, "Global", "wc2.1_30s_elev.tif")

unlink(file.path(temp_path, "srtm_30sec.tif"))
outfile <- file.path(temp_path,  "srtm_30sec.tif")
reproj_ras(infile, outfile, crs = crs(mask_30sec), ext = extent(mask_30sec), res = res(mask_30sec), method = "bilinear")

# slope and roughness
srtm <- raster(outfile)
slope <- terrain(srtm, "slope")
roughness <- terrain(srtm, "roughness")
writeRaster(slope, filename = file.path(temp_path, paste0("slope_30sec.tif")), driver = "GTiff", overwrite = TRUE)
writeRaster(roughness, filename = file.path(temp_path, paste0("roughness_30sec.tif")), driver = "GTiff", overwrite = TRUE)

# 2.e Bioclim variables
mask_30sec <- raster(file.path(temp_path, "mask_30sec.tif"))
biofiles <- list.files(file.path(raw_path, "Global", "wc21_30s_bio"), full.names = TRUE)
bionames <- sort(paste0("bio", 1:19))

for (i in 1:length(biofiles)){
  infile <- biofiles[i]
  outfile <- file.path(temp_path, paste0(bionames[i], "_30sec.tif"))
  reproj_ras(infile, outfile, crs = crs(mask_30sec), ext = extent(mask_30sec), res = res(mask_30sec), method = "bilinear")
}

# 2.f soils

#Full description: https://www.isric.org/explore/soilgrids/faq-soilgrids
mask_30sec <- raster(file.path(temp_path, "mask_30sec.tif"))

# Organic Carbon Density
infile <- file.path(raw_path, "Global", "soil_data", "ISRIC", "BLDFIE_M_sl3_1km_ll.tif")
outfile <- file.path(temp_path,  "ocdens_30sec.tif")
reproj_ras(infile, outfile, crs = crs(mask_30sec), ext = extent(mask_30sec), res = res(mask_30sec), method = "bilinear")

# Available Soil Water Capacity
infile <- file.path(raw_path, "Global", "soil_data", "ISRIC", "WWP_M_sl3_1km_ll.tif")
outfile <- file.path(temp_path,  "wwp_30sec.tif")
reproj_ras(infile, outfile, crs = crs(mask_30sec), ext = extent(mask_30sec), res = res(mask_30sec), method = "bilinear")

# pH Index measured in Water Solution
infile <- file.path(raw_path, "Global", "soil_data", "ISRIC", "PHIHOX_M_sl3_1km_ll.tif")
outfile <- file.path(temp_path,  "phihox_30sec.tif")
reproj_ras(infile, outfile, crs = crs(mask_30sec), ext = extent(mask_30sec), res = res(mask_30sec), method = "bilinear")

# Bulk density fine earth
infile <- file.path(raw_path, "Global", "soil_data", "ISRIC", "BLDFIE_M_sl3_1km_ll.tif")
outfile <- file.path(temp_path,  "bldfie_30sec.tif")
reproj_ras(infile, outfile, crs = crs(mask_30sec), ext = extent(mask_30sec), res = res(mask_30sec), method = "bilinear")

# cmip data (processed in parallel in separate script)
cmip5_path <- file.path("/Volumes", "external", "c3 processing", "gcm projections", "output quartiles")

q2_list <- list.files(cmip5_path, full.names = TRUE, pattern = "q2")
names <- sort(c(paste0("rcp45_", "bio", 1:19),
                paste0("rcp85_", "bio", 1:19)))

for (i in 1:length(q2_list)){
  infile <- q2_list[i]
  outfile <- file.path(temp_path,  paste0(names[i], "_30sec.tif"))
  reproj_ras(infile, outfile, crs = crs(mask_30sec), ext = extent(mask_30sec), res = res(mask_30sec), method = "bilinear")
}


#-------------------------------#
#### 3. Prepare land use data####
#-------------------------------#

mask_30sec <- raster(file.path(temp_path, "mask_30sec.tif"))
mask_30min <- raster(file.path(temp_path, "mask_30min.tif"))

# 3. a GLS Data

# Assign CRS and write back to disk
gls <- raster(file.path(raw_path, "Global", "GLS_data", "land_systems.asc"))
crs(gls) <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
writeRaster(gls, file.path(raw_path, "Global", "GLS_data", "land_systems.tif"), driver = "GTiff", overwrite = TRUE)

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

writeRaster(readAll(gls_layers), filename=file.path(temp_path, "gls_30min.tif"), options="INTERLEAVE=BAND", overwrite=TRUE)

# 3. b Harmonized land use data

# load data
harmonized <- list.files(file.path(raw_path, "Global", "harmonised land use downscaled"), pattern = ".bil", recursive = TRUE, full.name = TRUE)

# Get rid of ice class (no existe en Aus)

# harmonized <- harmonized[-which(grepl("ICE", harmonized))]
lu_names <- c("cropland", "pasture", "primary", "secondary", "urban")

# reproject harmonised to match gls and for downscaling
for (i in 1:length(harmonized)){
  
  lu_path <- harmonized[i]
  
  out <-  file.path(temp_path, paste0(lu_names[i], "_lu_30sec.tif"))
  reproj_ras(lu_path, out, crs = crs(mask_30sec), res = res(mask_30sec), method = "near", ext = extent(mask_30sec))
  
  out <- file.path(temp_path, paste0(lu_names[i], "_lu_30min.tif"))
  reproj_ras(lu_path, out, crs = crs(mask_30min), res = res(mask_30min), method = "near", ext = extent(mask_30min))
  removeTmpFiles(h=0)
}

# 3. c PREDICTS land use at 0.5 degree

# Make conversion table
conversion_table <- data.frame("GLS_class" = paste0("gls", 0:29))
conversion_table$cropland <- c(rep("minimal", 3), rep("light", 3), rep("intense", 5), "minimal", "light", rep("intense", 2), "minimal", "light", "intense", rep(NA, 12))
conversion_table$pasture <- c("light", rep("intense", 2), "light", rep("intense", 2), "light", rep("intense", 4), rep("light", 3), "intense", rep("light", 3), NA, "light", "intense", rep(NA,3), "light", "intense", NA, "light", rep(NA, 2))
conversion_table$primary <- c(rep(NA, 18), c("minimal", "light", "intense", rep("minimal", 3)), rep(NA, 2), "minimal", rep(NA, 3))
conversion_table$secondary <- conversion_table$primary
conversion_table$urban <- c(rep(NA, 28), "minimal", "intense")

write.csv(conversion_table, file = file.path("data", "gls_conversion_table.csv"))

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

gls_layers <- stack("/Volumes/external/c3 processing/temp/gls_30min.tif")
names(gls_layers) <- gls_classes

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
lu_30min <- stack(list.files(temp_path, pattern = "lu_30min", full.names = TRUE))

final_lu30min <- list()
for(i in 1:length(landuses)){
  r1 <- nl2[[which(grepl(landuses[i], names(nl2)))]]
  r2 <- lu_30min[[which(grepl(landuses[i], names(lu_30min)))]]
  r <- r1 * r2
  names(r) <- names(nl2)[which(grepl(landuses[i], names(nl2)))]
  final_lu30min[[i]] <- r
}

final_lu30min <- stack(final_lu30min)
names <- names(final_lu30min)
final_lu30min <- final_lu30min/sum(final_lu30min) # Here we are rescaling so cells sum to 1. Sometimes the harmonized classes in a cell aren't represented by any of the gls classes. In those areas, we assume that the gls class is more precise and assume no cover of the according harmoinzed class. Unallocated areas are distributed proportionately into existing harmonized classes.

names(final_lu30min) <- names

writeRaster(final_lu30min, filename=file.path(temp_path, paste0(names, "_predicts_30min.tif")), overwrite=TRUE, bylayer =TRUE)


# synch NA

# 30 min
mask_30min <- raster(file.path(temp_path, "mask_30min.tif"))
files_30min <- list.files(temp_path, pattern = "30min", full.names = TRUE)
names_30min <- list.files(temp_path, pattern = "30min")

for(i in 1:length(files_30min)){
  r <- raster(files_30min[[i]])
  mask_30min[is.na(r[])] <- NA
  print(i)
}

writeRaster(mask_30min, file.path(processed_path, "mask_30min.tif"), overwrite = TRUE)

for(i in 1:length(files_30min)){
  r <- raster(files_30min[[i]])
  r[is.na(mask_30min[])] <- NA
  writeRaster(r, file.path(processed_path, names_30min[[i]]), format = "GTiff", overwrite = TRUE)
  print(i)
}

# 30 sec 

files_30sec <- list.files(temp_path, pattern = "30sec", full.names = TRUE)
names_30sec <- list.files(temp_path, pattern = "30sec")
mask_30sec <- readAll(raster(file.path(temp_path, "mask_30sec.tif")))

for(i in 1:length(files_30sec)){
  r <- raster(files_30sec[[i]])
  mask_30sec <- readAll(mask(mask_30sec, r))
  print(i)
  raster::removeTmpFiles(h=0)
}
writeRaster(mask_30sec, file.path(processed_path, "mask_30sec.tif"), format = "GTiff", overwrite = TRUE)


mask_30sec <- raster(file.path(processed_path, "mask_30sec.tif"))
files_30sec <- list.files(temp_path, pattern = "30sec", full.names = TRUE)
names_30sec <- list.files(temp_path, pattern = "30sec")

for(i in 1:length(files_30sec)){
  system(paste0("gdal_calc.py -A '", file.path(processed_path, "mask_30sec.tif"),"'", " -B '", files_30sec[i],  "' --outfile='", file.path(temp_path, "temp.tif"), "'"," --calc='-99999*(A!=1) + B*(A==1)' --NoDataValue=-99999"))
  print(i)
}


