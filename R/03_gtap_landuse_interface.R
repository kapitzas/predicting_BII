#----------------------------#
#### 0. Packages, folders ####
#----------------------------#

rm(list = ls())
require(rgdal)
require(raster)

data_path <- file.path(".", "data")
temp_path <- file.path(data_path, "temp")
out_path <- file.path(".", "output")
raw_path <-  file.path("/Volumes", "external", "OneDrive - The University of Melbourne", "PhD - Large Files", "raw data")
processed_path <- file.path(raw_path, "Global", "processed rasters")

#-----------------------------------------------------#
#### 1. Match UN countries to our GTAP aggregation ####
#-----------------------------------------------------#

# Load data
mask_5min <- raster(file.path(processed_path, "mask_5min.tif"))

gtap_regions <- read.csv(file.path(data_path, "gtap_aggregation.csv"))
world_borders <- readOGR(file.path(raw_path, "Global", "TM_WORLD_BORDERS_SIMPL-0", "TM_WORLD_BORDERS_SIMPL-0.3.shp"))
regions_table <- world_borders@data
regions_table$ISO3 <- tolower(regions_table$ISO3)


# match countries to our GTAP aggregation, depending on available info

matching_list<- list(
  "by_ISO2" = list(
    "e25" = c("AT", "EE", "IT", "PT", "BE", "FI", "LV", "FR", "LT", "SK", "DE", "LU", 'SI', 'CY', 'EL', 'MT', 'ES', 'CZ', 'HU', 'NL', 'SE', 'DK', 'IE', 'PL'),
    "gbr" = "UK"
  ),
  
  "by_UN" = list(
    "xsu" = c(031, 051, 112, 233, 268, 398, 417, 428, 440, 498, 762, 795, 804, 860)
  ),
  
  "by_subregion" = list(
    "xoc" = c(54,57,61),
    "xea" = 30,
    "xse" = 35,
    "xsa" = 34,
    "xsm" = c(5, 13, 29),
    'xer' = c(39, 151, 154, 830, 155),
    'xws' = 145,
    'xnf' = 15,
    'ssa' = c(11, 14, 17, 18)
  )
)

regions_table['GTAP'] <- NA


# Allocate countries we retained in GTAP aggregation
existing <- which(regions_table$ISO3%in%gtap_regions$Country.code)
regions_table$GTAP[existing] <- regions_table$ISO3[existing]

# Match EU25 countries (ISO2 info)
matches <- matching_list[[1]]
for(i in 1:length(matches)){
  regions_table$GTAP[which(regions_table$ISO2%in%matches[[i]] & is.na(regions_table$GTAP))] <- names(matches[i])
}

# Match former soviet union countries (UN code info)
matches <- matching_list[[2]]

for(i in 1:length(matches)){
  regions_table$GTAP[which(regions_table$UN%in%matches[[i]] & is.na(regions_table$GTAP))] <- names(matches[i])
}

# Match other countries (UN SUBREGION code info)
matches <- matching_list[[3]]
for(i in 1:length(matches)){
  regions_table$GTAP[which(regions_table$SUBREGION%in%matches[[i]] & is.na(regions_table$GTAP))] <- names(matches[i])
}

# Match the unallocated rest
regions_table$GTAP[is.na(regions_table$GTAP)] <- "xtw" # This includes Taiwan.

# add our own GTAP code column (so we can make a numeric world raster with our regions)
regions_table$GTAP_code <- gtap_regions$code[match(regions_table$GTAP, gtap_regions$Country.code)]
world_borders@data <- regions_table

# Rasterize
writeOGR(world_borders, dsn = file.path(out_path), layer = "GTAP_aggregation_borders", driver = "ESRI Shapefile", overwrite_layer = TRUE)
output <- file.path(out_path, "GTAP_aggregation_borders_5min.tif")
input <- file.path(out_path, "GTAP_aggregation_borders.shp")
gdaltools::rasterize_shp(input, output, res = 0.083, c(-180, 180, -90, 90), attribute = "GTAP_code")

# Mask
gtap_aggregation <- raster::raster(output)
gtap_aggregation <- mask(gtap_aggregation, mask_5min)

# Save
writeRaster(readAll(gtap_aggregation), file.path(processed_path, "gtap_aggregation_5min.tif"))
