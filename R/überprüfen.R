# analysis
require(tidyverse)
raw_path <-  file.path("/Volumes", "external", "OneDrive - The University of Melbourne", "PhD - Large Files", "raw data")
data_path <- file.path(".", "data")
temp_path <- file.path(data_path, "temp")
out_path <- file.path("output")
int_path <- file.path(data_path, "lu_intensity")
processed_path <- file.path(raw_path, "Global", "processed rasters")

# checking demands

dmd <- readRDS("/Users/simon/ownCloud/PhD/chapter3/output/final_demands.rds")
dmd_or <- readRDS("/Users/simon/ownCloud/PhD/chapter3/output/demand_RCP85.rds")


dmd_85 <- dmd[[2]][c(1, 8)]
ma <- matrix(NA, nrow = 5, ncol = 30)
i <- 1
for(i in 1:30){
  dmd2 <- dmd_85[[2]][[i]][2,]
  dmd1 <- dmd_85[[1]][[i]][1,]
  ma[,i] <- (dmd2-dmd1)/dmd1 * 100
}

country_code <- read.csv(file.path(data_path, "GTAP_regions.csv"))
country_code <- distinct(country_code[,c(3,4)])

dmd_or <- readRDS(file.path(out_path, "demand_RCP85.rds"))
dmd_or <- dmd_or[,match(country_code$GTAP_aggregation, colnames(dmd_or))]


# checking that land use classes match with demands

files <- list.files(processed_path, paste0("rcp85"), full.names= TRUE)
files <- stack(files[which(grepl("lu_5min", files))])
fut <- as.data.frame(files)

list.files(processed_path)
files <- list.files(processed_path, paste0("lu_5min"), full.names= TRUE)
files <- stack(files[-which(grepl("rcp", files))])
pres <- as.data.frame(files)
cmeans <- colMeans(files, na.rm = TRUE)

means_pres <- colMeans(pres, na.rm = TRUE)
means_fut <- colMeans(fut, na.rm = TRUE)
(means_fut - means_pres)/means_pres * 100

rowMeans(ma)
rowMeans(dmd_or)

# intensity maps

files <- list.files(processed_path, paste0("rcp85"), full.names= TRUE)
files <- files[which(grepl(paste0(c("light", "intense", "minimal"), collapse = "|"), files))]
fut <- stack(files)

files <- list.files(processed_path, paste0("pres"), full.names= TRUE)
files <- files[which(grepl(paste0(c("light", "intense", "minimal"), collapse = "|"), files))]
pres <- stack(files)
diff_list <- list()
for(i in 1:nlayers(pres)){
  diff_list[[i]] <- fut[[i]] - pres[[i]]
}

st <- stack(diff_list)
names(st) <- names(fut)
spplot(st)
