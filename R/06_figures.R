source("./R/99_helperfunctions.R")

#detach("package:flutes", unload=TRUE)
#devtools::install_github("kapitzas/flutes", force = TRUE)

require(raster)
require(rasterVis)
require(viridis)

raw_path <-  file.path("/Volumes", "external", "OneDrive - The University of Melbourne", "PhD - Large Files", "raw data")
processed_path <- file.path(raw_path, "Global", "processed rasters")
temp_path <- file.path("/Volumes", "external", "c3 processing", "temp")
out_path <- file.path(".", "output")
data_path <- file.path(".", "data")

figure_thesis <- "/Users/simon/ownCloud/PhD/writing/thesis/chapters/figures/chapter4/"
figure_chapter <- "/Users/simon/ownCloud/PhD/chapter3/figures"


single_width <- 3.11 #(79mm, single column)
medium_width <- 4.33 #110mm, 1.5 x page)
full_width <- 6.61 #168 mm, two column width

# Figure 1
gtap <- raster(file.path(processed_path, "gtap_aggregation_5min.tif"))
classes <- read.csv("/Users/simon/ownCloud/PhD/chapter3/data/GTAP_regions.csv")
classes <- unique(classes[,c(4,3)])
classes <- classes[order(classes$GTAP_code),]
colnames(classes)[1] <- "ID"
tar <- levels(gtap)[[1]]
tar <- classes
levels(gtap) <- as.data.frame(tar)

                
l <- levelplot(gtap, margin = F, colorkey = list(space='right'),
          xlab="", ylab="",
          par.settings=list(
            strip.border=list(col='transparent'),
            strip.background=list(col='transparent'),
            axis.line=list(col='transparent')),
          scales=list(draw=FALSE),
          col.regions = magma(100),
          maxpixels = 1e7
          )

fig1_paths <- c(file.path(figure_thesis, "fig1.pdf"),
                file.path(figure_chapter, "fig1.pdf"))
multi_pdf(x = l, paths = fig1_paths, height = 3.5, width = full_width)

harvested <- read.csv("/Users/simon/ownCloud/PhD/chapter3/data/demand calculation/GTAP_harvestedbyregion.csv")
require(tidyverse)
regions <- harvested %>% select(GTAP_country) %>% unique()
regions[i]
head(harvested)
i <- 3
out <- matrix(nrow = 30, ncol = 9)
colnames(out) <- c("region", "c_b", "gro", "ocr", "osd", "pdr", "pfb", "v_f", "wht")
out <- as.data.frame(out)
for(i in 1:nrow(regions)){
  harv_coun <- 
  harvested %>% 
    filter(GTAP_country == regions[i,]) %>%  
    select(c(GTAP, Value))
  out[i,1] <- regions[i,]
  out[i, match(harv_coun$GTAP, colnames(out))] <- harv_coun$Value
}
out[,-1] <- round(out[,-1]/rowSums(out[,-1], na.rm = TRUE)  * 100, 2)
write.csv(out, "/Users/simon/ownCloud/PhD/writing/thesis/tables/ch4_fao_estimates.csv")
