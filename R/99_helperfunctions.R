
# https://stackoverflow.com/questions/48502802/r-raster-extent-conditional-on-cell-value
ftrim <- function(r, regcode){
  xy <- rasterToPoints(r, function(x){ x == regcode })  
  e <- extent(xy[,1:2])
  e <- alignExtent(e, r, snap='out')
  v <- crop(r, e)
  v[v[] != regcode] <- NA
  v[!is.na(v)] <- regcode
  v
}

# from Adriana
getJacAbSym <- function(s1, s2, data){
  
  # get the list of species that are present in site 1 (i.e., their abundance was greater than 0)
  s1species <- data %>%
    
    # filter out the SSBS that matches s1
    filter(SSBS == s1) %>%
    
    # filter out the species where the Measurement (abundance) is greater than 0
    filter(Measurement > 0) %>%
    
    # get the unique species from this dataset
    distinct(Taxon_name_entered) %>%
    
    # pull out the column into a vector
    pull
  
  # for site 2, get the total abundance of species that are also present in site 1
  
  s2abundance_s1species <- data %>%
    
    # filter out the SSBS that matches s2
    filter(SSBS == s2) %>%
    
    # filter out the species that are also present in site 1
    filter(Taxon_name_entered %in% s1species) %>%
    
    # pull out the Measurement into a vector
    pull(Measurement) %>%
    
    # calculate the sum
    sum()
  
  # calculate the total abundance of all species in site 2
  s2_sum <- data %>%
    
    # filter out the SSBS that matches s2
    filter(SSBS == s2) %>%
    
    # pull out the measurement column (the abundance)
    pull(Measurement) %>%
    
    # calculate the sum
    sum() 
  
  
  # Now calculate the compositional similarity
  # this is the number of individuals of species also found in s1, divided by the total abundance in s2 
  # so that equates to the proportion of individuals in s2 that are of species also found in s1
  
  sor <- s2abundance_s1species / s2_sum
  
  
  # if there are no taxa in common, then sor = 0
  # if abundances of all taxa are zero, then similarity becomes NaN.
  return(sor)
  
}

multi_pdf <- function(x, paths, height, width){
  for(i in 1:length(paths)){
  pdf(paths[i], height = height, width = width)
    print(x)
    dev.off()
  }
}
