
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
