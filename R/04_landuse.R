
# 1) Predict land use model for australia only (later for all of Ha's subregions)

devtools::install_github("kapitzas/flutes")

# 2) Disintegrate globaal data into Ha's subregions and run all model steps in chunks
#   - download global data sets (from last lu study)
#   - subset to eacah region (sequentially, overwriting files to avoid memory issues)
#   - make neighbourhood rasters
#   - do corr analysis and kick out predictors (take note of which aare kept in each chunk)
#   - run suit model
#   - make allocaations
#   - store future mappings to be combined later


