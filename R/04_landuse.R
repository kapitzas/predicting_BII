
# 1) link CGE sectors to land use type intensity classes
#   - this probably differs between Ha's subregions, think about it
#   - produce demand trajectories for each land use and intensity class (mean)

# 2) Disintegrate globaal data into Ha's subregions and run all model steps in chunks
#   - download global data sets (from last lu study)
#   - subset to eacah region (sequentially, overwriting files to avoid memory issues)
#   - make neighbourhood rasters
#   - do corr analysis and kick out predictors (take note of which aare kept in each chunk)
#   - run suit model
#   - make allocaations
#   - store future mappings to be combined later


