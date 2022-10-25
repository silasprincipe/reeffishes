#### Modelling of reef fishes ####
## Silas C. Principe - silascprincipe@gmail.com - 2022 ##

######## Integrated Species Distribution Model ###########

## Load environmental data ##

library(raster)

env <- stack("data/env/crop_layers/tempmax.tif",
             "data/env/crop_layers/salinitymean.tif",
             "data/env/crop_layers/chlomean.tif",
             "data/env/crop_layers/silicatemax.tif",
             "data/env/crop_layers/ph.tif",
             "data/env/crop_layers/windspeed.tif",
             "data/env/crop_layers/distcoast.tif")

env$chlomean <- log(env$chlomean)

env <- scale(env)



# m <- cellStats(env, "mean")
# sde <- cellStats(env, "sd")


# Load future layers ----
# ssp3 <- stack("data/env/fut_layers/BO21_tempmax_ss.tif",
#               "data/env/fut_layers/BO21_salinitymean_ss.tif",
#               "data/env/fut_layers/BO21_chlomean_ss.tif",
#               "data/env/fut_layers/BO_ph.tif")
# 
# ssp3$BO21_chlomean_ss <- log(ssp3$BO21_chlomean_ss)
# 
# ssp3 <- (ssp3 - m)/sde
###