#### Modelling of reef fishes ####
## Silas C. Principe - silascprincipe@gmail.com - 2022 ##

######## Preparing future layers ###########

# Note: see readme and files on the folder "cmip6" to understand
# how downscaled change factor files were produced

# Load packages ----
library(raster)
# Avoid GDAL to save xml aux file
rgdal::setCPLConfigOption("GDAL_PAM_ENABLED", "FALSE")

# Load original (present) Bio-ORACLE layers ----
env <- stack("data/env/crop_layers/tempmax.tif",
             "data/env/crop_layers/tempmean.tif",
             "data/env/crop_layers/salinitymean.tif",
             "data/env/crop_layers/chlomean.tif",
             "data/env/crop_layers/silicatemax.tif",
             "data/env/crop_layers/ph.tif",
             "data/env/crop_layers/windspeed.tif")

# Establish codes
codes <- c("tos", "tos", "sos", "chlos", "sios", "phos", "sfcWind")

# Change working dir
setwd("data/env/fut_layers")

# For each file, add the relative change for each scenario and save
for (i in 1:nlayers(env)) {
  
  cat("Running", codes[i], "\n")
  
  sel.cd <- codes[i]
  
  for (j in paste0("ssp", c(126, 245, 370, 585))) {
    cf <- list.files(paste0("change_factor/", sel.cd),
                     full.names = T)
    
    if (i <= 2) {
      cf <- cf[grep(ifelse(
        "tempmax" %in% names(env)[i],
        "max", "mean"
      ), cf)]
    }
    
    cf <- cf[grep(j, cf)]
    
    cat("Merging", cf, "and", names(env[[i]]), "\n")
    
    cf <- raster(cf)
    
    cf <- extend(cf, env)
    
    new <- cf * env[[i]]
    
    cat("Mean:", cellStats(new, 'mean'),
        "Max:", cellStats(new, 'max'),
        "Min:", cellStats(new, 'min'), "\n")
    
    final <- env[[i]] + new
    
    par(mfrow = c(1,2))
    plot(env[[i]], main = names(env)[i]);plot(final, main = j)
    
    writeRaster(final, paste0(names(env[[i]]), "_", j, ".tif"),
                overwrite = T)
    
  }
  
}

setwd("../../..")

# END