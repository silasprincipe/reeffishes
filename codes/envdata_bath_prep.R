#### Modelling of reef fishes ####
## Silas C. Principe - silascprincipe@gmail.com - 2021/2022

#### Bathymetry Layer preparation #####

### Important Note: to run this code you first need to download the
# GEBCO (General Bathymetric Chart of the Oceans) from www.gebco.net
# and place it in the data/env/bath_layers/ folder.
# Also note that you don't need to run it, as we already provide
# the croped layer.

# Load libraries ----
library(raster)

# Load data and define study parameters ----
# Define study extent
ext <- extent(-99, -29, -42.5, 42.5)

# Load GEBCO file
bath <- raster("data/env/bath_layers/gebco_america_africa.tif")

# Load shapefile
sel.area <- shapefile("data/gis/crop_shape.shp")

# Prepare the bath layer according to parameters ----
# Crop for the extension
bath.crop <- crop(bath, ext)

plot(bath.crop)

# Restrict according to depth
bath.2_300 <-
  calc(
    bath.crop,
    fun = function(x) {
      x[x > 2 | x < -300] <- NA
      return(x)
    }
  )

# Adjust resolution
bath.2_300 <- aggregate(bath.2_300, fact = 20, fun = mean)

# Mask in the selected area
m2_300 <- mask(bath.2_300, sel.area)

# Plot to verify
plot(m2_300)

# Convert areas above 0 to 0 (these areas are included just to avoid any
# difference between the Bio-ORACLE and the bath layer)
m2_300 <- calc(
  m2_300,
  fun = function(x){x[x >= 0] <- 0;x}
)

# See details to ensure everything is fine
m2_300

# Write final raster
writeRaster(m2_300, filename = "data/env/bath_layers/bath_2_300.tif",
            overwrite = T)

###END