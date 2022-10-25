#### Modelling of reef fishes ####
## Silas C. Principe - silascprincipe@gmail.com - 2021/2022

#### Preparing environmental layers #####

### Note: Wind Layer have to be prepared previously using wind_layer_prep.R
### Here in this data package, the edited wind layer is already provided.


library(sdmpredictors)
library(dplyr)
library(raster)
library(sp)
library(usdm)
# Avoid GDAL to save xml aux file
rgdal::setCPLConfigOption("GDAL_PAM_ENABLED", "FALSE")

#Define study extent
ext <- extent(-99, -29, -42.5, 42.5)

# Load bathymetry layer
bath <- raster("data/env/bath_layers/bath_2_300.tif")

#Load layer codes

#codes <- read.table("data/env/layercodes.txt")
codes <- c("BO22_tempmean_ss",
           "BO22_tempmax_ss",
           "BO22_salinitymean_ss",
           "BO22_salinitymin_ss",
           "BO_ph",
           "BO22_chlomean_ss",
           "BO22_chlomax_ss",
           "BO22_silicatemax_ss",
           "BO22_silicatemean_ss",
           "BO22_dissoxmean_ss",
           "BO22_dissoxmin_ss",
           "BO22_curvelmax_ss",
           "BO22_damean",
           "BO22_damax")

# Load environmental layers

#layer.codes <- as.vector(codes$V1)
options(timeout=300)

env <- load_layers(
  layercodes = codes,
  equalarea = FALSE,
  rasterstack = FALSE,
  datadir = "data/env/env_layers"
)

# Crop environmental layers to study extent

for (i in 1:length(env)) {
  env[[i]] <- crop(env[[i]], ext)
}

env <- stack(env)


# Mask layers to bathymetry and to selected area

env <- mask(env, bath)

plot(env$BO22_tempmean_ss)
plot(bath)

names(env) <- gsub("_ss", "", names(env))
names(env) <- gsub("BO.\\d_|BO_", "", names(env))



#### Prepare Wind Speed Layer
# Open layer
wind <- brick("data/env/env_layers/CERSAT-GLO-REP_WIND_L4-OBS_FULL_TIME_SERIE_1665001093187.nc")

# Get mean of monthly data
wind <- mean(wind)

# Crop
wind <- crop(wind, ext)

# Increase resolution by disagregating
wind <- disaggregate(wind, 3)

# There are some cells that are valid in the Bio-ORACLE data, but not in the
# wind layer. So, we fill those using IDW.
not.valid <- rasterToPoints(env[[1]])
not.valid <- data.frame(cbind(not.valid[,1:2], extract(wind, not.valid[,1:2])))

library(gstat)

val <- not.valid[!is.na(not.valid[,3]),]
not.valid <- not.valid[is.na(not.valid[,3]),]

coordinates(val) <- ~ x+y
coordinates(not.valid) <- ~x+y

mod <- gstat(formula = as.formula(paste(names(val),"~ 1")),
             data = val, nmax = 12)

pred <- predict(mod, not.valid)


wind[cellFromXY(wind, coordinates(pred))] <- pred$var1.pred

# Mask
wind <- mask(wind, bath)
wind <- mask(wind, env[[1]])


#### Prepare roughness layer
bath.original <- raster("data/env/bath_layers/gebco_america_africa.tif")
bath.original <- crop(bath.original, ext + c(0,10,0,0))

bath.original <- aggregate(bath.original, 20)

rough <- terrain(bath.original, v = "roughness")

rough <- crop(rough, ext)
rough <- mask(rough, bath)


#### Stack new layers
names(rough) <- "roughness"
names(wind) <- "windspeed"

env <- stack(env, rough, wind)



#### Colinearity verification
vifstep.env <- vifstep(env, th = 10)

vifstep.env


### NEW COLINEARITY VERIFICATION - INCLUDING LAYERS THAT WE DECIDED ARE IMPORTANT

#Exclude based on vifstep
env.2 <- exclude(env, vifstep.env)

names(env.2)

###Now exclude variations of the same variable based on the ones
### that have stronger biological connection

env.3 <- dropLayer(env.2,
                   c(
                     "BO21_salinitymean_ss"
                     ))

names(env.3)

#We can make another vifstep verification
vifstep.env3 <- vifstep(env.3, th = 10)
vifstep.env3

#Write final list of layers
write.table(names(env.3), "data/env/env_layers.txt", col.names = F)

#Save Vifstep outputs
capture.output(vifstep.env, file = "data/env/vifstep_result.txt")
capture.output(vifstep.env@corMatrix, file = "data/env/vifstep_matrix.txt")
capture.output(vifstep.env3, file = "data/env/vif_result_afterexcluding.txt")


##### Separate rasters
names(env)

if (dir.exists("data/env/crop_layers") == F) {
  dir.create("data/env/crop_layers", recursive = T)
}

writeRaster(
  env,
  filename = paste("data/env/crop_layers/", names(env), sep = ""),
  format = "GTiff",
  bylayer = TRUE,
  overwrite = T
)

write.table(names(env), "data/env/layers_names_full.txt", col.names = F)



### Produce a simplified outer bound polygon for the barrier model
ob <- raster("data/env/env_layers/BO2_tempmean_ss_lonlat.tif")

ob <- crop(ob, env)

# Load shapefile
sel.area <- shapefile("data/env/crop_shape.shp")

lines(sel.area)

ob <- mask(ob, sel.area)

innerpol <- env$BO21_tempmean_ss
innerpol[!is.na(innerpol)] <- 1

rmpol <- SpatialPolygons(list(
  Polygons(list(
    Polygon(data.frame(x = c(-58.71294, -55.02481, -34.61716,
                             -11.01313, -17.40589, -58.71294),
                       y = c(43.126527, 23.456503,  6.245232,
                             5.261731, 44.847654, 43.126527)))
  ), ID = 1)
))

innerpol <- mask(innerpol, rmpol, inverse = T)


innerpol <- rasterToPolygons(innerpol, dissolve = T)

simp.inner <- buffer(innerpol, 1)

ob <- mask(ob, simp.inner)
ob[!is.na(ob)] <- 1
ob <- rasterToPolygons(ob, dissolve = T)

shapefile(innerpol, "gis/starea.shp")
shapefile(simp.inner, "gis/starea_buffered.shp")
shapefile(ob, "gis/starea_simp.shp")


#### END

library(rnaturalearth)
library(rnaturalearthdata)

# Load coast shapefile
coast <- ne_coastline(scale = "large", returnclass = "sp")
coast <- terra::vect(coast)

# Load a Bio-ORACLE base file
base <- load_layers("BO22_tempmax_ss", datadir = "data/env/env_layers")
base <- terra::rast(base)

# base.coast <- base
# base.coast[] <- 1
# 
# coast <- crop(coast, raster::extent(-100, -20, -45, 45))
# 
# system.time(coast.rast <- rasterize(coast, base.coast, mask = T))

#coast.rast <- rast(coast.rast)

# # Get a raster with coastlines
base.coast <- base
base.coast <- terra::rasterize(coast, base.coast, touches = T)
base.coast[base.coast == 1] <- 5
# # Takes approximately 18 minutes

coast.rast <- crop(base, ext(-100, -20, -45, 45))
base.coast <- crop(base.coast, ext(-100, -20, -45, 45))

coast.rast[!is.na(coast.rast)] <- 1

plot(coast.rast)

base.coast[is.na(base.coast)] <- 0

coast.rast <- coast.rast + base.coast

plot(coast.rast)

coast.rast[is.na(coast.rast)] <- 7777777


base.vect <- terra::rast("data/env/crop_layers/tempmax.tif")

base.vect[!is.na(base.vect)] <- 1

base.vect <- terra::as.polygons(base.vect)

base.vect <- terra::buffer(base.vect, width = 700000)

plot(base.vect)

crs(base.vect) <- crs(coast.rast)

plot(coast.rast);lines(base.vect)

# coast.rast <- mask(coast.rast, coast.rast, inverse = T, updatevalue = 9999999)
# 
# coast.rast[is.na(coast.rast)] <- 1
# 
# plot(coast.rast)

coast.rast <- terra::mask(coast.rast, base.vect, inverse = F, updatevalue = 7777777)

plot(coast.rast)

# coast.rast[coast.rast == 9999999] <- NA
# 
# plot(coast.rast)

coast.rast[coast.rast == 1] <- NA

system.time(dist.coast <- terra::distance(coast.rast, exclude = 7777777)) # approx 3h

model <- terra::rast("data/env/crop_layers/tempmax.tif")

dist.coast.final <- crop(dist.coast, model)

dist.coast.final <- mask(dist.coast.final, model)

plot(dist.coast.final)

writeRaster(dist.coast.final, "data/env/crop_layers/distcoast.tif")