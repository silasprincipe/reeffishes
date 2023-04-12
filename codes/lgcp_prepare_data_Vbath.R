#### Modelling of reef fishes ####
## Silas C. Principe - silascprincipe@gmail.com - 2022 ##

######## Integrated Species Distribution Model ###########

## Prepare data to use with models ##
# - Prepare study area polygon
# - Prepare INLA mesh and save
# - Get integration points and save
# - Prepare full covariate dataset and save

# Load packages and define settings ----
# Modeling
library(INLA)
library(inlabru)
# Spatial data
library(raster)
library(rgeos)
library(gstat)
# Plotting and utilities
library(ggplot2)
# Settings
set.seed(2932) # Replicability of sampling
spatt <- theme_classic() # Theme for better ploting

# Prepare study area polygon ----
# Load base raster and convert to polygon
base <- raster("data/env/crop_layers/tempmean.tif")

starea <- calc(base, function(x){x[!is.na(x)] <- 1; x})
starea <- rasterToPolygons(starea, dissolve = T)

# Simplify study area
starea <- gSimplify(gBuffer(starea, width=0.1), tol=0.1);plot(starea)

# Get an outer bound that cover the whole study area
starea.pts <- rasterToPoints(base)
outer.bound <- inla.nonconvex.hull(starea.pts,
                                   convex = -0.05)
lines(outer.bound) # Plot to see


# Prepare INLA mesh ----
mesh <- inla.mesh.2d(
  boundary = list(starea, outer.bound),
  max.edge = c(0.4, 8),
  cutoff = 0.1,
  crs = crs(starea),
  offset = c(0.1)
)

ggplot()+gg(mesh)+coord_equal()+spatt;mesh$n



# Get integration points ----
ips <- ipoints(starea, mesh)

ggplot()+gg(mesh)+gg(ips)+coord_fixed()+spatt



# Prepare full covariate dataset ----
env <- stack("data/env/crop_layers/tempmean.tif",
             "data/env/crop_layers/tempmax.tif",
             "data/env/crop_layers/salinitymean.tif",
             "data/env/crop_layers/chlomean.tif",
             "data/env/crop_layers/silicatemax.tif",
             "data/env/crop_layers/ph.tif",
             "data/env/crop_layers/windspeed.tif",
             # In this file we include both the distance to coast
             # AND the bathymetry.
             # Note that, if we use distance to coast in the modeling
             # which is a certain proxy for bathymetry, than the model do fine
             # However, distance to coast does not incorporate all information
             # reagrding depth, so using bathymetry would be the ideal.
             "data/env/crop_layers/distcoast.tif",
             "data/env/bath_layers/bath_2_300.tif")

env$chlomean <- log(env$chlomean)

names(env)[8] <- "distcoast"

env <- mask(env, env[[1]]) # Just to ensure all layers are equal

env <- scale(env)

getd <- function(rast, ip){
  rast <- extend(rast, extent(-105, -20, -50, 50))
  
  epts <- extract(rast, ip)
  epts <- data.frame(epts, coordinates(ip))
  
  tofill <- epts[is.na(epts[,1]),]
  
  # get adjacent
  # adj <- adjacent(rast, cellFromXY(rast, tofill[, c("x", "y")]), pairs = F)
  # 
  # tofill <- rbind(
  #   tofill,
  #   cbind(extract(rast, adj), xyFromCell(rast, adj))
  # )
  
  tofill <- tofill[is.na(tofill[,1]),]
  
  epts <- data.frame(rasterToPoints(rast))
  
  coordinates(epts) <- ~x+y
  coordinates(tofill) <- ~x+y
  
  crs(epts) <- crs(tofill) <- crs(rast)
  
  for (i in 1:nlayers(rast)) {
    
    mod <- gstat(formula = as.formula(paste(names(epts)[i],"~ 1")),
                 data = epts, nmax = 12)
    
    pred <- predict(mod, tofill)
    
    rast[[i]][cellFromXY(rast[[i]], coordinates(tofill))] <- pred$var1.pred
  }
  
  rast
}

env.e <- getd(env, ip = ips)



# Save data ----
lgcp.data <- list(
  # Mesh
  mesh = mesh,
  # Study area shape
  starea = starea
)

saveRDS(lgcp.data, file = "data/lgcp_data.rds")

if (!dir.exists("data/env/ready_layers/")) {
  dir.create("data/env/ready_layers/")
}

# Convert to terra::rast to avoid problems when saving/reading in new versions
env.e <- terra::rast(env.e)

terra::writeRaster(env.e,
                   filename = paste0("data/env/ready_layers/",
                                     names(env.e), ".tif"),
                   overwrite = T)
