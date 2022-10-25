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
starea <- raster("data/env/crop_layers/tempmean.tif")
starea <- calc(starea, function(x){x[!is.na(x)] <- 1; x})
starea <- rasterToPolygons(starea, dissolve = T)

# Simplify study area
starea <- gSimplify(gBuffer(starea, width=0.5), tol=0.5);plot(starea)



# Prepare INLA mesh ----
mesh <- inla.mesh.2d(
  boundary = starea,
  max.edge = c(2, 4),
  cutoff = 0.8,
  crs = crs(starea),
  offset = c(0.1, 6)
)

ggplot()+gg(mesh)+coord_equal()+spatt;mesh$n



# Get integration points ----
ips <- ipoints(starea, mesh)

ggplot()+gg(mesh)+gg(ips)+coord_fixed()+spatt



# Prepare full covariate dataset ----
env <- stack("data/env/crop_layers/tempmax.tif",
             "data/env/crop_layers/salinitymean.tif",
             "data/env/crop_layers/chlomean.tif",
             "data/env/crop_layers/silicatemax.tif",
             "data/env/crop_layers/ph.tif",
             "data/env/crop_layers/windspeed.tif",
             "data/env/crop_layers/distcoast.tif")

env$chlomean <- log(env$chlomean)

env <- scale(env)

getd <- function(rast, ip){
  rast <- extend(rast, (extent(min(mesh$loc[,1]),
                               max(mesh$loc[,1]),
                               min(mesh$loc[,2]),
                               max(mesh$loc[,2])))+
                   c(-10, 10, -10, 10))
  
  epts <- extract(rast, ip)
  epts <- data.frame(epts, coordinates(ip))
  
  tofill <- epts[is.na(epts[,1]),]
  
  # get adjacent
  adj <- adjacent(rast, cellFromXY(rast, tofill[, c("x", "y")]), pairs = F)
  
  tofill <- rbind(
    tofill,
    cbind(extract(rast, adj), xyFromCell(rast, adj))
  )
  
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
  # Environmental layers expanded data
  env.e = env.e,
  # Mesh
  mesh = mesh,
  # Study area shape
  starea = starea
)

saveRDS(lgcp.data, file = "data/lgcp_data.rds")
