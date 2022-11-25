library(mapview)
source("functions/auxiliary_lgcp.R")
env <- raster("data/env/crop_layers/tempmax.tif")


sp <- 'acch'
pa.pts <- read.csv(paste0("data/", sp, "/", sp, "_pa.csv"))
pa.pts <- SpatialPointsDataFrame(pa.pts[, c("decimalLongitude", "decimalLatitude")],
                                 data.frame("presence" = pa.pts[,3]),
                                 proj4string = crs(env))

#pa.pts <- remove.duplicates(pa.pts, zero = 5)
pa.pts <- dup.cells(pa.pts, rast(env))

# Just to ensure all are falling inside study area
pa.pts <- pa.pts[!is.na(extract(env, pa.pts)),]

pa.pts <- pa.pts[pa.pts$presence == 1, ]


po.pts <- read.csv(paste0("data/", sp, "/", sp, "_final.csv"))
po.pts <- SpatialPoints(po.pts[, c("decimalLongitude", "decimalLatitude")],
                                                                                                   
                      proj4string = crs(env))
po.pts <- remove.duplicates(po.pts, zero = 5)


preds <- list.files(paste0("results/",sp,"/predictions/"), full.names = T, pattern = "q0.5")
preds <- preds[grep("_pa", preds)]
preds <- stack(preds)

mapview(po.pts, legend = F)+mapview(pa.pts, color = "red", legend = F)+
  mapview(preds, na.color = "transparent", legend = F, alpha.regions = 1)

