#### Modelling of coral reef herbivorous invertebrates ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2022 ##

# Calculate area change

# Load needed packages ----
library(raster)
library(terra)
source("functions/auxiliary_lgcp.R")

# Load species data ----
# Load a layer just as sample
env <- raster("data/env/crop_layers/tempmax.tif")

get.area <- function(sp, type){
  
  ##### Presence only dataset ----
  po.pts <- read.csv(paste0("data/", sp, "/", sp, "_final.csv"))
  po.pts <- SpatialPoints(po.pts[, c("decimalLongitude", "decimalLatitude")],
                          proj4string = crs(env))
  
  po.pts <- remove.duplicates(po.pts, zero = 5)
  
  ##### Presence-absence dataset ----
  pa.pts <- read.csv(paste0("data/", sp, "/", sp, "_pa.csv"))
  pa.pts <- SpatialPointsDataFrame(pa.pts[, c("decimalLongitude", "decimalLatitude")],
                                   data.frame("presence" = pa.pts[,3]),
                                   proj4string = crs(env))
  
  #pa.pts <- remove.duplicates(pa.pts, zero = 5)
  pa.pts <- dup.cells(pa.pts, rast(env))
  
  # Just to ensure all are falling inside study area
  pa.pts <- pa.pts[!is.na(extract(env, pa.pts)),]
  
  
  
  # Load results ----
  results <- list.files(paste0("results/", sp, "/predictions"),
                        pattern = type, full.names = T)
  results <- results[grep("q0.5", results)]
  
  # Load rasters generated before
  curr <- raster(results[1])
  ssp1 <- raster(results[2])
  ssp2 <- raster(results[3])
  ssp3 <- raster(results[4])
  ssp5 <- raster(results[5])
  
  # Convert to data.frame
  get.val <- function(x){
    temp <- as(x, "SpatialPixelsDataFrame")
    temp <- as.data.frame(temp)
    colnames(temp) <- c("val", "x", "y")
    return(temp)
  }
  
  curr.v <- get.val(curr)
  ssp1.v <- get.val(ssp1)
  ssp2.v <- get.val(ssp2)
  ssp3.v <- get.val(ssp3)
  ssp5.v <- get.val(ssp5)
  
  # Get thresholded contour
  if (type == "int" | type == "cont") {
    mtp <- min(extract(curr, po.pts))
    p10 <- rev(sort(extract(curr, po.pts)))[ceiling(length(po.pts) * 0.9)]
  } else {
    pres <- pa.pts[pa.pts$presence == 1,]
    mtp <- min(extract(curr, pres))
    p10 <- rev(sort(extract(curr, pres)))[ceiling(length(pres) * 0.9)]
  }
  
  get.cont <- function(rast, th){
    r <- rast
    r[r < th] <- NA
    r[r >= th] <- 1
    r <- rasterToPolygons(r, dissolve = T)
    #r <- st_as_sf(r)
    r
  }
  
  pols <- lapply(list(curr, ssp1, ssp2, ssp3, ssp5), get.cont, th = p10)
  
  df <- data.frame(
    species = sp,
    mode = ifelse(mod == "quant", "quant", NA),
    type = type,
    scenario = c("current", paste0("ssp", c(1,2,3,5))),
    area = NA
  )
  
  df$area <- sapply(1:5, function(x){raster::area(pols[[x]])/1e6})
  
  df$delta <- df$area-df$area[1]
  df$delta_perc <- (df$area*100)/df$area[1]
  
  return(df)
}

# Define species ----
sp.list <- c("acch", "scze", "spam", "mybo", "lujo")

# Get area ----

# Available types:
# [int = relative occurrence rate; cont = contrast (ROR without spatial);
# _pa = probability of occurrence (binomial likelihood)]
# NOTE: use _pa with spam, or the file will be wrong!
# If type is "int" or "cont" it's possible to use
# the quantile mode

int.area <- lapply(sp.list, get.area, type = "int", mod = "quant")
pa.area <- lapply(sp.list, get.area, type = "_pa", mod = "normal")
cont.area <- lapply(sp.list, get.area, type = "cont", mod = "quant")
int.area2 <- lapply(sp.list, get.area, type = "int", mod = "not")



# Overlap metrics ----
##### Between species ----
res <- list()

for (i in 1:length(sp.list)) {
  sp <- sp.list[i]
  
  results <- list.files(paste0("results/", sp, "/predictions"),
                        pattern = "_pa", full.names = T)
  results <- results[grep("q0.5", results)]
  
  # Load rasters generated before
  res[[i]] <- raster(results[1])
}

res <- stack(res)

names(res) <- sp.list

comb <- combn(sp.list, 2)

bspecies <- data.frame(species_a = NA, species_b = NA,
                       D = rep(NA, dim(comb)[2]), I = NA)

for (z in 1:nrow(bspecies)) {
  bspecies[z,"species_a"] <- comb[1,z];bspecies[z,"species_b"] <- comb[2,z]
  bspecies[z,"I"] <- dismo::nicheOverlap(res[[comb[1,z]]], res[[comb[2,z]]], stat='I', mask = F)
  bspecies[z,"D"] <- dismo::nicheOverlap(res[[comb[1,z]]], res[[comb[2,z]]], stat='D', mask = F)
}


##### Within species ----
wsp.res <- list()

for (i in 1:length(sp.list)) {
  sp <- sp.list[i]
  
  results <- list.files(paste0("results/", sp, "/predictions"),
                        pattern = "_pa", full.names = T)
  results <- results[grep("q0.5", results)]
  
  # Load rasters generated before
  res <- stack(results)
  names(res) <- c("current", paste0("ssp", c(1, 2, 3, 5)))
  
  comb <- combn(c("current", paste0("ssp", c(1, 2, 3, 5))), 2)
  
  wspecies <- data.frame(scen_a = NA, scen_b = NA,
                         D = rep(NA, dim(comb)[2]), I = NA)
  
  for (z in 1:nrow(bspecies)) {
    wspecies[z,"scen_a"] <- comb[1,z];wspecies[z,"scen_b"] <- comb[2,z]
    wspecies[z,"I"] <- dismo::nicheOverlap(res[[comb[1,z]]], res[[comb[2,z]]], stat='I', mask = F)
    wspecies[z,"D"] <- dismo::nicheOverlap(res[[comb[1,z]]], res[[comb[2,z]]], stat='D', mask = F)
  }
  
  wspecies$species <- sp
  
  wsp.res[[i]] <- wspecies
}

wsp.res <- do.call("rbind", wsp.res)



##### Within species Jaccard ----
get.jacc <- function(sp, type, crop.south = F){
  
  ##### Presence only dataset ----
  po.pts <- read.csv(paste0("data/", sp, "/", sp, "_final.csv"))
  po.pts <- SpatialPoints(po.pts[, c("decimalLongitude", "decimalLatitude")],
                          proj4string = crs(env))
  
  po.pts <- remove.duplicates(po.pts, zero = 5)
  
  ##### Presence-absence dataset ----
  pa.pts <- read.csv(paste0("data/", sp, "/", sp, "_pa.csv"))
  pa.pts <- SpatialPointsDataFrame(pa.pts[, c("decimalLongitude", "decimalLatitude")],
                                   data.frame("presence" = pa.pts[,3]),
                                   proj4string = crs(env))
  
  #pa.pts <- remove.duplicates(pa.pts, zero = 5)
  pa.pts <- dup.cells(pa.pts, rast(env))
  
  # Just to ensure all are falling inside study area
  pa.pts <- pa.pts[!is.na(extract(env, pa.pts)),]
  
  
  
  # Load results ----
  results <- list.files(paste0("results/", sp, "/predictions"),
                        pattern = type, full.names = T)
  results <- results[grep("q0.5", results)]
  
  # Load rasters generated before
  res <- stack(results)
  names(res) <- c("current", paste0("ssp", c(1, 2, 3, 5)))
  
  if (sp == "spam" | sp == "scze") {
    if (crop.south) {
      res <- crop(res, extent(-99, -29, -42.5, 3))
    }
  }
  
  # Get thresholded contour
  if (type == "int" | type == "cont") {
    mtp <- min(extract(res[[1]], po.pts))
    p10 <- rev(sort(extract(res[[1]], po.pts)))[ceiling(length(po.pts) * 0.9)]
  } else {
    pres <- pa.pts[pa.pts$presence == 1,]
    mtp <- min(extract(res[[1]], pres))
    p10 <- rev(sort(extract(res[[1]], pres)))[ceiling(length(pres) * 0.9)]
  }
  
  res <- res >= p10 # could also be mtp or other threshold
  
  comb <- combn(c("current", paste0("ssp", c(1, 2, 3, 5))), 2)
  
  wspecies <- data.frame(scen_a = NA, scen_b = NA, jacc = rep(NA, dim(comb)[2]))
  
  for (z in 1:nrow(bspecies)) {
    wspecies[z,"scen_a"] <- comb[1,z];wspecies[z,"scen_b"] <- comb[2,z]
    combination <- res[[comb[1,z]]] + res[[comb[2,z]]]
    intersection <- combination == 2
    union <- combination >= 1
    
    wspecies$jacc[z] <- freq(intersection)[2,2] / freq(union)[2,2]
  }
  
  wspecies$species <- sp
  
  return(wspecies)
}

int.jacc <- lapply(sp.list, get.jacc, type = "int")
pa.jacc <- lapply(sp.list, get.jacc, type = "_pa")

int.jacc <- do.call("rbind", int.jacc)
pa.jacc <- do.call("rbind", pa.jacc)

write.csv(int.jacc, "results/jaccard_int_lgcp.csv", row.names = F)
write.csv(pa.jacc, "results/jaccard_pa_lgcp.csv", row.names = F)


# Difference metrics ----
get.diffs <- function(sp, type = "_pa", crop.south = F){
  # Load species data ----
  # Load a layer just as sample
  env <- raster("data/env/crop_layers/tempmax.tif")
  
  ##### Presence only dataset ----
  po.pts <- read.csv(paste0("data/", sp, "/", sp, "_final.csv"))
  po.pts <- SpatialPoints(po.pts[, c("decimalLongitude", "decimalLatitude")],
                          proj4string = crs(env))
  
  po.pts <- remove.duplicates(po.pts, zero = 5)
  
  ##### Presence-absence dataset ----
  pa.pts <- read.csv(paste0("data/", sp, "/", sp, "_pa.csv"))
  pa.pts <- SpatialPointsDataFrame(pa.pts[, c("decimalLongitude", "decimalLatitude")],
                                   data.frame("presence" = pa.pts[,3]),
                                   proj4string = crs(env))
  
  #pa.pts <- remove.duplicates(pa.pts, zero = 5)
  pa.pts <- dup.cells(pa.pts, rast(env))
  
  # Just to ensure all are falling inside study area
  pa.pts <- pa.pts[!is.na(extract(env, pa.pts)),]
  
  
  
  # Load results ----
  results <- list.files(paste0("results/", sp, "/predictions"),
                        pattern = type, full.names = T)
  results <- results[grep("q0.5", results)]
  
  # Load rasters generated before
  curr <- raster(results[1])
  ssp1 <- raster(results[2])
  ssp2 <- raster(results[3])
  ssp3 <- raster(results[4])
  ssp5 <- raster(results[5])
  
  # Get thresholded rasters
  if (type == "int" | type == "cont") {
    mtp <- min(extract(curr, po.pts))
    p10 <- rev(sort(extract(curr, po.pts)))[ceiling(length(po.pts) * 0.9)]
  } else {
    pres <- pa.pts[pa.pts$presence == 1,]
    mtp <- min(extract(curr, pres))
    p10 <- rev(sort(extract(curr, pres)))[ceiling(length(pres) * 0.9)]
  }
  
  curr <- curr >= p10
  ssp1 <- ssp1 >= p10
  ssp2 <- ssp2 >= p10
  ssp3 <- ssp3 >= p10
  ssp5 <- ssp5 >= p10
  
  ssp1 <- (ssp1 * 2) - curr
  ssp2 <- (ssp2 * 2) - curr
  ssp3 <- (ssp3 * 2) - curr
  ssp5 <- (ssp5 * 2) - curr
  
  ar <- raster::area(curr)
  
  all.difs <- stack(ssp1, ssp2, ssp3, ssp5)
  
  if (sp == "spam" | sp == "scze") {
    if (crop.south) {
      all.difs <- crop(all.difs, extent(-99, -29, -42.5, 3))
      curr <- crop(curr, all.difs); ar <- crop(ar, all.difs)
    }
  }

  all.difs <- data.frame(rasterToPoints(all.difs))
  
  all.difs$area <- extract(ar, all.difs[,1:2])
  
  res <- data.frame(
    scenario = paste0("ssp", c(1,2,3,5)),
    lost = NA, gain = NA, kept = NA
  )
  
  for (i in 1:4) {
    res$lost[i] <- sum(all.difs[all.difs[,i+2] == -1, "area"]) 
    res$kept[i] <- sum(all.difs[all.difs[,i+2] == 1, "area"]) 
    res$gain[i] <- sum(all.difs[all.difs[,i+2] == 2, "area"]) 
  }
  
  curr <- data.frame(rasterToPoints(curr))
  res[5,1] <- c("current");res[5,2:3] <- 0
  res[5,4] <- sum(extract(ar, curr[curr$layer == 1,1:2]))
  
  res$total <- res$gain + res$kept
  
  res$species <- sp
  
  return(res)
  
}

sp.list <- c("acch", "scze", "spam", "mybo", "lujo")

pa.diffs <- lapply(sp.list, get.diffs)

pa.diffs <- do.call("rbind", pa.diffs)

write.csv(pa.diffs, "results/thresh_diffs_pa_lgcp.csv", row.names = F)

pa.diffs.cr <- lapply(sp.list, get.diffs, crop.south = T)

pa.diffs.cr <- do.call("rbind", pa.diffs.cr)

write.csv(pa.diffs.cr, "results/thresh_diffs_south_pa_lgcp.csv", row.names = F)
