#### Modelling of reef fishes ####
## Silas C. Principe - silascprincipe@gmail.com - 2022 ##

######## Integrated Species Distribution Model ###########

## Run integrated models ##

# Load packages and define settings ----
# Modeling
library(INLA)
library(inlabru)
# Spatial data
library(terra)
library(sp)
library(raster)
# Plotting and utilities
library(ggplot2)
library(patchwork)
library(fs)
source("functions/auxiliary_lgcp.R")
source("functions/response_curves.R")
# Settings
set.seed(2932) # Replicability of sampling
spatt <- theme_classic() # Theme for better ploting
nsamp <- 2000 # Define number of sampling for the final predictions
nsampcv <- 1000 # Define number of sampling in the CV predictions
itnumb <- 40 # Define number of maximum inlabru iterations
itnumbcv <- 10 # Define number of maximum inlabru iterations in cross-validation
intest <- "auto" # Integration strategy

# Species (each one is modeled separately)
sp <- "lujo"



# Load prepared data ----
# Data prepared with the 'lgcp_prepare_data.R' code
dat <- readRDS("data/lgcp_data.rds")
mesh <- dat$mesh
starea <- dat$starea
rm(dat)

##### Load environmental data ----
source("codes/lgcp_load_env_data.R")

# Get integration points
ips <- ipoints(starea, mesh)

# Plot mesh and integration points
ggplot()+gg(mesh)+gg(ips)+coord_fixed()+spatt



# Load species data ----
##### Presence only dataset ----
po.pts <- read.csv(paste0("data/", sp, "/", sp, "_final.csv"))
po.pts <- SpatialPoints(po.pts[, c("decimalLongitude", "decimalLatitude")],
                        proj4string = CRS(crs(env.e, proj = T)))

po.pts <- remove.duplicates(po.pts, zero = 5)

# Ensure coordnames are correct
coordnames(po.pts) <- c("x", "y")

##### Presence-absence dataset ----
pa.pts <- read.csv(paste0("data/", sp, "/", sp, "_pa.csv"))
pa.pts <- SpatialPointsDataFrame(pa.pts[, c("decimalLongitude", "decimalLatitude")],
                                 data.frame("presence" = pa.pts[,3]),
                                 proj4string = CRS(crs(env.e, proj = T)))

#pa.pts <- remove.duplicates(pa.pts, zero = 5)
pa.pts <- dup.cells(pa.pts, env[[1]])

# Just to ensure all are falling inside study area
pa.pts <- pa.pts[!is.na(extract(env[[1]], coordinates(pa.pts))[,1]),]



# Plot all to see
ggplot()+
  gg(po.pts, color = "blue", alpha = .5) +
  gg(pa.pts, aes(color = as.factor(presence)), alpha = .5, size = 4, shape = 2) +
  scale_color_manual(values = c("orange", "black")) +
  coord_equal() + spatt



# Prepare the barrier model ----
# Based on Bakka, H., J. Vanhatalo, J. Illian, D. Simpson, and H. Rue. 2016.
# “Accounting for Physical Barriers in Species Distribution Modeling with
# Non-Stationary Spatial Random Effects.” ArXiv preprint arXiv:1608.03787.
# Norwegian University of Science; Technology, Trondheim, Norway.
# https://haakonbakkagit.github.io/btopic128.html

# Get number of triangles of the mesh
tl <- length(mesh$graph$tv[,1])

# Create a matrix with the central coordinates of each triangle
tri.cord <- matrix(0, tl, 2)

for(i in 1:tl){
  # Take the vertex of triangles
  temp <- mesh$loc[mesh$graph$tv[i, ], ]
  
  # Compute center of each triangle
  tri.cord[i,] <- colMeans(temp)[c(1,2)]
}

# Convert to SpatialPoints
tri.cord <- SpatialPoints(tri.cord)

# Get intersection between mesh points and study area
crs(tri.cord) <- raster::crs(starea)
intersec <- over(starea, tri.cord, returnList = T)

# Remove the study area triangles from all to obtain the barrier ones
intersec <- unlist(intersec)
barrier.triangles <- setdiff(1:tl, intersec)

# Create a barrier polygon, i.e. all the polygons composing the "islands"
poly.barrier <- inla.barrier.polygon(mesh, barrier.triangles)



# Prepare matérn model for the spatial field ----
size <- diff(bbox(starea)[1,])
range0 <- as.numeric(round(size/10))
sigma <- 0.1

# Stationary model for comparison (if wanted)
# spde <- inla.spde2.pcmatern(mesh,
#                             prior.range = c(range0, 0.1),
#                             prior.sigma = c(sigma, 0.01))

# Barrier model
b.model <- inla.barrier.pcmatern(mesh,
                                 barrier.triangles = barrier.triangles,
                                 prior.range = c(range0, 0.01),
                                 prior.sigma = c(sigma, 0.01))

source("functions/barrier_model_plot.R")
# par(mfrow = c(1,2), mar = c(5, 5, 4, 5) + 0.1)
plot.bmodel(po.pts@coords[100,1:2], mesh = mesh, spde = b.model,
            areapol = poly.barrier, crs = CRS(proj), range = range0, msd = sigma)
# To plot the stationary version:
# plot.bmodel(po.pts@coords[100,1:2], mesh = mesh, spde = spde,
#             areapol = poly.barrier, spmode = T, range = range0)



# Prepare 1D SPDE for SST component ----
# The 1D SPDE acts similar to a GAM spline, fitting a [possible] non-linear
# relation. It's similar also to a "rw2" model (smoother than the "rw1"),
# but have as advantage that you don't need to group the values prior to fitting.
knots.st <- seq((minmax(env.e$tempmax)[1]-0.01), (minmax(ssp5$tempmax)[2]+0.01),
                length = 25)
d1mesh.st <- inla.mesh.1d(knots.st, degree = 2,
                          boundary = "free")

d1spde.st <- inla.spde2.pcmatern(d1mesh.st,
                                 prior.range = c(diff(c(min(knots.st), max(knots.st))), NA),
                                 prior.sigma = c(1, 0.01),
                                 constr = T)

# If you want to compare to the Random-walk of order 2 just uncoment:
# env.e$tempmax <- inla.group(env.e$tempmax, n = 30)
# sst.prior <- list(prior = "pcprec", param = c(1, 0.01))
# And use the following line in the components, instead of the d1spde.st
# tempmax(env.e, model = "rw1", hyper = list(prec = sst.prior), main_layer = "tempmax") +



# Prepare model components and formulas ----
# We will test 4 distinct models
# Model 1 - Essential covariates
# SST + SALmean + pH + spatial
# Model 2 - Mixing/hydrodynamics
# SST + SALmean + pH + wind + spatial
# Model 3 - Productivity
# SST + SALmean + pH + silicate + chlomean + spatial
# Model 4 - Full
# SST + SALmean + pH + wind + silicate + chlomean + spatial

# Components
cmp <- list(
  ~ tempmean(env.e, model = d1spde.st, main_layer = "tempmean") +
    salinitymean(env.e, model = "linear", mean.linear = 0, prec.linear = 0.01, main_layer = "salinitymean") +
    ph(env.e, model = "linear", mean.linear = 0, prec.linear = 0.01, main_layer = "ph") +
    spatial(coordinates, model = b.model, mapper = bru_mapper(mesh)) +
    spatial_pa(coordinates, copy = "spatial", fixed = FALSE) +
    Intercept(1)+
    intercept_pa(1)
)

cmp[[2]] <- update(cmp[[1]], ~ . + 
                     windspeed(env.e, model = "linear", mean.linear = 0, prec.linear = 0.01, main_layer = "windspeed"))
cmp[[3]] <- update(cmp[[1]], ~ . + 
                     silicatemax(env.e, model = "linear", mean.linear = 0, prec.linear = 0.01, main_layer = "silicatemax") +
                     chlomean(env.e, model = "linear", mean.linear = 0, prec.linear = 0.01, main_layer = "chlomean"))
cmp[[4]] <- update(cmp[[3]], ~ . + 
                     windspeed(env.e, model = "linear", mean.linear = 0, prec.linear = 0.01, main_layer = "windspeed"))

# Formulas
forms <- list(
  ~ tempmean + salinitymean + ph
)

forms[[2]] <- update(forms[[1]], ~ . + windspeed)
forms[[3]] <- update(forms[[1]], ~ . + silicatemax + chlomean)
forms[[4]] <- update(forms[[3]], ~ . + windspeed)



# Fit models ----
m <- list()

for (i in 1:length(forms)) {
  
  cat("\nRunning:", as.character(forms[[i]]), "\n")
  
  # Presence only (Log-Gaussian Cox Process)
  lik.lgcp <- like(data = po.pts,
                   family = "cp",
                   ips = ips,
                   domain = list(coordinates = mesh),
                   formula = update(forms[[i]], coordinates ~ . + Intercept + spatial))
  
  # Presence-absence (Binomial with cloglog link)
  lik.pa <- like(data =  pa.pts,
                 family = "binomial",
                 formula = update(forms[[i]], presence ~ . + intercept_pa + spatial_pa),
                 control.family = list(link = "cloglog"))
  
  # Fit
  m[[i]] <- bru(cmp[[i]],
                lik.lgcp,
                lik.pa,
                options = list(
                  bru_max_iter = itnumb, # Number of iterations for optimization
                  control.inla = list(int.strategy = intest),
                  #control.compute = list(waic = TRUE, cpo = TRUE, dic = TRUE, config = TRUE),
                  inla.mode = "classic"
                  # For verbosity, uncoment those lines:
                  , bru_verbose = TRUE
                  #, verbose = TRUE
                ))
}



# Verify results ----
tm <- 4 # <- change number to see each summary

# Summary and plots
summary(m[[tm]])

# See spatial effect --- Theta1 = Sigma, Theta2 = Range
exp(m[[tm]]$summary.hyperpar[2:3,])

# Predict to view
pxl <- as(raster::raster("data/env/crop_layers/tempmean.tif"),
          "SpatialPixelsDataFrame")

pred.f <- forms[[tm]] 

pred.fit <- predict(m[[tm]], pxl,
                    ~ list(
                      full = eval(parse(text = paste0("exp(",
                                                      update.formula(pred.f, ~ . + Intercept + spatial)[2],
                                                      ")"))),
                      lin = eval(parse(text = paste0(update.formula(pred.f, ~ . + Intercept + spatial)[2]))),
                      pa = eval(parse(text = paste0("1-exp(-exp(",
                                                      update.formula(pred.f, ~ . + intercept_pa + spatial_pa)[2],
                                                      "))"))),
                      spatial = spatial,
                      sst = tempmean,
                      sal = salinitymean
                      # Other variables can be added...
                    ))

plot.res("full", pred.fit)
plot.res("pa", pred.fit)
plot.res("spatial", pred.fit) + plot.res("sst", pred.fit) + plot.res("sal", pred.fit)

# Plot 1D SPDE (SST component spline)
plot.temp(m[[tm]])

# Get and plot response curves
resp.curves <- get.resp.curves(m[[tm]], forms[[tm]], mode = "exp")
plot(resp.curves)

# Get WAIC
deltaIC(m[[1]], m[[2]], m[[3]], m[[4]], criterion = "WAIC")

# Estimate abundance
pred.lamb <- predict(m[[tm]], ips,
                     eval(parse(text = paste0("~ sum(weight * exp(",
                                              update.formula(pred.f, ~ . + Intercept + spatial)[2],
                                              "))"))));pred.lamb



# Cross-validation ----
# Now we proceed cross-validation. We can CV all models or 
# just a subset to save time.

# Chose the subset to be CV
cv.m <- 1:4

# Prepare CV data
# Get spatial blocks
latpol <- raster::raster(raster::extent(starea),
                         crs = raster::crs(starea), nrows = 30, ncols = 1)
latpol <- raster::rasterToPolygons(latpol)

latpol$ID <- rep(1:5, 6)
latpol <- intersect(latpol, starea)
plot(latpol, col = as.factor(latpol$ID));lines(starea)

# Get the integration points of each block (we need to generate the IPS again
# with the argument blocks, so we correctly get the integration points involved 
# with the block underlying field)
ips.blocks <- ipoints(samplers = latpol, domain = mesh,
                      name = c("x", "y"), group = "ID")

# Get the presence-only blocks
po.blocks <- over(po.pts, latpol)[,2]

table(ips.blocks$ID)
table(po.blocks)

# Get blocks for the presence-absence data
pa.blocks <- over(pa.pts, latpol)[,2]
table(pa.pts$presence, pa.blocks)

# Run cross-validation
# Create a list to hold CVs
cv.res <- list()

for (z in cv.m) {
  
  cat("\n=== Running CV for model", z, "===\n")
  
  # Create data frames to hold results
  cv.df.pa <- data.frame(
    orig = pa.pts$presence,
    pred_pa_mean = NA,
    pred_pa_q025 = NA,
    pred_pa_q975 = NA
  )
  
  cv.vals <- data.frame(
    tvals = length(po.pts),
    pvals = rep(NA, max(po.blocks)),
    block_tvals = as.numeric(table(po.blocks)),
    block_pvals = NA
  )
  
  # Produce cross-validated plots (optional)
  #plot.list <- list()
  
  # Run blocks
  for (k in 1:max(po.blocks)) {
    
    cat("Running block", k, "\n")
    
    train.po <- po.pts[po.blocks != k,]

    train.pa <- pa.pts[pa.blocks != k,]
    
    # cat("-- Integration points:", length(train.ips), "\n")
    cat("-- Presence-only points:", length(train.po), "\n")
    cat("-- Presence:", sum(train.pa$presence == 1), "Absence:", 
        sum(train.pa$presence == 0), "\n")
    
    # Presence-only
    lik.lgcp <- like(data = train.po,
                     family = "cp",
                     ips = ips.blocks[ips.blocks$ID != k,],
                     domain = list(coordinates = mesh),
                     formula = update(forms[[z]], coordinates ~ . + Intercept + spatial))
    
    # Presence-absence (Binomial with cloglog link)
    lik.pa <- like(data =  train.pa,
                   family = "binomial",
                   formula = update(forms[[z]], presence ~ . + intercept_pa + spatial_pa),
                   control.family = list(link = "cloglog"))
    
    # Run CV model
    blockm <- bru(cmp[[z]],
                  lik.lgcp,
                  lik.pa,
                  #lik.abund,
                  options = list(
                    bru_max_iter = itnumbcv,
                    control.inla = list(int.strategy = intest),
                    control.mode = list(restart = T, theta = m[[z]]$mode$theta),
                    control.compute = list(waic = TRUE, cpo = TRUE, dic = TRUE, config = TRUE),
                    inla.mode = "classic"
                    # For verbosity, uncoment those lines:
                    , bru_verbose = TRUE
                    #, verbose = TRUE
                  ))
    
    # Predict
    # Presence-absence component
    cv.pred.pa <- predict(blockm, pa.pts[pa.blocks == k,],
                          as.formula(paste0(
                            "~1-",
                            update(forms[[z]], ~ exp(-exp(. + intercept_pa + spatial_pa)))[2]
                          )), n.samples = nsampcv)
    
    cv.df.pa$pred_pa_mean[pa.blocks == k] <- cv.pred.pa$mean
    cv.df.pa$pred_pa_q025[pa.blocks == k] <- cv.pred.pa$q0.025
    cv.df.pa$pred_pa_q975[pa.blocks == k] <- cv.pred.pa$q0.975
    
    # Integrated 
    pred.int <- predict(blockm, ips.blocks,
                        update(forms[[z]], ~ sum(weight * exp(. + Intercept + spatial))),
                        n.samples = nsampcv)
    
    pred.int.block <- predict(blockm, ips.blocks[ips.blocks$ID == k,],
                        update(forms[[z]], ~ sum(weight * exp(. + Intercept + spatial))),
                        n.samples = nsampcv)
    
    # Produce cross-validated maps (optional)
    # plot.list[[k]] <- predict(blockm, pxl,
    #                           update(forms[[z]], ~ exp(. + Intercept + spatial)),
    #                           n.samples = nsampcv)
    
    cv.vals$pvals[k] <- pred.int$mean
    cv.vals$block_pvals[k] <- pred.int.block$mean
    cv.vals$pvals_sd[k] <- pred.int$sd
    cv.vals$block_pvals_sd[k] <- pred.int.block$sd
    
    rm(pred.int, pred.int.block, cv.pred.pa)
    
  }
  
  cat("Getting full model predictions... \n")
  
  # Predict full model
  pred.pa <- predict(m[[z]], pa.pts,
                        as.formula(paste0(
                          "~1-",
                          update(forms[[z]], ~ exp(-exp(. + intercept_pa + spatial_pa)))[2]
                        )), n.samples = nsampcv)
  
  pred.int <- predict(m[[z]], ips,
                      update(forms[[z]], ~ sum(weight * exp(. + Intercept + spatial))),
                      n.samples = nsampcv)
  
  
  # Presence-absence metrics
  tdat <- cbind(
    ID = 1,
    cv.df.pa[,c("orig", "pred_pa_mean")],
    pred.pa$mean)
  colnames(tdat) <- c("ID", "y", "cv", "full")
  
  pa.metrics <- data.frame(
    auc = c(PresenceAbsence::auc(tdat, st.dev=F), 
            PresenceAbsence::auc(tdat, st.dev=F, which.model = 2))
  )
  
  pa.metrics <- cbind(ID = z, mod = c("cv", "full"), pa.metrics)
  
  cv.res[[z]] <- list(
    pa = pa.metrics,
    b_vals = cv.vals,
    b_full = pred.int$mean
  )
  
  rm(pred.pa, pred.int)
  
}

# See results
pa.cv <- lapply(cv.res, function(x){
  pa <- x$pa
  data.frame(full_auc = pa$auc[2],
             cv_auc = pa$auc[1])
})
metrics.cv <- lapply(cv.res, function(x){
  cv <- x$b_vals
  post.e <- cv$pvals
  post.e.b <- cv$block_pvals
  post.var <- cv$pvals + cv$pvals_sd^2
  post.var.b <- cv$block_pvals + cv$block_pvals_sd^2
  # Dawid-Sebastiani 
  dss <- (cv$tvals - post.e)^2 / post.var + log(post.var)
  dss.b <- (cv$block_tvals - post.e.b)^2 / post.var.b + log(post.var.b)
  # RMSE
  f <- sqrt(mean((cv$tvals - cv$pvals)^2))
  b <- sqrt(mean((cv$block_tvals - cv$block_pvals)^2))
  data.frame(dss = mean(dss), dss_sd = sd(dss),
             block_dss = mean(dss.b), block_dss_sd = sd(dss.b),
             rmse = f, block_rmse = b,
             full_pred_dev = cv$tvals[1] - x$b_full)
  })

pa.cv <- do.call("rbind", pa.cv)
metrics.cv <- do.call("rbind", metrics.cv)
pa.cv$model <- metrics.cv$model <- cv.m

pa.cv
metrics.cv

delta.metrics(pa.cv, all = T)
delta.metrics(metrics.cv, all = c("dss", "block_dss", "block_rmse"))



# Predict models ----
# After analyzing models and the cross-validation metrics we chose model:
smodel <- 4
# Now we get predictions for this model

# Remove unused objects to ensure predictions go smoothly
rm(pred.fit)

##### Save individual components effects for the chosen model ----
# This was already done after running the models as a exploratory tool
# Now we run again and save the results

# Create a function to save results
dir <- paste("results", sp, "effects", sep = "/")
dir_create(dir)

# Avoid GDAL to save xml aux file
rgdal::setCPLConfigOption("GDAL_PAM_ENABLED", "FALSE")

save.rast <- function(x, model){
  r <- as(pred.comp, "RasterStack")
  to.remove <- names(r)[!names(r) %in% c("mean", "sd", "q0.025", "q0.5", "q0.975")]
  r <- dropLayer(r, to.remove)
  for (i in 1:nlayers(r)) {
    writeRaster(r[[i]], paste0(dir, "/", sp, "_m", model, "_",
                               x, "_", names(r)[i], "_effect.tif"),
                overwrite = T, format="GTiff")
  }
  return(paste(x, "saved."))
}

# Predict model components
# Because predicting all effects at the same time produces a big Spatial object
# this can make the computation slower. So we predict each one separately
# and save to more efficiency:
topred <- c(
  m[[smodel]]$names.fixed[-grep("ntercept", m[[smodel]]$names.fixed)],
  names(m[[smodel]]$summary.random)
)

for (i in 1:length(topred)) {
  cat(topred[i], "\n")
  pred.comp <- predict(m[[smodel]], data = pxl,
                       formula = as.formula(paste0("~", topred[i])),
                       n.samples = nsamp, seed = 2932)
  
  # save
  save.rast(x = topred[i], model = smodel)
}

##### Predictions ----
# inlabru gets the environmental layer from the R environment
# thus, to predict to different scenarios we need to change the environmental
# layers object (i.e. 'env.e') to the new object for which we want the prediction

# We generate 3 predictions: 
# 1 with the spatial component, 1 without it (contrast) and 1 of the
# presence-absence component

dir <- paste("results", sp, "predictions", sep = "/")
dir_create(dir)

for (i in 1:3) {
  
  # Get formula
  pred.f <- list(
    update.formula(forms[[smodel]], ~exp(. + Intercept + spatial)),
    update.formula(forms[[smodel]], ~exp(. + Intercept)),
    as.formula(paste0(
      "~1-",
      update(forms[[smodel]], ~ exp(-exp(. + intercept_pa + spatial_pa)))[2]
    ))
  )[[i]]
  
  # Save name
  snam <- paste0("_m", smodel, "_", c("int", "cont", "pa")[i])
  
  cat("Running", snam, "\n")
  print(pred.f)
  
  # Predict to current scenario
  cat("Predicting for current scenario... \n")
  env.e <- env
  
  pred.cur <- predict(m[[smodel]],
                      data = pxl,
                      formula = pred.f,
                      n.samples = nsamp, seed = 2932)
  
  (plot.res(obj = pred.cur))
  
  # Predict to SSP 1
  cat("Predicting for SSP1 scenario... \n")
  env.e <- ssp1
  
  pred.ssp1 <- predict(m[[smodel]],
                      data = pxl,
                      formula = pred.f,
                      n.samples = nsamp, seed = 2932)
  
  (plot.res(obj = pred.ssp1))
  
  # Predict to SSP 2
  cat("Predicting for SSP2 scenario... \n")
  env.e <- ssp2
  
  pred.ssp2 <- predict(m[[smodel]],
                      data = pxl,
                      formula = pred.f,
                      n.samples = nsamp, seed = 2932)
  
  (plot.res(obj = pred.ssp2))
  
  # Predict to SSP 3
  cat("Predicting for SSP3 scenario... \n")
  env.e <- ssp3
  
  pred.ssp3 <- predict(m[[smodel]],
                      data = pxl,
                      formula = pred.f,
                      n.samples = nsamp, seed = 2932)
  
  (plot.res(obj = pred.ssp3))
  
  # Predict to SSP 5
  cat("Predicting for SSP5 scenario... \n")
  env.e <- ssp5
  
  pred.ssp5 <- predict(m[[smodel]],
                      data = pxl,
                      formula = pred.f,
                      n.samples = nsamp, seed = 2932)
  
  (plot.res(obj = pred.ssp5))
  
  
  # Convert all to raster (easier handling)
  cat("Saving files... \n")
  pred.cur <- as(pred.cur, "RasterStack")
  pred.ssp1 <- as(pred.ssp1, "RasterStack")
  pred.ssp2 <- as(pred.ssp2, "RasterStack")
  pred.ssp3 <- as(pred.ssp3, "RasterStack")
  pred.ssp5 <- as(pred.ssp5, "RasterStack")
  
  to.remove <- names(pred.cur)[!names(pred.cur) %in% 
                                 c("mean", "sd", "q0.025", "q0.5", "q0.975")]
  
  pred.cur <- dropLayer(pred.cur, to.remove)
  pred.ssp1 <- dropLayer(pred.ssp1, to.remove)
  pred.ssp2 <- dropLayer(pred.ssp2, to.remove)
  pred.ssp3 <- dropLayer(pred.ssp3, to.remove)
  pred.ssp5 <- dropLayer(pred.ssp5, to.remove)
  
  # Save results
  dir <- paste("results", sp, "predictions", sep = "/")
  dir_create(dir)
  
  # Raster was not being saved by layers for some reason
  # Thus we just use a small function to get the work done
  save.rast <- function(x, fnames){
    for (z in 1:nlayers(x)) {
      writeRaster(x[[z]], fnames[z], overwrite = T)
    }
    return(invisible(NULL))
  }
  
  save.rast(pred.cur, paste0(dir, "/", sp, "_",
                             names(pred.cur), snam,
                             "_current.tif"))
  
  save.rast(pred.ssp1, paste0(dir, "/", sp, "_",
                             names(pred.ssp1), snam,
                             "_ssp1.tif"))
  
  save.rast(pred.ssp2, paste0(dir, "/", sp, "_",
                             names(pred.ssp2), snam,
                             "_ssp2.tif"))
  
  save.rast(pred.ssp3, paste0(dir, "/", sp, "_",
                             names(pred.ssp3), snam,
                             "_ssp3.tif"))
  
  save.rast(pred.ssp5, paste0(dir, "/", sp, "_",
                             names(pred.ssp5), snam,
                             "_ssp5.tif"))
  
  cat("Done! \n")
}



# Save metrics and summaries ----
dir <- paste0("results/", sp)

# Extract summaries of the chosen model
model.summ <- m[[smodel]]$summary.fixed

sigma.mar <- inla.tmarginal(function(x) exp(x),
                            m[[smodel]]$marginals.hyperpar[[2]])
sigma.mar <- inla.zmarginal(sigma.mar)

range.mar <- inla.tmarginal(function(x) exp(x),
                            m[[smodel]]$marginals.hyperpar[[3]])
range.mar <- inla.zmarginal(range.mar)

sst.sd.mar <- inla.zmarginal(m[[smodel]]$marginals.hyperpar[[1]])

hyper <- rbind(cbind(as.data.frame(sigma.mar), mode = NA, kld = NA),
               cbind(as.data.frame(range.mar), mode = NA, kld = NA),
               cbind(as.data.frame(sst.sd.mar), mode = NA, kld = NA))
row.names(hyper) <- c("sigma", "range", "temp_sd")
hyper <- hyper[,-c(4, 6)]
names(hyper) <- names(model.summ)

model.summ <- rbind(model.summ, hyper)

# Save summaries
write.csv(model.summ, paste0(dir, "/", sp, "_model_summary.csv"))
write.csv(data.frame(
  model = 1:4,
  waic = sapply(1:4, function(x){m[[x]]$waic$waic})
), paste0(dir, "/", sp, "_model_waic.csv"))

# Save CV metrics
write.csv(metrics.cv, paste0(dir, "/", sp, "_cv_metrics.csv"))
write.csv(pa.cv, paste0(dir, "/", sp, "_pa_cv_metrics.csv"))

# Get plots of probability density
plots <- list()

# x label titles
xl <- c(
  Intercept = "Intercept", intercept_pa = "Intercept Presence-Absence",
  ph = "Mean pH", tempmean = "Mean SST", chlomean = "Mean log(Chl-a)",
  silicatemax = "Maximum silicate", windspeed = "Mean windspeed",
  salinitymean = "Mean salinity"
)
xl <- xl[grep(paste0(m[[smodel]]$names.fixed, collapse = "|"), names(xl))]

for (i in 1:length(m[[smodel]]$names.fixed)) {
  plots[[i]] <- plot(m[[smodel]], m[[smodel]]$names.fixed[i])+
    ylab(ifelse(i == 1, "Density", ""))+
    xlab(xl[grep(paste0(m[[smodel]]$names.fixed[i], collapse = "|"), names(xl))])+
    theme_classic()
}

pl <- eval(parse(text = paste("plots[[", 1:length(plots), "]]", collapse = "+")))
ggsave(paste0(dir, "/effects_density.jpg"), pl, quality = 100)

# Get plots of probability density for the SPDE parameters
sigp <- ggplot(data.frame(
  inla.smarginal(
    inla.tmarginal(function(x) exp(x),
                   m[[smodel]]$marginals.hyperpar$`Theta1 for spatial`))),
  aes(x, y)) + geom_line() + xlab("Sigma") + ylab("Density") + theme_classic()

rangep <- ggplot(data.frame(
  inla.smarginal(
    inla.tmarginal(function(x) exp(x),
                   m[[smodel]]$marginals.hyperpar$`Theta2 for spatial`))),
  aes(x, y)) + geom_line() + xlab("Range") + ylab("") + theme_classic()

sstp <- ggplot(data.frame(
  inla.smarginal(m[[smodel]]$marginals.hyperpar[[1]])),
  aes(x, y)) + geom_line() + xlab("SD(SST)") + ylab("") + theme_classic()

sigp + rangep + sstp
ggsave(paste0(dir, "/hyperpar_density.jpg"), height = 3, quality = 100)

# Save plot of the temperature component
ssteval <- predict(m[[smodel]],
                   NULL,
                   ~ tempmean_eval(seq(minmax(env$tempmean)[1], minmax(env$tempmean)[2], by = 0.1)),
                   nsamples = 1000)

ssteval$x <- seq(minmax(env$tempmean)[1], minmax(env$tempmean)[2], by = 0.1)

(p <- ggplot(ssteval) +
  geom_line(aes(x = x, y = mean)) +
  geom_line(aes(x = x, y = q0.5), linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_ribbon(aes(x = x, ymin = q0.025, ymax = q0.975), alpha = .4) +
  scale_x_continuous(expand = c(0,0))+
  spatt + ylab("Mean effect") + xlab("Mean SST"))
ggsave(paste0(dir, "/sst_effect.jpg"), quality = 100)

extvals <- c(
  seq(minmax(env$tempmean)[1], minmax(env$tempmean)[2], length.out = 100),
  seq(minmax(env$tempmean)[2], minmax(ssp5$tempmean)[2], by = 0.1)
)

sstext <- predict(m[[smodel]],
                  NULL,
                   ~ tempmean_eval(extvals), nsamples = 1000)
sstext$class <- c(rep("True values", 100), rep("Extrapolation", nrow(sstext)-100))
sstext$x <- extvals

(p <- ggplot(sstext) +
    geom_line(aes(x = x, y = mean, color = class)) +
    geom_line(aes(x = x, y = q0.5, color = class), linetype = "dotted") +
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = max(sstext$x[sstext$class == "True values"]),
               color = "red", linetype = "dashed")+
    geom_ribbon(aes(x = x, ymin = q0.025, ymax = q0.975, fill = class), alpha = .4) +
    scale_x_continuous(expand = c(0,0))+
    spatt + ylab("Mean effect") + xlab("Mean SST")) + theme(legend.title = element_blank())
ggsave(paste0(dir, "/sst_effect_extrapolation.jpg"), quality = 100)

# Save response curves
resp.curves <- get.resp.curves(m[[smodel]], forms[[smodel]], mode = NULL, samp = 1000)
plot(resp.curves)
ggsave(paste0(dir, "/resp_curves_lin.jpg"), quality = 100)

resp.curves <- get.resp.curves(m[[smodel]], forms[[smodel]], mode = "exp", samp = 1000)
plot(resp.curves)
ggsave(paste0(dir, "/resp_curves_exp.jpg"), quality = 100)

resp.curves <- get.resp.curves(m[[smodel]], forms[[smodel]], mode = "cloglog", samp = 1000)
plot(resp.curves)
ggsave(paste0(dir, "/resp_curves_pa.jpg"), quality = 100)



# Save session info ----
dir.info <- paste0(dir, "/", sp, "_sessioninfo.txt")
cat("\n ============= \n \n Research info \n \n
    Author: Silas C. Principe | silasprincipe@usp.br \n
    Co-authors: Tito M.C. Lotufo | André L. Acosta \n
    Modelling of reef fishes \n
    This work is part of a PhD project being held at the Oceanographic Institute - USP \n \n",
    "Settings used in the modeling: \n",
    "   Integration strategy:", intest, "\n",
    "   Number of samplings prediction:", nsamp, "\n",
    "   Number of samplings CV:", nsampcv, "\n",
    "   Number of inlabru iterations:", itnumb, "\n",
    "   Number of inlabru iterations on CV:", itnumbcv, "\n",
    file = dir.info)
cat("\n ============= \n \n Session info \n \n", file = dir.info, append = T)
write(capture.output(sessionInfo()), file = dir.info, append = T)
cat("\n ============= \n \n INLA info \n \n", file = dir.info, append = T)
write(capture.output(inla.version()), file = dir.info, append = T)

## END OF CODE
