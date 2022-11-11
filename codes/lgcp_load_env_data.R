#### Modelling of reef fishes ####
## Silas C. Principe - silascprincipe@gmail.com - 2022 ##

######## Integrated Species Distribution Model ###########

## Load environmental data ##

# Load packages ----
library(terra)

# Load current period data extended (for models) ----
env.e <- rast(c("data/env/ready_layers/tempmean.tif",
              "data/env/ready_layers/tempmax.tif",
              "data/env/ready_layers/salinitymean.tif",
              "data/env/ready_layers/chlomean.tif",
              "data/env/ready_layers/silicatemax.tif",
              "data/env/ready_layers/ph.tif",
              "data/env/ready_layers/windspeed.tif",
              "data/env/ready_layers/distcoast.tif"))


# Load current period data ----
env <- rast(c("data/env/crop_layers/tempmean.tif",
             "data/env/crop_layers/tempmax.tif",
             "data/env/crop_layers/salinitymean.tif",
             "data/env/crop_layers/chlomean.tif",
             "data/env/crop_layers/silicatemax.tif",
             "data/env/crop_layers/ph.tif",
             "data/env/crop_layers/windspeed.tif",
             "data/env/crop_layers/distcoast.tif"))

# Put chl-a in log
env$chlomean <- log(env$chlomean)

# Get mean and SD to scale future layers
m <- global(env, fun = mean, na.rm = T)[,1]
sde <- global(env, fun = sd, na.rm = T)[,1]

# Scale data
env <- scale(env)

# Correct names
names(env.e)[8] <- names(env)[8] <- "distcoast"



# Load future layers ----
lays <- c("data/env/fut_layers/tempmean.tif",
          "data/env/fut_layers/tempmax.tif",
          "data/env/fut_layers/salinitymean.tif",
          "data/env/fut_layers/chlomean.tif",
          "data/env/fut_layers/silicatemax.tif",
          "data/env/fut_layers/ph.tif",
          "data/env/fut_layers/windspeed.tif")

# SSP1
ssp1 <- rast(c(gsub(".tif", "_ssp126.tif", lays),
              "data/env/crop_layers/distcoast.tif"))

names(ssp1) <- names(env)

ssp1$chlomean <- log(ssp1$chlomean)

ssp1 <- (ssp1 - m)/sde

# SSP2
ssp2 <- rast(c(gsub(".tif", "_ssp245.tif", lays),
              "data/env/crop_layers/distcoast.tif"))

names(ssp2) <- names(env)

ssp2$chlomean <- log(ssp2$chlomean)

ssp2 <- (ssp2 - m)/sde

# SSP3
ssp3 <- rast(c(gsub(".tif", "_ssp370.tif", lays),
              "data/env/crop_layers/distcoast.tif"))

names(ssp3) <- names(env)

ssp3$chlomean <- log(ssp3$chlomean)

ssp3 <- (ssp3 - m)/sde

# SSP5
ssp5 <- rast(c(gsub(".tif", "_ssp585.tif", lays),
              "data/env/crop_layers/distcoast.tif"))

names(ssp5) <- names(env)

ssp5$chlomean <- log(ssp5$chlomean)

ssp5 <- (ssp5 - m)/sde

# Verify differences:
# delta1 <- ssp1 - env;plot(delta1)
# delta2 <- ssp2 - env;plot(delta2)
# delta3 <- ssp3 - env;plot(delta3)
# delta5 <- ssp5 - env;plot(delta5)

rm(m, sde, lays)

# END