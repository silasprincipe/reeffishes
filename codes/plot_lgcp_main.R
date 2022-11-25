#### Modelling of coral reef herbivorous invertebrates ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2022 ##

# Plot of LGCP results - all species
# Note: to save individual maps, just uncomment the lines with 'ggsave'

# Load needed packages ----
library(ggplot2)
library(sf)
library(raster)
library(terra)
library(patchwork)
source("functions/auxiliary_lgcp.R")




# Load base shapefiles ----
base <- shapefile("data/gis/basemaps/ne_110m_land_edited.shp")
# Crop to the extent
base <- buffer(base, 0)
base <- crop(base, extent(-120, -10, -55, 65))
# Convert to SF
base <- st_as_sf(base)

base.pac <- shapefile("data/gis/basemaps/ne_110m_land_no_oc.shp")
# Crop to the extent
base.pac <- buffer(base.pac, 0)
base.pac <- crop(base.pac, extent(-120, -10, -55, 65))
# Convert to SF
base.pac <- st_as_sf(base.pac)



# Define species
sp <- "acch"

# Define type 
# [int = relative occurrence rate; cont = contrast (ROR without spatial);
# _pa = probability of occurrence (binomial likelihood)]
# NOTE: use _pa with spam, or the file will be wrong!
type <- "int"

# If type is "int" or "cont" it's possible to use
# the quantile mode for ploting
mod <- "quant"

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


# Get thresholded contour
if (type == "int" | type == "cont") {
  mtp <- min(extract(curr, po.pts))
  p10 <- rev(sort(extract(curr, po.pts)))[ceiling(length(po.pts) * 0.9)]
} else {
  pres <- pa.pts[pa.pts$presence == 1,]
  mtp <- min(extract(curr, pres))
  p10 <- rev(sort(extract(curr, pres)))[ceiling(length(pres) * 0.9)]
}


get.cont <- function(th, rast){
  r <- rast
  r[r < th] <- NA
  r[r >= th] <- 1
  r <- rasterToPolygons(r, dissolve = T)
  r <- st_as_sf(r)
  r
}

curr.pol <- get.cont(p10, curr)
ssp1.pol <- get.cont(p10, ssp1)
ssp2.pol <- get.cont(p10, ssp2)
ssp3.pol <- get.cont(p10, ssp3)
ssp5.pol <- get.cont(p10, ssp5)


if (type == "int" & mod == "quant" | type == "cont" & mod == "quant") {
  qq <- quantile(curr, seq(0,1, by = .1))
  qq[1] <- 0 # Set as 0 because min may be slightly more
  
  fint <- function(x){ # Function to get the interval
    findInterval(x, qq, rightmost.closed = T, all.inside = T)
  }
  # Classify
  curr <- calc(curr, fint)
  ssp1 <- calc(ssp1, fint)
  ssp2 <- calc(ssp2, fint)
  ssp3 <- calc(ssp3, fint)
  ssp5 <- calc(ssp5, fint)
}

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

# Just if type == "int" or "cont"
# Highlight (set to NA) values higher than the maximum in current period
# and also normalize to 0-1 to easier handling of legends
if (type == "int" & mod != "quant" | type == "cont" & mod != "quant") {
  ssp1.v$val <- ifelse(ssp1.v$val > max(curr.v$val), NA, ssp1.v$val)
  ssp2.v$val <- ifelse(ssp2.v$val > max(curr.v$val), NA, ssp2.v$val)
  ssp3.v$val <- ifelse(ssp3.v$val > max(curr.v$val), NA, ssp3.v$val)
  ssp5.v$val <- ifelse(ssp5.v$val > max(curr.v$val), NA, ssp5.v$val)
  
  ssp1.v$val <- (ssp1.v$val - min(curr.v$val))/(max(curr.v$val) - min(curr.v$val))
  ssp2.v$val <- (ssp2.v$val - min(curr.v$val))/(max(curr.v$val) - min(curr.v$val))
  ssp3.v$val <- (ssp3.v$val - min(curr.v$val))/(max(curr.v$val) - min(curr.v$val))
  ssp5.v$val <- (ssp5.v$val - min(curr.v$val))/(max(curr.v$val) - min(curr.v$val))
  curr.v$val <- (curr.v$val - min(curr.v$val))/(max(curr.v$val) - min(curr.v$val))
}



# Prepare theme/ploting stuff ----
# Guide for legend
if (type == "int" | type == "cont") {
  step.guide <- guide_colorbar(title = "Relative occurrence rate",
                               show.limits = TRUE,
                               barheight = unit(0.12, "in"),
                               barwidth = unit(3.5, "in"),
                               ticks = F,
                               ticks.colour = "grey20",
                               frame.colour = "grey20",
                               title.position = "top")
  
  if (mod != "quant") {
    # Scale of values
    sca <- scale_fill_gradientn(
      colours = rev(c("#7d1500", "#da4325", "#eca24e", "#e7e2bc", "#5cc3af", "#0a6265")),
      limits = c(0,1),
      guide = step.guide,
      labels = c("Low", "", "", "", "High"),
      na.value = "#2b0700")
  } else {
    # Scale of values
    sca <- scale_fill_stepsn(
      colours = rev(c("#7d1500", "#da4325", "#eca24e", "#e7e2bc", "#5cc3af", "#0a6265")),
      breaks = 0:10,
      show.limits = T,
      limits = c(0,10),
      guide = step.guide,
      labels = c("Low", "", "","","", "", "", "", "", "", "High"),
      na.value = "#2b0700")
  }
  
} else {
  step.guide <- guide_colorbar(title = "Probability of occurrence",
                               show.limits = TRUE,
                               barheight = unit(0.12, "in"),
                               barwidth = unit(3.5, "in"),
                               ticks = T,
                               ticks.colour = "grey20",
                               frame.colour = "grey20",
                               title.position = "top")
  
  # Scale of values
  sca <- scale_fill_gradientn(
    colours = rev(c("#7d1500", "#da4325", "#eca24e", "#e7e2bc", "#5cc3af", "#0a6265")),
    limits = c(0,1),
    guide = step.guide,
    labels = c("0", "", "0.5", "", "1"),
    na.value = "#2b0700")
}


# Themes
nlt <- theme_classic()+
        theme(axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.title.x.top = element_text(vjust = 2, size = 16),
              axis.title.y.left = element_text(size = 16),
              legend.position="none",
              panel.grid.major = element_line(linetype = 'dashed', 
                                              colour = "grey70",
                                              size = .05),
              panel.background = element_blank()#,
              #plot.background = element_rect(color = "grey30", fill = 'white')
        )

wlt <- theme_classic()+
        theme(axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_text(colour = "grey60", size = 14),
              axis.title.x.top = element_text(vjust = 2, size = 16),
              axis.title.y.left = element_text(size = 16),
              panel.background = element_blank(),
              panel.border = element_rect(fill = NA, color = "grey60"),
              panel.grid.major = element_line(linetype = 'dashed', 
                                              colour = "grey70",
                                              size = .1),
              legend.position="bottom",
              legend.title.align=0.5,
              legend.text = element_text(size = 16),
              legend.title = element_text(size = 16),
              legend.background = element_rect(fill = "white"),
              
        )

int <- theme_classic()+
        theme(axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              legend.position="none",
              panel.grid.major = element_blank(),
              panel.background = element_blank(),
              plot.background = element_rect(color = "grey30", fill = 'white'),
              plot.margin = margin(0.5, 0.5, 0, 0, "pt")
        )



# Generate plots ----
##### Current ----
# Rectangles to highlight
rects <- list(data.frame(x1 = -85, x2 = -72, y1 = 20, y2= 30, label = "a"),
              data.frame(x1 = -40, x2 = -32, y1 = -21, y2= -15, label = "b"),
              data.frame(x1 = -37, x2 = -32, y1 = -8, y2= -3, label = "c"))

(pc <- ggplot()+
         # Base maps
         geom_sf(data = base.pac,
                 color = NA,
                 fill = c("#D9D9D9"),
                 alpha = 0.5) +
         geom_sf(data = base,
                 color = c("#CFCFCF"),
                 size = 0.6,
                 fill = "white") +
         # Get the raster
         geom_tile(data = curr.v, aes(x = x, y = y, fill = val)) +
         sca + # Add color scale
         # Add MTP polygon
         geom_sf(data = curr.pol,
                 color = "gray30",
                 alpha = 0.4,
                 size = 0.2,
                 fill = NA)+
         # Add occurrence points
         # geom_point(data = data.frame(pres@coords),
         #            aes(x = x, y = y), size = 1.2,
         #            alpha = .6, shape = 16, color = "yellow")+
         # Establish area
         coord_sf(xlim = c(-99, -29),
                  ylim = c(-42.5, 42.5),
                  datum = st_crs(base),
                  expand = F,
                  label_axes = list(
                          top = "E",
                          left = "N",
                          top = ""
                  )) +
         # Draw rectangles
         geom_rect(data = rects[[1]],
                   aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
                   fill = NA, color = "grey20", size = .3)+
         # geom_text(data = rects[[1]],
         #            aes(x = (x2 - 100), y = (y2 - 100), label = label))+
         geom_rect(data = rects[[2]],
                   aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
                   fill = NA, color = "grey20", size = .3)+
         # geom_text(data = rects[[1]],
         #            aes(x = (x2 - 100), y = (y2 - 100), label = label))+
         geom_rect(data = rects[[3]],
                   aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
                   fill = NA, color = "grey20", size = .3)+
         # geom_text(data = rects[[1]],
         #            aes(x = (x2 - 100), y = (y2 - 100), label = label))+
         # Add title
         geom_label(aes(label = "A", x = -35, y = 38),
                    size = 10, fontface = "bold", color = "grey30",
                    label.size = 0)+
         # Add grid
         scale_x_discrete(position = "top", breaks = c(-95, -75, -55, -35)) +
         # Remove axis labels and add theme
         xlab("") + ylab("") + wlt
         
)

# Prepare insets
pc.i1 <- ggplot()+
                # Base maps
                geom_sf(data = base,
                        color = c("#CFCFCF"),
                        size = 0.6,
                        fill = "white") +
                # Get the raster
                geom_tile(data = curr.v, aes(x = x, y = y, fill = val)) +
                sca + # Add color scale
                # Add MTP polygon
                geom_sf(data = curr.pol,
                        color = "gray30",
                        alpha = 0.4,
                        size = 0.2,
                        fill = NA)+
                # Establish area
                coord_sf(xlim = c(rects[[1]]$x1, rects[[1]]$x2),
                         ylim = c(rects[[1]]$y1, rects[[1]]$y2),
                         datum = st_crs(base),
                         expand = F) + int

pc.i2 <- ggplot()+
                # Base maps
                geom_sf(data = base,
                        color = c("#CFCFCF"),
                        size = 0.6,
                        fill = "white") +
                # Get the raster
                geom_tile(data = curr.v, aes(x = x, y = y, fill = val)) +
                sca + # Add color scale
                # Add MTP polygon
                geom_sf(data = curr.pol,
                        color = "gray30",
                        alpha = 0.4,
                        size = 0.2,
                        fill = NA)+
                # Establish area
                coord_sf(xlim = c(rects[[2]]$x1, rects[[2]]$x2),
                         ylim = c(rects[[2]]$y1, rects[[2]]$y2),
                         datum = st_crs(base),
                         expand = F) + int

pc.i3 <- ggplot()+
                # Base maps
                geom_sf(data = base,
                        color = c("#CFCFCF"),
                        size = 0.6,
                        fill = "white") +
                # Get the raster
                geom_tile(data = curr.v, aes(x = x, y = y, fill = val)) +
                sca + # Add color scale
                # Add MTP polygon
                geom_sf(data = curr.pol,
                        color = "gray30",
                        alpha = 0.4,
                        size = 0.2,
                        fill = NA)+
                # Establish area
                coord_sf(xlim = c(rects[[3]]$x1, rects[[3]]$x2),
                         ylim = c(rects[[3]]$y1, rects[[3]]$y2),
                         datum = st_crs(base),
                         expand = F) + int

(pcf <- pc +
        annotate(
                "segment",
                x = c(rects[[1]]$x2, -68, -62),
                y = c((rects[[1]]$y1+rects[[1]]$y2)/2,
                      (rects[[2]]$y1+rects[[2]]$y2)/2,
                      (rects[[3]]$y1+rects[[3]]$y2)/2),
                xend = c(-58, rects[[2]]$x1, rects[[3]]$x1),
                yend = c((rects[[1]]$y1+rects[[1]]$y2)/2,
                         (rects[[2]]$y1+rects[[2]]$y2)/2,
                         (rects[[3]]$y1+rects[[3]]$y2)/2),
                lineend = "round",
                colour = "grey20",
                size = 0.3
        ) +
        annotation_custom(
                grob = ggplotGrob(pc.i1),
                xmin = -63,
                xmax = -40,
                ymin = 21,
                ymax = 39) +
        annotation_custom(
                grob = ggplotGrob(pc.i2),
                xmin = -93,
                xmax = -68,
                ymin = -39,
                ymax = -14) +
        annotation_custom(
                grob = ggplotGrob(pc.i3),
                xmin = -75,
                xmax = -55,
                ymin = -12,
                ymax = 3))


# ggsave(paste0("figures/", sp, "_lgcp_current_m", m, "_", type, ".jpg"), pcf,
#        width = 16, height = 18, units = "cm", quality = 100)



##### SSP1 ----
(ps1 <- ggplot()+
   # Base maps
   geom_sf(data = base.pac,
           color = NA,
           fill = c("#D9D9D9"),
           alpha = 0.5) +
   geom_sf(data = base,
           color = c("#CFCFCF"),
           size = 0.6,
           fill = "white") +
   # Get the raster
   geom_tile(data = ssp1.v, aes(x = x, y = y, fill = val)) +
   sca + # Add color scale
   # Add MTP polygon
   geom_sf(data = ssp1.pol,
           color = "gray30",
           alpha = 0.4,
           size = 0.2,
           fill = NA)+
   # Add occurrence points
   # geom_point(data = data.frame(pres@coords),
   #            aes(x = x, y = y), size = 1.2,
   #            alpha = .6, shape = 16, color = "yellow")+
   # Establish area
   coord_sf(xlim = c(-99, -29),
            ylim = c(-42.5, 42.5),
            datum = st_crs(base),
            expand = F,
            label_axes = list(
              top = "E",
              left = "",
              top = ""
            )) +
   # Draw rectangles
   geom_rect(data = rects[[1]],
             aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
             fill = NA, color = "grey20", size = .3)+
   # geom_text(data = rects[[1]],
   #            aes(x = (x2 - 100), y = (y2 - 100), label = label))+
   geom_rect(data = rects[[2]],
             aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
             fill = NA, color = "grey20", size = .3)+
   # geom_text(data = rects[[1]],
   #            aes(x = (x2 - 100), y = (y2 - 100), label = label))+
   geom_rect(data = rects[[3]],
             aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
             fill = NA, color = "grey20", size = .3)+
   # geom_text(data = rects[[1]],
   #            aes(x = (x2 - 100), y = (y2 - 100), label = label))+
   # Add title
   geom_label(aes(label = "B", x = -35, y = 38),
              size = 10, fontface = "bold", color = "grey30",
              label.size = 0)+
   # Add grid
   scale_x_discrete(position = "top", breaks = c(-95, -75, -55, -35)) +
   # Remove axis labels and add theme
   xlab("") + ylab("") + wlt
 
)

# Prepare insets
ps1.i1 <- ggplot()+
  # Base maps
  geom_sf(data = base,
          color = c("#CFCFCF"),
          size = 0.6,
          fill = "white") +
  # Get the raster
  geom_tile(data = ssp1.v, aes(x = x, y = y, fill = val)) +
  sca + # Add color scale
  # Add MTP polygon
  geom_sf(data = ssp1.pol,
          color = "gray30",
          alpha = 0.4,
          size = 0.2,
          fill = NA)+
  # Establish area
  coord_sf(xlim = c(rects[[1]]$x1, rects[[1]]$x2),
           ylim = c(rects[[1]]$y1, rects[[1]]$y2),
           datum = st_crs(base),
           expand = F) + int

ps1.i2 <- ggplot()+
  # Base maps
  geom_sf(data = base,
          color = c("#CFCFCF"),
          size = 0.6,
          fill = "white") +
  # Get the raster
  geom_tile(data = ssp1.v, aes(x = x, y = y, fill = val)) +
  sca + # Add color scale
  # Add MTP polygon
  geom_sf(data = ssp1.pol,
          color = "gray30",
          alpha = 0.4,
          size = 0.2,
          fill = NA)+
  # Establish area
  coord_sf(xlim = c(rects[[2]]$x1, rects[[2]]$x2),
           ylim = c(rects[[2]]$y1, rects[[2]]$y2),
           datum = st_crs(base),
           expand = F) + int

ps1.i3 <- ggplot()+
  # Base maps
  geom_sf(data = base,
          color = c("#CFCFCF"),
          size = 0.6,
          fill = "white") +
  # Get the raster
  geom_tile(data = ssp1.v, aes(x = x, y = y, fill = val)) +
  sca + # Add color scale
  # Add MTP polygon
  geom_sf(data = ssp1.pol,
          color = "gray30",
          alpha = 0.4,
          size = 0.2,
          fill = NA)+
  # Establish area
  coord_sf(xlim = c(rects[[3]]$x1, rects[[3]]$x2),
           ylim = c(rects[[3]]$y1, rects[[3]]$y2),
           datum = st_crs(base),
           expand = F) + int

(ps1f <- ps1 +
    annotate(
      "segment",
      x = c(rects[[1]]$x2, -68, -62),
      y = c((rects[[1]]$y1+rects[[1]]$y2)/2,
            (rects[[2]]$y1+rects[[2]]$y2)/2,
            (rects[[3]]$y1+rects[[3]]$y2)/2),
      xend = c(-58, rects[[2]]$x1, rects[[3]]$x1),
      yend = c((rects[[1]]$y1+rects[[1]]$y2)/2,
               (rects[[2]]$y1+rects[[2]]$y2)/2,
               (rects[[3]]$y1+rects[[3]]$y2)/2),
      lineend = "round",
      colour = "grey20",
      size = 0.3
    ) +
    annotation_custom(
      grob = ggplotGrob(ps1.i1),
      xmin = -63,
      xmax = -40,
      ymin = 21,
      ymax = 39) +
    annotation_custom(
      grob = ggplotGrob(ps1.i2),
      xmin = -93,
      xmax = -68,
      ymin = -39,
      ymax = -14) +
    annotation_custom(
      grob = ggplotGrob(ps1.i3),
      xmin = -75,
      xmax = -55,
      ymin = -12,
      ymax = 3))

ps1f.s <- ps1f + theme(legend.position = "none")

# ggsave(paste0("figures/", sp, "_lgcp_ssp1_m", m, "_", type, ".jpg"), ps1f.s,
#        width = 16, height = 18, units = "cm", quality = 100)




##### SSP2 ----
(ps2 <- ggplot()+
   # Base maps
   geom_sf(data = base.pac,
           color = NA,
           fill = c("#D9D9D9"),
           alpha = 0.5) +
   geom_sf(data = base,
           color = c("#CFCFCF"),
           size = 0.6,
           fill = "white") +
   # Get the raster
   geom_tile(data = ssp2.v, aes(x = x, y = y, fill = val)) +
   sca + # Add color scale
   # Add MTP polygon
   geom_sf(data = ssp2.pol,
           color = "gray30",
           alpha = 0.4,
           size = 0.2,
           fill = NA)+
   # Add occurrence points
   # geom_point(data = data.frame(pres@coords),
   #            aes(x = x, y = y), size = 1.2,
   #            alpha = .6, shape = 16, color = "yellow")+
   # Establish area
   coord_sf(xlim = c(-99, -29),
            ylim = c(-42.5, 42.5),
            datum = st_crs(base),
            expand = F,
            label_axes = list(
              top = "E",
              left = "",
              top = ""
            )) +
   # Draw rectangles
   geom_rect(data = rects[[1]],
             aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
             fill = NA, color = "grey20", size = .3)+
   # geom_text(data = rects[[1]],
   #            aes(x = (x2 - 100), y = (y2 - 100), label = label))+
   geom_rect(data = rects[[2]],
             aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
             fill = NA, color = "grey20", size = .3)+
   # geom_text(data = rects[[1]],
   #            aes(x = (x2 - 100), y = (y2 - 100), label = label))+
   geom_rect(data = rects[[3]],
             aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
             fill = NA, color = "grey20", size = .3)+
   # geom_text(data = rects[[1]],
   #            aes(x = (x2 - 100), y = (y2 - 100), label = label))+
   # Add title
   geom_label(aes(label = "C", x = -35, y = 38),
              size = 10, fontface = "bold", color = "grey30",
              label.size = 0)+
   # Add grid
   scale_x_discrete(position = "top", breaks = c(-95, -75, -55, -35)) +
   # Remove axis labels and add theme
   xlab("") + ylab("") + wlt
 
)

# Prepare insets
ps2.i1 <- ggplot()+
  # Base maps
  geom_sf(data = base,
          color = c("#CFCFCF"),
          size = 0.6,
          fill = "white") +
  # Get the raster
  geom_tile(data = ssp2.v, aes(x = x, y = y, fill = val)) +
  sca + # Add color scale
  # Add MTP polygon
  geom_sf(data = ssp2.pol,
          color = "gray30",
          alpha = 0.4,
          size = 0.2,
          fill = NA)+
  # Establish area
  coord_sf(xlim = c(rects[[1]]$x1, rects[[1]]$x2),
           ylim = c(rects[[1]]$y1, rects[[1]]$y2),
           datum = st_crs(base),
           expand = F) + int

ps2.i2 <- ggplot()+
  # Base maps
  geom_sf(data = base,
          color = c("#CFCFCF"),
          size = 0.6,
          fill = "white") +
  # Get the raster
  geom_tile(data = ssp2.v, aes(x = x, y = y, fill = val)) +
  sca + # Add color scale
  # Add MTP polygon
  geom_sf(data = ssp2.pol,
          color = "gray30",
          alpha = 0.4,
          size = 0.2,
          fill = NA)+
  # Establish area
  coord_sf(xlim = c(rects[[2]]$x1, rects[[2]]$x2),
           ylim = c(rects[[2]]$y1, rects[[2]]$y2),
           datum = st_crs(base),
           expand = F) + int

ps2.i3 <- ggplot()+
  # Base maps
  geom_sf(data = base,
          color = c("#CFCFCF"),
          size = 0.6,
          fill = "white") +
  # Get the raster
  geom_tile(data = ssp2.v, aes(x = x, y = y, fill = val)) +
  sca + # Add color scale
  # Add MTP polygon
  geom_sf(data = ssp2.pol,
          color = "gray30",
          alpha = 0.4,
          size = 0.2,
          fill = NA)+
  # Establish area
  coord_sf(xlim = c(rects[[3]]$x1, rects[[3]]$x2),
           ylim = c(rects[[3]]$y1, rects[[3]]$y2),
           datum = st_crs(base),
           expand = F) + int

(ps2f <- ps2 +
    annotate(
      "segment",
      x = c(rects[[1]]$x2, -68, -62),
      y = c((rects[[1]]$y1+rects[[1]]$y2)/2,
            (rects[[2]]$y1+rects[[2]]$y2)/2,
            (rects[[3]]$y1+rects[[3]]$y2)/2),
      xend = c(-58, rects[[2]]$x1, rects[[3]]$x1),
      yend = c((rects[[1]]$y1+rects[[1]]$y2)/2,
               (rects[[2]]$y1+rects[[2]]$y2)/2,
               (rects[[3]]$y1+rects[[3]]$y2)/2),
      lineend = "round",
      colour = "grey20",
      size = 0.3
    ) +
    annotation_custom(
      grob = ggplotGrob(ps2.i1),
      xmin = -63,
      xmax = -40,
      ymin = 21,
      ymax = 39) +
    annotation_custom(
      grob = ggplotGrob(ps2.i2),
      xmin = -93,
      xmax = -68,
      ymin = -39,
      ymax = -14) +
    annotation_custom(
      grob = ggplotGrob(ps2.i3),
      xmin = -75,
      xmax = -55,
      ymin = -12,
      ymax = 3))

ps2f.s <- ps2f + theme(legend.position = "none")

# ggsave(paste0("figures/", sp, "_lgcp_ssp2_m", m, "_", type, ".jpg"), ps2f.s,
#        width = 16, height = 18, units = "cm", quality = 100)




##### SSP3 ----
(ps3 <- ggplot()+
   # Base maps
   geom_sf(data = base.pac,
           color = NA,
           fill = c("#D9D9D9"),
           alpha = 0.5) +
   geom_sf(data = base,
           color = c("#CFCFCF"),
           size = 0.6,
           fill = "white") +
   # Get the raster
   geom_tile(data = ssp3.v, aes(x = x, y = y, fill = val)) +
   sca + # Add color scale
   # Add MTP polygon
   geom_sf(data = ssp3.pol,
           color = "gray30",
           alpha = 0.4,
           size = 0.2,
           fill = NA)+
   # Add occurrence points
   # geom_point(data = data.frame(pres@coords),
   #            aes(x = x, y = y), size = 1.2,
   #            alpha = .6, shape = 16, color = "yellow")+
   # Establish area
   coord_sf(xlim = c(-99, -29),
            ylim = c(-42.5, 42.5),
            datum = st_crs(base),
            expand = F,
            label_axes = list(
              top = "E",
              left = "",
              top = ""
            )) +
   # Draw rectangles
   geom_rect(data = rects[[1]],
             aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
             fill = NA, color = "grey20", size = .3)+
   # geom_text(data = rects[[1]],
   #            aes(x = (x2 - 100), y = (y2 - 100), label = label))+
   geom_rect(data = rects[[2]],
             aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
             fill = NA, color = "grey20", size = .3)+
   # geom_text(data = rects[[1]],
   #            aes(x = (x2 - 100), y = (y2 - 100), label = label))+
   geom_rect(data = rects[[3]],
             aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
             fill = NA, color = "grey20", size = .3)+
   # geom_text(data = rects[[1]],
   #            aes(x = (x2 - 100), y = (y2 - 100), label = label))+
   # Add title
   geom_label(aes(label = "D", x = -35, y = 38),
              size = 10, fontface = "bold", color = "grey30",
              label.size = 0)+
   # Add grid
   scale_x_discrete(position = "top", breaks = c(-95, -75, -55, -35)) +
   # Remove axis labels and add theme
   xlab("") + ylab("") + wlt
 
)

# Prepare insets
ps3.i1 <- ggplot()+
  # Base maps
  geom_sf(data = base,
          color = c("#CFCFCF"),
          size = 0.6,
          fill = "white") +
  # Get the raster
  geom_tile(data = ssp3.v, aes(x = x, y = y, fill = val)) +
  sca + # Add color scale
  # Add MTP polygon
  geom_sf(data = ssp3.pol,
          color = "gray30",
          alpha = 0.4,
          size = 0.2,
          fill = NA)+
  # Establish area
  coord_sf(xlim = c(rects[[1]]$x1, rects[[1]]$x2),
           ylim = c(rects[[1]]$y1, rects[[1]]$y2),
           datum = st_crs(base),
           expand = F) + int

ps3.i2 <- ggplot()+
  # Base maps
  geom_sf(data = base,
          color = c("#CFCFCF"),
          size = 0.6,
          fill = "white") +
  # Get the raster
  geom_tile(data = ssp3.v, aes(x = x, y = y, fill = val)) +
  sca + # Add color scale
  # Add MTP polygon
  geom_sf(data = ssp3.pol,
          color = "gray30",
          alpha = 0.4,
          size = 0.2,
          fill = NA)+
  # Establish area
  coord_sf(xlim = c(rects[[2]]$x1, rects[[2]]$x2),
           ylim = c(rects[[2]]$y1, rects[[2]]$y2),
           datum = st_crs(base),
           expand = F) + int

ps3.i3 <- ggplot()+
  # Base maps
  geom_sf(data = base,
          color = c("#CFCFCF"),
          size = 0.6,
          fill = "white") +
  # Get the raster
  geom_tile(data = ssp3.v, aes(x = x, y = y, fill = val)) +
  sca + # Add color scale
  # Add MTP polygon
  geom_sf(data = ssp3.pol,
          color = "gray30",
          alpha = 0.4,
          size = 0.2,
          fill = NA)+
  # Establish area
  coord_sf(xlim = c(rects[[3]]$x1, rects[[3]]$x2),
           ylim = c(rects[[3]]$y1, rects[[3]]$y2),
           datum = st_crs(base),
           expand = F) + int

(ps3f <- ps3 +
    annotate(
      "segment",
      x = c(rects[[1]]$x2, -68, -62),
      y = c((rects[[1]]$y1+rects[[1]]$y2)/2,
            (rects[[2]]$y1+rects[[2]]$y2)/2,
            (rects[[3]]$y1+rects[[3]]$y2)/2),
      xend = c(-58, rects[[2]]$x1, rects[[3]]$x1),
      yend = c((rects[[1]]$y1+rects[[1]]$y2)/2,
               (rects[[2]]$y1+rects[[2]]$y2)/2,
               (rects[[3]]$y1+rects[[3]]$y2)/2),
      lineend = "round",
      colour = "grey20",
      size = 0.3
    ) +
    annotation_custom(
      grob = ggplotGrob(ps3.i1),
      xmin = -63,
      xmax = -40,
      ymin = 21,
      ymax = 39) +
    annotation_custom(
      grob = ggplotGrob(ps3.i2),
      xmin = -93,
      xmax = -68,
      ymin = -39,
      ymax = -14) +
    annotation_custom(
      grob = ggplotGrob(ps3.i3),
      xmin = -75,
      xmax = -55,
      ymin = -12,
      ymax = 3))

ps3f.s <- ps3f + theme(legend.position = "none")

# ggsave(paste0("figures/", sp, "_lgcp_ssp3_m", m, "_", type, ".jpg"), ps3f.s,
#        width = 16, height = 18, units = "cm", quality = 100)



##### SSP5 ----
(ps5 <- ggplot()+
   # Base maps
   geom_sf(data = base.pac,
           color = NA,
           fill = c("#D9D9D9"),
           alpha = 0.5) +
   geom_sf(data = base,
           color = c("#CFCFCF"),
           size = 0.6,
           fill = "white") +
   # Get the raster
   geom_tile(data = ssp5.v, aes(x = x, y = y, fill = val)) +
   sca + # Add color scale
   # Add MTP polygon
   geom_sf(data = ssp5.pol,
           color = "gray30",
           alpha = 0.4,
           size = 0.2,
           fill = NA)+
   # Add occurrence points
   # geom_point(data = data.frame(pres@coords),
   #            aes(x = x, y = y), size = 1.2,
   #            alpha = .6, shape = 16, color = "yellow")+
   # Establish area
   coord_sf(xlim = c(-99, -29),
            ylim = c(-42.5, 42.5),
            datum = st_crs(base),
            expand = F,
            label_axes = list(
              top = "E",
              left = "N",
              top = ""
            )) +
   # Draw rectangles
   geom_rect(data = rects[[1]],
             aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
             fill = NA, color = "grey20", size = .3)+
   # geom_text(data = rects[[1]],
   #            aes(x = (x2 - 100), y = (y2 - 100), label = label))+
   geom_rect(data = rects[[2]],
             aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
             fill = NA, color = "grey20", size = .3)+
   # geom_text(data = rects[[1]],
   #            aes(x = (x2 - 100), y = (y2 - 100), label = label))+
   geom_rect(data = rects[[3]],
             aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
             fill = NA, color = "grey20", size = .3)+
   # geom_text(data = rects[[1]],
   #            aes(x = (x2 - 100), y = (y2 - 100), label = label))+
   # Add title
   # geom_label(aes(label = "C", x = -35, y = 38),
   #            size = 10, fontface = "bold", color = "grey30",
   #            label.size = 0)+
   # Add grid
   scale_x_discrete(position = "top", breaks = c(-95, -75, -55, -35)) +
   # Remove axis labels and add theme
   xlab("") + ylab("") + wlt
 
)

# Prepare insets
ps5.i1 <- ggplot()+
  # Base maps
  geom_sf(data = base,
          color = c("#CFCFCF"),
          size = 0.6,
          fill = "white") +
  # Get the raster
  geom_tile(data = ssp5.v, aes(x = x, y = y, fill = val)) +
  sca + # Add color scale
  # Add MTP polygon
  geom_sf(data = ssp5.pol,
          color = "gray30",
          alpha = 0.4,
          size = 0.2,
          fill = NA)+
  # Establish area
  coord_sf(xlim = c(rects[[1]]$x1, rects[[1]]$x2),
           ylim = c(rects[[1]]$y1, rects[[1]]$y2),
           datum = st_crs(base),
           expand = F) + int

ps5.i2 <- ggplot()+
  # Base maps
  geom_sf(data = base,
          color = c("#CFCFCF"),
          size = 0.6,
          fill = "white") +
  # Get the raster
  geom_tile(data = ssp5.v, aes(x = x, y = y, fill = val)) +
  sca + # Add color scale
  # Add MTP polygon
  geom_sf(data = ssp5.pol,
          color = "gray30",
          alpha = 0.4,
          size = 0.2,
          fill = NA)+
  # Establish area
  coord_sf(xlim = c(rects[[2]]$x1, rects[[2]]$x2),
           ylim = c(rects[[2]]$y1, rects[[2]]$y2),
           datum = st_crs(base),
           expand = F) + int

ps5.i3 <- ggplot()+
  # Base maps
  geom_sf(data = base,
          color = c("#CFCFCF"),
          size = 0.6,
          fill = "white") +
  # Get the raster
  geom_tile(data = ssp5.v, aes(x = x, y = y, fill = val)) +
  sca + # Add color scale
  # Add MTP polygon
  geom_sf(data = ssp5.pol,
          color = "gray30",
          alpha = 0.4,
          size = 0.2,
          fill = NA)+
  # Establish area
  coord_sf(xlim = c(rects[[3]]$x1, rects[[3]]$x2),
           ylim = c(rects[[3]]$y1, rects[[3]]$y2),
           datum = st_crs(base),
           expand = F) + int

(ps5f <- ps5 +
    annotate(
      "segment",
      x = c(rects[[1]]$x2, -68, -62),
      y = c((rects[[1]]$y1+rects[[1]]$y2)/2,
            (rects[[2]]$y1+rects[[2]]$y2)/2,
            (rects[[3]]$y1+rects[[3]]$y2)/2),
      xend = c(-58, rects[[2]]$x1, rects[[3]]$x1),
      yend = c((rects[[1]]$y1+rects[[1]]$y2)/2,
               (rects[[2]]$y1+rects[[2]]$y2)/2,
               (rects[[3]]$y1+rects[[3]]$y2)/2),
      lineend = "round",
      colour = "grey20",
      size = 0.3
    ) +
    annotation_custom(
      grob = ggplotGrob(ps5.i1),
      xmin = -63,
      xmax = -40,
      ymin = 21,
      ymax = 39) +
    annotation_custom(
      grob = ggplotGrob(ps5.i2),
      xmin = -93,
      xmax = -68,
      ymin = -39,
      ymax = -14) +
    annotation_custom(
      grob = ggplotGrob(ps5.i3),
      xmin = -75,
      xmax = -55,
      ymin = -12,
      ymax = 3))

# ps5f.s <- ps5f + theme(legend.position = "none")

if (type == "_pa") {
  type <- "pa"
}

# Save SSP5 separately
ggsave(paste0("figures/", sp, "_lgcp_ssp5_", type, ".jpg"), ps5f,
       width = 16, height = 18, units = "cm", quality = 100)



# Save composite ----
final <- pcf + ps1f + ps2f + ps3f + plot_layout(nrow = 1, guides = "collect") &
        theme(legend.position='bottom')

ggsave(paste0("figures/", sp, "_lgcp_comp_", type, ".jpg"), final,
       width = 50, height = 18, units = "cm", quality = 100)
