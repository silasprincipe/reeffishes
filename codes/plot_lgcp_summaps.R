#### Modelling of coral reef herbivorous invertebrates ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2022 ##

# Plot of the summed map (i.e. aggregating the results of all species)

# Load needed packages ----
library(ggplot2)
library(sf)
library(raster)
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



# Calculate results ----
group <- "carn"

type = "int"

if (group == "herb") {
  species <- c("acch", "scze", "spam")
} else {
  if (group == "carn") {
    species <- c("lujo", "mybo")
  } else {
    species <- c("acch", "scze", "spam", "lujo", "mybo")
  }
}

# Load raster [new version]
scen.rasts <- lapply(species, function(x){
  
  all.rast <- lapply(c("current", paste0("ssp", c(1,2,3,5))), function(z){
    # Open files
    self <- list.files(paste0("results/", x, "/predictions/"),
                       pattern = "q0.5", full.names = T)
    self <- self[grep(paste0(type, "_", z), self)]
    r <- raster(self)
    return(r)
  })
  
  all.rast <- stack(all.rast)
  
  env <- raster("data/env/crop_layers/tempmax.tif")
  
  # Load occurrence data
  if (type == "_pa") {
    pa.pts <- read.csv(paste0("data/", x, "/", x, "_pa.csv"))
    pa.pts <- SpatialPointsDataFrame(pa.pts[, c("decimalLongitude", "decimalLatitude")],
                                     data.frame("presence" = pa.pts[,3]),
                                     proj4string = crs(env))
    
    #pa.pts <- remove.duplicates(pa.pts, zero = 5)
    pa.pts <- dup.cells(pa.pts, rast(env))
    
    # Just to ensure all are falling inside study area
    pa.pts <- pa.pts[!is.na(extract(env, pa.pts)),]
    
    pa.pts <- pa.pts[pa.pts$presence == 1, ]
    
    occ <- pa.pts
  } else {
    po.pts <- read.csv(paste0("data/", x, "/", x, "_final.csv"))
    po.pts <- SpatialPoints(po.pts[, c("decimalLongitude", "decimalLatitude")],
                            proj4string = crs(env))
    
    po.pts <- remove.duplicates(po.pts, zero = 5)
    
    occ <- po.pts
  }
  
  
  ### Threhsolding
  get.thresh <- function(res, pts){
    vals <- raster::extract(res, pts)
    p10 <- ceiling(length(vals) * 0.9)
    thresh <- rev(sort(vals))[p10]
    
    res[res < thresh] <- NA
    
    secvals <- raster::extract(res, pts)
    secvals <- secvals[!is.na(secvals)]
    secthresh <- quantile(secvals, 0.5)
    
    return(c(thresh, secthresh))
  }
  
  threshs <- get.thresh(all.rast[[1]], occ)
  
  classify <- function(scen, threshold){
    scen[scen < threshold] <- 0
    scen[scen >= threshold] <- 1
    return(scen)
  }
  
  for (i in 1:5) {
    all.rast[[i]] <- classify(all.rast[[i]], threshs[1])
  }
  
  names(all.rast) <- paste(x, c("current", paste0("ssp", c(1,2,3,5))), sep = "_")
  
  return(all.rast)
  
})



sum.rasts <- eval(parse(
  text = paste0("scen.rasts[[", 1:length(scen.rasts), "]]", collapse = "+")
))

curr <- sum.rasts[[1]]
ssp1 <- sum.rasts[[2]]
ssp2 <- sum.rasts[[3]]
ssp3 <- sum.rasts[[4]]
ssp5 <- sum.rasts[[5]]

# Convert to data.frame
get.val <- function(x){
  temp <- as(x, "SpatialPixelsDataFrame")
  temp <- as.data.frame(temp)
  colnames(temp) <- c("val", "x", "y")
  temp$val <- as.factor(temp$val)
  return(temp)
}

curr.v <- get.val(curr)
ssp1.v <- get.val(ssp1)
ssp2.v <- get.val(ssp2)
ssp3.v <- get.val(ssp3)
ssp5.v <- get.val(ssp5)


# Prepare theme/ploting stuff ----
# Guide for legend
step.guide <- guide_legend(title = "Number of species",
                             show.limits = TRUE,
                             keyheight = unit(0.2, "in"),
                             keywidth = unit(0.7, "in"),
                             nrow = 1,
                             ticks = T,
                             ticks.colour = "grey20",
                             frame.colour = "grey20",
                             title.position = "top",
                             title.hjust = 0.5,
                             label.position = "bottom",
                             label.hjust = 0.5)

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
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.background = element_rect(fill = "white"),
        legend.spacing.x = unit(0.001, 'cm')
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

sca <- scale_fill_manual(values = c("#E9E7E7",
                                    if (group != "total") {
                                      RColorBrewer::brewer.pal(n = 4, "YlGnBu")[2:4]
                                    }else{
                                      RColorBrewer::brewer.pal(n = 6, "YlGnBu")[2:6]
                                    }),
                         guide = step.guide)


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
    # Add occurrence points
    # geom_point(data = data.frame(occ@coords),
    #            aes(x = decimalLongitude, y = decimalLatitude), size = 1,
    #            alpha = .2, shape = 16)+
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

# ggsave(paste0("figures/summap_lgcp_current.jpg"), pcf,
#        width = 16, height = 18, units = "cm")

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
   # Add occurrence points
   # geom_point(data = data.frame(occ@coords),
   #            aes(x = decimalLongitude, y = decimalLatitude), size = 1,
   #            alpha = .2, shape = 16)+
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

# ggsave(paste0("figures/summap_lgcp_ssp1.jpg"), ps1f.s,
#        width = 16, height = 18, units = "cm")




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
   # Add occurrence points
   # geom_point(data = data.frame(occ@coords),
   #            aes(x = decimalLongitude, y = decimalLatitude), size = 1,
   #            alpha = .2, shape = 16)+
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


# ggsave(paste0("figures/summap_lgcp_ssp2.jpg"), ps2f.s,
#        width = 16, height = 18, units = "cm")




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
   # Add occurrence points
   # geom_point(data = data.frame(occ@coords),
   #            aes(x = decimalLongitude, y = decimalLatitude), size = 1,
   #            alpha = .2, shape = 16)+
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

# ggsave(paste0("figures/summap_lgcp_ssp3.jpg"), ps3f.s,
#        width = 16, height = 18, units = "cm")


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
   # Add occurrence points
   # geom_point(data = data.frame(occ@coords),
   #            aes(x = decimalLongitude, y = decimalLatitude), size = 1,
   #            alpha = .2, shape = 16)+
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
   geom_label(aes(label = "E", x = -35, y = 38),
              size = 10, fontface = "bold", color = "grey30",
              label.size = 0)+
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

ps5f.s <- ps5f + theme(legend.position = "none")

# ggsave(paste0("figures/summap_lgcp_ssp3.jpg"), ps3f.s,
#        width = 16, height = 18, units = "cm")

# Save composite ----
final <- pcf + ps1f + ps2f + ps3f + ps5f + plot_layout(nrow = 1, guides = "collect") &
  theme(legend.position='bottom')

ggsave(paste0("figures/summap_lgcp_composite_", group,"_", ifelse(
  type == "_pa", "pa", type
), ".jpg"), final,
       width = 62, height = 18, units = "cm", quality = 100)




## Q0.5 quantile version ----

# Load raster [new version]
scen.rasts <- lapply(species, function(x){
  
  all.rast <- lapply(c("current", paste0("ssp", c(1,2,3,5))), function(z){
    # Open files
    self <- list.files(paste0("results/", x, "/predictions/"),
                       pattern = "q0.5", full.names = T)
    self <- self[grep(paste0(type, "_", z), self)]
    r <- raster(self)
    return(r)
  })
  
  all.rast <- stack(all.rast)
  
  env <- raster("data/env/crop_layers/tempmax.tif")
  
  # Load occurrence data
  if (type == "_pa") {
    pa.pts <- read.csv(paste0("data/", x, "/", x, "_pa.csv"))
    pa.pts <- SpatialPointsDataFrame(pa.pts[, c("decimalLongitude", "decimalLatitude")],
                                     data.frame("presence" = pa.pts[,3]),
                                     proj4string = crs(env))
    
    #pa.pts <- remove.duplicates(pa.pts, zero = 5)
    pa.pts <- dup.cells(pa.pts, rast(env))
    
    # Just to ensure all are falling inside study area
    pa.pts <- pa.pts[!is.na(extract(env, pa.pts)),]
    
    pa.pts <- pa.pts[pa.pts$presence == 1, ]
    
    occ <- pa.pts
  } else {
    po.pts <- read.csv(paste0("data/", x, "/", x, "_final.csv"))
    po.pts <- SpatialPoints(po.pts[, c("decimalLongitude", "decimalLatitude")],
                            proj4string = crs(env))
    
    po.pts <- remove.duplicates(po.pts, zero = 5)
    
    occ <- po.pts
  }
  
  
  ### Threhsolding
  get.thresh <- function(res, pts){
    vals <- raster::extract(res, pts)
    p10 <- ceiling(length(vals) * 0.9)
    thresh <- rev(sort(vals))[p10]
    
    res[res < thresh] <- NA
    
    secvals <- raster::extract(res, pts)
    secvals <- secvals[!is.na(secvals)]
    secthresh <- quantile(secvals, 0.5)
    
    return(c(thresh, secthresh))
  }
  
  threshs <- get.thresh(all.rast[[1]], occ)
  
  classify <- function(scen, threshold){
    scen[scen < threshold] <- 0
    scen[scen >= threshold] <- 1
    return(scen)
  }
  
  for (i in 1:5) {
    all.rast[[i]] <- classify(all.rast[[i]], threshs[2])
  }
  
  names(all.rast) <- paste(x, c("current", paste0("ssp", c(1,2,3,5))), sep = "_")
  
  return(all.rast)
  
})



sum.rasts <- eval(parse(
  text = paste0("scen.rasts[[", 1:length(scen.rasts), "]]", collapse = "+")
))

curr <- sum.rasts[[1]]
ssp1 <- sum.rasts[[2]]
ssp2 <- sum.rasts[[3]]
ssp3 <- sum.rasts[[4]]
ssp5 <- sum.rasts[[5]]

# Convert to data.frame
get.val <- function(x){
  temp <- as(x, "SpatialPixelsDataFrame")
  temp <- as.data.frame(temp)
  colnames(temp) <- c("val", "x", "y")
  temp$val <- as.factor(temp$val)
  return(temp)
}

curr.v <- get.val(curr)
ssp1.v <- get.val(ssp1)
ssp2.v <- get.val(ssp2)
ssp3.v <- get.val(ssp3)
ssp5.v <- get.val(ssp5)

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
    # Add occurrence points
    # geom_point(data = data.frame(occ@coords),
    #            aes(x = decimalLongitude, y = decimalLatitude), size = 1,
    #            alpha = .2, shape = 16)+
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

# ggsave(paste0("figures/summap_lgcp_current.jpg"), pcf,
#        width = 16, height = 18, units = "cm")

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
   # Add occurrence points
   # geom_point(data = data.frame(occ@coords),
   #            aes(x = decimalLongitude, y = decimalLatitude), size = 1,
   #            alpha = .2, shape = 16)+
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

# ggsave(paste0("figures/summap_lgcp_ssp1.jpg"), ps1f.s,
#        width = 16, height = 18, units = "cm")




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
   # Add occurrence points
   # geom_point(data = data.frame(occ@coords),
   #            aes(x = decimalLongitude, y = decimalLatitude), size = 1,
   #            alpha = .2, shape = 16)+
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


# ggsave(paste0("figures/summap_lgcp_ssp2.jpg"), ps2f.s,
#        width = 16, height = 18, units = "cm")




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
   # Add occurrence points
   # geom_point(data = data.frame(occ@coords),
   #            aes(x = decimalLongitude, y = decimalLatitude), size = 1,
   #            alpha = .2, shape = 16)+
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

# ggsave(paste0("figures/summap_lgcp_ssp3.jpg"), ps3f.s,
#        width = 16, height = 18, units = "cm")


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
   # Add occurrence points
   # geom_point(data = data.frame(occ@coords),
   #            aes(x = decimalLongitude, y = decimalLatitude), size = 1,
   #            alpha = .2, shape = 16)+
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
   geom_label(aes(label = "E", x = -35, y = 38),
              size = 10, fontface = "bold", color = "grey30",
              label.size = 0)+
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

ps5f.s <- ps5f + theme(legend.position = "none")

# ggsave(paste0("figures/summap_lgcp_ssp3.jpg"), ps3f.s,
#        width = 16, height = 18, units = "cm")


# Save composite ----
final <- pcf + ps1f + ps2f + ps3f + ps5f + plot_layout(nrow = 1, guides = "collect") &
  theme(legend.position='bottom')

ggsave(paste0("figures/summap_lgcp_composite_", group,"_", ifelse(
  type == "_pa", "pa", type
), "_q5v.jpg"), final,
       width = 62, height = 18, units = "cm", quality = 100)

### END