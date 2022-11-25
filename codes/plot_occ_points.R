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



# Load species data ----
sp.list <- c("acch", "scze", "spam", "lujo", "mybo")

# Load a layer just as sample
env <- raster("data/env/crop_layers/tempmax.tif")

# Occurrence datasets
sp.data <- lapply(sp.list, function(sp){
  po.pts <- read.csv(paste0("data/", sp, "/", sp, "_final.csv"))
  po.pts <- SpatialPoints(po.pts[, c("decimalLongitude", "decimalLatitude")],
                          proj4string = crs(env))
  
  po.pts <- remove.duplicates(po.pts, zero = 5)
  
  pa.pts <- read.csv(paste0("data/", sp, "/", sp, "_pa.csv"))
  pa.pts <- SpatialPointsDataFrame(pa.pts[, c("decimalLongitude", "decimalLatitude")],
                                   data.frame("presence" = pa.pts[,3]),
                                   proj4string = crs(env))
  
  #pa.pts <- remove.duplicates(pa.pts, zero = 5)
  pa.pts <- dup.cells(pa.pts, rast(env))
  
  # Just to ensure all are falling inside study area
  pa.pts <- pa.pts[!is.na(extract(env, pa.pts)),]
  
  fdat <- data.frame(
    rbind(
      cbind(coordinates(pa.pts),
            type = ifelse(pa.pts$presence == 1,
                          "Presence", "Absence"),
            symb = 2),
      cbind(coordinates(po.pts),
            type = "Presence-only",
            symb = 1)
    )
  )
  
  fdat$x <- as.numeric(fdat$x)
  fdat$y <- as.numeric(fdat$y)
  fdat$symb <- as.numeric(fdat$symb)
  
  return(fdat)
  
})

names(sp.data) <- sp.list

sp.names <- c(
  "Acanthurus chirurgus",
  "Scarus zelindae",
  "Sparisoma amplum",
  "Lutjanus jocu",
  "Mycteroperca bonaci"
)

for (i in 1:length(sp.data)) {
  sp.data[[i]]$species <- sp.names[i]
  if (i == 1) {
    sp.data.f <- sp.data[[i]]
  }else{
    sp.data.f <- rbind(sp.data.f,
                       sp.data[[i]])
  }
}



# Load environmental layers ----
bath <- raster("data/env/bath_layers/bath_2_300.tif")
bath <- bath * -1
bath <- mask(bath, env)

# Convert to data.frame
get.val <- function(x){
        temp <- as(x, "SpatialPixelsDataFrame")
        temp <- as.data.frame(temp)
        colnames(temp) <- c("val", "x", "y")
        return(temp)
}

bath.v <- get.val(bath)



# Prepare theme/ploting stuff ----
# Guide for legend
step.guide <- guide_colorbar(title = "Depth",
                             show.limits = TRUE,
                             barheight = unit(0.12, "in"),
                             barwidth = unit(3.5, "in"),
                             ticks = T,
                             ticks.colour = "grey20",
                             frame.colour = "grey20",
                             title.position = "top")

# Scale of values
sca <- scale_fill_gradientn(
  colours = rev(RColorBrewer::brewer.pal(9, "YlGnBu")),
  guide = step.guide)

dt.guide <- guide_legend(title.position = "top",
                         override.aes = list(size=3))


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
ggplot()+
  # Base maps
  geom_sf(data = base.pac,
          color = NA,
          fill = c("#D9D9D9"),
          alpha = 0.5) +
  geom_sf(data = base,
          color = c("#CFCFCF"),
          size = 0.6,
          fill = "white") +
  # # Get the raster
  geom_raster(data = bath.v, aes(x = x, y = y), fill = "gray70") +
  # sca + # Add color scale
  # ggnewscale::new_scale_color()+
  # Add occurrence points
  geom_point(data = sp.data.f,
             aes(x = x, y = y, shape = type, color = type, size = symb),
             alpha = 0.8)+
  scale_color_manual(name = "Data type",
                     values = c("#FF7300", "#A515E8", "gray10"),#"#28C7A7"
                     guide = dt.guide)+
  scale_shape_manual(name = "Data type",
                     values = c(17, 17, 1),
                     guide = dt.guide)+
  scale_size(guide = "none", range = c(1.5, 3))+
  # # Establish area
  coord_sf(xlim = c(-99, -27),
           ylim = c(-42.5, 42.5),
           datum = st_crs(base),
           expand = F) +
  # # Add title
  # geom_label(aes(label = "A", x = -35, y = 38),
  #            size = 10, fontface = "bold", color = "grey30",
  #            label.size = 0)+
  # Add grid
  scale_x_discrete(position = "top", breaks = c(-95, -75, -55, -35)) +
  # Remove axis labels and add theme
  xlab("") + ylab("") + wlt + 
  theme(
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0, size = 14, face = "italic")
  ) +
  facet_wrap(~species, nrow = 1)


ggsave("figures/occurrence_points.jpg",
       width = 45, height = 15, units = "cm", quality = 100)
