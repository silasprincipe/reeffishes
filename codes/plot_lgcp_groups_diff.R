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

# Set themes/color scales ----
wlt <- theme_classic()+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(colour = "grey60", size = 14),
        axis.title.x.top = element_text(vjust = 2, size = 16),
        axis.title.y.left = element_text(size = 16),
        panel.background = element_rect(fill = "#e8e8e8"),
        panel.border = element_rect(fill = NA, color = "grey60"),
        panel.grid.major = element_line(linetype = 'dashed', 
                                        colour = "grey70",
                                        linewidth = .1),
        legend.position="bottom",
        legend.title.align=0.5,
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.background = element_rect(fill = "white"),
        strip.background = element_blank(),
        strip.text.y = element_text(face = "italic"),
        strip.text = element_text(size = 14))


# Define groups ----
groups <- list(
  carnivores = list(
    acronym = c("lujo", "mybo"),
    full_name = c("Lutjanus jocu", "Mycteroperca bonaci")
  ),
  invertivores = list(
    acronym = c("boru", "hara"),
    full_name = c("Bodianus rufus", "Halichoeres radiatus")
  ),
  herbivores = list(
    acronym = c("acch", "scze", "spam"),
    full_name = c("Acanthurus chirurgus", "Scarus zelindae", "Sparisoma amplum")
  )
)


# Define metric(s) to plot ----
metrics <- c("q0.5", "mean")

types <- c("pa_clamp", "pa_noclamp")

# Get the plots ----
for (j in 1:length(types)) {
  
  type <- types[j]
  
  for (m in metrics) {
    
    if (m == "sd") {
      # Set color scale
      step.guide <- guide_colorbar(title = "SD(Difference in the probability of occurrence)",
                                   show.limits = TRUE,
                                   barheight = unit(0.12, "in"),
                                   barwidth = unit(3.5, "in"),
                                   ticks = T,
                                   ticks.colour = "grey20",
                                   frame.colour = "grey20",
                                   title.position = "top")
      
      # Scale of values
      sca <- scale_fill_distiller(palette = "YlGnBu",
                                  guide = step.guide,
                                  limits = c(-1,1),
                                  labels = c("0", "", "0.5", "", "1"))
    } else {
      # Set color scale
      step.guide <- guide_colorbar(title = "Difference (probability of occurrence)",
                                   show.limits = TRUE,
                                   barheight = unit(0.12, "in"),
                                   barwidth = unit(3.5, "in"),
                                   ticks = T,
                                   ticks.colour = "grey20",
                                   frame.colour = "grey20",
                                   title.position = "top")
      
      # Scale of values
      cs <- RColorBrewer::brewer.pal(11, "PuOr")
      #cs[6] <- "#bababa"
      
      sca <- scale_fill_gradientn(
        colours = cs,
        limits = c(-1,1),
        guide = step.guide,
        labels = c("-1", "", "0", "", "+1"),
        na.value = "#2b0700")
    }
    
    for (g in 1:length(groups)) {
      species <- groups[[g]]$acronym
      group.name <- names(groups)[g]
      
      for (i in 1:length(species)) {
        sp <- species[i]
        
        preds <- list.files(paste0("results/", sp, "/predictions"),
                            pattern = type, full.names = T)
        
        preds <- preds[grep(m, preds)]
        
        preds <- rast(preds)
        
        names(preds) <- c("Current", paste0("SSP", c(1,2,3,5)))
        
        delta.preds <- preds[[2:5]] - preds[[1]]
        
        delta.preds <- as.data.frame(delta.preds, xy = T)
        
        delta.preds <- pivot_longer(delta.preds,
                              cols = 3:6, names_to = "scenario", values_to = "value")
        
        delta.preds$species <- sp
        
        if (i == 1) {
          pred.data <- delta.preds
        } else {
          pred.data <- rbind(pred.data, delta.preds)
        }
      }
      
      pred.data$species <- factor(pred.data$species,
                                  labels = groups[[g]]$full_name)
      
      (p <- ggplot() +
          geom_sf(data = base,
                  color = c("#CFCFCF"),
                  size = 0.6,
                  fill = "white") +
          coord_sf(xlim = c(-99, -29),
                   ylim = c(-42.5, 42.5),
                   datum = st_crs(base),
                   expand = F,
                   label_axes = list(
                     top = "",
                     left = "N",
                     bottom = "E"
                   )) +
          scale_x_discrete(position = "top", breaks = c(-85, -60, -35)) +
          xlab("") + ylab("") +
          geom_raster(data = pred.data, aes(x = x, y = y, fill = value)) +
          sca +
          wlt +
          facet_grid(species ~ scenario))
      
      # Define plot height and save
      ph <- 9 * length(species)
      
      ggsave(plot = p,
             filename = paste0("figures/", group.name, "_", gsub("\\.", "-", m), "_", type, "_diff.jpg"),
             width = 30, height = ph, units = "cm", quality = 100)
    }
    
  }
}

# END