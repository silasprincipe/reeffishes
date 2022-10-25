#### Modelling of reef fishes ####
## Silas C. Principe - silascprincipe@gmail.com - 2021/2022

### Species data conversion to 1 per cell ###

# Load packages ----
library(raster)

# Load base layer ----
base <- raster("data/env/crop_layers/BO21_tempmean_ss.tif")

# Convert species points to one per cell ----

# List species
sp.codes <- c("acch", "scze", "spam", "lujo", "mybo")

# Create a function to convert to 1 point per cell
to.cell <- function(species){
        
        sp <- read.csv(paste0("data/", species, "/", species, "_final.csv"))
        
        #Remove points out of the area/bathymetry
        out.p <- raster::extract(base, sp[,2:3])
        sp <- sp[!is.na(out.p),]
        
        #Create a clean raster
        r <- base
        r[] <- NA
        
        #Put 1 in presence cells
        r[cellFromXY(r, sp[,2:3])] <- 1
        
        #Convert raster to dataframe
        data <- rasterToPoints(r)
        
        #Change column names
        colnames(data) <- c("decimalLongitude", "decimalLatitude",
                            species)
        
        data
}

# Apply to datasets
cell.data <- lapply(sp.codes, to.cell)

# Plot to see
# plot(base, col = 'grey', legend = F)
# points(cell.data[[5]][,1:2], pch = 20, cex = 0.5, col = "blue")

# Save datasets ----
for (i in 1:length(cell.data)) {
        
        write.csv(cell.data[[i]], paste0("data/",
                                          sp.codes[i],
                                          "/", sp.codes[i],
                                          "_cell.csv"),
                  row.names = F)
}

#END of code