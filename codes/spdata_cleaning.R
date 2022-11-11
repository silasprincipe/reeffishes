#### Modelling of reef fishes ####
## Silas C. Principe - silascprincipe@gmail.com - 2021/2022

### Reef fishes data download and cleaning ###
# We also aggregate the data from bibliographic sources

# Download in 2021-10-19

# Load packages ----
library(tidyverse)
library(robis)
library(rgbif)
library(obistools)
library(fuzzySim)
library(raster)



# Download data | OBIS ----

# List species that will be downloaded
sp.list <- c("Acanthurus chirurgus",
             "Scarus zelindae",
             "Sparisoma amplum",
             "Lutjanus jocu",
             "Mycteroperca bonaci")

# Create a unique code for each species
sp.codes <- tolower(spCodes(sp.list, nchar.gen = 2, nchar.sp = 2))

# Add synonims for OBIS search
syn <- c("Acanthurus chirurgicus", # 6:8 Acanthurus chirurgus
         "Acanthurus phlebotomus",
         "Chaetodon chirurgus",
         "Scarus amplus", # 9 Sparisoma amplum
         "Anthias jocu", # 10:11 Lutjanus jocu
         "Mesoprion litura", 
         "Serranus arara", # 12:18 Mycteroperca bonaci
         "Serranus bonaci",
         "Serranus brunneus",
         "Serranus cyclopomatus",
         "Serranus decimalis",
         "Serranus latepictus",
         "Trisotropis aguaji") 

# Download from OBIS
obis.data <- lapply(c(sp.list, syn), occurrence)

# Bind synonims results
obis.data[[1]] <- bind_rows(obis.data[[1]], obis.data[6:8]) # A. chirurgus
obis.data[[3]] <- bind_rows(obis.data[[3]], obis.data[[9]]) # S. amplum
obis.data[[4]] <- bind_rows(obis.data[[4]], obis.data[10:11]) # L. jocu
obis.data[[5]] <- bind_rows(obis.data[[5]], obis.data[12:18]) # M. bonaci

obis.data <- obis.data[1:5]

lapply(obis.data, nrow) # Verify number of records

# Add codes of species to the list to easy handling
names(obis.data) <- sp.codes




# Download data | GBIF ----

# Download from GBIF

gbif.data <- list()

for (i in 1:length(sp.list)) {
        
        cat("Downloading data for", sp.list[i], "\n")      
        
        gbif <- occ_data(scientificName = sp.list[i], hasCoordinate = T,
                         limit = 50000)
        gbif.data[[i]] <- as.data.frame(gbif$data)
        rm(gbif)
}

lapply(gbif.data, nrow) # Verify number of records

# Add codes of species to the list to easy handling
names(gbif.data) <- sp.codes




# Data cleaning | OBIS ----

# Define year, record type and area of work
years <- as.character(c(1950:2022))
record <- c("HumanObservation", "Occurrence", "PreservedSpecimen")
area <- c(-42.5, 42.5, -99, -29) #Lat 1, Lat 2, Long 1, Long 2

# Remove points that don't meet requirements
obis.data <- lapply(obis.data, function(x){
        x <- x %>%
                dplyr::select(-flags) %>%
                filter(year %in% years) %>%
                filter(basisOfRecord %in% record) %>%
                filter(decimalLatitude >= area[1] & 
                               decimalLatitude <= area[2]) %>%
                filter(decimalLongitude >= area[3] & 
                               decimalLongitude <= area[4])
                
})

lapply(obis.data, nrow)


# Acanthurus chirurgus
unique(obis.data$acch$institutionCode)
#We use this piece of code to get the institutions codes

obis.data$acch <- obis.data$acch %>% 
        filter(!is.na(institutionCode)) %>%
        filter(institutionCode != "Diveboard")
                                

# Scarus zelindae
unique(obis.data$scze$institutionCode)

obis.data$scze <- obis.data$scze %>% 
        filter(!is.na(institutionCode))


# Sparisoma amplum
unique(obis.data$spam$institutionCode)

obis.data$spam <- obis.data$spam %>% 
        filter(!is.na(institutionCode)) %>%
        filter(institutionCode != "Diveboard") #Remove Diveboard data


# Lutjanus jocu
unique(obis.data$lujo$institutionCode)

obis.data$lujo <- obis.data$lujo %>% 
        filter(!is.na(institutionCode)) %>%
        filter(institutionCode != "Diveboard") #Remove Diveboard data


# Mycteroperca bonaci
unique(obis.data$mybo$institutionCode)

obis.data$mybo <- obis.data$mybo %>% 
        filter(!is.na(institutionCode)) %>%
        filter(institutionCode != "Diveboard") #Remove Diveboard data


# Optional: plot each data to verify
plot_map_leaflet(obis.data$mybo)


# Ensure all points are on water / inside study area
# First we load a base layer
base <- raster("data/env/crop_layers/tempmax.tif")

# Then we create a function that will see if points fall inside the area
# and return only those
get.onarea <- function(data){
  ov <- raster::extract(base, data[,c("decimalLongitude", "decimalLatitude")])
  data[!is.na(ov), ]
}

obis.data <- lapply(obis.data, get.onarea)

# Final check
lapply(obis.data, nrow)

# Optional: plot each data to verify
plot(base)
points(obis.data$mybo[,c("decimalLongitude", "decimalLatitude")],
       pch = 16, cex = 0.5)



# Data cleaning | GBIF ----

# Define year, record type and area of work
years <- as.character(c(1950:2022))
record <- c("HUMAN_OBSERVATION", 
            "PRESERVED_SPECIMEN", 
            "MACHINE_OBSERVATION",
            "LIVING_SPECIMEN")
area <- c(-42.5, 42.5, -99, -29) #Lat 1, Lat 2, Long 1, Long 2

# Remove points that don't meet requirements
gbif.data <- lapply(gbif.data, function(x){
        x <- x %>%
                filter(year %in% years) %>%
                filter(basisOfRecord %in% record) %>%
                filter(decimalLatitude >= area[1] & 
                               decimalLatitude <= area[2]) %>%
                filter(decimalLongitude >= area[3] & 
                               decimalLongitude <= area[4]) %>%
                dplyr::select(-networkKeys)
})

lapply(gbif.data, nrow)


# Check each dataset for source (this one is done individually)

# A. chirurgus
unique(gbif.data$acch$institutionCode)

gbif.data$acch <- gbif.data$acch %>% 
        filter(!is.na(institutionCode)) %>%
        filter(!institutionCode %in% c("iNaturalist", "Diveboard",
                                       "naturgucker", "NO DISPONIBLE"))

# S. zelindae
unique(gbif.data$scze$institutionCode)

gbif.data$scze <- gbif.data$scze %>% 
        filter(!is.na(institutionCode)) %>%
        filter(!institutionCode == "iNaturalist")


# S. amplum
unique(gbif.data$spam$institutionCode)

gbif.data$spam <- gbif.data$spam %>% 
        filter(!is.na(institutionCode)) %>%
        filter(!institutionCode %in% c("iNaturalist", "Diveboard"))


# L. jocu
unique(gbif.data$lujo$institutionCode)

gbif.data$lujo <- gbif.data$lujo %>% 
        filter(!is.na(institutionCode)) %>%
        filter(!institutionCode %in% c("iNaturalist", "Diveboard",
                                       "NO DISPONIBLE", "naturgucker"))

# M. bonaci
unique(gbif.data$mybo$institutionCode)

gbif.data$mybo <- gbif.data$mybo %>% 
        filter(!is.na(institutionCode)) %>%
        filter(!institutionCode %in% c("iNaturalist", "Diveboard",
                                       "NO DISPONIBLE", "SEAK"))



# Optional: plot each data to verify
plot_map_leaflet(gbif.data$lujo)


# Ensure all points are on water / inside study area
gbif.data <- lapply(gbif.data, get.onarea)

# Final check
lapply(gbif.data, nrow)

# Optional: plot each data to verify
plot(base)
points(gbif.data$lujo[,c("decimalLongitude", "decimalLatitude")],
       pch = 16, cex = 0.5)

## End of data cleaning 




# Save full datasets ----
for (i in 1:length(obis.data)) {
        
        if (dir.exists(paste0("data/", sp.codes[i])) == F) {
                dir.create(paste0("data/", sp.codes[i]), recursive = T)
        }
        
        
        write.csv(obis.data[[i]], paste0("data/",
                                         sp.codes[i],
                                         "/", sp.codes[i],
                                         "_obis_fulldata.csv"),
                  row.names = F)
        
        write.csv(gbif.data[[i]], paste0("data/",
                                         sp.codes[i],
                                         "/", sp.codes[i],
                                         "_gbif_fulldata.csv"),
                  row.names = F)
}





# Bind data from bibliographic sources and databases ----

# Open data from bibliographic sources
biblio.data <- lapply(paste0("data/", sp.codes, "/", sp.codes, "_biblio.csv"), read.csv)

# Add a column with the species name (acronym)
for (i in 1:length(sp.codes)) {
        biblio.data[[i]]$scientificName <- sp.codes[i]
}

# Select relevant collumns
biblio.data <- lapply(biblio.data, function(x){
        x <- x %>% dplyr::select(scientificName, decimalLongitude, decimalLatitude)
})

# Select relevant collumns of OBIS and GBIF data
obis.data <- lapply(obis.data, function(x){
        x <- x %>% dplyr::select(scientificName, decimalLongitude, decimalLatitude)
})

gbif.data <- lapply(gbif.data, function(x){
        x <- x %>% dplyr::select(scientificName, decimalLongitude, decimalLatitude)
})

# Bind datasets
final.data <- list() # Create an empty list to bind datasets

#Change the names of OBIS and GBIF data to the acronyms and bind all data
for (i in 1:length(sp.codes)) {
        obis.data[[i]]$scientificName <- sp.codes[i]
        gbif.data[[i]]$scientificName <- sp.codes[i]
        final.data[[i]] <- rbind(obis.data[[i]],
                                 gbif.data[[i]],
                                 biblio.data[[i]])
}


# Remove points on land (we do this again because we added the bibliographic
# data and its better to have a second check!).
final.data <- lapply(final.data, get.onarea)

lapply(final.data, nrow)


# Save final datasets ----
for (i in 1:length(final.data)) {
        
        write.csv(final.data[[i]], paste0("data/",
                                         sp.codes[i],
                                         "/", sp.codes[i],
                                         "_final.csv"),
                  row.names = F)
}

# See points
lapply(final.data, robis::map_ggplot)



# Abundance datasets preparing ----

# Open GASPAR dataset
gaspar <- read.csv("data/databases/gaspar_dataset.csv")

# Select only species from the study and also standardize abundance per 100m^2
gaspar <- gaspar %>%
  filter(Species %in% sp.list) %>%
  mutate(Abundance = (100 * Abundance)/Search_Area_m2) %>%
  # Then we group by Lat and Long, species to get the median
  # abundance over the sampled years.
  group_by(Latitude, Longitude, Species) %>%
  summarise(m_abundance = median(Abundance))

# Open RLS dataset
rls <- lapply(list.files("data/databases", pattern = "rls", full.names = T),
              read.csv)
rls <- bind_rows(rls)

# Select only species from the study and also standardize abundance per 100m^2
rls <- rls %>%
  filter(Taxon %in% sp.list) %>%
  mutate(Total = (100 * Total)/250) %>%
  # Then we group by Lat and Long, species to get the median
  # abundance over the sampled years.
  group_by(SiteLat, SiteLong, Taxon) %>%
  summarise(m_abundance = median(Total))


# Standardize colnames
colnames(rls) <- colnames(gaspar) <- c("decimalLatitude", "decimalLongitude",
                                       "scientificName", "abundance")
# Put on x/y format
gaspar <- gaspar[,c(2,1,3,4)]
rls <- rls[,c(2,1,3,4)]

# Create a key to know the dataset origin and bind
abundance <- data.frame(
  rbind(
    cbind(gaspar, dataset = "gaspar"),
    cbind(rls, dataset = "rls")
  )
)

# Round values to the next integer to enable the use of Poisson (avoiding xpoisson)
# TO VERIFY TO SEE IF ITS NOT POSSIBLE TO USE GAUSSIAN
abundance$abundance <- ceiling(abundance$abundance)

# Plot to verify
ggplot(abundance) +
  geom_point(aes(x = decimalLongitude, y = decimalLatitude,
                 color = abundance)) + coord_equal() +
  facet_wrap(~scientificName)

# Save by species
for (i in 1:length(sp.list)) {
  sabund <- abundance[abundance$scientificName == sp.list[i], ]
  sabund$scientificName <- sp.codes[i]
  
  write.csv(sabund, paste0("data/", sp.codes[i], "/",
                           sp.codes[i], "_abundance.csv"), row.names = F)
  rm(sabund)
}



# Presence/absence datasets ----
# Open GASPAR dataset
gaspar <- read.csv("data/databases/gaspar_dataset.csv")

# Open RLS dataset
rls <- lapply(list.files("data/databases", pattern = "rls", full.names = T),
              read.csv)
rls <- bind_rows(rls)

# Open literature dataset
lit <- read.csv("data/databases/pa_dataset_literature.csv")

# Open Quimbayo et al. dataset
quim <- read.csv("data/databases/quimbayo_etal_2021.csv")

# Select only relevant collumns
gaspar <- gaspar %>%
  dplyr::select(Longitude, Latitude, Species, Abundance)

rls <- rls %>%
  dplyr::select(SiteLong, SiteLat, Taxon, Total)

quim <- quim %>%
  dplyr::select(lon, lat, species, abun)

# We start working with the RLS and Gaspar datasets
# Standardize colnames
colnames(rls) <- colnames(gaspar)

# Bind data
padata <- bind_rows(rls, gaspar)

# Select area of study
area <- c(-42.5, 42.5, -99, -29)

# Filter and get data for presence and absence
padata <- padata %>%
  mutate(
    acch = ifelse(Species == "Acanthurus chirurgus", 1, 0),
    scze = ifelse(Species == "Scarus zelindae", 1, 0),
    spam = ifelse(Species == "Sparisoma amplum", 1, 0),
    lujo = ifelse(Species == "Lutjanus jocu", 1, 0),
    mybo = ifelse(Species == "Mycteroperca bonaci", 1, 0)
  ) %>%
  group_by(Longitude, Latitude) %>%
  summarise(acch = sum(acch),
            scze = sum(scze),
            spam = sum(spam),
            lujo = sum(lujo),
            mybo = sum(mybo)) %>%
  filter(Latitude >= area[1] & Latitude <= area[2]) %>%
  filter(Longitude >= area[3] & Longitude <= area[4])

# Assign only presence/absence
padata[,3:7] <- apply(padata[,3:7], 2, function(x){ifelse(x > 0, 1, 0)})

# Now we work with the Quimbayo dataset
quim <- quim %>%
  mutate(
    acch = ifelse(species == "acanthurus_chirurgus", 1, 0),
    scze = ifelse(species == "scarus_zelindae", 1, 0),
    spam = ifelse(species == "sparisoma_amplum", 1, 0),
    lujo = ifelse(species == "lutjanus_jocu", 1, 0),
    mybo = ifelse(species == "mycteroperca_bonaci", 1, 0)
  ) %>%
  group_by(lon, lat) %>%
  summarise(acch = sum(acch),
            scze = sum(scze),
            spam = sum(spam),
            lujo = sum(lujo),
            mybo = sum(mybo)) %>%
  filter(lat >= area[1] & lat <= area[2]) %>%
  filter(lon >= area[3] & lon <= area[4])

# Assign only presence/absence
quim[,3:7] <- apply(quim[,3:7], 2, function(x){ifelse(x > 0, 1, 0)})

# Now we change colnames to standardize
colnames(quim) <- colnames(padata)

# Include dataset with data from the bibliographic sources
lit <- lit %>%
  dplyr::select(decimalLongitude, decimalLatitude,
                acch, scze, spam, lujo, mybo) %>%
  filter(decimalLatitude >= area[1] & decimalLatitude <= area[2]) %>%
  filter(decimalLongitude >= area[3] & decimalLongitude <= area[4])


colnames(lit) <- colnames(padata)

# Merge datasets
padata <- bind_rows(padata, quim, lit)

# If wanted, we can plot the data
plotdat <- padata %>% pivot_longer(3:7, names_to = "species", values_to = "pa")

ggplot(plotdat) +
  geom_point(aes(x = Longitude, y = Latitude, color = pa, size = pa))+
  #geom_point(data = abundance[abundance$scientificName == "Acanthurus chirurgus",], 
  #           aes(x = decimalLongitude, y = decimalLatitude), size = 2, alpha = .1, color = "white")+
  coord_equal() +
  facet_wrap(~species)


# Save by species
for (i in 1:length(sp.codes)) {
  pasel <- padata[ , c("Longitude", "Latitude", sp.codes[i])]
  nr <- nrow(pasel)
  # We remove lines with NA (data was unavailable or with a kind of error
  # and we prefered to assign NA instead of assuming as absence)
  pasel <- pasel[!is.na(pasel[,3]),]
  cat(nr-nrow(pasel), "NAs removed from", sp.codes[i], "\n")
  
  colnames(pasel)[1:2] <- c("decimalLongitude", "decimalLatitude")
  
  write.csv(pasel, paste0("data/", sp.codes[i], "/",
                           sp.codes[i], "_pa.csv"), row.names = F)
  rm(pasel)
}


#END of code