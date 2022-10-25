library(readxl)
library(tidyverse)

scopus <- read.csv("data/databases/scopus_15082022.csv")
wos <- read_xls("data/databases/wos_15082022.xls")

scopus <- scopus %>% select(Authors, Title, DOI, Abstract)
wos <- wos %>% select(Authors, Title = `Article Title`, DOI, Abstract)

total <- bind_rows(scopus, wos)

total <- total %>%
  mutate(Title = tolower(Title)) %>%
  distinct(Title, DOI, .keep_all = T)

write.csv(total, "scopus_wos_merged.csv", row.names = F)


#WOS = ALL=(atlantic AND "reef* fish*" AND (checklist*  OR  abundance  OR  count   OR  census))
#Scopus = ( TITLE-ABS-KEY ( atlantic  "reef* fish*" ) )  AND  ( abundance  OR  count  OR  census  OR  checklist ) 