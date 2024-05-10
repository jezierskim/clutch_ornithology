#### DATA PREPARATION CODE
#### Manuscript 'Birds that nest exclusively on islands have smaller clutches'. By M.T.Jezierski.
#### Versioned 16th January 2024.

#### Required packages.
library(tidyverse)
library(sf)

##### SET WORKING DIRECTORY.
#setwd('path/to/this/directory')

##### LOAD IN THE CLUTCH DATSET.
clutch <- read_csv('clutch_dataset.csv')

#### PREPARE DISTRIBUTION DATA.
## Need to obtain latitude and breeding area estimates using BirdLife distribution maps.
## To that end the original distributions need to be 1) filtered to breeding range and unionised
## 2) limited to land only; 3) latitude and area need to be calculated.
sf_use_s2(F)
maps <- st_read('BOTW.gpkg') ### BOTW file as BOTW.gdb can be obtained from BirdLife International upon research request. 
maps <- maps %>% filter(binomial %in% clutch$binomial) %>%
  filter(presence == 1) %>% ## Filter to include extant, native, breeding/year round only.
  filter(origin == 1) %>%
  filter(seasonal == 1 | seasonal == 2) %>%
  mutate(migratory = ifelse(seasonal == 2, 'migratory', 'resident')) # extract migrant information based on presence of seasonal breeding range.
migratory <- maps %>% st_drop_geometry() %>% select(binomial, migratory)
migrant_sp <- migratory %>% filter(migratory == 'migratory') %>% distinct() %>% write_csv('migrant_sp.csv')
resident_sp <- migratory %>% filter(migratory == 'resident') %>% distinct() %>% anti_join(migrant_sp, by = 'binomial') %>% write_csv('resident_sp.csv')
maps <- maps %>% group_by(binomial) %>% summarise(geom = st_union(geom))
land <- st_read('ne_50m_land/ne_50m_land.shp') %>% ## ne_50m_land.shp is available from NaturalEarth data 
  summarise(geometry = st_union(geometry))
maps <- st_intersection(maps, land) ## intersect with land at 50m resolution 

maps %>% group_by(binomial) %>% ## calculate 'latitude' as range centroid; calculate total area of the range used. 
  mutate(latitude = as.data.frame(st_coordinates(st_centroid(geom)))[,2]) %>%
  mutate(area = st_area(geom)) %>% st_drop_geometry() %>% drop_na() %>%
  right_join(clutch, by = 'binomial') %>% write_csv('clutch_dataset_geo.csv')

#### ADD LIFE HISTORY INFORMATION
# Merge the geographic data dataset with information and developmental mode, habitat and migratory status.
# Also calculate the geometric mean to analyse. 
clutch <- read_csv('clutch_dataset_geo.csv')
migrant_sp <- read_csv('migrant_sp.csv') # species that are migratory.
dev <- read_csv('dev_sp.csv') # bonus datasets with bird families categorised according to development and being land/seabirds.
geom_mean <- clutch %>% dplyr::select(clumax, clumin)
geom_mean <- exp(rowMeans(log(geom_mean))) # calculate geometric mean of clutch size.
mass_mean <- clutch %>% dplyr::select(bodmax, bodmin)
mass_mean <- exp(rowMeans(log(mass_mean))) # calculate geometric mean of mass.

clutch <- clutch %>% left_join(migrant_sp, by = 'binomial') %>% left_join(dev, by = 'family') %>% # add all datasets togegether to form the analysis file.
  mutate(clu_geom = geom_mean, mass_geom = mass_mean) %>% write_csv('analysis_dataset.csv')
