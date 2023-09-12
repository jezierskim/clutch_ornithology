#####
#### PREPARING THE DATA FOR THE ANALYSIS IN CLUTCH SIZE PAPER
#####
library(tidyverse)
library(sf)

##### LOAD IN THE CLUTCH DATSET
clutch <- read_csv('~/Desktop/clutch/clutch_dataset.csv') %>% 
  dplyr::select(order, family, ebird_name, treelabels, binomial, clumin, clumax, bodmax, is) %>%
  drop_na()
write_csv(clutch, '~/Desktop/clutch/NEW_clutch_dataset.csv')

#### PREPARE DISTRIBUTION DATA
## Need to obtain latitude and breeding area estimates using BirdLife distribution maps.
## To that end the original distributions need to be 1) filtered to breeding range and unionised
## 2) limited to land only; 3) latitude and area need to be calculated.
sf_use_s2(F)
maps <- st_read('~/Desktop/BOTW.gpkg')
maps <- maps %>% filter(binomial %in% clutch$binomial) %>%
  filter(presence == 1) %>% ## Filter to include extant, native, breeding/year round only.
  filter(origin == 1) %>%
  filter(seasonal == 1 | seasonal == 2) %>%
  group_by(binomial) %>% summarise(geom = st_union(geom))
#st_write(maps, '~/Desktop/FILTERED.gpkg') 
maps <- st_read('~/Desktop/clutch/FILTERED.gpkg')
land <- st_read('~/Desktop/clutch/ne_50m_land/ne_50m_land.shp') %>%
  summarise(geometry = st_union(geometry))
maps <- st_intersection(maps, land) ## intersect with land at 50m resolution 

maps %>% group_by(binomial) %>% 
  mutate(latitude = as.data.frame(st_coordinates(st_centroid(geom)))[,2]) %>%
  mutate(area = st_area(geom)) %>% st_drop_geometry() %>% drop_na() %>%
  right_join(clutch, by = 'binomial') %>% write_csv('~/Desktop/clutch/clutch_dataset_geo.csv')
