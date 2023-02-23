#####
#### PREPARING THE DATA FOR THE ANALYSIS IN CLUTCH SIZE PAPER
#####
library(tidyverse)
library(sf)

##### THE ENTIRE CLUTCH DATASET
clutch <- read_csv('/Users/user/Library/Mobile Documents/com~apple~CloudDocs/Desktop/island clutch size/clutch_dataset.csv') %>% 
  select(order, family, ebird_name, treelabels, binomial, clumin, clumax, bodmax, is) %>%
  drop_na()
write_csv(clutch, 'NEW_clutch_dataset.csv')

#### DISTRIBUTION DATA - FILTER, OBTAIN CENTROID.
sf_use_s2(F)
maps <- st_read('~/Desktop/BOTW.gpkg')
maps <- maps %>% filter(binomial %in% clutch$binomial) %>%
  filter(presence == 1) %>% 
  filter(origin == 1) %>%
  filter(seasonal == 1 | seasonal == 2) %>%
  group_by(binomial) %>% summarise(geom = st_union(geom))
st_write(maps, '~/Desktop/FILTERED.gpkg')

maps %>% group_by(binomial) %>% 
  mutate(latitude = as.data.frame(st_coordinates(st_centroid(geom)))[,2]) %>%
  mutate(area = st_area(geom)) %>% st_drop_geometry() %>% drop_na() %>%
  right_join(clutch, by = binomial) %>% write_csv('~/Desktop/clutch_dataset_geo.csv')
