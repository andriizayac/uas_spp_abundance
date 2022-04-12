library(sf)
library(tidyverse)

# === helper fn from Andrie: https://stackoverflow.com/questions/7963898/extracting-the-last-n-characters-from-a-string-in-r
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
# ===

f <- read.csv("D:/UAV_analysis/demographic_data/CRP/regen3/crp_regen3_2021_recr.csv")

fpts <- f %>% 
  filter(plant == 0) %>% 
  mutate(lid = substrRight(species, 3)) %>% 
  st_as_sf(coords = c("long", "lat")) %>% 
  st_set_crs(4326) 

fpts %>% 
  st_transform(32612) %>% 
  select(lid, elev) %>% 
  group_by(lid) %>% 
  summarize(m = mean(elev)) %>% 
  st_cast('LINESTRING') %>% 
  st_buffer(1, endCapStyle = 'FLAT') %>% 
  st_cast('POLYGON') %>% 
  ungroup() %>% 
  select(lid) %>% 
  mutate(id = paste0("L", lid)) %>% 
  st_transform(4326) %>%
  st_write("D:/Andrii/abundance_sodaBurn_2021/uas_spp_abundance/data/transect_shapes/crp_regen3_transects.shp", delete_layer = TRUE)



  