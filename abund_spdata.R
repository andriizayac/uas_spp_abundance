pkgs <- c("lidR", "raster","sf", "viridisLite", "stars", "tidyverse", "ggplot2")
sapply(pkgs, require, character.only = TRUE)

# this script combines the following datasets for the abundance anaylys
# 1. Field data
# 2. Additional hand-digitized large plants
# 3. Segmented output

# generates outputs
# 1. Plant data set
# 2. Cell data set

# Variables Plants: id, type, transect, area, height, lat, long, y, ht_true, z, cellid
# Variables Cells: elev,	cid,	lat,	long,	burn,	cht,	n,	dist




# === 0. bring transect polygons, clip shape (subset), CHM, and DTM
mask <- st_read("data/mask/crp_regen1_subset_clip.shp") %>% st_transform(32612)

tr <- st_read("data/transect_shapes/crp_regen1_transects.shp") %>% 
  st_transform(32612) %>% 
  select(-id) %>% 
  mutate(var = lengths(st_intersects(., mask))) %>% filter(var == 1) %>% select(-var)

chm <- raster("./../crp_regen1_20210804_CHM_WGS84UTM12.tif") %>% 
  crop(extent(mask))
dtm <- raster("./../crp_regen1_20210804_DTM_WGS84UTM12.tif") %>% 
 crop(chm@extent)

# === 1. Field data: points and polygons
fpts <- read.csv("data/digitized/crp_regen1_2021_recr.csv") %>% 
  filter(!is.na(long)) %>% 
  mutate(x = long, y = lat, z = elev,
         across(c(height, width1, width2), ~.x/100)) %>%
  st_as_sf(coords = c("x", "y"), crs = 4326) %>%
  filter(plant == 1, live == 1, species == "tridentata") 

fpol <- st_read("data/digitized/crp_regen1_2021_digitized.shp") %>% 
  st_transform(32612) %>% 
  mutate(point = ifelse(nchar(id) < 2, paste0("DOWN00", id), paste0("DOWN0", id))) %>% 
  select(-id, -notes) %>% 
  filter(point %in% fpts$point) %>% 
  left_join(as.data.frame(fpts) %>% select(-geometry, -date), by = "point") %>% 
  mutate(#height.f = fpts$height, 
         var = as.numeric(st_intersects(., mask))) %>% 
  filter(var == 1) %>% select(-var) %>% 
  st_transform(4326) %>% 
  mutate(area = as.numeric(st_area(st_transform(.,32612))), 
         height = raster::extract(chm, st_transform(., 32612), fun = max, method = "bilinear", na.rm = TRUE))

fdf <- fpts %>% 
  mutate(var = as.numeric(st_intersects(., st_transform(mask, 4326)))) %>% 
  filter(var == 1) %>% select(-var) %>% 
  mutate(id = point, type = "field") %>% 
  select(id, lat, long, type) %>% 
  mutate(transect = as.numeric(st_intersects(., st_buffer( st_transform(tr, 4326), 0.01)) )) %>% 
  st_join(fpol %>% select(area, height)) %>% 
  mutate(transect = ifelse(is.na(transect), 0, transect)) %>% 
  as.data.frame() %>% select(-geometry)

# --- end field data


# === 2. Additional hand-digitized plants
adf <- st_read("data/digitized/crp_regen1_2021_additional_crowns.shp") %>% 
  st_transform(4326) %>% 
  mutate(id = paste0("DOWND", id), type = "digit") %>% 
  mutate(area = as.numeric(st_area(st_transform(.,32612))), 
         transect = as.numeric(st_intersects(., st_buffer( st_transform(tr, 4326), 0.01)) ))  %>% 
  select(-notes, -spp) %>% 
  mutate(var = lengths(st_intersects(., st_transform(mask, 4326)))) %>% 
  filter(var > 0) %>% select(-var) %>% 
  mutate(transect = ifelse(is.na(transect), 0, transect), 
         height = raster::extract(chm, st_transform(., 32612), fun = max, method = "bilinear", na.rm = TRUE)) %>% 
  # mutate(height = extract(chm, st_buffer(st_transform(.,32611), .01), fun = max, na.rm = TRUE),
  #        elev = extract(dtm, st_buffer(st_transform(.,32611), .01), fun = max, na.rm = TRUE)) %>%
  st_point_on_surface() %>% 
  mutate(lat = as.numeric(st_coordinates(.)[,2]),
         long = as.numeric(st_coordinates(.)[,1])) %>% 
  as.data.frame() %>% select(-geometry)
# --- end hand digitized


# === 3. Segmented output (note the area)
spts <- st_read("data/segmentation/crp_regen1_ttops.shp") %>% 
  rename(id = treeID, height = Z) %>% 
  st_transform(4326)

spol <- st_read("data/segmentation/crp_regen1_crowns.shp") %>% 
  st_transform(32612) %>% 
  rename(id = 1) %>%
  mutate(id = spts$id) %>% 
  st_buffer(0) %>% 
  mutate(area = as.numeric(st_area(.)) ) %>% 
  mutate(transect = as.numeric(st_intersects(., st_buffer(tr, 0.01)) )) %>% 
  mutate(transect = ifelse(is.na(transect), 0, transect)) %>% 
  st_transform(4326)


sdf <- spts %>% 
  st_join(spol) %>% select(-id.y) %>% rename(id = id.x) %>% 
  mutate(id = paste0("DOWNS", id), 
         lat = as.numeric(st_coordinates(.)[,2]),
         long = as.numeric(st_coordinates(.)[,1]), 
         type = "segment") %>% 
  as.data.frame() %>% select(-geometry)
# --- end segmented data


# === bind all data together

df <- rbind(fdf, adf, sdf) %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>% 
  mutate(lat = as.numeric(st_coordinates(.)[,2]),
         long = as.numeric(st_coordinates(.)[,1])) 

ggplot(df) + geom_point(aes(x = long, y = lat, color = height))
  
# --- export
df %>% 
  st_write("data/model_input/shapes/spdat_crp_regen1.shp", delete_layer = FALSE)

df %>% 
  st_write("data/model_input/spdat_crp_regen1.csv", delete_layer = FALSE)

# === recast the data into unique id observations with individual covariates
df <- st_read("data/model_input/shapes/spdat_crp_regen1.shp")

detectedID <- spol %>% 
  mutate(id = paste0("DOWN", id)) %>% 
  st_intersects(., st_buffer( filter(df, type %in% c("field", "digit")) , 0.05) ) %>% 
  unlist() 

detectedID <- ifelse(nchar(detectedID) == 1, paste0("DOWN00", detectedID), paste0("DOWN0", detectedID))

tmp <- filter(df, type %in% c("field", "digit")) %>% 
  mutate(z = 1, y = ifelse(id %in% detectedID, 1, 0)) %>% 
  st_join(fpts %>% 
              rename(id = point, ht_true = height) %>% 
              select(id, ht_true)) %>% 
  rename(id = id.x) %>% select(-id.y)


tmp %>% 
  ggplot() + geom_point(aes(y = height, x = ht_true)) + 
  facet_wrap(.~as.factor(y)) + labs(x = "True height, m", y = "CHM height", title = "Detected 0/1")

out <- df %>% 
  filter(type == "segment") %>% 
  mutate(y = 1, ht_true = NA) %>% 
  mutate(var = lengths(st_intersects(spol, st_buffer(tmp, .01))) ) %>% 
  mutate(z = var) %>% select(-var) %>%  # true detection
  mutate(z = ifelse(transect == 0 & z == 0, 101, z)) %>% 
  filter(z != 1) %>% 
  bind_rows(tmp) 

st_write(out, "data/model_input/shapes/spdat_joint_crp_regen1.shp", delete_layer = TRUE)
st_write(out, "data/model_input/spdat_joint_crp_regen1.csv", delete_layer = TRUE)

out %>% 
  ggplot() + geom_point(aes(y = height, x = ht_true)) + 
  facet_wrap(.~ as.factor(y))

out %>% 
  mutate(detected = fct_recode(as.factor(y), detected = '1', not_detected = '0'), 
         true_presence = fct_recode(as.factor(z), present = '1', absent = '0', unknown = '101')) %>% 
  ggplot() + geom_histogram(aes(x = height)) + 
  geom_vline(xintercept = .1, colour = "red", size = 1.5) + 
  labs(x = "CHM height, m", y = "Count", title = "Size distributions") +
  facet_wrap(.~detected + true_presence, scales = "free")


# === Step 2: the cell data frame
# --- generate cell ID and assign to obs ID
rout <- aggregate(dtm, 500, expand = FALSE)
chmagg <- aggregate(chm , 500, expand = FALSE)

# create cell IDs
routid <- rout
values(routid) <- 1:length(routid)

# load in the binary index for burnt/not burnt
bi <- st_read("data/mask/crp_regen1_subset_clip_burnIndex.shp") 

# extract cell ID to the main point df and export out
df <- st_read("data/model_input/shapes/spdat_joint_crp_regen1.shp") %>% 
  mutate(cellid = raster::extract(routid, st_transform(., 32612), method = "simple", na.rm = TRUE)) %>% 
  filter(is.na(cellid) == FALSE)

st_write(df, "data/model_input/shapes/spdat_joint_crp_regen1.shp", delete_layer = TRUE)
st_write(df, "data/model_input/spdat_joint_crp_regen1_v1.csv", delete_layer = TRUE)

rpol <- rasterToPolygons(rout) %>%
  st_as_sf() %>%
  mutate(cid = 1:n()) %>%
  rename(elev = 1)

celldf <- rpol %>% 
  st_transform(4326) %>% 
  st_point_on_surface() %>% 
  mutate(lat = as.numeric(st_coordinates(.)[,2]),
         long = as.numeric(st_coordinates(.)[,1]),
         burn = lengths(st_intersects(st_transform(., 32612), st_transform(bi[1,], 32612))), 
         cht = raster::extract(chmagg, st_transform(., 32612), method = "simple", na.rm = TRUE)) %>% 
  left_join(df %>% 
              count(cellid) %>% 
              as.data.frame() %>% select(-geometry), by = c("cid" = "cellid")) %>% 
  mutate(n = replace_na(n, 0)) %>% 
  as.data.frame() %>% select(-geometry)

# add distance
dmat <- celldf %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  st_transform(32612) %>% 
  st_distance(by_element = FALSE)

idx <- which(celldf$burn == 1)
idxb <- dmat[-idx,]

celldf$dist <- 0
for(i in idx){
  celldf$dist[i] = min(idxb[,i])
}
# --- end distance

# write.csv(celldf, "data/model_input/celldat_crp_regen1.csv", row.names = FALSE)
celldf <- read.csv("data/model_input/celldat_crp_regen1.csv")

# --- intensity raster
routint <- reclassify(rout, rcl = matrix(c(-Inf, Inf, 0), byrow = TRUE))
values(routint)[celldf$cid] <- celldf$n
writeRaster(routint, "data/model_input/raster_agg_intensity_crp_regen1.tif", overwrite = FALSE) 


