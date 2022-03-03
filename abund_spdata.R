pkgs <- c("lidR", "raster", "gstat", "sf", "sp", "viridisLite", "stars", "tidyverse", "ggplot2")
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
tr <- st_read("data/transect_shapes/Soda_regen_transects.shp") %>% 
  select(-id)

mask <- st_read("data/mask/sodaBurn_subset_clip.shp")

chm <- raster("D:/Andrii/abundance_sodaBurn_2021/SodaBurn_20210615_CHM_WGS84UTM11_subset.tif")
dtm <- raster("D:/Andrii/abundance_sodaBurn_2021/SodaBurn_20210615_DTM_WGS84UTM11.tif") %>% 
 crop(chm@extent)

# === 1. Field data: points and polygons
fpts <- read.csv("data/digitized/Soda_regen_2021_recr_august.csv") %>% 
  mutate(x = long, y = lat, z = elev) %>% 
  st_as_sf(coords = c("x", "y"), crs = 4326) %>%
  filter(plant == 1) 

fpol <- st_read("data/digitized/soda_regen_2021_digitized.shp") %>% 
  mutate(height = fpts$height, 
         var = as.numeric(st_intersects(., mask))) %>% 
  filter(var == 1) %>% select(-var) %>% 
  st_transform(4326) %>% 
  mutate(area = as.numeric(st_area(st_transform(.,32611))), 
         height = raster::extract(chm, st_transform(., 32611), fun = max, method = "bilinear", na.rm = TRUE)) %>% 
  select(area, height) 

fdf <- fpts %>% 
  mutate(var = as.numeric(st_intersects(., mask))) %>% 
  filter(var == 1) %>% select(-var) %>% 
  filter(plant == 1) %>% 
  mutate(id = point, type = "field") %>% 
  select(id, lat, long, type) %>% 
  mutate(transect = as.numeric(st_intersects(., st_buffer(tr, 0.01)) )) %>% 
  st_join(fpol) %>% 
  mutate(transect = ifelse(is.na(transect), 0, transect)) %>% 
  as.data.frame() %>% select(-geometry)

# --- end field data


# === 2. Additional hand-digitized plants
adf <- st_read("data/digitized/SodaBurn_2021_additional_crowns.shp") %>% 
  mutate(id = paste0("SND", id), type = "digit") %>% 
  mutate(area = as.numeric(st_area(st_transform(.,32611)))) %>% 
  mutate(transect = as.numeric(st_intersects(., st_buffer(tr, 0.01)) )) %>% 
  select(-note) %>% 
  mutate(transect = ifelse(is.na(transect), 0, transect), 
         height = raster::extract(chm, st_transform(., 32611), fun = max, method = "bilinear", na.rm = TRUE)) %>% 
  # mutate(height = extract(chm, st_buffer(st_transform(.,32611), .01), fun = max, na.rm = TRUE),
  #        elev = extract(dtm, st_buffer(st_transform(.,32611), .01), fun = max, na.rm = TRUE)) %>% 
  st_point_on_surface() %>% 
  mutate(lat = as.numeric(st_coordinates(.)[,2]),
         long = as.numeric(st_coordinates(.)[,1])) %>% 
  as.data.frame() %>% select(-geometry)
# --- end hand digitized


# === 3. Segmented output (note the area)
spts <- st_read("data/segmentation/ttops.shp") %>% 
  rename(id = treeID, height = Z) %>% 
  st_transform(4326)

spol <- st_read("data/segmentation/crowns.shp") %>% 
  rename(id = SB_2021) %>%
  mutate(id = spts$id) %>% 
  st_buffer(0) %>% 
  st_transform(32611) %>% 
  mutate(area = as.numeric(st_area(.))) %>% 
  st_transform(4326) %>% 
  mutate(transect = as.numeric(st_intersects(., st_buffer(tr, 0.01)) )) %>% 
  mutate(transect = ifelse(is.na(transect), 0, transect))


sdf <- spts %>% 
  mutate(id = paste0("SNS", id)) %>% 
  st_join(spol) %>% select(-id.y) %>% rename(id = id.x) %>% 
  mutate(lat = st_coordinates(.)[,2],
         long = st_coordinates(.)[,1], type = "segment") %>% 
  mutate(lat = as.numeric(lat), long = as.numeric(long)) %>% 
  as.data.frame() %>% select(-geometry)
# --- end segmented data


# === bind all data together

df <- rbind(fdf, adf, sdf) %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>% 
  mutate(lat = st_coordinates(.)[,2],
         long = st_coordinates(.)[,1]) %>% 
  mutate(lat = as.numeric(lat), long = as.numeric(long))

ggplot(df) + geom_point(aes(x = long, y = lat, color = height))
  
# --- export
df %>% 
  st_write("data/model_input/shapes/spdat.shp", delete_layer = FALSE)

df %>% 
  st_write("data/model_input/spdat.csv", delete_layer = FALSE)

# === recast the data into unique id observations with individual covariates
df <- st_read("data/model_input/shapes/spdat.shp")

detectedID <- unlist(st_intersects(spol, st_buffer( filter(df, type %in% c("field", "digit")) , 0.05)))

tmp <- filter(df, type %in% c("field", "digit")) %>% 
  mutate(z = 1, 
         y = ifelse(1:n() %in% detectedID, 1, 0), 
         ht_true = c(fpts %>% slice(3:n()) %>% 
                       mutate(height = height/100) %>%
                       pull(height), rep(NA, 2)) )
tmp %>% 
  ggplot() + geom_point(aes(y = height, x = ht_true)) + 
  facet_wrap(.~as.factor(y)) + labs(x = "True heught, m", y = "CHM height", title = "Detected 0/1")

out <- df %>% 
  filter(type == "segment") %>% 
  mutate(y = 1, ht_true = NA) %>% 
  mutate(var = lengths(st_intersects(spol, st_buffer(tmp, .05))) ) %>% 
  mutate(z = var) %>% select(-var) %>%  # true detection
  mutate(z = ifelse(transect == 0 & z == 0, 101, z)) %>% 
  filter(z != 1) %>% 
  bind_rows(tmp) 

st_write(out, "data/model_input/shapes/spdat_joint.shp", delete_layer = FALSE)
st_write(out, "data/model_input/spdat_joint.csv", delete_layer = FALSE)

out %>% 
  ggplot() + geom_point(aes(y = height, x = ht_true)) + 
  facet_wrap(.~as.factor(y))

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
bi <- st_read("data/mask/sodaBurn_subset_clip_burnIndex.shp")

# extract cell ID to the main point df and export out
df <- st_read("data/model_input/shapes/spdat_joint.shp") %>% 
  mutate(cellid = raster::extract(routid, st_transform(., 32611), method = "simple", na.rm = TRUE)) %>% 
  filter(is.na(cellid) == FALSE)

st_write(df, "data/model_input/shapes/spdat_joint.shp", delete_layer = TRUE)
st_write(df, "data/model_input/spdat_joint0.csv", delete_layer = TRUE)

rpol <- rasterToPolygons(rout) %>%
  st_as_sf() %>%
  mutate(cid = 1:n()) %>%
  rename(elev = 1)

celldf <- rpol %>% 
  st_transform(4326) %>% 
  st_point_on_surface() %>% 
  mutate(lat = as.numeric(st_coordinates(.)[,2]),
         long = as.numeric(st_coordinates(.)[,1]),
         burn = lengths(st_intersects(st_transform(., 32611), bi[1,])), 
         cht = raster::extract(chmagg, st_transform(., 32611), method = "simple", na.rm = TRUE)) %>% 
  left_join(df %>% 
              count(cellid) %>% 
              as.data.frame() %>% select(-geometry), by = c("cid" = "cellid")) %>% 
  mutate(n = replace_na(n, 0)) %>% 
  as.data.frame() %>% select(-geometry)

# add distance
dmat <- celldf %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  st_transform(32611) %>% 
  st_distance(by_element = FALSE)

idx <- which(celldf$burn == 1)
idxb <- dmat[-idx,]

celldf$dist <- 0
for(i in idx){
  celldf$dist[i] = min(idxb[,i])
}
# --- end distance

# write.csv(celldf, "data/model_input/celldat.csv", row.names = FALSE)
celldf <- read.csv("data/model_input/celldat.csv")

# --- intensity raster
routint <- reclassify(rout, rcl = matrix(c(-Inf, Inf, 0), byrow = TRUE))
values(routint)[celldf$cid] <- celldf$n
writeRaster(routint, "data/model_input/raster_agg_intensity.tif", overwrite = FALSE) 


