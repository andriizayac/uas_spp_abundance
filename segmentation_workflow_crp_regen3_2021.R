pkgs <- c("lidR", "raster", "gstat", "sf", "sp", "viridisLite", "stars", "dplyr")
sapply(pkgs, require, character.only = TRUE)

# NOTE: this script operates on large files that are currently not available from github
# outputs:
# - 1 - tree tops as a point layer: data/segmentation/ttops.shp
# - 2 - crown segments as a polygon layer: data/segmentation/crowns.shp
# [for a subset of of the field]

# === helper stuff
cols <- gray.colors(255)
colsid <- sample(viridis(50, alpha = 1, direction = 1, begin = 0))

# === Clouds NOT USED at this point
{
  # ---read in the las catalog
  path <- "D:/UAV_analysis/SodaBurn2021/"
  # --- select a subset
  tiles <- readLAScatalog(paste0(path, "tiles/")) 
  ts <- catalog_select(tiles)
  # --- retile if needed
  # opt_output_files(ts) <- "C:/tmp/dtm_workflow_example/crp_regen3_2021_{ID}"
  # opt_chunk_buffer(ts) <- 0
  # opt_chunk_size(ts) <- 8
  # newts <- catalog_retile(ts)
}

path <- "D:/UAV_analysis/SodaBurn2021/"

# bring dtm 
dtm <- raster(list.files(path, pattern = "20210615_CHM.*\\.tif", full.names = TRUE))

# normalize las, applies for cloud-based segmentations
las <- readLAS(ts, select = "xyz")
# lasn <- normalize_height(las, dtm)
# lasn <- filter_poi(lasn, Z > -.2)

# import and crop the CHM
path <- "./../"
clist <- list.files(path, pattern = "regen3.*2021.*_CHM_.*.tif", full.names = TRUE)
chm <- raster(clist)


# - reclassify pixels below zero if needed
chm <- reclassify(chm, rcl = matrix(c(-Inf, 0, 0), byrow = TRUE))
# --[not used at this point]-- create roi (buffer = -.9m) to avoid half-captured plants 
roi <- chm %>% 
  st_bbox() %>% 
  st_as_sfc() %>% 
  st_buffer(0)

fn <- function(x) { 0.2 + 2 * tanh(2*x)^3 }
# find trees (each year separately)
ttops <- find_trees(chm, lmf(ws = fn, hmin = .1, shape = "square"))

# segment trees in CHM
cro <- silva2016(chm, ttops, max_cr_factor = .5, exclusion = 0.05)()

# === visualize
plot(chm, col = cols)
plot(cro, add = TRUE, col = colsid)
plot(ttops, add = TRUE, pch = 19, cex = .5, col = "black")

# --- subset to roi
tout <- ttops %>% 
  st_as_sf() %>% 
  st_crop(st_bbox(roi)) 

cout <- cro %>% 
  st_as_stars() %>% 
  st_as_sf(as_points = FALSE, merge = TRUE) %>% 
  mutate(var = lengths(st_intersects(., st_buffer(tout, 0.01)) )) %>% 
  filter(var > 0) %>% 
  select(-var)

# tout <- tout %>% 
#  mutate(var = lengths(st_intersects(., st_buffer(cout, 0.01)) ))

# ---
orphans <- cout %>% 
  filter(var == 1) %>% 
  select(-var) %>% 
  st_buffer(0.01) 
mutate(mer = unlist(st_intersects(., filter(cout, var == 0))))
# finish up - what to do with orphaned segments - have id's but not points 

# === export
pathout <- "D:/Andrii/abundance_sodaBurn_2021/uas_spp_abundance/data/segmentation/"
tout %>% 
  st_write(paste0(pathout, "crp_regen3_ttops.shp"), delete_layer = TRUE)

cout %>% 
  st_write(paste0(pathout, "crp_regen3_crowns.shp"), delete_layer = TRUE)



