# =================== packages ===================
library(sf)
library(dplyr)
library(dbscan)
library(ggplot2)

# =================== files ===============
filters <-  read.csv("~/Documents/Uganda_DRC/Data_UG_filters.csv")      # has X_id and Watercontactsite(s)
occ_raw     <- read.csv("~/OneDrive - KU Leuven/UG_DRC_CHAD/Data/occurrence_UG_DRC.csv") ## This is the output from GBIF
# https://www.gbif.org/dataset/1b001cca-29dd-4446-a0cc-f81b7c939f69 
# Running from line 90 it's recommended

# =================== CRS in metres ==============
# Uganda -> UTM 36N; change per country later
crs_m <- 32636
pick_site_col <- function(df){
  cand <- c("Watercontactsite")
  site_col <- cand[cand %in% names(df)][1]
  if (is.na(site_col)) stop("Couldn't find a Watercontactsite column in Data_UG_filters.csv")
  site_col
}

# =================== STEP A: MAD (median absolut deviation) radius =========

occ_raw = occ_raw[ , colSums(is.na(occ_raw)) == 0 ]
occ_raw = occ_raw %>% filter(countryCode=='UG')

# Identify the site column name in filters
site_col <- pick_site_col(filters)

# Build a SHORT occurrence table (one row per eventID) and join site labels via X_id
# This is to not have all the occurrences making things difficult
occ_short <- occ_raw %>%
  transmute(
    eventID          = trimws(as.character(eventID)),
    decimalLongitude = as.numeric(decimalLongitude),
    decimalLatitude  = as.numeric(decimalLatitude),
    eventDate        = as.Date(eventDate)
  ) %>%
  filter(!is.na(eventID), !is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
  group_by(eventID) %>% slice(1) %>% ungroup()

# Join Watercontactsite from filters (X_id -> eventID)

occ_short$eventID = as.integer(occ_short$eventID)

occ_mad <- occ_short %>%
  left_join(filters %>% select(X_id, !!site_col), by = c("eventID" = "X_id"))

# Keep rows that have a Watercontactsite label
occ_mad <- occ_mad %>% filter(!is.na(.data[[site_col]]))

# To sf (metres)
occ_mad_sf <- st_as_sf(occ_mad,
                       coords = c("decimalLongitude","decimalLatitude"),
                       crs = 4326) %>%
  st_transform(crs_m)

# One centre per Watercontactsite (centroid of MULTIPOINT)
site_centres <- occ_mad_sf %>%
  group_by(site_id = .data[[site_col]]) %>%
  summarise(.groups = "drop") %>%
  st_centroid()

# Distances event -> its site centre (aligned vectors)
idx <- match(occ_mad_sf[[site_col]], site_centres$site_id)
centre_vec <- site_centres$geometry[idx]
d_m <- as.numeric(st_distance(st_geometry(occ_mad_sf), centre_vec, by_element = TRUE))

dist_df <- tibble::tibble(
  site_id = occ_mad_sf[[site_col]],
  dist_m  = d_m
)

# MAD clean within each Watercontactsite; ignore sites with <5 points in the radius calc
site_clean <- dist_df %>%
  group_by(site_id) %>%
  mutate(n_site = dplyr::n(),
         med = median(dist_m, na.rm=TRUE),
         mad = mad(dist_m, constant = 1.4826, na.rm=TRUE),
         keep = dist_m <= (med + 3*mad)) %>%
  ungroup() %>%
  filter(n_site >= 5)

# Robust radius = 95th percentile of kept distances × small safety factor
radius_m <- quantile(site_clean$dist_m[site_clean$keep], 0.95, na.rm=TRUE) * 1.10
radius_m <- as.numeric(radius_m)
message(sprintf("Chosen robust radius (MAD on Watercontactsite): %.1f m", radius_m))

rm(list=ls())

### REPLICATE FROM HERE

# =================== STEP B: DBSCAN by recordedBy (occurrence only) ===========
# Build a SHORT occurrence table with recorder

# radius : 191.7 m

occ_raw     <- read.csv("~/OneDrive - KU Leuven/UG_DRC_CHAD/Data/occurrence_UG_DRC.csv")
occ_raw = occ_raw[ , colSums(is.na(occ_raw)) == 0 ]
occ_raw = occ_raw %>% filter(countryCode=='UG')

# just events

occ_short2 <- occ_raw %>%
  transmute(
    eventID          = as.character(eventID),
    decimalLongitude = as.numeric(decimalLongitude),
    decimalLatitude  = as.numeric(decimalLatitude),
    eventDate        = as.Date(eventDate),
    recordedBy       = recordedBy
  ) %>%
  group_by(eventID) %>% slice(1) %>% ungroup()

crs_m <- 32636

occ_sf <- st_as_sf(occ_short2,
                   coords = c("decimalLongitude","decimalLatitude"),
                   crs = 4326) %>%
  st_transform(crs_m)

# DBSCAN per recorder (so different citizens don't merge)
minPts_user <- 5  # 

radius_m= 191.7

groups <- split(occ_sf, occ_sf$recordedBy)
res_list <- vector("list", length(groups)); gi <- 0L
for (nm in names(groups)){
  gi <- gi + 1L
  g  <- groups[[nm]]
  cl <- dbscan(st_coordinates(g), eps = radius_m, minPts = minPts_user)$cluster
  lab <- ifelse(cl == 0, NA_character_, paste0("CS", gi, "_C", sprintf("%03d", cl)))
  g$site_id_dbscan <- lab
  res_list[[gi]] <- g
}
occ_db <- do.call(rbind, res_list)

# Inspect cluster sizes
cluster_counts <- occ_db %>%
  st_drop_geometry() %>%
  dplyr::count(recordedBy, site_id_dbscan, name = "n_points") %>%
  arrange(desc(n_points))
print(head(cluster_counts, 20))
message("Total clusters (non-NA): ",
        length(unique(na.omit(occ_db$site_id_dbscan))))

# Filter clusters with at least 10 points 
kept_clusters <- occ_db %>%
  st_drop_geometry() %>%
  filter(!is.na(site_id_dbscan)) %>%                  # ignore noise
  count(recordedBy, site_id_dbscan, name = "n_points") %>%
  filter(n_points >= 10) %>%
  arrange(recordedBy, site_id_dbscan) %>%
  mutate(site_id_renum = sprintf("S%04d", row_number()))   # neat new IDs

# Tag every point with n_points, keep flag, and the renumbered ID
occ_db_tag <- occ_db %>%
  left_join(kept_clusters, by = c("recordedBy","site_id_dbscan")) %>%
  mutate(keep = !is.na(site_id_renum),
         n_points = ifelse(is.na(n_points), 0L, n_points))

# Filter to kept points (>=10) for discs
occ_db_filt <- occ_db_tag %>% filter(keep)

# Build equal-area discs for the kept clusters (uses your radius_m)
site_centres_db_10 <- occ_db_filt %>%
  group_by(site_id_renum) %>%
  summarise(.groups = "drop") %>%
  st_centroid()

site_discs_db_10 <- st_buffer(site_centres_db_10, dist = radius_m)

# Join back to the original occurrence table (by eventID)
#    Make sure eventID  match
occ_join_cols <- occ_db_tag %>%
  st_drop_geometry() %>%
  transmute(
    eventID,
    recordedBy,
    site_id_dbscan,
    site_id_renum,
    n_points,
    keep
  )

occ_raw2 <- occ_raw %>%
  mutate(eventID = trimws(as.character(eventID))) %>%
  left_join(occ_join_cols %>%
              mutate(eventID = trimws(as.character(eventID))),
            by = "eventID")

occ_raw2 = occ_raw2 %>% filter(n_points >= 10)

# save
write.csv(occ_raw2, "~/OneDrive - KU Leuven/UG_DRC_CHAD/Data/occurrence_with_clusters_UG.csv", row.names = FALSE)
st_write(occ_db_filt, "~/OneDrive - KU Leuven/UG_DRC_CHAD/Data/events_points_clustered_UG.gpkg", layer = "points", delete_layer = TRUE)

## But we need centroids to extract Temp / Precipitation (Centroids of the new created circles - This is just for the big variables)

# Build centroids for kept clusters (one point per cluster)
site_centres_db_10 <- occ_db_filt %>%
  group_by(site_id_renum) %>%                  # your tidy cluster ID
  summarise(.groups = "drop") %>%              # MULTIPOINT per cluster
  st_centroid()                                # -> POINT per cluster

# Add some metadata (cs + cluster size)
centroid_meta <- occ_db_filt %>%
  st_drop_geometry() %>%
  group_by(site_id_renum) %>%
  summarise(
    n_points   = first(n_points),              # from add_count earlier
    recordedBy = paste(sort(unique(recordedBy)), collapse = "; "),
    .groups = "drop"
  )

site_centroids <- site_centres_db_10 %>%
  left_join(centroid_meta, by = "site_id_renum")

# Write GeoPackage (recommended)
gpkg <- "~/OneDrive - KU Leuven/UG_DRC_CHAD/Data/dbscan_sites_UG.gpkg" 
st_write(site_centroids,   gpkg, layer = "site_centroids", delete_layer = TRUE)  
st_write(site_discs_db_10, gpkg, layer = "site_discs",     delete_layer = TRUE)## this is to add to the previous layer

# Write CSV in UTM 36N (metres)
centroids_utm <- site_centroids %>%
  mutate(x_m = st_coordinates(.)[,1],
         y_m = st_coordinates(.)[,2]) %>%
  st_drop_geometry()
write.csv(centroids_utm, "~/OneDrive - KU Leuven/UG_DRC_CHAD/Data/dbscan_centroids_utm36n_UG.csv", row.names = FALSE)

# Just in case, in WGS84 lon/lat 
centroids_wgs84 <- st_transform(site_centroids, 4326) %>%
  mutate(lon = st_coordinates(.)[,1],
         lat = st_coordinates(.)[,2]) %>%
  st_drop_geometry()
write.csv(centroids_wgs84, "~/OneDrive - KU Leuven/UG_DRC_CHAD/Data/dbscan_centroids_wgs84_UG.csv", row.names = FALSE)

# Map: grey = dropped (<10 or noise), color = kept clusters; discs outlined
ggplot() +
  geom_sf(data = occ_db_tag %>% filter(!keep), color = "grey70", size = 0.8, alpha = 0.6) +
  geom_sf(data = site_discs_db_10, fill = NA, color = "black") +
  geom_sf(data = occ_db_filt, aes(color = site_id_renum), size = 1.1) +
  guides(color = guide_legend(title = "Kept cluster")) +
  coord_sf() + theme_minimal() +
  labs(title = "DBSCAN clusters (kept ≥ 10 points)",
       subtitle = sprintf("eps ≈ %.0f m; noise/ small clusters shown in grey", radius_m))
