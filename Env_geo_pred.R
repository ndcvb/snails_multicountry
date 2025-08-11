# All common predictors filtered P3

rm(list=ls())

library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)

# Read the output of the GEE Precipitation script 

df_pp <- read.csv("~/OneDrive - KU Leuven/UG_DRC_CHAD/Data/Uganda/Predictors/CHIRPS_nearest_UG.csv", sep = ",", header = TRUE)

# Find CHIRPS band columns: XYYYYMMDD_precipitation
band_cols <- grep("^X\\d{8}_precipitation$", names(df_pp), value = TRUE)

# Only needed data + nice headers

chirps = as.data.frame(cbind(site_id = df_pp$site_id_renum, df_pp[band_cols]))
names(chirps) <- gsub("^X", "", names(chirps))
names(chirps) <- gsub("_precipitation", "", names(chirps))

# CS data
cs_data = read.csv("~/OneDrive - KU Leuven/UG_DRC_CHAD/Data/Uganda/Snails/events_points_clustered_UG.csv", sep = ",", header = TRUE)

# Ensure proper date formats
cs_data$eventDate <- as.Date(cs_data$eventDate)
new_dates_chirps <- as.Date(names(chirps)[-1], format = "%Y%m%d")

head(new_dates_chirps)

# Create a new column in cs_data for precipitation

cs_data$oneweek = cs_data$eventDate - 7

cs_data$pp_oneweek <- NA  # Initialize with NA

for (i in 1:nrow(cs_data)) {
  site <- cs_data$site_id_renum[i]
  sta_date <- cs_data$oneweek[i]
  end_date <- cs_data$eventDate[i]
  
  site_temps <- chirps %>% filter(site_id == site)
  
  if (nrow(site_temps) == 0) {
    warning(paste("No site match in new_df for:", site, "at row", i))
    next
  }
  
  period <- seq(sta_date, end_date, by = "days")
  period <- gsub("-", "", period)
  
  matching_period <- colnames(site_temps) %in% period
  
  if (sum(matching_period) == 0) {
    warning(paste("No matching date columns for site:", site, "at row", i))
    next
  }
  
  cs_data$pp_oneweek[i] <- rowSums(site_temps[1, matching_period], na.rm = TRUE)
}

summary(cs_data$pp_oneweek)

# Create density plot of precipitation
ggplot(cs_data, aes(x = pp_oneweek)) +
  geom_density(fill = "blue", color = "darkblue", alpha = 0.5) +
  labs(
    x = "Precipitation",
    y = "Density") + theme_bw()

# Temperature 

rm(list = setdiff(ls(), c('cs_data')))

df_st <- read.csv("~/OneDrive - KU Leuven/UG_DRC_CHAD/Data/Uganda/Predictors/MODIS_nearest_UG.csv", sep = ",", header = TRUE)

# Find CHIRPS band columns: XYYYYMMDD_precipitation
band_cols <- grep("^X", names(df_st), value = TRUE)

# Fix format for the GEE output

modis = as.data.frame(cbind(site = df_st$site_id_renum, df_st[band_cols]))
names(modis) <- gsub("^X", "", names(modis))
names(modis) <- gsub("_", "", names(modis))
names(modis) <- gsub("LSTDay1km", "", names(modis))

# Transpose data to start with the interpolation 
modis_transposed <- as.data.frame(t(modis))
colnames(modis_transposed) <- modis_transposed[1, ]
modis_transposed <- modis_transposed[-1, ]

# Convert in numeric
modis_transposed[] <- lapply(modis_transposed, function(x) as.numeric(as.character(x)))
modis_transposed$date <- rownames(modis_transposed)
modis_transposed <- modis_transposed[, c("date", setdiff(names(modis_transposed), "date"))]
rownames(modis_transposed) <- NULL

# Interpolation 
library(imputeTS)
library(zoo)

# Convert the date column to Date format
modis_transposed$date <- as.Date(as.character(modis_transposed$date), format = "%Y%m%d")

summary(modis_transposed)

# Create a sequence of all dates from the minimum to the maximum date
alldates <- seq.Date(from = min(modis_transposed$date), to = max(modis_transposed$date), by = 1)

# Create a data frame with all dates and NA values
NAframes <- data.frame(date = alldates)

# Function to calculate Mean Absolute Error (MAE)
calculate_mae <- function(true_values, predicted_values) {
  mean(abs(true_values - predicted_values), na.rm = TRUE)
}

# Initialize vectors to store the error for each interpolation method
mae_linear <- numeric(length(colnames(modis_transposed)[-1]))
mae_cubic <- numeric(length(colnames(modis_transposed)[-1]))

# Loop through each site column (except the date column)
for (i in seq_along(colnames(modis_transposed)[-1])) {
  
  site <- colnames(modis_transposed)[i + 1]
  
  # Extract data for the current site
  site_data <- data.frame(date = modis_transposed$date, temp = modis_transposed[[site]])
  
  # Remove 50% of the data randomly
  set.seed(123)  # For reproducibility
  missing_indices <- sample(1:nrow(site_data), size = floor(0.5 * nrow(site_data)), replace = FALSE)
  
  # Create a copy of site_data with 50% of the data removed
  site_data_missing <- site_data
  site_data_missing$temp[missing_indices] <- NA
  
  # Merge with NAframes to ensure all dates are present
  merged_data <- merge(NAframes, site_data_missing, by = "date", all.x = TRUE)
  
  # Perform linear interpolation
  interpolated_ts_linear <- na_interpolation(merged_data$temp, option = "linear")
  
  # Perform cubic spline interpolation
  interpolated_ts_cubic <- na.spline(merged_data$temp)
  
  # Compare the interpolated values with the original values
  mae_linear[i] <- calculate_mae(site_data$temp[missing_indices], interpolated_ts_linear[missing_indices])
  mae_cubic[i] <- calculate_mae(site_data$temp[missing_indices], interpolated_ts_cubic[missing_indices])
  
  # Optionally, print or store the results for each site
  cat("Site:", site, "\n")
  cat("MAE Linear:", mae_linear[i], "\n")
  cat("MAE Cubic:", mae_cubic[i], "\n\n")
}

# Summary of the results
results <- data.frame(
  Sites = colnames(modis_transposed)[-1],
  MAE_Linear = mae_linear,
  MAE_Cubic = mae_cubic
)

print(results)

# Determine which method is better for each site
results$Better_Method <- ifelse(results$MAE_Linear < results$MAE_Cubic, "Linear", "Cubic Spline")

print(results)

# Proceed with linear interpolation 

# Initialize empty data frames to store the interpolated results
interpolated_df_linear <- data.frame(date = alldates)
#interpolated_df_cubic <- data.frame(date = alldates)

# Loop through each site column (except the date column)
for (site in colnames(modis_transposed)[-1]) {
  
  # Extract data for the current site
  site_data <- data.frame(date = modis_transposed$date, temp = modis_transposed[[site]])
  
  # Merge with NAframes to ensure all dates are present
  merged_data <- merge(NAframes, site_data, by = "date", all.x = TRUE)
  
  # Perform linear interpolation
  interpolated_ts_linear <- na_interpolation(merged_data$temp, option = "linear")
  
  # Perform cubic spline interpolation
  interpolated_ts_cubic <- na.spline(merged_data$temp)
  
  # Add the interpolated data to the respective data frames
  interpolated_df_linear[[site]] <- interpolated_ts_linear
  interpolated_df_cubic[[site]] <- interpolated_ts_cubic
}

# Save the linear interpolated data to a CSV file
write.csv(interpolated_df_linear, "~/OneDrive - KU Leuven/UG_DRC_CHAD/Data/Uganda/Predictors/temperature_interpolated.csv", row.names = FALSE)

# Function to add source data to df_excel
add_source_data <- function(main_df, source_df, source_prefix) {
  # Initialize the new column
  new_column <- paste0(source_prefix)
  main_df[[new_column]] <- NA
  
  # Iterate through each row of main_df
  for (i in seq_len(nrow(main_df))) {
    current_date <- main_df$eventDate[i]
    current_site <- main_df$site_id_renum[i]
    
    # Match the date in source_df
    matched_row <- source_df %>%
      filter(date == current_date)
    
    if (nrow(matched_row) > 0) {
      # Extract the column name from the current site
      col_name <- as.character(current_site)
      
      # Check if the column exists in the matched row
      if (col_name %in% colnames(matched_row)) {
        # Assign the value from the column
        main_df[[new_column]][i] <- matched_row[[col_name]]
      } else {
        # If the column doesn't exist, assign NA
        main_df[[new_column]][i] <- NA
      }
    }
  }
  
  return(main_df)
}

cs_data <- add_source_data(cs_data, interpolated_df_linear, "temp_modis")

# Create density plot of temperature
ggplot(cs_data, aes(x = temp_modis)) +
  geom_density(fill = "red", color = "darkred", alpha = 0.5) +
  labs(
    x = "Temperature",
    y = "Density") + theme_bw()

rm(list = setdiff(ls(), c('cs_data')))

cs_data <- cs_data %>% dplyr::select(-oneweek)

write.csv(cs_data, "~/OneDrive - KU Leuven/UG_DRC_CHAD/Data/Uganda/UG_cs_pp_temp.csv", row.names = FALSE)

rm(list=ls())

## Geomorphological variables

library(st)
library(sf)

cs_data <- read.csv("~/OneDrive - KU Leuven/UG_DRC_CHAD/Data/Uganda/UG_cs_pp_temp.csv", sep = ",", header = TRUE) 

# Take geometry from the gpkg file (already in UTM)

st_layers("~/OneDrive - KU Leuven/UG_DRC_CHAD/Data/Uganda/Snails/events_points_clustered_UG.gpkg")   # see layer names
pts <- st_read("~/OneDrive - KU Leuven/UG_DRC_CHAD/Data/Uganda/Snails/events_points_clustered_UG.gpkg", layer = "points")  # points or centroids

# Group by Watercontactsite and create convex hulls
polygons_sf <- pts %>%
  group_by(site_id_renum) %>%
  summarize(geom = st_combine(geom)) %>%
  st_convex_hull()  # Create the convex hull

saveRDS(polygons_sf,"~/OneDrive - KU Leuven/UG_DRC_CHAD/Data/Uganda/Snails/Polygons_UG.RDS")

# If you want to visualize the polygons
plot(st_geometry(polygons_sf), col = 'lightblue', border = 'darkblue')

## Geo variables

rm(list = setdiff(ls(), c('cs_data')))

# Extract DEM

library(raster)

# Nothing of this will change for every site so let's calculate for the simple dataframe 

# Use the utm one 

simple_df <- read.csv("~/OneDrive - KU Leuven/UG_DRC_CHAD/Data/Uganda/Snails/dbscan_centroids_utm36n_UG.csv", sep = ",", header = TRUE) 

# Polygons should be in UTM (normally they are)

polygons = readRDS("~/OneDrive - KU Leuven/UG_DRC_CHAD/Data/Uganda/Snails/Polygons_UG.RDS")

st_crs(polygons) # verify 

# Elevation 

image_file <- "~/Library/CloudStorage/OneDrive-KULeuven/Uganda_Congo/Data/Uganda/Clean/D_Copernicus_tiles_reproj_clip_upd.tif"

image_raster <- raster(image_file)  

# Function to compute Q90
q90 <- function(x, na.rm = TRUE) {
  quantile(x, probs = 0.90, na.rm = na.rm)
}

# Function to extract Q90 value for each polygon
extract_value <- function(raster, polygon_sp) {
  value <- tryCatch({
    raster::extract(raster, polygon_sp, fun = q90, na.rm = TRUE, df = TRUE)[, 2]
  }, error = function(e) {
    NA
  })
  return(value)
}


# Create an empty column 'elev' in the Excel file

simple_df$elev <- NA

# Loop through each polygon in polygons_latlon
for (j in 1:nrow(polygons)) {
  # Extract the current polygon
  polygon <- polygons[j, ]
  
  # Convert the current polygon to Spatial object
  polygon_sp <- as(polygon, "Spatial")
  
  # Extract the location from the polygon object
  location <- polygon$site_id_renum
  
  # Extract Q90 values for the current polygon
  elev_value <- extract_value(image_raster, polygon_sp)
  
  # Find the corresponding row in simple_df and update the 'elev' column
  simple_df[simple_df$site_id_renum == location, "elev"] <- elev_value
}

# Flow accumulation 

image_file <- "~/Library/CloudStorage/OneDrive-KULeuven/Uganda_Congo/Data/Uganda/Clean/D_flow_acc_up.tif"

image_raster <- raster(image_file)  

# Create an empty column 'flow_acc' in the Excel file

simple_df$flow_acc <- NA

# Loop through each polygon in polygons_latlon
for (j in 1:nrow(polygons)) {
  # Extract the current polygon
  polygon <- polygons[j, ]
  
  # Convert the current polygon to Spatial object
  polygon_sp <- as(polygon, "Spatial")
  
  # Extract the location from the polygon object
  location <- polygon$site_id_renum
  
  # Extract Q90 values for the current polygon
  flow_value <- extract_value(image_raster, polygon_sp)
  
  # Find the corresponding row in simple_df and update the 'flow_acc' column
  simple_df[simple_df$site_id_renum == location, "flow_acc"] <- flow_value
}

# TPI 

image_file <- "~/Library/CloudStorage/OneDrive-KULeuven/Uganda_Congo/Data/Uganda/Clean/D_TPI_up.tif"

image_raster <- raster(image_file)  

# Create an empty column 'TPI' in the Excel file

simple_df$TPI <- NA

# Loop through each polygon in polygons_latlon
for (j in 1:nrow(polygons)) {
  # Extract the current polygon
  polygon <- polygons[j, ]
  
  # Convert the current polygon to Spatial object
  polygon_sp <- as(polygon, "Spatial")
  
  # Extract the location from the polygon object
  location <- polygon$site_id_renum
  
  # Extract Q90 values for the current polygon
  TPI_value <- extract_value(image_raster, polygon_sp)
  
  # Find the corresponding row in simple_df and update the 'flow_acc' column
  simple_df[simple_df$site_id_renum == location, "TPI"] <- TPI_value
}

# SS Copernicus

image_file <- "~/Library/CloudStorage/OneDrive-KULeuven/Uganda_Congo/Data/Uganda/Clean/D_SS_up.tif"

image_raster <- raster(image_file)  

# Create an empty column 'SS' in the Excel file

simple_df$SS <- NA

# Loop through each polygon in polygons_latlon
for (j in 1:nrow(polygons)) {
  # Extract the current polygon
  polygon <- polygons[j, ]
  
  # Convert the current polygon to Spatial object
  polygon_sp <- as(polygon, "Spatial")
  
  # Extract the location from the polygon object
  location <- polygon$site_id_renum
  
  # Extract Q90 values for the current polygon
  SS_value <- extract_value(image_raster, polygon_sp)
  
  # Find the corresponding row in simple_df and update the 'flow_acc' column
  simple_df[simple_df$site_id_renum == location, "SS"] <- SS_value
}


# Perform the merge
cs_data <- merge(cs_data, simple_df[, c("site_id_renum", "elev", "flow_acc", "TPI", "SS")],
                 by = "site_id_renum", 
                 all.x = TRUE)

write.csv(cs_data, "~/OneDrive - KU Leuven/UG_DRC_CHAD/Data/Uganda/UG_cs_pp_temp_geo.csv", row.names = FALSE)

rm(list = setdiff(ls(), c('cs_data', 'polygons')))

# this to not rename everything again 

polygons_sf = polygons 

#### NDVI

## Polygons already created, then they need to extract the ndvi values 
# Load NDVI raster directory 
ndvi_dir <- "/Volumes/T7 Shield/NDVI_mosaic_output"

# List all NDVI images
ndvi_files <- list.files(ndvi_dir, pattern = "\\.tif$", full.names = TRUE)

# Custom function to extract the 90th percentile [q90 for upland]
q90 <- function(x, na.rm = TRUE) {
  quantile(x, probs = 0.90, na.rm = na.rm)
}

# Function to extract Q90 NDVI value for each polygon
extract_ndvi <- function(ndvi_raster, polygon_sp) {
  # Extract Q90 values
  ndvi_value <- tryCatch({
    raster::extract(ndvi_raster, polygon_sp, fun = q90, na.rm = TRUE, df = TRUE)[, 2]
  }, error = function(e) {
    NA
  })
  return(ndvi_value)
}

# Initialize an empty data frame to store results
ndvi_q90_t <- data.frame()
# Loop through each NDVI image
for (ndvi_file in ndvi_files) {
  # Extract the date from the file name
  date_str <- sub("new_mosaic_by_date_(\\d{8})_NDVI\\.tif", "\\1", basename(ndvi_file))
  
  # Load the NDVI raster
  ndvi_raster <- raster(ndvi_file)
  
  # Loop through each polygon
  for (i in 1:nrow(polygons_sf)) {
    # Extract the current polygon
    polygon <- polygons_sf[i, ]
    
    # Convert the current polygon to Spatial object
    polygon_sp <- as(polygon, "Spatial")
    
    # Extract the location from the polygon object
    location <- polygon$site_id_renum
    
    # Extract NDVI values for the current polygon
    ndvi_value <- extract_ndvi(ndvi_raster, polygon_sp)
    
    # Create a data frame for the current result
    ndvi_data <- data.frame(
      date = date_str,
      location = location,
      ndvi_value = ndvi_value
    )
    
    # Append the result to the all_ndvi_data data frame
    ndvi_q90_t <- rbind(ndvi_q90_t, ndvi_data)
  }
}
# Display combined NDVI data
print(ndvi_q90_t)

ndvi_q90_t_w = ndvi_q90_t %>% pivot_wider(names_from = location,  values_from = ndvi_value)

#write.csv(ndvi_q90_t_w, "~/Library/CloudStorage/OneDrive-KULeuven/Paper3/Data/D_NDVI_q90.csv", row.names = FALSE )

# Convert the date column to Date format
ndvi_q90_t_w$date <- as.Date(as.character(ndvi_q90_t_w$date), format = "%Y%m%d")

# Create a sequence of all dates from the minimum to the maximum date
alldates <- seq.Date(from = min(ndvi_q90_t_w$date), to = max(ndvi_q90_t_w$date), by = 1)

# Create a data frame with all dates and NA values
NAframes <- data.frame(date = alldates)

# Function to calculate Mean Absolute Error (MAE)
calculate_mae <- function(true_values, predicted_values) {
  mean(abs(true_values - predicted_values), na.rm = TRUE)
}

# Initialize vectors to store the error for each interpolation method
mae_linear <- numeric(length(colnames(ndvi_q90_t_w)[-1]))
mae_cubic <- numeric(length(colnames(ndvi_q90_t_w)[-1]))

# Loop through each Watercontactsite column (except the date column)
for (i in seq_along(colnames(ndvi_q90_t_w)[-1])) {
  
  site <- colnames(ndvi_q90_t_w)[i + 1]
  
  # Extract data for the current site
  site_data <- data.frame(date = ndvi_q90_t_w$date, NDVI = ndvi_q90_t_w[[site]])
  
  # Remove 50% of the data randomly
  set.seed(123)  # For reproducibility
  missing_indices <- sample(1:nrow(site_data), size = floor(0.5 * nrow(site_data)), replace = FALSE)
  
  # Create a copy of site_data with 50% of the data removed
  site_data_missing <- site_data
  site_data_missing$NDVI[missing_indices] <- NA
  
  # Merge with NAframes to ensure all dates are present
  merged_data <- merge(NAframes, site_data_missing, by = "date", all.x = TRUE)
  
  # Perform linear interpolation
  interpolated_ts_linear <- na_interpolation(merged_data$NDVI, option = "linear")
  
  # Perform cubic spline interpolation
  interpolated_ts_cubic <- na.spline(merged_data$NDVI)
  
  # Compare the interpolated values with the original values
  mae_linear[i] <- calculate_mae(site_data$NDVI[missing_indices], interpolated_ts_linear[missing_indices])
  mae_cubic[i] <- calculate_mae(site_data$NDVI[missing_indices], interpolated_ts_cubic[missing_indices])
  
  # Optionally, print or store the results for each site
  cat("Site:", site, "\n")
  cat("MAE Linear:", mae_linear[i], "\n")
  cat("MAE Cubic:", mae_cubic[i], "\n\n")
}

# Summary of the results
results <- data.frame(
  Watercontactsite = colnames(ndvi_q90_t_w)[-1],
  MAE_Linear = mae_linear,
  MAE_Cubic = mae_cubic
)

print(results)

# Determine which method is better for each site
results$Better_Method <- ifelse(results$MAE_Linear < results$MAE_Cubic, "Linear", "Cubic Spline")

print(results)

# Initialize empty data frames to store the interpolated results
interpolated_df_linear <- data.frame(date = alldates)
interpolated_df_cubic <- data.frame(date = alldates)

# Loop through each Watercontactsite column (except the date column)
for (site in colnames(ndvi_q90_t_w)[-1]) {
  
  # Extract data for the current site
  site_data <- data.frame(date = ndvi_q90_t_w$date, NDVI = ndvi_q90_t_w[[site]])
  
  # Merge with NAframes to ensure all dates are present
  merged_data <- merge(NAframes, site_data, by = "date", all.x = TRUE)
  
  # Perform linear interpolation
  interpolated_ts_linear <- na_interpolation(merged_data$NDVI, option = "linear")
  
  # Perform cubic spline interpolation
  interpolated_ts_cubic <- na.spline(merged_data$NDVI)
  
  # Add the interpolated data to the respective data frames
  interpolated_df_linear[[site]] <- interpolated_ts_linear
  interpolated_df_cubic[[site]] <- interpolated_ts_cubic
}

write.csv(interpolated_df_linear, "~/OneDrive - KU Leuven/UG_DRC_CHAD/Data/Uganda/Predictors/UG_NDVI_interpolated.csv", row.names = FALSE)

NDVI_q90_t = interpolated_df_linear

rm(list = setdiff(ls(), c('cs_data','NDVI_q90_t')))

## Add terrestrial NDVI to the datasets 

# Function to add source data to df_excel
add_source_data <- function(main_df, source_df, source_prefix) {
  # Initialize the new column
  new_column <- paste0(source_prefix)
  main_df[[new_column]] <- NA
  
  # Iterate through each row of main_df
  for (i in seq_len(nrow(main_df))) {
    current_date <- main_df$eventDate[i]
    current_site <- main_df$site_id_renum[i]
    
    # Match the date in source_df
    matched_row <- source_df %>%
      filter(date == current_date)
    
    if (nrow(matched_row) > 0) {
      # Extract the column name from the current site
      col_name <- as.character(current_site)
      
      # Check if the column exists in the matched row
      if (col_name %in% colnames(matched_row)) {
        # Assign the value from the column
        main_df[[new_column]][i] <- matched_row[[col_name]]
      } else {
        # If the column doesn't exist, assign NA
        main_df[[new_column]][i] <- NA
      }
    }
  }
  
  return(main_df)
}

# Add columns for each source dataset
cs_data <- add_source_data(cs_data, NDVI_q90_t, "NDVI_q90_t")

# Save the terrestrial NDVI interpolated data to a CSV file
write.csv(cs_data, "~/OneDrive - KU Leuven/UG_DRC_CHAD/Data/Uganda/UG_geo_env.csv", row.names = FALSE)



