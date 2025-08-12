
# Anthropogenic variables 

rm(list=ls())

library(dplyr)
library(tidyr)
library(sf)

cs_data = read.csv("~/OneDrive - KU Leuven/UG_DRC_CHAD/Data/Uganda/UG_geo_env.csv")

# 1. Load both shapefiles (with extracted values per buffer)
shp_2000 <- st_read("~/OneDrive - KU Leuven/UG_DRC_CHAD/Data/Uganda/Predictors/UG_2000_buffer_pop.shp")
shp_2020 <- st_read("~/OneDrive - KU Leuven/UG_DRC_CHAD/Data/Uganda/Predictors/UG_2020_buffer_pop.shp")


# 2. Keep only relevant columns (site and value)
df_2000 <- shp_2000 %>% st_drop_geometry() %>% dplyr::select(site_id_re, value_2000 = X_sum)
df_2020 <- shp_2020 %>% st_drop_geometry() %>% dplyr::select(site_id_re, value_2020 = X_sum)

# 3. Merge into a single dataframe
comparison_df <- left_join(df_2000, df_2020, by = "site_id_re")

# 4. Add a difference column
comparison_df <- comparison_df %>%
  mutate(pop_rate_change= (value_2020 - value_2000) / value_2000
  )

# 5. View the result
print(comparison_df)

colnames(comparison_df)[1] = 'site_id_renum'

cs_data <- merge(cs_data, comparison_df[, c("site_id_renum", "pop_rate_change")],
                 by = "site_id_renum", 
                 all.x = TRUE)

## Dynamic World - Land use 

rm(list = setdiff(ls(), c('cs_data')))

# Load your CSV
dw <- read.csv("~/OneDrive - KU Leuven/UG_DRC_CHAD/Data/Uganda/Predictors/Dynamic_World_UG.csv")

# Convert to character if needed
dw$histogram <- as.character(dw$histogram)

# 1: Replace curly braces and split into rows
dw_long <- dw %>%
  mutate(histogram = str_remove_all(histogram, "[{}]")) %>%
  separate_rows(histogram, sep = ",\\s*") %>%                   # One class=value per row
  separate(histogram, into = c("class", "value"), sep = "=", convert = TRUE) %>%
  mutate(class = paste0("class_", class))                       # Format class names

# 2: Pivot wider — fill missing with 0
dw_wide <- dw_long %>%
  pivot_wider(names_from = class, values_from = value, values_fill = 0)

# 3.Rename class columns to readable land cover names
dw_named <- dw_wide %>%
  rename(
    water = class_0,
    trees = class_1,
    grass = class_2,
    flooded_vegetation = class_3,
    crops = class_4,
    shrub_and_scrub = class_5,
    built = class_6,
    bare = class_7,
    snow_and_ice = class_8
  )

# 4.Calculate total and proportions
dw_props <- dw_named %>%
  mutate(total_pixels = water + trees + grass + flooded_vegetation +
           crops + shrub_and_scrub + built + bare + snow_and_ice) %>%
  mutate(across(
    c(water, trees, grass, flooded_vegetation, crops, shrub_and_scrub, built, bare, snow_and_ice),
    ~ .x / total_pixels,
    .names = "prop_{.col}"
  ))

# 5. Convert date if needed
dw_props <- dw_props %>%
  mutate(date = as.Date(end_date))  # this is just for the plot 

# Select two sites 
sites_to_plot <- c("S0001", "S0002")

dw_long <- dw_props %>%
  filter(site_id_re %in% sites_to_plot) %>%
  dplyr::select(site_id_re, date, starts_with("prop_")) %>%
  pivot_longer(cols = starts_with("prop_"), names_to = "class", values_to = "proportion") %>%
  mutate(class = str_remove(class, "prop_"))

dw_palette <- c(
  water = "#419bdf",
  trees = "#397d49",
  grass = "#88b053",
  flooded_vegetation = "#7a87c6",
  crops = "#e49635",
  shrub_and_scrub = "#dfc35a",
  built = "#c4281b",
  bare = "#a59b8f",
  snow_and_ice = "#b39fe1"
)

ggplot(dw_long, aes(x = date, y = proportion, fill = class)) +
  geom_area(position = "stack") +
  facet_wrap(~site_id_re, scales = "free_x") +
  scale_fill_manual(values = dw_palette) +
  labs(    x = "Date",
           y = "Proportion of pixels",
           fill = "Land cover class"
  ) +
  theme_minimal()

dw_slim <- dw_props %>%
  dplyr::select(site_id_re, start_date, end_date, starts_with("prop_"))

# Make sure types line up
cs_data <- cs_data %>%
  mutate(site_id_renum = as.character(site_id_renum),
         eventDate     = as.Date(eventDate)) # this is the actual obs day 

dw_slim <- dw_slim %>%
  mutate(site_id_re = as.character(site_id_re),
         start_date = as.Date(start_date),
         end_date   = as.Date(end_date))

# eventDate is between [start_date, end_date] ## Please check the amount of rows 
cs_with_dw <- cs_data %>%
  left_join(
    dw_slim,
    by = join_by(
      site_id_renum == site_id_re,
      between(eventDate, start_date, end_date)
    )
  )
# Save to file
write.csv(cs_with_dw, "~/OneDrive - KU Leuven/UG_DRC_CHAD/Data/Uganda/UG_complete_event.csv", row.names = FALSE)

