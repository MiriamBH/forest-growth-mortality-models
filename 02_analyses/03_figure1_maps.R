##############################################################################
##
## Script name: Fig 1 maps
## Purpose: Create productivity components maps and BL-NL proportion map (Figure1)
##
## Input data: species_list.rds
## Output data: none 
## Output figures: Maps in Figure 1 
##
## Date: 2026-02-25
##
##############################################################################

# 1. Libraries and read data -------------------------------------------------
library(tidyverse)
library(here)
library(rnaturalearth)
library(sf)
library(ggmap)
library(tidyterra)
library(terra)
library(ggspatial)
library(viridis)

species_list <- here("01_data", "species_list.rds") |> read_rds()

# 2. Create BL-NL proportion map ---------------------------------------------

data <- species_list |> list_rbind() |> 
  group_by(id, longitude, latitude) |> 
  summarise(conifer_prop = mean(conifer_prop))

#Europe's base map 
europe <- ne_countries(
  scale = "medium", returnclass = "sf", continent = "Europe")

europe <- st_transform(europe, crs = 4326)

st_is_valid(europe)

# filter all values over the 95% quantile
q_high <- quantile(data$conifer_prop, probs = 0.95, na.rm = T)

# sp_data2 <- data_1[!(data_1[var_name] < q_low),]
data2 <- data[(data$conifer_prop < q_high),]

# verctor con datos de sp_data
dat_df<- terra::vect(
  x = data2,
  geom = c("longitude", "latitude"),
  crs = "EPSG:4326", 
  keep = T
)

figure1_bl_nl <- ggplot() +
  geom_spatvector(
    data = europe,
    fill = "grey80",
    color = "grey80"
  ) +
  stat_summary_hex(
    data = dat_df,
    aes(x = longitude, y = latitude, z = conifer_prop),
    bins = 50, # or whatever number of bins you want
    fun = mean, # or sum, median, etc.
    alpha = 1
  ) +
  labs(
    title = "Proportion of needle-leaved species",
  ) +
  scale_fill_viridis(option = "D", direction = 1) +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.title = element_blank(),
    axis.text = element_text(size = 8), 
    title = element_text(size = 10), 
    legend.position = c(0.05, 1), 
    legend.justification = c(0, 1),  
    legend.text = element_text(size = 10),
    legend.title = element_blank(),
    legend.key.size = unit(0.8, "cm"), 
    legend.direction = "vertical"
  ) +
  coord_sf(
    xlim = c(-15, 37.5),
    ylim = c(35, 75),
    expand = FALSE
  ) +
  annotation_scale(
    location = "br",         # bottom left
    width_hint = 0.2         
  ) +
  annotation_north_arrow(
    location = "br",
    which_north = "true",
    pad_x = unit(0.2, "in"),
    pad_y = unit(0.4, "in"),
    style = north_arrow_fancy_orienteering
  )

figure1_bl_nl

ggsave(
  path = here("03_results"),
  filename = paste0("Figure1_nl_proportion.png"),
  plot = figure1_bl_nl, 
  width = 130,
  height = 160,
  units = "mm",
  dpi = 300
)


# 3. Create productivity components maps ----------------------------------------

data <- species_list |> 
  list_rbind() |> 
  group_by(id) |> 
  summarise(growth = sum(growth_annual), 
            moccurrence = sum(mort_occ_annual), 
            mintensity = sum(mort_int_annual), 
            productivity = growth + moccurrence * mintensity, 
            longitude = mean(longitude), 
            latitude = mean(latitude)) |> 
  mutate(moccurrence = case_when(moccurrence > 0 ~ 1, T ~ 0))

table(data$moccurrence)

#Europe's base map 
europe <- ne_countries(
  scale = "medium", returnclass = "sf", continent = "Europe")

europe <- st_transform(europe, crs = 4326)

st_is_valid(europe)

create_productivity_components_maps <- function(var_name, var_title) {

  
  if (var_name == "mintensity") {
    data <- data |> filter(mintensity > 0)} else {
    data <- data
  }
  
  if (var_name == "moccurrence") {
    dat_df <- terra::vect(
      x = data,
      geom = c("longitude", "latitude"),
      crs = "EPSG:4326",
      keep = T
    )}else {
    
    ## Filter data taking only data lower than the 95% quantile
    q_high <- quantile(data[var_name], probs = 0.95, na.rm = T)
    
    # sp_data2 <- data_1[!(data_1[var_name] < q_low),]
    data2 <- data[(data[var_name] < q_high),]
    
    dat_df<- terra::vect(
      x = data2,
      geom = c("longitude", "latitude"),
      crs = "EPSG:4326",
      keep = T
    )
    
  }
  
  if (var_name == "productivity") {option = "D"} else {option = "B"}

  if (var_name == "moccurrence"){
    
    plot <- ggplot() +
      geom_spatvector(
        data = europe,
        fill = "grey80",
        color = "grey80"
      ) +
      stat_summary_hex(
        data = dat_df,
        aes(x = longitude, y = latitude, 
            z = .data[[var_name]]),
        bins = 40, 
        fun = mean, # or sum, median, etc.
        alpha = 1, 
      ) +
      scale_fill_viridis(option = option, direction = 1) +
      labs(
        title = var_title,
      ) +
      theme(
        panel.background = element_rect(fill = "white"),
        axis.title = element_blank(),
        axis.text = element_text(size = 8), 
        title = element_text(size = 10), 
        legend.position = c(0.05, 1),  
        legend.justification = c(0, 1),  
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "cm"), 
        legend.direction = "vertical"
      ) +
      coord_sf(
        xlim = c(-15, 37.5),
        ylim = c(35, 75),
        expand = FALSE
      ) +
      annotation_scale(
        location = "br",         # bottom left
        width_hint = 0.2         
      ) +
      annotation_north_arrow(
        location = "br",
        which_north = "true",
        pad_x = unit(0.2, "in"),
        pad_y = unit(0.4, "in"),
        style = north_arrow_fancy_orienteering
      )
    
  } else {
    
    plot <- ggplot() +
      geom_spatvector(
        data = europe,
        fill = "grey80",
        color = "grey80"
      ) +
      stat_summary_hex(
        data = dat_df,
        aes(x = longitude, y = latitude, z = .data[[var_name]]),
        bins = 50, 
        fun = median, # or sum, median, etc.
        alpha = 1
      ) +
      scale_fill_viridis(option = option, direction = 1) +
      labs(
        title = var_title,
      ) +
      theme(
        panel.background = element_rect(fill = "white"),
        axis.title = element_blank(),
        axis.text = element_text(size = 8), 
        title = element_text(size = 10), 
        legend.position = c(0.05, 1),  
        legend.justification = c(0, 1),  
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "cm"), 
        legend.direction = "vertical"
      ) +
      coord_sf(
        xlim = c(-15, 37.5),
        ylim = c(35, 75),
        expand = FALSE
      ) +
      annotation_scale(
        location = "br",         # bottom left
        width_hint = 0.2         
      ) +
      annotation_north_arrow(
        location = "br",
        which_north = "true",
        pad_x = unit(0.2, "in"),
        pad_y = unit(0.4, "in"),
        style = north_arrow_fancy_orienteering
      )
  }
  
  plot
  
  ggsave(
    path = here("03_results"),
    filename = paste0("Figure1_", var_name, ".png"),
    plot = plot,
    width = 130,
    height = 160,
    units = "mm",
    dpi = 300
  )
  
  return(plot)
  
}

var_name <- c("growth", "moccurrence", "mintensity", "productivity")
var_title <- c("Growth", "Mortality occurrence", "Mortality intensity", "Productivity")

map2(var_name, var_title, \(x, y) create_productivity_components_maps(x, y))


















