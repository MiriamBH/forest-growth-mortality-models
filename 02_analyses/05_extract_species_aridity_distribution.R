##############################################################################
##
## Script name: Extract species aridity distribution
## Purpose: Extract species aridity distribution to make species the predictions
##
## Input data: species_list.rds, cmi data from Chelsa, Species global distribution from 
##             Caudullo et al. (2019) database
## Output data: cmi_sp_distribution.rds 
## Output figures: none
##
## Date: 2026-02-18
##
## THIS SCRIPT CAN BE SKIPPED, THE OUTPUT DATA IS AVAILABLE in 01_data FOLDER
##
##############################################################################

# 1. Libraries and read data -------------------------------------------------

library(tidyverse)
library(here)
library(patchwork)
library(sf)
library(ggmap)
library(tidyterra)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)

species_list <- here("01_data", "species_list.rds") |> read_rds()

# 2. Obtain species-specific distribution limits -------------------------
### Species global distribution shapefiles
### (Caudullo et al. 2019): https://data.mendeley.com/datasets/hr5h2hcgg4/6
### Download and include in "01_data" folder to run the script

# New object with the names with a "_" to read the .shp archives
names <- names(species_list)
files <- tidytable::map_chr(names(species_list), \(x) sub(" ", "_", x))
files
length(files)

# Modify names to match the Caudullo files 
fs <- which(files == "Fagus_sylvatica")
files[fs] <- "Fagus_sylvatica_sylvatica"

Pn <- which(files == "Pinus_nigra")
files[Pn] <- "Pinus_nigra_salzmannii" # Olsson, S., Grivet, D., Cattonaro, F. et al. Evolutionary relevance of lineages in the European black pine (Pinus nigra) in the transcriptomic era. Tree Genetics & Genomes 16, 30 (2020). https://doi.org/10.1007/s11295-020-1424-8

Qi <- which(files == "Quercus_ilex")
files[Qi] <- "Quercus_ilex_rotundifolia" 

Ldd <- which(files == "Larix_decidua")
files[Ldd] <- "Larix_decidua_decidua"

# Creating a second line for Q.ilex, then we will join both again later
names[22] <- "Quercus ilex"
files[22] <- "Quercus_ilex_ilex" 

names[23] <- "Larix decidua"
files[23] <- "Larix_decidua_carpatica"

names
files

### Function to read the files 
read_sp_shp <- function(i){
  st_read(
    paste0("01_data/SpeciesDistributionEurope_CMD/", names[i], "/shapefiles/", files[i], "_plg.shp"))
}

### Shapefiles with the species distribution 
species_shp <- map(seq_along(files), \(x) read_sp_shp(x))

## 2.1. Extract coordinates from .shp of the limits for the different species #### 
### Species X Y coordinates 
species_coordinates <- map(species_shp, \(x) as.data.frame(st_coordinates(x)))

## 2.2. Restrict species distribution limits to our NFI longitudinal limits -------------------------------------------------------------
### We will take into account the latitudinal limits of the species in the world, but as 
### longitudinal limits we will just take the ones inside our database. 

### Latitudinal global distribution limits for each species 
min_lat <- map(species_coordinates, \(x) min(x$Y))
max_lat <- map(species_coordinates, \(x) max(x$Y))

### Longitudinal limits for each species restricted by the longitudinal limits of our database
lon <- map(species_coordinates, \(x) filter(x, X > -10 & X < 25))
min_lon <- map(lon, \(x) min(x$X))
max_lon <- map(lon, \(x) max(x$X))

### Creating a database with the maximum and minimum latitude and longitude for each species 
sp_limits <- data.frame(min_lon = unlist(min_lon), 
                        max_lon = unlist(max_lon), 
                        min_lat = unlist(min_lat), 
                        max_lat = unlist(max_lat)
) |> 
  mutate(species = names) |> 
  group_by(species) |> 
  summarise(
    min_lon = min(min_lon),
    max_lon = max(max_lon),
    min_lat = min(min_lat), 
    max_lat = max(max_lat)
  ) 

### List with the min and max coordinates for each species
list_sp_limits <- split(sp_limits, seq(nrow(sp_limits)))
names(list_sp_limits) <- names(species_list)

list_sp_limits[[1]]

## 2.3. Generate a grid of regularly distributed points considering all species limits -----------------
### Coordinates of the points to be generated inside the latitude and longitude limits and 
### how many points to generate considering the global limits of all species taken together
latitude_coords <- seq(min(sp_limits$min_lat), max(sp_limits$max_lat), length.out = 400) # lenght.out is the number of points generated, in the latitude axe in this case
longitude_coords <- seq(min(sp_limits$min_lon), max(sp_limits$max_lon), length.out = 400)

### Grid with all the points generated 
points <- expand.grid(longitude_coords, latitude_coords) 

### Create a variable to name each point (id)
points <- mutate(points, id = seq(1:nrow(points))) 
colnames(points) <- c("longitude", "latitude", "id")

### Create an sf object to visualize the points in a map 
points_sf <- st_as_sf(points, coords = c("longitude", "latitude"), 
                      crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

## 2.4. Cut the grid with the European coast ------------------------------------------------
europe <- ne_countries(scale = "medium", returnclass = "sf", continent = "Europe")

### Make sure they are in the same coordinate system
points_sf <- st_transform(points_sf, st_crs(europe))

### Europe map with the latitudinal and longitudinal limits of the generated grid 
points_sf_cut <- st_intersection(points_sf, europe)

# 3. Extract aridity (CMI) data for all the points of the grid -----------------------------------------
### Read names of the url of the CMI maps 
url <- read.table("01_data/Chelsa/cmi.txt")

# ### Download cmi .tif files from Chelsa database (Karger et al. 2021): https://envicloud.wsl.ch/#/?bucket=https%3A%2F%2Fos.unil.cloud.switch.ch%2Fchelsa02%2F&prefix=chelsa%2Fglobal%2Fmonthly%2F
# download_timeseries <- function(i){
#   options(timeout = 100000)
#   download.file(url = url[i, 1],
#                 destfile = paste0("01_data/Chelsa/downloads/cmi/", files[i]),
#                 mode = "wb"
#   )
# }
# 
# plan(multisession, workers = 4)
# 
# future_walk(seq_along(files), \(x) download_timeseries(i = x))

### Cut the url names to keep just the name of the variable, the year, month, etc. 
files <- str_split(url$V1, pattern = "_") |>
  map(\(x) paste0("CHELSA_", paste(x[3:6], collapse = "_"))) |>
  unlist()

### Function to read the .tif files 
read_raster_timeseries <- function(i){
  rast(paste0("01_data/Chelsa/downloads/cmi/", files[i]))
}

### .tif files 
timeseries <- map(seq_along(files), \(x) read_raster_timeseries(x))

### Coordinates of our points 
data_coords <- terra::vect(points_sf_cut)

### Function to extract the CMI value for each of our point coordinates 
extract_timeseries <- function(i){
  options(timeout = 100000)
  out <- terra::extract(timeseries[[i]], data_coords)
  outdf <- data.frame(value = out[, 2],
                      var = timeseries[[i]]@ptr[["names"]], ## cpp o ptr o pntr
                      id = points_sf_cut$id)
}

### Database with the CMI value for each of our points and for each month and year between 1980 and 2018
data_timeseries <- map_df(seq_along(files), \(x) extract_timeseries(x))

### Function to average the CMI values to obtain one mean value per point
get_climatologies <- function(mean_var) {
  data_timeseries |> 
    group_by(id) |> 
    summarise(
      !! mean_var := mean(value, na.rm = T)
    )
}

### CMI mean values for each coordinate 
data_climatology <- get_climatologies(mean_var = "cmi")

### Join the CMI values with the points database to keep the latitude and longitude variables  
cmi_data <- inner_join(data_climatology, points)

### save data
write_rds(x = cmi_data, file = "02_data/cmi_data_distribution.rds")

## 3.1. Create species specific CMI datasets ---------------------------------------------------------------------
### Cut the cmi_data for each species inside its range to then calculate the species specific quantiles
create_species_cmi_data <- function(sp_names) {
  
  sp_lim <- list_sp_limits[[sp_names]]
  
  db <- cmi_data |> 
    filter(
      longitude >= sp_lim$min_lon & 
        longitude <= sp_lim$max_lon & 
        latitude >= sp_lim$min_lat & 
        latitude <= sp_lim$max_lat
    )
  
  return(db)
}

cmi_sp_distribution <- map(names(species_list), \(x) create_species_cmi_data(x))
names(cmi_sp_distribution) <- names(species_list)

# 4. Save aridity species-specific distribution data -------------------------------------------------------------------------------------------------------
write_rds(x = cmi_sp_distribution, file = "01_data/cmi_sp_distribution.rds")

