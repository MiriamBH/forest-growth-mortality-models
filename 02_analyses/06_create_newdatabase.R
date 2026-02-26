##############################################################################
##
## Script name: Create species predictions
## Purpose: Create species-specific predictions in each aridity and stand 
##          development level across the species VPD anomaly gradient 
##
## Input data: species_list.rds, cmi_sp_distribution.rds, models.rds
## Output data: newdatabase_sp_list.rds
## Output figures: none
##
## Date: 2026-02-25
##
##############################################################################

# 1. Libraries and read data -------------------------------------------------

library(tidyverse)
library(here)

species_list <- here("01_data", "species_list.rds") |> read_rds()
cmi_sp_distribution <- here("01_data", "cmi_sp_distribution.rds") |> read_rds()
models <- here("01_data", "models.rds") |> read_rds()

# 2. Create new database with fixed values for the species-specific predictions -------

create_sp_newdatabase <- function(sp_name) {
  
  sp <- species_list[[sp_name]] 

  # ARIDITY ####
  ## Values specific for the species based on their european distribution, 
  ## restricted for the longitudinal limits of our database (see script 04)
  q_arid <- quantile(cmi_sp_distribution[[sp_name]]$cmi, probs = 0.15, na.rm = T)
  q_mild <- quantile(cmi_sp_distribution[[sp_name]]$cmi, probs = 0.5, na.rm = T)
  q_wet <- quantile(cmi_sp_distribution[[sp_name]]$cmi, probs = 0.85, na.rm = T)
  
  # STAND DEVELOPMENT ####
  ## Calulate the quantiles for each species
  q_early <- quantile(sp$sdevelopment, probs = 0.25, na.rm = T)
  q_late <- quantile(sp$sdevelopment, probs = 0.75, na.rm = T)
  
  # VPD ANOMALY #### 
  ## For each species we calculate its minimum and maximum vpd values and space at regular intervals, to have the 
  ## same values for all species. Different species will have different ranges (larger or smaller databases), 
  ## but the same values so then 
  vpd_min <- round(min(sp$vpd_anomaly), 3)
  vpd_max <- round(max(sp$vpd_anomaly), 3)
  length <- length(seq(vpd_min, vpd_max, by = 0.001))
  length # this value is the one that controls the length of the database for each species. 
  
  ## OTHER VARIABLES ####
  ba_mean <- mean(sp$ba_sp)
  dens_mean <- mean(sp$density)
  harv_mean <- mean(sp$harvest_reg)
  
  ## COUNTRY ####  
  ## Different country per species, because for species specific models we can not include a country 
  ## in which the species is not present, we select the country in which the species was most abundant

  ## Number of plots in which the species is present in each country
  sp_count <- sp |> 
    count(country) 
  
  # Number of plots in which the species is present in each country for each aridity group
  # First group goes from lower aridity values to the percentile 5
  n_arid <- sp |> 
    filter(cmi < q_arid) |> 
    count(country) |> 
    rename(q_arid = n)
  
  # Second group goes from the percentile 5 to the mean (percentile 50)
  n_mild <- sp |> 
    filter(cmi > q_arid) |> 
    filter(cmi < q_wet) |> 
    count(country) |> 
    rename(q_mild = n)
  
  # Third group goes from the percentile 50 to the percentile 95
  n_wet <- sp |> 
    filter(cmi > q_wet) |> 
    count(country) |> 
    rename(q_wet = n) 
  
  # Countries with higher number of plots in which the species is present for each quantile. 
  # If the species is not present in one of the quantiles, then it just gives the other two. 
  ncountries <- full_join(sp_count, n_arid, by = "country") |> 
    full_join(n_mild) |> 
    full_join(n_wet) |> 
    mutate(across(where(is.integer), ~replace_na(., 0))) |> 
    pivot_longer(cols = 3:5, names_to = "quantile") |> 
    mutate(quantile_perc = value*100/n) |> 
    group_by(quantile) |> 
    filter(value == case_when(max(value) == 0 ~ NA, T ~ max(value))) |> 
    summarise(country = first(country), # because in the cases where value is the same for two countries, it will pick four countries instead of three (i.e. Quercus petraea)
              n = first(n), 
              value = first(value), 
              quantile_perc = first(quantile_perc)) |> 
    arrange(quantile) 
  
  ncountries
  
  # Names of the countries selected 
  vect_countries <- ncountries$country
  vect_countries
  
  third_country <- if (length(vect_countries) == 2) rep(vect_countries[2], 1) else vect_countries
  
  vect_countries2 <- if (length(vect_countries) == 2) append(vect_countries, third_country) else vect_countries
  
  ## New database
  # Differentent outputs (nrows) depending if 2 or 3 countries got selected for the species, 
  # which depends on the species being present on all three aridity quantiles or not 
  df_newdata <- data.frame(species = names(species_list[sp_name]),
                           funct_group = unique(species_list[[sp_name]]$funct_group),
                           ba_sp = rep(ba_mean, times = length*3), # 3 times because we have 3 aridiry quantiles
                           density = rep(dens_mean, times = length*3),
                           harvest_reg = rep(harv_mean, times = length*3),
                           vpd_anomaly = rep(seq(vpd_min, vpd_max, by = 0.001), times = 3),
                           country = rep(vect_countries2, each = length), 
                           cmi = rep(c(q_arid, q_mild, q_wet), each = length), 
                           # variable to see which quantiles of aridity have been picked for each species
                           quantiles = rep(c("q_arid", "q_mild", "q_wet"), each = length) 
  ) |>  
    mutate(quantiles = 
             case_when(quantiles == "5%" ~ "q_arid", T ~ case_when(
               quantiles == "50%" ~ "q_mild", T ~ case_when(
                 quantiles == "95%" ~ "q_wet", T ~ quantiles))), 
           aridity_factor = 
             case_when(quantiles == "q_arid" ~ "Arid", T ~ case_when(
               quantiles == "q_mild" ~ "Mild", T ~ "Wet")))
  
  # Join sdevelopment variable and create a categorical variable
  df_newdata_early <- df_newdata |>
    mutate(sdevelopment = unname(q_early), 
           sdev_factor = "Early") 
  
  df_newdata_late <- df_newdata |>
    mutate(sdevelopment = unname(q_late), 
           sdev_factor = "Late") 
  
  df_newdata_final <- rbind(df_newdata_early, df_newdata_late) 
  
  return(df_newdata_final)
  
}

sp_names <- names(species_list) # vector with the species names
newdatabase_sp_list <- map(sp_names, \(x) create_sp_newdatabase(x)) 
names(newdatabase_sp_list) <- names(species_list)

# 3. Save data -------------------------------------------------------------------
write_rds(newdatabase_sp_list, "01_data/newdatabase_sp_list.rds")

