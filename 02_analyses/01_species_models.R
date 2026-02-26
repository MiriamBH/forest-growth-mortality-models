##############################################################################
##
## Script name: species_models
## Purpose: Create and save the models 
##
## Input data: species_list.rds
## Output data: models.rds
##
## Date: 2026-02-12
##
##############################################################################

# 1. Libraries and read data -------------------------------------------------
library(tidyverse)
library(mgcv)
library(here)

## Read the data - a list with 21 dataframes (1 per species) with all data
## necessary for the models and figures
species_list <- here("01_data", "species_list.rds") |> read_rds()

# 2. Function to create specific species models ------------------------------
create_species_models <- function(sp_names) {
  
  sp <- species_list[[sp_names]] # dataframe with the data of 1 of the target species
  
  # if the species is present in only 1 country, we do not include country in the model
  if (length(unique(sp$country)) < 2) {
    
    mg <- gam(growth_annual ~
                log(ba_sp) + log(density) + harvest_reg +
                te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
              family = Gamma(link = "log"),
              data = sp,
              method = "REML")
    
    mo <- gam(mort_occ_annual ~ 
                log(ba_sp) + log(density) + harvest_reg +
                te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
              family = "binomial",
              data = sp,
              method = "REML")
    
    mi <- gam(mort_int_annual ~
                log(ba_sp) + log(density) + harvest_reg +
                te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
              family = Gamma(link = "log"),
              data = sp |> filter(mort_int_annual > 0), 
              method = "REML")
    
  } else {
    
    mg <- gam(growth_annual ~
                log(ba_sp) + log(density) + harvest_reg + country +
                te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
              family = Gamma(link = "log"),
              data = sp,
              method = "REML")
    
    mo <- gam(mort_occ_annual ~ 
                log(ba_sp) + log(density) + harvest_reg + country +
                te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
              family = "binomial",
              data = sp,
              method = "REML")
    
    mi <- gam(mort_int_annual ~
                log(ba_sp) + log(density) + harvest_reg + country +
                te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
              family = Gamma(link = "log"),
              data = sp |> filter(mort_int_annual > 0),
              method = "REML")
    
  }
  
  mg_s <- summary(mg)
  mo_s <- summary(mo)
  mi_s <- summary(mi) F0
  .
  
  return(list(growth = mg, 
              growth_s = mg_s,
              mort_occ = mo, 
              mort_occ_s = mo_s, 
              mort_int = mi,
              mort_int_s = mi_s))
  
}

## Run function 
sp_names <- names(species_list)

models <- purrr::map(sp_names, \(x) create_species_models(x))

names(models) <- names(species_list)

# 3. Save models -------------------------------------------------------------
write_rds(x = models, file = "01_data/models.rds")

