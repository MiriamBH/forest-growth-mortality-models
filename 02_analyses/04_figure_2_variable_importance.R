##############################################################################
##
## Script name: Fig 2 variable importance
## Purpose: Create Variable importance figure 
##
## Input data: species_list.rds, CV_growth.rds, CV_mort_occ.rds, CV_mort_int.rds
## Output data: none 
## Output figures: Variable importance figures for growth, mortality 
##                 occurrence and mortality intensity (Figure 2 in the ms)
##
## Date: 2026-02-13
##
##############################################################################

# 1. Libraries and read data -------------------------------------------------
library(tidyverse)
library(here)
library(patchwork)

species_list <- here("01_data", "species_list.rds") |> read_rds()

crossg <- readRDS("01_data/CV_growth.rds") |> 
  rename(species = sp, 
         g.r2 = r2_mean, 
         g.nrmse = nrmse_mean, 
         g.nmae = nmae_mean) |> 
  select(species, g.r2, g.nrmse, g.nmae)

crossmo <- readRDS("01_data/CV_mort_occ.rds")|> 
  rename(species = sp, 
         mo.auc = auc_mean, 
         mo.sesp = SESP_mean) |> 
  select(species, mo.auc, mo.sesp)

crossmi <- readRDS("01_data/CV_mort_int.rds")|> 
  rename(species = sp,
         mi.r2 = r2_mean, 
         mi.nrmse = nrmse_mean, 
         mi.nmae = nmae_mean) |> 
  select(species, mi.r2, mi.nrmse, mi.nmae)

# 2. Create Explained deviance databases -----------------------------------------------
## Funtions to get the table for deviance explained by each variable in the model by the 
# leave-one-covariate-out (LOCO) approach 
calculate_de_growth <- function(sp_names) {
  
  sp <- species_list[[sp_names]] 

  ## growth models
  ### total model
  if (length(unique(sp$country)) < 2) {
    mg <- gam(growth_annual ~
                log(ba_sp) + log(density) + harvest_reg +
                te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
              family = Gamma(link = "log"),
              data = sp,
              method = "REML")
  } else {
    mg <- gam(growth_annual ~
                log(ba_sp) + log(density) + harvest_reg + country +
                te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
              family = Gamma(link = "log"),
              data = sp,
              method = "REML")
  }
  
  s <- summary(mg)
  de <- s$dev.expl
  
  ### model without stand density 
  if (length(unique(sp$country)) < 2) {
    mg_d <- gam(growth_annual ~
                  log(ba_sp) + harvest_reg +
                  te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
                family = Gamma(link = "log"),
                data = sp,
                method = "REML")
  } else {
    mg_d <- gam(growth_annual ~
                  log(ba_sp) + harvest_reg + country +
                  te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
                family = Gamma(link = "log"),
                data = sp,
                method = "REML")
  }
  
  s_d <- summary(mg_d)
  de_d <- s_d$dev.expl

  ### model without stand development (sdevelopment)
  if (length(unique(sp$country)) < 2) {
    mg_stdv <- gam(growth_annual ~
                    log(ba_sp) + log(density) + harvest_reg +
                    te(cmi, vpd_anomaly, bs = "tp"),
                  family = Gamma(link = "log"),
                  data = sp,
                  method = "REML")
  } else {
    mg_stdv <- gam(growth_annual ~
                    log(ba_sp) + log(density) + harvest_reg + country +
                    te(cmi, vpd_anomaly, bs = "tp"),
                  family = Gamma(link = "log"),
                  data = sp,
                  method = "REML")
  }
  
  s_stdv <- summary(mg_stdv)
  de_stdv <- s_stdv$dev.expl

  ### model without initial basal area 
  if (length(unique(sp$country)) < 2) {
    mg_ba <- gam(growth_annual ~
                   log(density) + harvest_reg +
                   te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
                 family = Gamma(link = "log"),
                 data = sp,
                 method = "REML")
  } else {
    mg_ba <- gam(growth_annual ~
                   log(density) + harvest_reg + country +
                   te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
                 family = Gamma(link = "log"),
                 data = sp,
                 method = "REML")
  }
  
  s_ba <- summary(mg_ba)
  de_ba <- s_ba$dev.expl

  ### model without mean climate (CMI)
  if (length(unique(sp$country)) < 2) {
    mg_cmi <- gam(growth_annual ~
                    log(ba_sp) + log(density) + harvest_reg + 
                    te(vpd_anomaly, sdevelopment, bs = "tp"),
                  family = Gamma(link = "log"),
                  data = sp,
                  method = "REML")
  } else {
    mg_cmi <- gam(growth_annual ~
                    log(ba_sp) + log(density) + harvest_reg + country +
                    te(vpd_anomaly, sdevelopment, bs = "tp"),
                  family = Gamma(link = "log"),
                  data = sp,
                  method = "REML")
  }
  
  s_cmi <- summary(mg_cmi)
  de_cmi <- s_cmi$dev.expl

  ### model without climate change (VPD anomaly)
  if (length(unique(sp$country)) < 2) {
    mg_cc <- gam(growth_annual ~
                   log(ba_sp) + log(density) + harvest_reg + 
                   te(cmi, sdevelopment, bs = "tp"),
                 family = Gamma(link = "log"),
                 data = sp,
                 method = "REML")
  } else {
    mg_cc <- gam(growth_annual ~
                   log(ba_sp) + log(density) + harvest_reg + country +
                   te(cmi, sdevelopment, bs = "tp"),
                 family = Gamma(link = "log"),
                 data = sp,
                 method = "REML")
  }
  
  s_cc<- summary(mg_cc)
  de_cc <- s_cc$dev.expl

  ### model without harvest regime
  if (length(unique(sp$country)) < 2) {
    mg_h <- gam(growth_annual ~
                  log(ba_sp) + log(density) + 
                  te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
                family = Gamma(link = "log"),
                data = sp,
                method = "REML")
  } else {
    mg_h <- gam(growth_annual ~
                  log(ba_sp) + log(density) + country +
                  te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
                family = Gamma(link = "log"),
                data = sp,
                method = "REML")
  }
  
  s_h <- summary(mg_h)
  de_h <- s_h$dev.expl

  # matrix with explained deviances
  dev_mat <- matrix(c(names(species_list[sp_names]), de, de_stdv, de_d, de_ba, de_cmi, de_cc, de_h), ncol = 8)
  colnames(dev_mat) = c("species", "de", "de_stdv", "de_d", "de_ba", "de_cmi", "de_cc", "de_h")
  
  # as dataframe
  table <- as.data.frame(dev_mat) |> 
    mutate(nrow = nrow(sp))
  
  return(table)
}

calculate_de_moccurrence <- function(sp_names) {
  
  sp <- species_list[[sp_names]] 
  
  ## mortality occurrence models
  ### total model
  if (length(unique(sp$country)) < 2) {
    mo <- gam(mort_occ_annual ~ 
                log(ba_sp) + log(density) + harvest_reg +
                te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
              family = "binomial",
              data = sp,
              method = "REML")
  } else {
    mo <- gam(mort_occ_annual ~ 
                log(ba_sp) + log(density) + harvest_reg + country +
                te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
              family = "binomial",
              data = sp,
              method = "REML")
  }
  
  s <- summary(mo)
  de <- s$dev.expl
  
  ### model without mean_dbh_large
  if (length(unique(sp$country)) < 2) {
    mo_stdv <- gam(mort_occ_annual ~ 
                    log(ba_sp) + log(density) + harvest_reg +
                    te(cmi, vpd_anomaly, bs = "tp"),
                  family = "binomial",
                  data = sp,
                  method = "REML")
  } else {
    mo_stdv <- gam(mort_occ_annual ~ 
                    log(ba_sp) + log(density) + harvest_reg + country +
                    te(cmi, vpd_anomaly, bs = "tp"),
                  family = "binomial",
                  data = sp,
                  method = "REML")
  }
  
  s_stdv <- summary(mo_stdv)
  de_stdv <- s_stdv$dev.expl
  
  ### without density 
  if (length(unique(sp$country)) < 2) {
    mo_d <- gam(mort_occ_annual ~ 
                  log(ba_sp) + harvest_reg +
                  te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
                family = "binomial",
                data = sp,
                method = "REML")
  } else {
    mo_d <- gam(mort_occ_annual ~ 
                  log(ba_sp) + harvest_reg + country +
                  te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
                family = "binomial",
                data = sp,
                method = "REML")
  }
  
  s_d <- summary(mo_d)
  de_d <- s_d$dev.expl
  
  ### model without initial basal area 
  if (length(unique(sp$country)) < 2) {
    mo_ba <- gam(mort_occ_annual ~ 
                   log(density) + harvest_reg +
                   te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
                 family = "binomial",
                 data = sp,
                 method = "REML")
  } else {
    mo_ba <- gam(mort_occ_annual ~ 
                   log(density) + harvest_reg + country +
                   te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
                 family = "binomial",
                 data = sp,
                 method = "REML")
  }
  
  s_ba <- summary(mo_ba)
  de_ba <- s_ba$dev.expl
  
  ### without mean climate (CMI)
  if (length(unique(sp$country)) < 2) {
    mo_cmi <- gam(mort_occ_annual ~ 
                    log(ba_sp) + log(density) + harvest_reg +
                    te(vpd_anomaly, sdevelopment, bs = "tp"),
                  family = "binomial",
                  data = sp,
                  method = "REML")
  } else {
    mo_cmi <- gam(mort_occ_annual ~ 
                    log(ba_sp) + log(density) + harvest_reg + country +
                    te(vpd_anomaly, sdevelopment, bs = "tp"),
                  family = "binomial",
                  data = sp,
                  method = "REML")
  }
  
  s_cmi <- summary(mo_cmi)
  de_cmi <- s_cmi$dev.expl
  
  ### model without climate change
  if (length(unique(sp$country)) < 2) {
    mo_cc <- gam(mort_occ_annual ~ 
                   log(ba_sp) + log(density) + harvest_reg +
                   te(cmi, sdevelopment, bs = "tp"),
                 family = "binomial",
                 data = sp,
                 method = "REML")
  } else {
    mo_cc <- gam(mort_occ_annual ~ 
                   log(ba_sp) + log(density) + harvest_reg + country +
                   te(cmi, sdevelopment, bs = "tp"),
                 family = "binomial",
                 data = sp,
                 method = "REML")
  }
  
  s_cc <- summary(mo_cc)
  de_cc <- s_cc$dev.expl
  
  ### model without harvest
  if (length(unique(sp$country)) < 2) {
    mo_h <- gam(mort_occ_annual ~ 
                  log(ba_sp) + log(density) +
                  te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
                family = "binomial",
                data = sp,
                method = "REML")
  } else {
    mo_h <- gam(mort_occ_annual ~ 
                  log(ba_sp) + log(density) + country +
                  te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
                family = "binomial",
                data = sp,
                method = "REML")
  }
  
  s_h <- summary(mo_h)
  de_h <- s_h$dev.expl
  
  # long format
  dev_mat <- matrix(c(names(species_list[sp_names]), de, de_stdv, de_d, de_ba, de_cmi, de_cc, de_h), ncol = 8)
  colnames(dev_mat) = c("species", "de", "de_stdv", "de_d", "de_ba", "de_cmi", "de_cc", "de_h")
  
  # as dataframe
  table <- as.data.frame(dev_mat) |> 
    mutate(nrow = nrow(sp))
  
  return(table)
}

calculate_de_mintensity <- function(sp_names) {
  
  sp <- species_list[[sp_names]] 

  ## mortality intensity models
  ### full model
  if (length(unique(sp$country)) < 2) {
    mi <- gam(mort_int_annual ~
                log(ba_sp) + log(density) + harvest_reg +
                te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
              family = Gamma(link = "log"),
              data = sp |> filter(mort_int_annual > 0),
              method = "REML")
  } else {
    mi <- gam(mort_int_annual ~
                log(ba_sp) + log(density) + harvest_reg + country +
                te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
              family = Gamma(link = "log"),
              data = sp |> filter(mort_int_annual > 0),
              method = "REML")
  }
  
  s <- summary(mi)
  de <- s$dev.expl

  ### model without stand development (sdevelopment)
  if (length(unique(sp$country)) < 2) {
    mi_stdv <- gam(mort_int_annual ~
                    log(ba_sp) + log(density) + harvest_reg +
                    te(cmi, vpd_anomaly, bs = "tp"),
                  family = Gamma(link = "log"),
                  data = sp |> filter(mort_int_annual > 0),
                  method = "REML")
  } else {
    mi_stdv <- gam(mort_int_annual ~
                    log(ba_sp) + log(density) + harvest_reg + country +
                    te(cmi, vpd_anomaly, bs = "tp"),
                  family = Gamma(link = "log"),
                  data = sp |> filter(mort_int_annual > 0),
                  method = "REML")
  }
  
  s_stdv <- summary(mi_stdv)
  de_stdv <- s_stdv$dev.expl

  ### without density
  if (length(unique(sp$country)) < 2) {
    mi_d <- gam(mort_int_annual ~
                  log(ba_sp) + harvest_reg +
                  te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
                family = Gamma(link = "log"),
                data = sp |> filter(mort_int_annual > 0),
                method = "REML")
  } else {
    mi_d <- gam(mort_int_annual ~
                  log(ba_sp) + harvest_reg + country +
                  te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
                family = Gamma(link = "log"),
                data = sp |> filter(mort_int_annual > 0),
                method = "REML")
  }
  
  s_d <- summary(mi_d)
  de_d <- s_d$dev.expl

  ### model without initial basal area 
  if (length(unique(sp$country)) < 2) {
    mi_ba <- gam(mort_int_annual ~
                   log(density) + harvest_reg +
                   te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
                 family = Gamma(link = "log"),
                 data = sp |> filter(mort_int_annual > 0),
                 method = "REML")
  } else {
    mi_ba <- gam(mort_int_annual ~
                   log(density) + harvest_reg + country +
                   te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
                 family = Gamma(link = "log"),
                 data = sp |> filter(mort_int_annual > 0),
                 method = "REML")
  }
  
  s_ba <- summary(mi_ba)
  de_ba <- s_ba$dev.expl

  ### model without mean climate (CMI)
  if (length(unique(sp$country)) < 2) {
    mi_cmi <- gam(mort_int_annual ~
                    log(ba_sp) + log(density) + harvest_reg +
                    te(vpd_anomaly, sdevelopment, bs = "tp"),
                  family = Gamma(link = "log"),
                  data = sp |> filter(mort_int_annual > 0),
                  method = "REML")
  } else {
    mi_cmi <- gam(mort_int_annual ~
                    log(ba_sp) + log(density) + harvest_reg + country +
                    te(vpd_anomaly, sdevelopment, bs = "tp"),
                  family = Gamma(link = "log"),
                  data = sp |> filter(mort_int_annual > 0),
                  method = "REML")
  }
  
  s_cmi <- summary(mi_cmi)
  de_cmi <- s_cmi$dev.expl

  ### model without climate change (VPD anomaly)
  if (length(unique(sp$country)) < 2) {
    mi_cc <- gam(mort_int_annual ~
                   log(ba_sp) + log(density) + harvest_reg +
                   te(cmi, sdevelopment, bs = "tp"),
                 family = Gamma(link = "log"),
                 data = sp |> filter(mort_int_annual > 0),
                 method = "REML")
  } else {
    mi_cc <- gam(mort_int_annual ~
                   log(ba_sp) + log(density) + harvest_reg + country +
                   te(cmi, sdevelopment, bs = "tp"),
                 family = Gamma(link = "log"),
                 data = sp |> filter(mort_int_annual > 0),
                 method = "REML")
  }
  
  s_cc <- summary(mi_cc)
  de_cc <- s_cc$dev.expl

  ### model without harvest 
  if (length(unique(sp$country)) < 2) {
    mi_h <- gam(mort_int_annual ~
                  log(ba_sp) + log(density) +
                  te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
                family = Gamma(link = "log"),
                data = sp |> filter(mort_int_annual > 0),
                method = "REML")
  } else {
    mi_h <- gam(mort_int_annual ~
                  log(ba_sp) + log(density) + country +
                  te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
                family = Gamma(link = "log"),
                data = sp |> filter(mort_int_annual > 0),
                method = "REML")
  }
  
  s_h <- summary(mi_h)
  de_h <- s_h$dev.expl

  # matrix with explained deviances
  dev_mat <- matrix(c(names(species_list[sp_names]), de, de_stdv, de_d, de_ba, de_cmi, de_cc, de_h), ncol = 8)
  colnames(dev_mat) = c("species", "de", "de_stdv", "de_d", "de_ba", "de_cmi", "de_cc", "de_h")

  table <- as.data.frame(dev_mat) |> 
    mutate(nrow = nrow(sp))

  return(table)
}

## Run functions 
sp_v <- sort(names(species_list)) # vector with the species names

table_de_growth <- purrr::map(sort(sp_v), \(x) calculate_de_growth(x)) |>  
  list_rbind()

table_de_moccurrence <- purrr::map(sort(sp_v), \(x) calculate_de_moccurrence(x)) |>  
  list_rbind()

table_de_mintensity <- purrr::map(sort(sp_v), \(x) calculate_de_mintensity(x)) |>  
  list_rbind()

## Save databases 
write_rds(x = table_de_growth, file = "01_data/table_de_growth.rds")
write_rds(x = table_de_moccurrence, file = "01_data/table_de_moccurrence.rds")
write_rds(x = table_de_mintensity, file = "01_data/table_de_mintensity.rds")

## Read databases 
# table_de_growth <- readRDS("01_data/table_de_growth.rds")
# table_de_moccurrence <- readRDS("01_data/table_de_moccurrence.rds")
# table_de_mintensity <- readRDS("01_data/table_de_mintensity.rds")

# 3. Function to calculate variable_importance variable using the dev tables ----------------------------
calculate_var_imp_dev <- function(table) {
  
  ## All variables 
  all_tables <- table |> 
    mutate(across(c(de, de_stdv, de_d, de_ba, de_cmi, de_cc, de_h), as.numeric),
           # deviance of the total model (de) - deviance of the model without the variable (de_stdv)
           # to see how much deviance is lost when the variable is taken out of the model.
           # The deviance lost is the variable importance 
           std = de - de_stdv, 
           d = de - de_d,
           ba = de - de_ba, 
           cmi= de - de_cmi,
           cc = de - de_cc,
           h = de - de_h) |>  
    select(!c(de_stdv, de_d, de_ba, de_cmi, de_cc, de_h)) |> 
    pivot_longer(cols = c("stdv", "d", "ba", "cmi", "cc", "h"), names_to = "model") |> 
    rename(dev_expl = value) |> 
    mutate(dev_expl_perc = dev_expl*100/de)
  
  return(all_tables)
  
}

## Run function 
growth_table <- calculate_var_imp_dev(table_de_growth) 
moccurrence_table <- calculate_var_imp_dev(table_de_moccurrence)
mintensity_table <- calculate_var_imp_dev(table_de_mintensity)

# 4. Join cross-val data for weighting ----------------------------------------------------
growth_table <- growth_table |> left_join(crossg)
moccurrence_table <- moccurrence_table |> left_join(crossmo)
mintensity_table <- mintensity_table |> left_join(crossmi)

# 5. Join functional group data -----------------------------------------------

### Assign a funct_group to each species 
species <- names(species_list) # a vector with the name of the species 

funct_group <- map_chr(species, \(x) unique(species_list[[x]]$funct_group)) # a vector with the type of the species 

funct_group_df <- data.frame(species = species, funct_group = funct_group)  # create a dataframe with both vectors 

### Join with the tables 
tables <- list(growth_table, moccurrence_table, mintensity_table)
names(tables) <- c("growth_table", "moccurrence_table", "mintensity_table")

new_tables <- map(tables, \(x) left_join(x, funct_group_df))

### Modify the tables to include type 
growth_table <- new_tables$growth_table 
moccurrence_table <- new_tables$moccurrence_table
mintensity_table <- new_tables$mintensity_table

# 6. Create and safe graphic -----------------------------------------------------------
create_figure_2 <- function(variable){
  
  tau <- 0.4
  model_levels = rev(c("ba", "d", "stdv", "h", "cmi", "cc"))
  legend_labels = rev(c("Species basal area", "Stand density", "Stand development","Harvest regime", "Aridity", "VPD anomaly"))
  colors = rev(c("darkolivegreen", "darkolivegreen1", "darkolivegreen3","#75801f","#6633CC","#9999FF"))
  
  if (variable == "g") {
    
    table = growth_table
    title = "Growth"
    # include the weighting by the predictive capacity of each species model
    table_w <- table |> 
      group_by(species) |> 
      mutate(w = exp(-g.nmae/tau)/sum(exp(-g.nmae/tau)), na.rm = T) |> 
      mutate(w = case_when(is.na(w) ~ 0, T ~ w)) |>
      ungroup() |> 
      group_by(model) |> 
      summarise(dev_expl_perc_p = mean(dev_expl_perc * w))|> 
      mutate(dev_expl_perc_p = case_when(dev_expl_perc_p < 0 ~ 0, T ~ dev_expl_perc_p))
    
  } else {
    
    if (variable == "mi") {
      
      table = mintensity_table
      title = "Mortality intensity"
      table_w <- table |> 
        group_by(species) |> 
        mutate(w = exp(-mi.nmae/tau)/sum(exp(-mi.nmae/tau)), na.rm = T) |> 
        mutate(w = case_when(is.na(w) ~ 0, T ~ w)) |>
        ungroup() |> 
        group_by(model) |> 
        summarise(dev_expl_perc_p = mean(dev_expl_perc * w)) |> 
        mutate(dev_expl_perc_p = case_when(dev_expl_perc_p < 0 ~ 0, T ~ dev_expl_perc_p))
      
    } else {
      
      table = moccurrence_table 
      title = "Mortality occurrence"
      table_w <- table |> 
        group_by(species) |> 
        mutate(w = exp(mo.auc/tau)/sum(exp(mo.auc/tau)), na.rm = T) |> 
        mutate(w = case_when(is.na(w) ~ 0, T ~ w)) |>
        ungroup() |> 
        group_by(model) |> 
        summarise(dev_expl_perc_p = mean(dev_expl_perc * w)) |> 
        mutate(dev_expl_perc_p = case_when(dev_expl_perc_p < 0 ~ 0, T ~ dev_expl_perc_p))
      
    }
    
  }
 
  gg_sp <- table |> 
    
    group_by(species) |> 
    
    filter(dev_expl_perc > 0) |>
    
    ungroup() |> 
    
    ggplot(aes(y = dev_expl_perc, 
               x = fct_relevel(species, rev(sort(unique(species)))), 
               fill = fct_relevel(model, model_levels))) +
    
    geom_bar(position = "fill", stat = "identity") +
    
    geom_abline(slope = 0, intercept = 100, col = "#999999", lty = 2) +
    
    labs(title = title,
         subtitle = "(a)",
         x = element_blank(), y = "Explained deviance (%)"
    ) +
    
    scale_fill_manual(labels = legend_labels, values = colors) +
    
    theme(plot.title = element_text(hjust = 0.5, size = 12, vjust = 0.1),
          plot.subtitle = element_text(hjust = 0, size = 8, vjust = 0.1),
          legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 10, face = "italic"),
          axis.text.x = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y = element_line(colour = "lightgrey", linewidth = 0.1),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(fill = "white", color = "black"),
          strip.text.y.left = element_text(angle = 90, 
                                           color = "grey20",
                                           face = "bold",
                                           size = 9),
    ) +
    
    guides(fill = guide_legend(title = element_blank(), reverse = TRUE)) +
    
    facet_grid(funct_group ~ ., scales = "free_y", space = "free_y", switch = "y",
               labeller = as_labeller(c("BL" = "Broad-leaved", "NL" = "Needle-leaved"))) +
    coord_flip()
  
  gg_sp
  
  gg_func <- table |>
    
    group_by(species) |>
    
    filter(dev_expl_perc > 0) |>
    
    mutate(type = case_when(funct_group == "BL" ~ "Broad-leaved", T ~ "Needle-leaved")) |>
    
    ggplot(aes(y = dev_expl_perc,
               x = fct_relevel(type, c("Needle-leaved", "Broad-leaved")),
               fill = fct_relevel(model, model_levels))) +
    
    geom_bar(position = "fill", stat = "identity") +
    
    geom_abline(slope = 0, intercept = 100, col = "#999999", lty = 2) +
    
    labs(title = element_blank(),
         x = element_blank(), y = "Variable importance"
    ) +
    
    scale_fill_manual(labels = legend_labels, values = colors) +
    
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size = 9),
      axis.title.x = element_text(size = 9),
      axis.text.y = element_text(size = 9,
                                 face = "bold",
                                 colour = "grey20"),
      axis.text.x = element_text(size = 8),
      panel.background = element_blank(),
      panel.grid.major.y = element_line(colour = "lightgrey", linewidth = 0.1),
      axis.ticks.y = element_blank(), 
      plot.title = element_text(hjust = 0, size = 8)
    ) +
    
    labs(title = "(b)") +
    
    guides(fill = guide_legend(title = element_blank(), reverse = TRUE)) +
    
    coord_flip()
  
  gg_func
  
  gg <- gg_sp / gg_func +
    theme(plot.margin = margin(0, 0, -5, 0),
          axis.text.y = element_blank()) +
    plot_layout(heights = c(7,0.75)) 
  
  
  ggsave(
    path = here("03_results"),
    filename = paste0("Figure2_vi_", variable, ".png"),
    plot = gg,
    width = 140,
    height = 130,
    units = "mm",
    dpi = 300
  )
  
  return(gg)
  
}

# Create and safe figures
variables <- c("g", "mi", "mo")
gg <- map(variables, \(x) create_figure_2(x))

