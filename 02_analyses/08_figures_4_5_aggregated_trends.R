##############################################################################
##
## Script name: Fig 4 and 5 aggregated trends
## Purpose: Create Figures 4 and 5 (aggregated species and functional 
##          group trends)
##
## Input data: species_list.rds, models.rds, newdatabase_sp_list.rds 
## Output data: none 
## Output figures: Aggregated species and functional group 
##                trends figures (Figures 4 and 5 in ms)
##
## Date: 2026-02-25
##
##############################################################################

# 1. Libraries and read data -------------------------------------------------

library(tidyverse)
library(here)
library(patchwork)
library(mgcv)

species_list <- here("01_data", "species_list.rds") |> read_rds()
newdatabase_sp_list <- here("01_data", "newdatabase_sp_list.rds") |> read_rds()
models <- here("01_data", "models.rds") |> read_rds()

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

# 2. Create species-specific predictions for the VPD anomaly range of the species -------------
create_sp_predictions <- function(sp_name) {
  
  data <- newdatabase_sp_list[[sp_name]]
  
  mg <- models[[sp_name]][["growth"]]
  mo <- models[[sp_name]][["mort_occ"]]
  mi <- models[[sp_name]][["mort_int"]]
  
  # Add prediction vairables (response gains variables predicted by the model for the explanatory variables we just created)
  pred_g <- predict(mg, newdata = data, type = "response", se.fit = TRUE)
  pred_mo <- predict(mo, newdata = data, type = "response", se.fit = TRUE)
  pred_mi <- predict(mi, newdata = data, type = "response", se.fit = TRUE)
  
  data$pred_g <- as.numeric(pred_g$fit)
  data$dev_g <- models[[sp_name]][["growth_s"]]$dev.exp
  data$se_g <- as.numeric(pred_g$se.fit)
  
  data$pred_mo <- as.numeric(pred_mo$fit)
  data$dev_mo <- models[[sp_name]][["mort_occ_s"]]$dev.exp
  data$se_mo <- as.numeric(pred_mo$se.fit)
  
  data$pred_mi <- as.numeric(pred_mi$fit)
  data$dev_mi <- models[[sp_name]][["mort_int_s"]]$dev.exp
  data$se_mi <- as.numeric(pred_mi$se.fit)
  
  
  return(data)
}

sp_names <- names(species_list) # vector with the species names
newdatabase_sp_predictions_list <- map(sp_names, \(x) create_sp_predictions(x)) 
names(newdatabase_sp_predictions_list) <- sp_names

# 3. Dataframe with all the species growth and mortality aggregated ----------------
newdata <- list_rbind(newdatabase_sp_predictions_list) |> 
  left_join(crossg) |> 
  left_join(crossmo) |> 
  left_join(crossmi) |> 
  mutate(across(c(species, country, sdev_factor), as.factor)) 

# 4. Set VPD range as the range shared in each given point for at least 75% of the species ---------------------------------
## All sp together ####
n_species_vpd <- newdata |> 
  filter(sdev_factor == "Early") |> 
  filter(aridity_factor == "Arid") |> 
  mutate(vpd_anomaly = round(vpd_anomaly, 3)) |> 
  group_by(vpd_anomaly) |> 
  summarise(n_species = n()) |> 
  ungroup() 

# 75% of 21 species is 15 species 
n_species_vpd |> filter(n_species == 15) # value max VPD = 0.455
n_species_vpd |> filter(n_species == 16) # values min VPD = 0.082 

newdata_vpd_reduced75 <- newdata |> 
  filter(vpd_anomaly >= 0.082) |> 
  filter(vpd_anomaly <= 0.455)

## Per functional group ####
n_species_vpd_fg <- newdata |> 
  filter(sdev_factor == "Early") |> 
  filter(aridity_factor == "Arid") |> 
  mutate(vpd_anomaly = round(vpd_anomaly, 3)) |> 
  group_by(vpd_anomaly, funct_group) |> 
  summarise(n_species = n()) |> 
  ungroup() 

## 75% of BL species is 10 sp
n_species_vpd_fg |> filter(funct_group == "BL") |> filter(n_species == 10) # value max VPD = 0.459
n_species_vpd_fg |> filter(funct_group == "BL") |> filter(n_species == 11) # value min VPD = 0.082

## 75% of NL species is 5 sp
n_species_vpd_fg |> filter(funct_group == "NL") |> filter(n_species == 5) |> 
  print(n = 30) # value min VPD = 0.066 and max VPD = 0.407

newdata_vpd_reduced75_fg <- newdata |>
  filter(
    (funct_group == "NL" & vpd_anomaly >= 0.066 & vpd_anomaly <= 0.407) |
      (funct_group != "NL" & vpd_anomaly >= 0.082 & vpd_anomaly <= 0.459)
  )

# 5. Mean weighted predictions for all the species together (with their confidence intervals) -------------------------------------

## generate mean predictions weighting each species by the predictive performance of its model, 
## to penalize by the worse models and calculate confidence intervals
tau <- 0.4

newdata_mean_sp <- newdata_vpd_reduced75 |> 
  mutate(aridity_factor = as.factor(aridity_factor), 
         vpd_anomaly_factor = as.factor(vpd_anomaly)) |> 
  select(aridity_factor, sdev_factor, vpd_anomaly_factor, vpd_anomaly, 
         pred_g, pred_mo, pred_mi,
         se_g, se_mo, se_mi, dev_g, dev_mo, dev_mi, 
         g.nmae, mi.nmae, mo.auc) |> 
  group_by(aridity_factor, sdev_factor, vpd_anomaly_factor) |> 
  mutate(
    w.g = exp(-g.nmae/tau)/sum(exp(-g.nmae/tau)), 
    w.mi = exp(-mi.nmae/tau)/sum(exp(-mi.nmae/tau)),
    w.mo = exp(mo.auc/tau)/sum(exp(mo.auc/tau))
  ) |> 
  summarise(
    wmean.pred.g = sum(pred_g * w.g, na.rm = T), 
    wmean.pred.mo = sum(pred_mo * w.mo, na.rm = T), 
    wmean.pred.mi = sum(pred_mi * w.mi, na.rm = T), 
    wse_g = sqrt(sum((w.g^2) * (se_g^2))),
    wse_mo = sqrt(sum((w.mo^2) * (se_mo^2))),
    wse_mi = sqrt(sum((w.mi^2) * (se_mi^2))),
    vpd_anomaly = round(first(vpd_anomaly), 3)
  ) |> 
  mutate(
    upper.ci.g = wmean.pred.g + 1.96 * wse_g, 
    upper.ci.mo = wmean.pred.mo + 1.96 * wse_mo, 
    upper.ci.mi = wmean.pred.mi + 1.96 * wse_mi, 
    lower.ci.g = wmean.pred.g - 1.96 * wse_g, 
    lower.ci.mo = wmean.pred.mo - 1.96 * wse_mo, 
    lower.ci.mi = wmean.pred.mi - 1.96 * wse_mi,
  ) |> 
  ungroup()

# Group by functional group (for Figure 5)
newdata_mean_fg <- newdata_vpd_reduced75_fg |> 
  mutate(aridity_factor = as.factor(aridity_factor), 
         vpd_anomaly_factor = as.factor(vpd_anomaly)) |> 
  select(aridity_factor, sdev_factor, vpd_anomaly_factor, vpd_anomaly, funct_group, 
         pred_g, pred_mo, pred_mi,
         se_g, se_mo, se_mi, dev_g, dev_mo, dev_mi, 
         g.nmae, mi.nmae, mo.auc) |> 
  group_by(aridity_factor, sdev_factor, vpd_anomaly_factor, funct_group) |> 
  mutate(
    w.g = exp(-g.nmae/tau)/sum(exp(-g.nmae/tau)), 
    w.mi = exp(-mi.nmae/tau)/sum(exp(-mi.nmae/tau)),
    w.mo = exp(mo.auc/tau)/sum(exp(mo.auc/tau))
  ) |> 
  summarise(
    wmean.pred.g = sum(pred_g * w.g, na.rm = T), 
    wmean.pred.mo = sum(pred_mo * w.mo, na.rm = T), 
    wmean.pred.mi = sum(pred_mi * w.mi, na.rm = T), 
    wse_g = sqrt(sum((w.g^2) * (se_g^2), na.rm = T)),
    wse_mo = sqrt(sum((w.mo^2) * (se_mo^2), na.rm = T)),
    wse_mi = sqrt(sum((w.mi^2) * (se_mi^2), na.rm = T)),
    vpd_anomaly = round(first(vpd_anomaly), 3)
  ) |> 
  mutate(
    upper.ci.g = wmean.pred.g + 1.96 * wse_g, 
    upper.ci.mo = wmean.pred.mo + 1.96 * wse_mo, 
    upper.ci.mi = wmean.pred.mi + 1.96 * wse_mi, 
    lower.ci.g = wmean.pred.g - 1.96 * wse_g, 
    lower.ci.mo = wmean.pred.mo - 1.96 * wse_mo, 
    lower.ci.mi = wmean.pred.mi - 1.96 * wse_mi,
  ) |> 
  ungroup()


# 6. Create and save figure 4, aggregating all species -------------------------------

## Generate smooth trends for visualization and better overlap between the trend and the confidence interval
smooth_data <- newdata_mean_sp |>
  group_by(aridity_factor, sdev_factor) |>
  mutate(
    wmean.pred.g.smooth = predict(smooth.spline(vpd_anomaly, wmean.pred.g, spar = 1), x = vpd_anomaly)$y,
    upper.ci.g.smooth = predict(smooth.spline(vpd_anomaly, upper.ci.g, spar = 1), x = vpd_anomaly)$y,
    lower.ci.g.smooth = predict(smooth.spline(vpd_anomaly, lower.ci.g, spar = 1), x = vpd_anomaly)$y,
    wmean.pred.mo.smooth = predict(smooth.spline(vpd_anomaly, wmean.pred.mo, spar = 1), x = vpd_anomaly)$y,
    upper.ci.mo.smooth = predict(smooth.spline(vpd_anomaly, upper.ci.mo, spar = 1), x = vpd_anomaly)$y,
    lower.ci.mo.smooth = predict(smooth.spline(vpd_anomaly, lower.ci.mo, spar = 1), x = vpd_anomaly)$y,
    wmean.pred.mi.smooth = predict(smooth.spline(vpd_anomaly, wmean.pred.mi, spar = 1), x = vpd_anomaly)$y,
    upper.ci.mi.smooth = predict(smooth.spline(vpd_anomaly, upper.ci.mi, spar = 1), x = vpd_anomaly)$y,
    lower.ci.mi.smooth = predict(smooth.spline(vpd_anomaly, lower.ci.mi, spar = 1), x = vpd_anomaly)$y
  ) |>
  ungroup() |>
  select(aridity_factor, sdev_factor, vpd_anomaly,
         wmean.pred.g.smooth, upper.ci.g.smooth, lower.ci.g.smooth,
         wmean.pred.mo.smooth, upper.ci.mo.smooth, lower.ci.mo.smooth,
         wmean.pred.mi.smooth, upper.ci.mi.smooth, lower.ci.mi.smooth) |>
  pivot_longer(
    cols = -c(aridity_factor, sdev_factor, vpd_anomaly),  # todo lo demás se transforma
    names_to = c(".value", "response"),             # .value crea columnas mean.pred, upper, lower
    names_pattern = "(wmean\\.pred|upper\\.ci|lower\\.ci)\\.(.+)\\.smooth"
    )
  
  # main graphic 
  gg_main <- ggplot(
    smooth_data |> 
      mutate(
        response = factor(
          response,
          levels = c("g", "mo", "mi"),
          labels = c(
            "Growth",
            "Mortality occurrence",
            "Mortality intensity"
          )
        ), 
        sdev_factor = factor(
          sdev_factor,
          levels = c("Early", "Late"), 
          labels = c("Early stand development", "Late stand development"))
      ) |> 
      filter(aridity_factor != "Mild"),
    aes(x = vpd_anomaly, y = wmean.pred, color = aridity_factor)
  ) +
    
    geom_line(linewidth = 0.7) +
    geom_ribbon(
      aes(ymin = lower.ci, ymax = upper.ci, fill = aridity_factor),
      alpha = 0.4,
      color = NA,
      show.legend = FALSE
    ) +
    
    scale_color_manual(values = c(
      "Wet"  = "skyblue3",
      "Arid" = "darkorange3"
    )) +
    scale_fill_manual(values = c(
      "Wet"  = "skyblue3",
      "Arid" = "darkorange3"
    )) +
    
    facet_grid(
      sdev_factor ~ response,
      scales = "free_y"
    ) +
    
    scale_y_continuous(
      name = expression(
        paste(
          "Growth or M. intensity (",
          m^2, " ", ha^-1, " ", yr^-1, " or Mortality occurrence (prob.))"
        )
      )
    ) +
    scale_x_continuous(limits = c(0.082, 0.455)) +
    
    theme_minimal() +
    theme(
      strip.background = element_blank(),
      strip.placement = "inside", 
      strip.text = element_text(size = 10),
      axis.text = element_text(color = "black", size = 8),
      axis.ticks = element_blank(),
      axis.line = element_line(color = "lightgrey"),
      panel.background = element_rect(
        fill = "white",
        colour = "lightgrey",
        linetype = "solid"
      ),
      legend.position = "bottom",
      legend.key = element_blank(),
      legend.background = element_blank(),
      legend.text = element_text(size = 10), 
      plot.background = element_rect("white")
    ) +
    
    labs(x = "VPD anomaly") +
    guides(color = guide_legend(title = NULL))
  
  gg_main
  
  ggsave(
    path = here("03_results"),
    filename = paste0("Figure4_gm_trends.png"),
    plot = gg_main,
    width = 135,
    height = 110,
    dpi = 300,
    units = "mm"
  )
  
  # Number of species bar 
  gg_tile <- ggplot(
    n_species_vpd |> filter(vpd_anomaly >= 0.082, vpd_anomaly <= 0.455),
    aes(x = vpd_anomaly, y = 0, fill = n_species)) +
    geom_tile(width = 0.04, height = 0.2) +
    scale_fill_gradient(low = "lightgrey", high = "black", name = "No. of species") +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "top", 
      legend.text = element_text(size = 10)
    )
  
  gg_tile
  
  # save tile (number of species bar)
  ggsave(
    path = here("03_results"),
    filename = paste0("Figure4_tile.png"),
    plot = gg_tile,
    width = 110,
    height = 50,
    dpi = 300,
    units = "mm"
  )
  

# 7. Create and save figure 5, aggregating by functional group ---------

## Generate smooth trends for visualization and better overlap between the trend and the confidence interval
smooth_data <- newdata_mean_fg |> 
  group_by(aridity_factor, sdev_factor, funct_group) |> 
  mutate(
    wmean.pred.g.smooth = predict(smooth.spline(vpd_anomaly, wmean.pred.g, spar = 1), x = vpd_anomaly)$y,
    upper.ci.g.smooth = predict(smooth.spline(vpd_anomaly, upper.ci.g, spar = 1), x = vpd_anomaly)$y,
    lower.ci.g.smooth = predict(smooth.spline(vpd_anomaly, lower.ci.g, spar = 1), x = vpd_anomaly)$y, 
    wmean.pred.mo.smooth = predict(smooth.spline(vpd_anomaly, wmean.pred.mo, spar = 1), x = vpd_anomaly)$y,
    upper.ci.mo.smooth = predict(smooth.spline(vpd_anomaly, upper.ci.mo, spar = 1), x = vpd_anomaly)$y,
    lower.ci.mo.smooth = predict(smooth.spline(vpd_anomaly, lower.ci.mo, spar = 1), x = vpd_anomaly)$y, 
    wmean.pred.mi.smooth = predict(smooth.spline(vpd_anomaly, wmean.pred.mi, spar = 1), x = vpd_anomaly)$y,
    upper.ci.mi.smooth = predict(smooth.spline(vpd_anomaly, upper.ci.mi, spar = 1), x = vpd_anomaly)$y,
    lower.ci.mi.smooth = predict(smooth.spline(vpd_anomaly, lower.ci.mi, spar = 1), x = vpd_anomaly)$y
  ) |> 
  ungroup() |> 
  select(aridity_factor, sdev_factor, vpd_anomaly, funct_group, 
         wmean.pred.g.smooth, upper.ci.g.smooth, lower.ci.g.smooth, 
         wmean.pred.mo.smooth, upper.ci.mo.smooth, lower.ci.mo.smooth, 
         wmean.pred.mi.smooth, upper.ci.mi.smooth, lower.ci.mi.smooth) |> 
  pivot_longer(
    cols = -c(aridity_factor, sdev_factor, vpd_anomaly, funct_group),  # todo lo demás se transforma
    names_to = c(".value", "response"),             # .value crea columnas mean.pred, upper, lower
    names_pattern = "(wmean\\.pred|upper\\.ci|lower\\.ci)\\.(.+)\\.smooth"
  ) 
  
create_figure_5_fg <- function(aridity_level) {
  
  if (aridity_level == "Arid") {
    color = "darkorange3"
  } else {
    color = "skyblue3"
  }
  
  gg <- ggplot(smooth_data |> 
                      mutate(response = factor(response, levels = c("g", "mo", "mi"), 
                                               labels = c("Growth", "Mortality occurrence", "Mortality intensity")), 
                             sdev_factor = factor(
                               sdev_factor,
                               levels = c("Early", "Late"), 
                               labels = c("Early stand development", "Late stand development"))) |> 
                      filter(aridity_factor == aridity_level), 
                    aes(x = vpd_anomaly, y = wmean.pred, group = funct_group, linetype = funct_group)) +
    
    geom_line(size = 0.7, color = color) +
    
    geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), 
                fill = color, 
                alpha = 0.3, color = NA, show.legend = FALSE) +
    
    scale_linetype_manual(
      values = c("BL" = "dashed", "NL" = "solid"),
      labels = c("BL" = "Broad-leaved", "NL" = "Needle-leaved")
    ) +
    
    facet_grid(sdev_factor ~ response) +
    
    labs(
      x = "VPD anomaly", 
      y = expression(
        paste(
          "Growth or M. intensity (",
          m^2, " ", ha^-1, " ", yr^-1, " or Mortality occurrence (prob.))"
        )
      )
    ) +
    
    scale_x_continuous(limits = c(0.082, 0.455)) +
    
    theme_minimal() +
    
    theme(
      strip.background = element_blank(),
      strip.placement = "outside",
      axis.text = element_text(color = "black", size = 8),
      strip.text = element_text(size = 10),
      axis.line = element_line(color = "lightgrey"),
      panel.background = element_rect(fill = "white", colour = "lightgrey", linetype = "solid"),
      legend.position = "bottom",
      legend.key = element_blank(),
      legend.background = element_blank(), 
      legend.title = element_blank(), 
      plot.title = element_text(hjust = 0.5, size = 12), 
      legend.text = element_text(size = 10), 
      plot.background = element_rect("white")
    ) +
    
    labs(x = "VPD anomaly") +
    
    guides(color = guide_legend(title = NULL)) +
    
    coord_cartesian(ylim = c(0, .8))
  
  gg
  
  ggsave(
    path = here("03_results"),
    filename = paste0("Figure5_trends", aridity_level, ".png"),
    plot = gg,
    width = 135,
    height = 110,
    dpi = 300,
    units = "mm"
  )
  
  return(gg)
  
}

aridity_levels <- c("Arid", "Wet")
map(aridity_levels, \(x) create_figure_5_fg(x))

  
create_figure_5_tiles <- function(functional_group) {
  
  ## TILES (number of species bars)
  gg_tile <- ggplot(
    n_species_vpd_fg |> 
      filter(funct_group == functional_group) |> 
      filter(vpd_anomaly >= 0.07, vpd_anomaly <= 0.41),
    aes(x = vpd_anomaly, y = 0, fill = n_species)) +
    geom_tile(width = 0.04, height = 0.2) +
    scale_fill_gradient(low = "lightgrey", high = "black", name = "No. of species") +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "top", 
      legend.text = element_text(size = 10)
    )
  
  # save tile (number of species bar)
  ggsave(
    path = here("03_results"),
    filename = paste0("Figure5_tile_", functional_group, ".png"),
    plot = gg_tile,
    width = 110,
    height = 50,
    dpi = 300,
    units = "mm"
  )
  
  return(gg_tile)
  
}

functional_groups <- c("BL", "NL")
map(functional_groups, \(x) create_figure_5_tiles(x))


  
