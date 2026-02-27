##############################################################################
##
## Script name: Fig 3 species histograms
## Purpose: Create Histograms figure (Figure 3 in the ms)
##
## Input data: species_list.rds, models.rds, newdatabase_sp_list.rds 
## Output data: none 
## Output figures: Histograms figures (Figure 3 in the ms)
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

# 2. Create species-specific predictions for species low and high VPD anomaly values ----------

newdatabse_reduced_list <- map(newdatabase_sp_list, \(x) x |> 
                     # filter the minimum and maximum value of the species
                      filter(vpd_anomaly == min(vpd_anomaly) | vpd_anomaly == max(vpd_anomaly)) |>
                     # create a new factor variable names vpd factor
                      mutate(vpd_factor = case_when(vpd_anomaly == min(vpd_anomaly) ~ "Low", T ~ "High")) |> 
  group_by(aridity_factor, sdev_factor, vpd_factor) |> 
  # keep the fixed values of the variables to predict 
  summarise(ba_sp = first(ba_sp), 
            density = first(density),
            country = first(country),
            harvest_reg = first(harvest_reg), 
            species = first(species), 
            funct_group = first(funct_group), 
            cmi = first(cmi), 
            sdevelopment = first(sdevelopment), 
            vpd_anomaly = first(vpd_anomaly)) |> 
  ungroup())

names(newdatabse_reduced_list) <- names(species_list)

# Function to create predictions 
create_sp_predictions <- function(sp_name) {
  
  data <- newdatabse_reduced_list[[sp_name]]
  
  mg <- models[[sp_name]][["growth"]]
  mo <- models[[sp_name]][["mort_occ"]]
  mi <- models[[sp_name]][["mort_int"]]
  
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

sp_names <- names(species_list)
histograms_pred_list <- map(sp_names, \(x) create_sp_predictions(x))
names(histograms_pred_list) <- sp_names 
glimpse(histograms_pred_list[[1]])

# 3. Calculate significance of the pattern ------------------------------------
## Does the CI interval of the difference in the prediction between high 
## VPD anomaly and low VPD anomaly includes cero? 

sp_table_hist_low_vpd <- list_rbind(histograms_pred_list) |> 
  filter(vpd_factor == "Low") |>
  select(!c(vpd_factor, vpd_anomaly)) |>
  rename(pred_g_low = pred_g, 
         pred_g_se_low = se_g, 
         pred_mo_low = pred_mo, 
         pred_mo_se_low = se_mo, 
         pred_mi_low = pred_mi, 
         pred_mi_se_low = se_mi) 

sp_table_hist_high_vpd <- list_rbind(histograms_pred_list) |> 
  filter(vpd_factor == "High") |>
  select(!c(vpd_factor, vpd_anomaly)) |>
  rename(pred_g_high = pred_g, 
         pred_g_se_high = se_g, 
         pred_mo_high = pred_mo, 
         pred_mo_se_high = se_mo, 
         pred_mi_high = pred_mi, 
         pred_mi_se_high = se_mi) 

sp_table_all <- full_join(sp_table_hist_low_vpd, sp_table_hist_high_vpd) 

sp_table_signif_95_ci <- sp_table_all |>
  mutate(
    ## 1. differences
    g_diff =
      pred_g_high - pred_g_low,
    mo_diff =
      pred_mo_high - pred_mo_low,
    mi_diff =
      pred_mi_high - pred_mi_low,
    g_pattern = case_when(g_diff > 0 ~ "Increase", T ~ "Decline"),
    mo_pattern = case_when(mo_diff > 0 ~ "Increase", T ~ "Decline"),
    mi_pattern = case_when(mi_diff > 0 ~ "Increase", T ~ "Decline"),
    
    
    ## 2. SE of differences (Gelman et al pag 61 (75 del pdf):https://books.google.es/books?hl=es&lr=&id=fILoDwAAQBAJ&oi=fnd&pg=PR11&dq=Regression+and+Other+Stories&ots=Wy-SryxqT9&sig=aAR5HNQCJ6IEJg0FBIByWFjljDk#v=onepage&q=Regression%20and%20Other%20Stories&f=false)
    g_diff_se = sqrt(
      pred_g_se_high^2 +
        pred_g_se_low^2
    ),
    mo_diff_se = sqrt(
      pred_mo_se_high^2 +
        pred_mo_se_low^2
    ),
    mi_diff_se = sqrt(
      pred_mi_se_high^2 +
        pred_mi_se_low^2
    ),
    
    ## 3. 95% CI for differences
    g_diff_ci_low =
      g_diff - 1.96 * g_diff_se,
    g_diff_ci_high =
      g_diff + 1.96 * g_diff_se,

    mo_diff_ci_low =
      mo_diff - 1.96 * mo_diff_se,
    mo_diff_ci_high =
      mo_diff + 1.96 * mo_diff_se,

    mi_diff_ci_low =
      mi_diff - 1.96 * mi_diff_se,
    mi_diff_ci_high =
      mi_diff + 1.96 * mi_diff_se,

    ## 4. significance
    g_sig =
      as.integer(!(g_diff_ci_low <= 0 & g_diff_ci_high >= 0)),
    
    mo_sig =
      as.integer(!(mo_diff_ci_low <= 0 & mo_diff_ci_high >= 0)),
    
    mi_sig =
      as.integer(!(mi_diff_ci_low <= 0 & mi_diff_ci_high >= 0))
    
  ) |> 
  select(c(aridity_factor, sdev_factor, species, funct_group, g_diff,
           mo_diff, mi_diff, g_pattern, mo_pattern, mi_pattern,
           g_sig, mo_sig, mi_sig)) |>
  mutate(sdev_factor = case_when(sdev_factor == "Early" ~ "Early stand development", T ~ "Late stand development")) |>
  pivot_longer(
    cols = matches("^(g|mo|mi)_(diff|pattern|sig)$"),
    names_to = c("variable", ".value"),
    names_pattern = "(g|mo|mi)_(.*)"
  ) 

# 4. Get significance numbers (numbers of the figures) --------------------------------
significant_species_95 <- sp_table_signif_95_ci |> filter(sig ==1)

significant_species_95 |> group_by(aridity, stdev, variable, pattern) |>
  summarise(n = n()) |> print(n = 36)

# Per type
significant_species_95 |> group_by(type, aridity, stdev, variable, pattern) |>
  summarise(n = n()) |> print(n = 30)

# 5. Create and save histogram figures ------------------------------------------------
create_figure_3_histograms <- function(aridity_level) {
  
  sp_table_signif <- sp_table_signif_95_ci |> 
    mutate(
      sig_level = case_when(sig == 1 ~ "sig",
                            TRUE       ~ "ns"
      )
    )

  
  if (aridity_level == "Arid") {
    color = "darkorange3"
    df <- sp_table_signif |> filter(aridity_factor == "Arid")
    title = "a)"

  } else {
    if (aridity_level == "Mild") {
      color = "darkgreen"
      df <- sp_table_signif |> filter(aridity_factor == "Mild")
      title = "b)"
      
    } else {
      color = "skyblue3"
      df <- sp_table_signif |> filter(aridity_factor == "Wet")
      title = "b)"
    }
  }
  
  create_figure <- function(var) {
    
    figure <- ggplot(
      
      df |> mutate(variable = factor(variable, levels = c("g", "mo", "mi"), 
                                     labels = c("Growth", 
                                                "Mortality occurrence", 
                                                "Mortality intensity"))) |> 
        filter(variable == var),
      
      aes(x = diff, fill = factor(sig_level))) +
      
      geom_histogram(color = "black", alpha = 0.7) +
      
      facet_grid(sdev_factor ~ variable) +  
      
      scale_fill_manual(values = c("ns" = "grey","sig" = color),
                        labels =  c("ns" = "Non-significant", "sig" = "significant 95%"),
                        name = "Statistical significance") +
      
      labs(
        x = "Sensitivity to VPD anomaly",
        y = "No. species",
        title = title
      ) +
      
      geom_vline(xintercept = 0, linetype = "dashed") +
      
      theme_minimal() +
      
      theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 10),
        strip.text.y.right = element_text(angle = -90),
        axis.text = element_text(color = "black", size = 10),
        axis.ticks = element_blank(),
        axis.line = element_line(color = "lightgrey"),
        panel.background = element_rect(fill = "white", colour = "lightgrey", linetype = "solid"),
        legend.position = "none",
        legend.key = element_blank(),
        legend.background = element_blank(), 
        plot.title = element_text(hjust = -0.1, vjust = .8, size = 8), 
        plot.background = element_rect(color = "white")
      ) +
      
      xlim(-2.5, 2.5) +
      
      ylim(0, 13)
    
    return(figure)
    
  }
  
  var = list("Growth", "Mortality occurrence", "Mortality intensity")
  
  figure <- map(var, \(x) create_figure(x))
  
  names(figure) <- c("g", "mo", "mi")
  
  g <- figure$g + theme(
    axis.title.x = element_blank(),
    strip.text.y.right = element_blank(),
  )
  
  mo <- figure$mo + theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), 
    strip.text.y.right = element_blank(),
    plot.title = element_blank()
  )
  
  mi <- figure$mi + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y.left = element_blank(),
    strip.background = element_blank(),
    plot.title = element_blank()
  )
  
  all <- (g + mo + mi)
  
  ggsave(
    path = here("03_results"),
    filename = paste0("Figure3_", aridity_level, ".png"),
    plot = all,
    width = 200,
    height = 120,
    dpi = 300,
    units = "mm"
  )
  
  return(all)
  
}

aridity_levels <- unique(sp_table_signif_95_ci$aridity_factor)
map(aridity_levels, \(x) create_figure_3_histograms(x))




