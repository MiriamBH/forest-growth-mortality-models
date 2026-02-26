##############################################################################
##
## Script name: cross_validation
## Purpose: Cross validation for every model
##
## Input data: species_list.rds
## Output data: CV_growth.rds, CV_mort_occ.rds, CV_mort_int.rds
##              CV_growth.word, CV_mort_occ.word, CV_mort_int.word
##
## Date: 2026-02-12
##
##############################################################################

# libraries ####
library(tidyverse)
library(here)
library(mgcv)
library(modelr)
library(yardstick)
library(caret)
library(ROCR)
library(parallel) 
library(future) 
library(furrr) 

# 1. Read data --------------------------------------------------------------
species_list <- here("01_data", "species_list.rds") |> read_rds()

# 2. Functions to select sp and repes of the cross validation ---------------
select_sp <- function(sp_name){
  
  ## Select species 
  sp <- species_list[[sp_name]]
  
  ## Select repes
  n_repes <- 500
  cv <- modelr::crossv_mc(sp, test = .2, n = n_repes)

  crossval <- function(x){ 
    
    # Select train (FIT) and test (PREDICT) 
    train <- as_tibble(cv$train[[x]])
    test <- as_tibble(cv$test[[x]])
    
    # Modelos
    if (length(unique(train$country)) < 2) {
      
      mg <- gam(growth_annual ~
                  log(ba_sp) + log(density) + harvest_reg +
                  te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
                family = Gamma(link = "log"),
                data = train,
                method = "REML")
      
      mo <- gam(mort_occ_annual ~ 
                  log(ba_sp) + log(density) + harvest_reg +
                  te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
                family = "binomial",
                data = train,
                method = "REML")
      
      mi <- gam(mort_int_annual ~
                  log(ba_sp) + log(density) + harvest_reg +
                  te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
                family = Gamma(link = "log"),
                data = train |> filter(mort_int_annual > 0), 
                method = "REML")
      
    } else {
      
      mg <- gam(growth_annual ~
                  log(ba_sp) + log(density) + harvest_reg + country +
                  te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
                family = Gamma(link = "log"),
                data = train,
                method = "REML")
      
      mo <- gam(mort_occ_annual ~ 
                  log(ba_sp) + log(density) + harvest_reg + country +
                  te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
                family = "binomial",
                data = train,
                method = "REML")
      
      mi <- gam(mort_int_annual ~
                  log(ba_sp) + log(density) + harvest_reg + country +
                  te(cmi, vpd_anomaly, sdevelopment, bs = "tp"),
                family = Gamma(link = "log"),
                data = train |> filter(mort_int_annual > 0),
                method = "REML")
      
    }
    
    # Growth ####
    pred_g <- tibble(
      obs  = test$growth_annual,
      pred = predict(mg, 
                     newdata = test, 
                     type = "response") |> as.numeric()
    )
    
    summary(pred_g)
    
    results_growth <- tibble(
      sp = sp_name, 
      mean_growth = mean(pred_g$obs),
      median_growth = median(pred_g$obs),
      range_growth = quantile(pred_g$obs, probs = .75) - quantile(pred_g$obs, probs = .25),
      # Higher, better:
      r2 = yardstick::rsq(pred_g, truth = obs, estimate = pred)$.estimate,
      # Lower, better:
      rmse = yardstick::rmse(pred_g, truth = obs, estimate = pred)$.estimate,
      nrmse = ((rmse - median_growth) / range_growth),
      # Lower, better:
      mae = yardstick::mae(pred_g, truth = obs, estimate = pred)$.estimate,
      nmae = ((mae - median_growth) / range_growth),
      # To see the consequences of the country filter
      test_row = nrow(test), 
      train_row = nrow(train)
      
    )
    
    
    # Mortality intensity ####
    train_mi <- train |> filter(mort_int_annual > 0)
    test_mi <-  test |> filter(mort_int_annual > 0)
    
    pred_mi <- tibble(
      obs  = test_mi$mort_int_annual,
      pred = predict(mi, 
                     newdata = test_mi, 
                     type = "response") |> as.numeric()
    )
    
    summary(pred_mi)
    
    results_mort_int <- tibble(
      sp = sp_name, 
      mean_mort_int = mean(pred_mi$obs),
      median_mort_int = median(pred_mi$obs),
      range_mort_int = quantile(pred_mi$obs, probs = .75) - quantile(pred_mi$obs, probs = .25),
      # Higher, better:
      r2 = yardstick::rsq(pred_mi, truth = obs, estimate = pred)$.estimate,
      # Lower, better:
      rmse = yardstick::rmse(pred_mi, truth = obs, estimate = pred)$.estimate,
      nrmse = ((rmse - median_mort_int) / range_mort_int),
      # Lower, better:
      mae = yardstick::mae(pred_mi, truth = obs, estimate = pred)$.estimate,
      nmae = ((mae - median_mort_int) / range_mort_int),
      # To see the consequences of the country filter
      test_row = nrow(test_mi), 
      train_row = nrow(train_mi)
      
    )
    
    
    # Mortality occurrence (binomial) ####
    ## Predictions with test 
    predictions <- predict(
      mo, type = "response", newdata = test) |> 
      as.numeric()
    
    procr <- prediction(predictions, test |> select(mort_occ_annual))
    
    ## Calculate sens y spec #
    temporal <- performance(procr,measure = "sens","spec")
    diferencia <- abs(temporal@x.values[[1]]- temporal@y.values[[1]])
    
    se <- temporal@y.values[[1]][which.min(diferencia)]
    sp <- temporal@x.values[[1]][which.min(diferencia)]
    SESPmil <- (se+sp)/2
    
    ## Youden_index for dynamic threshold #
    perf <- performance(procr, "tpr", "fpr")
    sens <- perf@y.values[[1]]   # vector de sensitivities
    fpr  <- perf@x.values[[1]]   # vector de false positive rate
    spec <- 1 - fpr              # vector de specificities
    
    ### Corresponding thresholds 
    thresh <- slot(performance(procr, "sens", "spec"), "alpha.values")[[1]]
    
    ### Get rid of infine values (ROCR includes at beginning/end)
    valid <- is.finite(thresh)
    sens <- sens[valid]
    spec <- spec[valid]
    thresh <- thresh[valid]
    
    ### Youden index
    youden <- sens + spec - 1
    best_idx <- which.max(youden)
    best_threshold <- thresh[best_idx]
    
    best_threshold
    
    ## Perfomrance with auc #
    auc <- performance(procr, measure = "auc")
    aucmil <- auc@y.values[[1]]
    
    # create a confusion matrix with thershold value = 0.25
    results_acc <- tibble(
      obs = test$mort_occ_annual,
      preds = case_when(
        predictions >= best_threshold ~ 1, 
        T ~ 0)) |> 
      mutate(across(everything(), as.factor))
    
    conf <- caret::confusionMatrix(
      data = results_acc$obs,
      reference = results_acc$preds,
      positive = "0")
    
    ## Results
    results_mort_occ <- tibble(
      sp = sp_name, 
      auc = aucmil, 
      SESP = SESPmil,
      Sensiv = conf[["byClass"]][["Sensitivity"]], 
      Specif = conf[["byClass"]][["Specificity"]],
      test_row = nrow(test), 
      train_row = nrow(train)
    )
    
    # Return
    sp_repes <- list(
      results_growth, results_mort_occ, results_mort_int)
    names(sp_repes) <- c("growth", "mort_occ", "mort_int")
    
    return(sp_repes)
    
  } # Run by rep
  
  # Repes 3 models ####
  safe_crossval <- purrr::possibly(crossval, otherwise = NULL)
  repes <- purrr::map(seq_len(n_repes), \(x) safe_crossval(x))
  
  
  # Unir repeticiones por tipo de modelo
  growth_table <- purrr::map_dfr(repes, \(x) x |> pluck("growth"))
  mort_occ_table <- purrr::map_dfr(repes, \(x) x |> pluck("mort_occ"))
  mort_int_table <- purrr::map_dfr(repes, \(x) x |> pluck("mort_int"))
  
  # Return
  all_models_bysp <- list(
    growth = growth_table,
    mort_occ  = mort_occ_table,
    mort_int  = mort_int_table
  )
  
  return(all_models_bysp)
  
} 

## 2.1. Run by species --------------------------------------------------------------
future::plan(multisession, workers = 2)
allrepes_cv <- furrr::future_map(
  names(species_list), \(x) select_sp(x), 
  .options = furrr::furrr_options(seed = TRUE)
); names(allrepes_cv) <- names(species_list)

write_rds(allrepes_cv, here("01_data", "crossval_allrepes.rds"))


## 2.2. Summary by species ----------------------------------------------------------
# allrepes_cv <- here("01_data", "crossval_allrepes.rds") |> read_rds()
models_name <- c("growth", "mort_occ", "mort_int")

all_models_bysp <- purrr::map(
  models_name, function(model) {
    
    allrepes_cv |> purrr::map(
      \(x) x[[model]]) |> list_rbind()
    
  }
); names(all_models_bysp) <- models_name

final_results <- all_models_bysp |> 
  purrr::map(
    \(x) x |> group_by(sp) |> 
      summarise(across(
        !c(train_row, test_row), 
        list(mean = mean, sd = sd), 
        na.rm = TRUE))
  )

# 3. Save final tables ----------------------------------------------------------
purrr::map2(
  final_results, names(final_results), 
  \(x, y){ 
    x |> write_rds(paste0("01_data/CV_", y, ".rds"))
    x |> sjPlot::tab_df(
      file = paste0("03_results/CV_", y, ".word"), 
      digits = 4)
  }
)


