# K Fold Cross Validation Script

# Note: need to run : model_design_outsample.R and create_expectations.R script and then
# devtools::load_all()
# this will allow K fold cross validation to run

cross_validation(
  project,
  mod.name="logit_c_mod1",
  zone.dat = "ZoneID",
  groups = "Observations",
  k = 5,
  time_var = NULL,
  use.scalers = FALSE,
  scaler.func = NULL
)

###################
# --- OVERWRITE THE BROKEN WRAPPER FUNCTION ---
model_design_outsample <- function(project, 
                                   mod_name, 
                                   outsample_mod_name = NULL, 
                                   CV = FALSE, 
                                   CV_dat = NULL, 
                                   use_scalers = FALSE, 
                                   scaler_func = NULL){
  
  flag <- 0
  if(!CV){ 
    tryCatch({
      suppressWarnings(outsample_dat <- readRDS(paste0(locoutput(project), project, "filtered_outsample.rds")))
    }, error = function(e) {flag <<- 1})
    if(flag == 1) stop('A filtered out-of-sample dataset is required.')
  } else { 
    if(is.null(CV_dat)) stop('One or more cross validation datasets are empty.')
    outsample_dat <- CV_dat 
  }
  
  if (table_exists(paste0(project, "ModelInputData"), project)) {
    mdf <- model_design_list(project)
  } else {
    stop('Model input table does not exist.')
  }
  
  mdf_n <- model_names(project)
  mdf <- mdf[[which(mdf_n == mod_name)]]
  
  if(length(mdf$expectcatchmodels) > 0){
    e_list <- unlist(mdf$expectcatchmodels)
    e_settings <- expected_catch_list(project)[paste0(e_list,"_settings")]
  } else {
    e_settings <- NULL
  }
  
  alt_insample <- unserialize_table(paste0(project,"AltMatrix"), project)
  if(!is.null(alt_insample$spatname) && alt_insample$spatname == "spat"){
    alt_insample$spatname <- NULL
  }
  
  create_alternative_choice(outsample_dat, project, 
                            occasion = alt_insample$occasion, 
                            occasion_var = alt_insample$occasion_var,
                            alt_var = alt_insample$alt_var, 
                            dist.unit = alt_insample$altChoiceUnits, 
                            min.haul = 0, 
                            zoneID = alt_insample$zoneID,
                            zone.cent.name = alt_insample$zone_cent_name, 
                            fish.cent.name = alt_insample$fish_cent_name,
                            spatname = alt_insample$spatname, 
                            spatID = alt_insample$spatID, 
                            outsample = TRUE)
  
  if(!is.null(e_settings)){
    for(i in 1:length(e_list)){
      tmp_settings <- e_settings[paste0(e_list[i], "_settings")][[1]]
      create_expectations(dat = outsample_dat, project = project, name = "TEMP_NAME",
                          catch = tmp_settings$catch, price = tmp_settings$price,
                          defineGroup = tmp_settings$defineGroup, temporal = tmp_settings$temporal,
                          calc_method = tmp_settings$calc_method, outsample = TRUE)
    }
  }
  
  if(is.null(outsample_mod_name) | is_empty(outsample_mod_name)){
    outsample_mod_name <- paste0(mod_name,"_outsample")  
  }
  
  # --- THIS IS THE CRITICAL FIX ---
  make_model_design(project = project, 
                    catchID = mdf$catchID, 
                    likelihood = mdf$likelihood, 
                    initparams = mdf$initparams,
                    optimOpt = mdf$optimOpt, 
                    methodname = mdf$methodname, 
                    mod.name = outsample_mod_name, # CHANGED FROM mod_name TO mod.name
                    vars1 = mdf$vars1, 
                    vars2 = mdf$vars2, 
                    priceCol = mdf$priceCol, 
                    expectcatchmodels = mdf$expectcatchmodels,
                    startloc = mdf$startloc, 
                    polyn = mdf$polyn, 
                    crs = mdf$crs, 
                    outsample = TRUE, 
                    CV_dat = CV_dat)
}
#####################

model_design_outsample <- function(project, 
                                   mod_name, 
                                   outsample_mod_name = NULL, 
                                   CV = FALSE, 
                                   CV_dat = NULL, 
                                   use_scalers = FALSE, 
                                   scaler_func = NULL){
  
  flag <- 0
  if(!CV){ 
    tryCatch({
      suppressWarnings(outsample_dat <- readRDS(paste0(locoutput(project), project, "filtered_outsample.rds")))
    }, error = function(e) {flag <<- 1})
    if(flag == 1) stop('A filtered out-of-sample dataset is required.')
  } else { 
    if(is.null(CV_dat)) stop('One or more cross validation datasets are empty.')
    outsample_dat <- CV_dat 
  }
  
  if (table_exists(paste0(project, "ModelInputData"), project)) {
    mdf <- model_design_list(project)
  } else {
    stop('Model input table does not exist.')
  }
  
  mdf_n <- model_names(project)
  mdf <- mdf[[which(mdf_n == mod_name)]]
  
  if(length(mdf$expectcatchmodels) > 0){
    e_list <- unlist(mdf$expectcatchmodels)
    e_settings <- expected_catch_list(project)[paste0(e_list,"_settings")]
  } else {
    e_settings <- NULL
  }
  
  alt_insample <- unserialize_table(paste0(project,"AltMatrix"), project)
  if(!is.null(alt_insample$spatname) && alt_insample$spatname == "spat"){
    alt_insample$spatname <- NULL
  }
  
  create_alternative_choice(outsample_dat, project, 
                            occasion = alt_insample$occasion, 
                            occasion_var = alt_insample$occasion_var,
                            alt_var = alt_insample$alt_var, 
                            dist.unit = alt_insample$altChoiceUnits, 
                            min.haul = 0, 
                            zoneID = alt_insample$zoneID,
                            zone.cent.name = alt_insample$zone_cent_name, 
                            fish.cent.name = alt_insample$fish_cent_name,
                            spatname = alt_insample$spatname, 
                            spatID = alt_insample$spatID, 
                            outsample = TRUE)
  
  if(!is.null(e_settings)){
    for(i in 1:length(e_list)){
      tmp_settings <- e_settings[paste0(e_list[i], "_settings")][[1]]
      t_var <- ifelse(!is.null(tmp_settings$temp_var), tmp_settings$temp_var, "date_time")
      
      create_expectations(dat = outsample_dat, project = project, name = "TEMP_NAME",
                          catch = tmp_settings$catch, price = tmp_settings$price,
                          defineGroup = tmp_settings$defineGroup, 
                          temp_var = t_var, 
                          temporal = tmp_settings$temporal,
                          calc_method = tmp_settings$calc_method, outsample = TRUE)
    }
  }
  
  if(is.null(outsample_mod_name) | is_empty(outsample_mod_name)){
    outsample_mod_name <- paste0(mod_name,"_outsample")  
  }
  
  # --- FIXING THE expectcatchmodels ARGUMENT ---
  make_model_design(project = project, 
                    catchID = mdf$catchID, 
                    likelihood = mdf$likelihood, 
                    initparams = mdf$initparams,
                    optimOpt = mdf$optimOpt, 
                    methodname = mdf$methodname, 
                    mod.name = outsample_mod_name, 
                    vars1 = mdf$vars1, 
                    vars2 = mdf$vars2, 
                    priceCol = mdf$priceCol, 
                    
                    # FORCING THIS TO 'individual' TO MATCH YOUR ORIGINAL DESIGN
                    expectcatchmodels = list('individual'), 
                    
                    startloc = mdf$startloc, 
                    polyn = mdf$polyn, 
                    crs = mdf$crs, 
                    outsample = TRUE, 
                    CV_dat = CV_dat)
}

# --- A. RE-DEFINE THE HELPER WITH THE CORRECT MODEL NAME ---
prepare_rtmb_matrices <- function(project, outsample_dat, core_vessels_list, low_freq_list, master_zones) {
  
  # 1. Generate FishSET matrices for the current subset
  mdf_out <- model_design_outsample(project = project, 
                                    mod_name = "logit_c_mod1", 
                                    CV = TRUE, 
                                    CV_dat = outsample_dat)
  
  # 2. Extract Matrices and FIX alignment
  # FishSET stores the haul_id as rownames in these matrices
  dist_mat <- mdf_out$distance
  rev_mat  <- mdf_out$gridVaryingVariables$exp1
  
  # Convert to long format using the ACTUAL rownames from the matrices
  dist_long <- as.data.frame(dist_mat) %>%
    mutate(haul_id = rownames(dist_mat)) %>% # USE ROWNAMES, NOT outsample_dat$haul_id
    pivot_longer(-haul_id, names_to = "ZoneID", values_to = "dist")
  
  rev_long <- as.data.frame(rev_mat) %>%
    mutate(haul_id = rownames(rev_mat)) %>% # USE ROWNAMES
    pivot_longer(-haul_id, names_to = "ZoneID", values_to = "rev")
  
  # 3. Reconstruct the long-form data
  df_out_long <- outsample_dat %>%
    dplyr::select(haul_id, vessel_id, first_haul_dummy, not_first_dummy, new_zoneID, price_per_km, mean_hauls) %>%
    inner_join(dist_long, by = "haul_id") %>%
    inner_join(rev_long, by = c("haul_id", "ZoneID")) %>%
    mutate(
      dist_km = as.numeric(dist) * 1.60934,
      fuel = (price_per_km * dist_km),
      weighted_fuel = fuel / mean_hauls,
      profit = rev - weighted_fuel
    ) %>%
    filter(ZoneID %in% master_zones)
  
  # 4. Pivot to Wide
  model_mats <- pivot_to_wide_matrices(
    data = df_out_long,
    id_col = "haul_id",
    names_from_col = "ZoneID",
    values_to_spread = c("selected", "profit", "first_haul_dummy", "not_first_dummy")
  )
  
  Y <- model_mats$selected
  profit_mat <- model_mats$profit
  
  # Re-sync vessel mapping with the rows that survived the join
  v_mapping <- df_out_long$vessel_id[match(rownames(Y), df_out_long$haul_id)]
  
  v_covs <- list()
  for(v_id in core_vessels_list) {
    v_mask <- as.numeric(v_mapping == v_id)
    v_covs[[paste0("beta_vessel_first_", v_id)]] <- v_mask * (profit_mat * model_mats$first_haul_dummy)
    v_covs[[paste0("beta_vessel_not_first_", v_id)]] <- v_mask * (profit_mat * model_mats$not_first_dummy)
  }
  
  low_mask <- as.numeric(v_mapping %in% low_list)
  v_covs[["beta_vessel_first_LOW"]] <- low_mask * (profit_mat * model_mats$first_haul_dummy)
  v_covs[["beta_vessel_not_first_LOW"]] <- low_mask * (profit_mat * model_mats$not_first_dummy)
  
  all_covs <- c(list(beta_profit_first = profit_mat * model_mats$first_haul_dummy, 
                     beta_profit_not_first = profit_mat * model_mats$not_first_dummy), v_covs)
  
  return(list(Y = Y, covs = all_covs))
}

# --- B. RUN THE CV LOOP ---
k <- 5
set.seed(42)
main_data$fold <- sample(rep(1:k, length.out = nrow(main_data)))
cv_metrics <- list()

for (i in 1:k) {
  message(paste0("\n>>> Processing Fold ", i, " of ", k, "..."))
  train_df <- main_data[main_data$fold != i, ]
  test_df  <- main_data[main_data$fold == i, ]
  
  # Use the string "NEWP" to ensure project variable is correctly passed
  train_assets <- prepare_rtmb_matrices("NEWP", train_df, core_vessels, low_freq_vessels, unique_zones)
  
  fit_fold <- cond_logit_model(train_assets$Y, train_assets$covs, all_params)
  
  test_assets <- prepare_rtmb_matrices("NEWP", test_df, core_vessels, low_freq_vessels, unique_zones)
  
  preds <- predict_choice_probs(fit_fold, test_assets$covs)
  
  actual_zone_idx <- max.col(test_assets$Y)
  prob_of_chosen <- preds$trip_probabilities[cbind(1:nrow(test_assets$Y), actual_zone_idx)]
  
  cv_metrics[[i]] <- data.frame(fold = i, 
                                avg_prob_chosen = mean(prob_of_chosen),
                                log_loss = -mean(log(prob_of_chosen + 1e-9)), 
                                obs_count = nrow(test_df))
}

# --- C. RESULTS ---
cv_summary <- bind_rows(cv_metrics)
print(cv_summary)
