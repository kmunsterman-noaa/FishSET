# =================================================================================================
# File: rtmb_gfbt.R
# Description: This script performs a full analysis workflow for a fishery choice model using RTMB.
#              It includes data preparation, model fitting, prediction of choice probabilities,
#              simulation of a zone closure policy, and a welfare analysis of the closure.
#              Finally, it saves model outputs to a FishSET project database.
# Version: 2.0
# =================================================================================================

# LOAD LIBRARIES ----------------------------------------------------------------------------------
library(FishSET)
library(RTMB)
library(tidyverse)

# FUNCTIONS ---------------------------------------------------------------------------------------
#' Fit a Conditional Logit Choice Model using RTMB
#'
#' @description This function constructs and fits a general conditional logit model. It dynamically
#'              builds the utility function based on a list of provided covariates.
#'
#' @param response_matrix A numeric matrix (trips x zones) where each row is an observation and 
#'                        contains a 1 for the chosen alternative and 0s otherwise.
#' @param covariate_list A named list of numeric matrices. Each matrix should have the same
#'                       dimensions as `response_matrix` and represent a covariate (e.g., catch, 
#'                       distance). The names of the list elements should correspond to the 
#'                       parameter names in `start_params`.
#' @param start_params A named list of starting values for the parameters to be estimated. The 
#'                     names must begin with "beta_" and match the names in `covariate_list`.
#'
#' @return A list containing the fitted model object (`fit`), the `sdreport` object (`sdr`), and a
#'         dataframe of the formatted results (`results`) with estimates, standard errors, 
#'         and p-values.
#'        
cond_logit_model <- function(response_matrix, covariate_list, start_params, map = list()) {
  
  nll_func <- function(par) {
    # 1. Calculate utility from standard covariates
    utility_components <- mapply(function(cov, beta) cov * beta,
                                 covariate_list, par, SIMPLIFY = FALSE)
    utility <- Reduce('+', utility_components)
    
    # 2. Log-likelihood calculation
    utility_of_chosen <- rowSums(response_matrix * utility)
    # logsumexp for stability
    log_sum_exp_utility <- log(rowSums(exp(utility)))
    nll <- -sum(utility_of_chosen - log_sum_exp_utility)
    
    return(nll)
  }
  
  # Prepare data
  data_list <- c(list(Y = response_matrix), covariate_list)
  
  # Create RTMB object with map support
  obj <- RTMB::MakeADFun(func = nll_func, 
                         data = data_list, 
                         parameters = start_params,
                         map = map,
                         silent = TRUE)
  
  fit <- nlminb(start = obj$par, 
                objective = obj$fn, 
                gradient = obj$gr,
                control = list(iter.max = 1000, eval.max = 1000))
  
  sdr <- sdreport(obj)
  summary_sdr <- summary(sdr)
  
  # z-score
  z_scores <- summary_sdr[,"Estimate"]/summary_sdr[,"Std. Error"]
  
  # log-likelihood
  ll_full <- -fit$objective 
  
  # log-likelihood of null model
  null_params <- rep(0, length(start_params))
  ll_null <- -obj$fn(null_params)
  
  # calculate McFadden's Pseudo-R2
  pseudo_r2 <- 1 - (ll_full / ll_null)
  
  # Return results
  results <- data.frame(
    Estimate = summary_sdr[,"Estimate"],
    Std_error = summary_sdr[,"Std. Error"],
    z_scores = z_scores,
    p_value = 2 * pnorm(-abs(z_scores)),
    pseudo_r2 = pseudo_r2,
    log_likelihood = ll_full
  )
  
  
  return(list(fit = fit, sdr = sdr, results = results, obj = obj))
}

#' Predict Choice Probabilities
#'
#' @description Calculates the predicted probability of choosing each alternative for each 
#'              observation, based on a fitted conditional logit model.
#'
#' @param model_fit The list object returned by `cond_logit_model()`.
#' @param covariate_list The same named list of covariate matrices used to fit the model.
#'
#' @return A list containing:
#'         - `trip_probabilities`: A matrix of choice probabilities for each trip (rows) and 
#'                                 zone (columns).
#'         - `zone_probabilities`: A data frame summarizing the average probability for each zone.
predict_choice_probs <- function(model_fit, covariate_list) {
  
  # Extract estimated parameters
  betas <- model_fit$results$Estimate
  
  # Re-calculate utility using estimated parameters
  utility_components <- mapply(function(cov, beta) cov * beta,
                               covariate_list, betas, SIMPLIFY = FALSE)
  utility_hat <- Reduce('+', utility_components)
  
  # Softmax function to convert utility to probability
  softmax <- function(x) exp(x) / sum(exp(x))

  # Calculate probabilities for each trip
  trip_probs <- t(apply(utility_hat, 1, softmax))

  # Calculate average probability for each zone
  zone_probs <- colMeans(trip_probs)
  zone_probs <- data.frame(
    ZoneID = names(zone_probs),
    Probability = unname(zone_probs)
  )
  
  return(list(trip_probabilities = trip_probs,
              zone_probabilities = zone_probs))
}

#' Predict Redistributed Probabilities After a Zone Closure
#'
#' @description Calculates how choice probabilities are redistributed among remaining alternatives
#'              after one or more zones are closed.
#'
#' @param model_fit The list object returned by `cond_logit_model()`.
#' @param covariate_list The same named list of covariate matrices used to fit the model.
#' @param closed_zones A character vector containing the names of the zones to be closed. These
#'                     names must match the column names of the covariate matrices.
#'
#' @return A list containing:
#'         - `trip_probabilities_redist`: A matrix of the ### choice probabilities for each trip 
#'                                        and zone.
#'         - `zone_probabilities_redist`: A data frame of the ### average probabilities for 
#'                                        each zone.

predict_redistributed_probs <- function(model_fit, covariate_list, closed_zones) {
  betas <- model_fit$results$Estimate
  
  utility_components <- mapply(function(cov, beta) cov * beta,
                               covariate_list, betas, SIMPLIFY = FALSE)
  utility_hat <- Reduce('+', utility_components)
  
  # Set utility of closed zones to -Inf to remove them from the choice set
  if(!all(closed_zones %in% colnames(utility_hat))) {
    stop("One or more 'closed_zones' not found in the data column names.")
  }
  utility_hat[, closed_zones] <- -Inf
  
  softmax <- function(x) exp(x) / sum(exp(x))
  
  # Calculate redistributed probabilities
  trip_probs_redist <- t(apply(utility_hat, 1, softmax))
  zone_probs_redist <- colMeans(trip_probs_redist)
  zone_probs_redist <- data.frame(
    ZoneID = names(zone_probs_redist),
    Probability = unname(zone_probs_redist)
  )
  return(list(trip_probabilities_redist = trip_probs_redist,
              zone_probabilities_redist = zone_probs_redist))
}


#' Calculate Welfare Change from a Zone Closure
#'
#' @description Estimates the change in economic welfare (compensating variation) resulting from a
#'              zone closure. This function incorporates parameter uncertainty by simulating draws
#'              from the estimated beta coefficients' multivariate normal distribution using the 
#'              MASS package.
#'
#' @param model_fit The list object returned by `cond_logit_model()`.
#' @param covariate_list The same named list of covariate matrices used to fit the model.
#' @param closed_zones A character vector containing the names of the zones to be closed.
#' @param cost_variable_index An integer specifying the position of the cost-related beta
#'                            coefficient in the `model_fit$results` data frame.
#' @param beta_samples A positive integer specifying the number of simulations to run for the
#'                     welfare error estimation.
#'
#' @return A data frame containing:
#'         - `mean_welfare_loss_per_trip`: The average welfare loss per trip, averaged across 
#'                                         all simulations.
#'         - `se_welfare_loss_per_trip`: The standard error of the welfare loss per trip.
#'         - `mean_total_welfare_loss`: The total welfare loss for the sample, averaged across 
#'                                      all simulations.
#'         - `se_total_welfare_loss`: The standard error of the total welfare loss.
#'         
#'
#'         
calculate_welfare_change <- function(model_fit, 
                                               covariate_list, 
                                               closed_zones, 
                                               first_haul_idx,  
                                               other_haul_idx,  
                                               is_first_haul_vec, 
                                               beta_samples = 1000) {
  
  beta_means <- model_fit$results$Estimate
  beta_vcov <- model_fit$sdr$cov.fixed
  beta_draws <- MASS::mvrnorm(n = beta_samples, mu = beta_means, Sigma = beta_vcov)
  
  num_obs <- nrow(covariate_list[[1]])
  welfare_change_per_obs <- matrix(NA, nrow = num_obs, ncol = beta_samples)
  
  for (i in 1:beta_samples) {
    betas_drawn <- beta_draws[i, ]
    
    # 1. Calculate utilities
    utility_components <- mapply(function(cov, beta) cov * beta, 
                                 covariate_list, betas_drawn, SIMPLIFY = FALSE)
    utility_before <- Reduce('+', utility_components)
    
    utility_after <- utility_before
    utility_after[, closed_zones] <- -Inf 
    
    # 2. Assign the correct Marginal Utility (Alpha) for each row
    # If is_first_haul_vec is TRUE, use first_haul_idx; else use other_haul_idx
    alpha_vec <- ifelse(is_first_haul_vec, 
                        abs(betas_drawn[first_haul_idx]), 
                        abs(betas_drawn[other_haul_idx]))
    
    # 3. Logsum Calculation
    ls_before <- log(rowSums(exp(utility_before)))
    ls_after <- log(rowSums(exp(utility_after)))
    
    # 4. Welfare change (delta LS / alpha)
    # We use alpha_vec so each row is divided by its specific Beta
    welfare_change_per_obs[, i] <- (ls_after - ls_before) / alpha_vec
  }
  
  return(welfare_change_per_obs)
}

#' Convert Long Data Frame to Multiple Wide Matrices
#'
#' @description A helper function that takes a long-format data frame and efficiently
#'              pivots multiple specified value columns into separate wide-format matrices.
#'
#' @param data A long-format data frame.
#' @param id_col A character string specifying the name of the column that identifies
#'               unique observations (e.g., "TRIPID"). This column's values will become the rows.
#' @param names_from_col A character string specifying the name of the column whose unique
#'                       values will become the ### column names (e.g., "ZoneID").
#' @param values_to_spread A character vector of column names whose values will be spread
#'                         into the ### wide matrices. A separate matrix will be created
#'                         for each column name provided.
#'
#' @return A named list of wide-format matrices. The names of the list elements correspond
#'         to the column names provided in `values_to_spread`.
pivot_to_wide_matrices <- function(data, id_col, names_from_col, values_to_spread) {
  
  # Use lapply to iterate over each variable that needs to be pivoted
  wide_matrices <- lapply(values_to_spread, function(value_col) {
    data %>%
      # Select only the necessary columns for the pivot
      dplyr::select(dplyr::all_of(c(id_col, names_from_col, value_col))) %>%
      # Pivot the data
      tidyr::pivot_wider(
        id_cols = dplyr::all_of(id_col),
        names_from = dplyr::all_of(names_from_col),
        values_from = dplyr::all_of(value_col)
      ) %>%
      # Remove the ID column to leave just the matrix data
      dplyr::select(-dplyr::all_of(id_col)) %>%
      # Convert the resulting data frame to a matrix
      as.matrix()
  })
  
  # Name the list elements for easy access later (e.g., list$selected, list$expected_catch)
  names(wide_matrices) <- values_to_spread
  
  return(wide_matrices)
}

# DEFINE VARIABLES --------------------------------------------------------------------------------
## Note: Run all FishSET function prior to this; up until discretefish_subroutine()
## Note: This is completed in rtmb_gfbt_prep.R

# Set the FishSET project name
project <- "NEWP"
update_folderpath()

# DATA PREPARATION --------------------------------------------------------------------------------

# Load data from the FishSET project
main_data <- table_view(paste0(project, "MainDataTable"), project)

#####################################################################
# Create dummy variable
main_data <- main_data %>%
  group_by(ftid) %>%
  mutate(first_haul_dummy = ifelse(haul_counter == min(haul_counter), 1, 0)) %>%
  ungroup()

main_data$not_first_dummy <- as.numeric(!main_data$first_haul_dummy)
######################################################################

altc_data <- unserialize_table(paste0(project,"AltMatrix"), project)
mdf <- model_design_list(project)[[1]]

# Filter data based on zones present in the alternative catch matrix
unique_zones <- unique(altc_data$greaterNZ)
main_data <- main_data %>%
  filter(new_zoneID %in% unique_zones)

# Create a long format dataframe for all possible choices
# 'haul_id' in this example is the unique row/observation ID

df <- data.frame(
  zones = rep(unique_zones, length(main_data$haul_id)),
  obsID = rep(main_data$haul_id, each = length(unique_zones))
)

# Identify the chosen zone for each haul
df$zone_obs <- paste0(df$zones, df$obsID)
selected <- paste0(main_data$new_zoneID, main_data$haul_id)
df$selected <- 0
df$selected[which(df$zone_obs %in% selected)] <- 1
df$zone_obs <- NULL

# Reshape covariate data (distance, revenue, and fuel) to long format and join

# distance
distance_long <- as.data.frame(mdf$distance) %>%
  mutate(haul_id = main_data$haul_id) %>%
  pivot_longer(
    cols = -c(haul_id),
    names_to = "zones",
    # change this to distance from haul
    values_to = "distance_from_port")

# expected revenue
exp_revenue_long <- as.data.frame(mdf$gridVaryingVariables$exp1) %>%
  mutate(haul_id = main_data$haul_id) %>%
  pivot_longer(
    cols = -c(haul_id),
    names_to = "zones",
    values_to = "expected_revenue")

############################

# Join all data together
df_long <- df %>%
  mutate(zones = as.character(zones)) %>%
  left_join(distance_long, by = c("zones" = "zones", "obsID" = "haul_id")) %>%
  left_join(exp_revenue_long, by = c("zones" = "zones", "obsID" = "haul_id")) %>%
  
  # Join vessel_id and dummies from main_data (haul-specific)
  left_join(main_data %>% dplyr::select(haul_id, vessel_id, first_haul_dummy, not_first_dummy),
            by = c("obsID" = "haul_id")) %>%
  rename(ZoneID = zones, haul_id = obsID) %>%
  dplyr::select(ZoneID, haul_id, vessel_id, selected, distance_from_port, 
                expected_revenue, first_haul_dummy, not_first_dummy) %>%
  mutate(distance_from_port = as.numeric(distance_from_port) * 1.60934)

# Final data cleaning: remove zones with NA distance values
zones_to_remove <- unique(df_long[which(is.na(df_long$distance_from_port)),]$ZoneID)
df_long <- df_long %>% 
  filter(!(ZoneID %in% zones_to_remove))

# Convert data from long to wide matrix format for the model
model_matrices <- pivot_to_wide_matrices(
  data = df_long,
  id_col = "haul_id",
  names_from_col = "ZoneID",
  values_to_spread = c("selected", "expected_revenue", "distance_from_port",
                       "first_haul_dummy", "not_first_dummy")
)

Y <- model_matrices$selected
revenue <- model_matrices$expected_revenue
distance <- model_matrices$distance_from_port

# Fuel Covariate

# fuel cost
price_per_km <- main_data$price_per_km
mean_hauls <- main_data$mean_hauls

fuel <- (price_per_km * distance)
weighted_fuel <- (price_per_km * distance)/mean_hauls

# profit
profit <- revenue - weighted_fuel

# set up dummy variables
profit_first <- profit * model_matrices$first_haul_dummy
profit_not_first <- profit * model_matrices$not_first_dummy

# --- Create Vessel Fixed Effects ---

# 1. Get the unique haul IDs in the exact order they appear in the wide matrix
# Usually, this is the unique 'id_col' from your df_long
ordered_haul_ids <- unique(df_long$haul_id)

# 2. Map vessel_id to the rows of Y
# We look up the vessel for each haul_id in our ordered list
vessel_mapping <- main_data$vessel_id[match(ordered_haul_ids, main_data$haul_id)]

unique_vessels <- unique(vessel_mapping)

# 3. Create the vessel-specific profit interactions
vessel_covs <- list()
vessel_starts <- list()

# Use 2:length to set a reference vessel (baseline)
for(i in 2:length(unique_vessels)) {
  v_id <- unique_vessels[i]
  
  # Binary mask: 1 if this row belongs to vessel i
  vessel_mask <- as.numeric(vessel_mapping == v_id)
  
  # Multiply profit matrix by the mask (broadcasting the vector to match the matrix)
  # This makes the matrix 'profit' for vessel i and '0' for all others
  v_profit_matrix <- vessel_mask * profit 
  
  v_name <- paste0("beta_vessel_", v_id)
  vessel_covs[[v_name]] <- v_profit_matrix
  vessel_starts[[v_name]] <- 0
}

# MODEL FITTING -----------------------------------------------------------------------------------

# Profit #

# Define inputs for conditional logit model
covariates_base <- list(profit_first = profit_first,
                        profit_not_first = profit_not_first)

# Combine lists
all_covariates <- c(covariates_base, vessel_covs)
all_params <- c(list(beta_profit_first = 0, beta_profit_not_first = 0), 
                vessel_starts)

# Fit the conditional logit model
results_profit_vessel <- cond_logit_model(response_matrix = Y,
                                   covariate_list = all_covariates,
                                   start_params = all_params)

print(results_profit_vessel$results)

###########################

# --- Save model results ---

master_results <- results_profit_vessel$results %>%
  # 1. Convert row names into an actual column called 'Variable'
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Variable") %>%
  
  mutate(
    # 2. Clean up the variable names
    variable = case_when(
      Variable == "beta_profit_first"     ~ "profit_first",
      Variable == "beta_profit_not_first" ~ "profit_not_first",
      grepl("beta_vessel_", Variable)     ~ gsub("beta_vessel_", "", Variable),
      TRUE                                ~ Variable
    ),
    
    # 2. Add your specific metadata
    IOPAC_PORT_GROUP = "NEWP",
    min_haul         = 1,
    obs              = nrow(main_data),
    unique_zones     = length(unique_zones),
    
    # 3. Calculate AIC (using the NLL/objective from the fit)
    AIC = (2 * length(results_profit_vessel$fit$par)) + (2 * results_profit_vessel$fit$objective)
  )

# Save the combined file
saveRDS(master_results, 
        file = here::here("data", "confidential", "FishSETfolder", "NEWP", "output", "complete_model_outputs_NEWP.rds"))

## check

# Count hauls per vessel in Newport
vessel_haul_counts <- main_data %>%
  group_by(vessel_id) %>%
  tally(name = "num_hauls") %>%
  arrange(num_hauls)

print(vessel_haul_counts)

# Summary of the distribution
summary(vessel_haul_counts$num_hauls)
nrow(vessel_haul_counts[vessel_haul_counts$num_hauls < 20, ])

# POLICY SIMULATION AND WELFARE ANALYSIS ---------------------------------------

# ----------------------------- PROFIT MODEL -----------------------------------

# --- Predict baseline probabilities ---
predicted_probabilities_profit <- predict_choice_probs(results_profit_vessel, all_covariates)

# Save data
predicted_probs_zones_profit <- predicted_probabilities_profit[[2]]
saveRDS(predicted_probs_zones_profit, file=here::here("data", "confidential", "FishSETfolder", "NEWP", "output", "predicted_probs_zones_profit_NEWP.rds"))

# --- Load and process the zone closure scenario ---

## SCENARIO 1

scenario_name <- "closure_1"

filename <- paste0(locoutput(project), scenario_name, "_closures.yaml")

scenario_1 <- readRDS("~/Documents/GitHub/FishSET/data/non-confidential/other/scen_1.rds")

yaml::write_yaml(scenario_1, filename)

closures <- yaml::read_yaml(filename)

closed_zones <- gsub("Zone_", "", closures$GRID5KM_ID[closures$scenario == scenario_name])

unique_zones <- unique(main_data$new_zoneID)
closed_zones_scen1 <- intersect(closed_zones, unique_zones)

# --- Predict redistributed probabilities under the closure ---

# First item in list is redistributed probabilities for each trip
# Second item is redistributed probabilities for each zone

redistributed_probabilities_scen1_profit <- predict_redistributed_probs(results_profit_vessel, all_covariates, closed_zones_scen1)

# Save data
redist_probs_zones_scen1_profit <- redistributed_probabilities_scen1_profit[[2]]
saveRDS(redist_probs_zones_scen1_profit, file=here::here("data", "confidential", "FishSETfolder", "NEWP", "output", "redist_probs_zones_scen1_profit_NEWP.rds"))

# --- Run welfare analysis ---
# The losses are reported as positive values in the output (thus, negative 
# values would indicate gains)

# ------ Marginal Utility of Income: Profit ------
# Create the logical vector: TRUE if it's the first haul, FALSE otherwise
main_data$is_first <- main_data$haul_counter == 1

if (length(closed_zones_scen1) > 0) {

## profit
welfare_output_scen1_profit <- calculate_welfare_change(model_fit = results_profit_vessel,
                                                        covariate_list = all_covariates,
                                                        closed_zones = closed_zones_scen1,
                                                        first_haul_idx = 1,           # Position of Profit (First Haul) Beta
                                                        other_haul_idx = 2,           # Position of Profit (Subsequent) Beta
                                                        is_first_haul_vec = main_data$is_first, 
                                                        beta_samples = 1000)

# Calculate mean welfare loss per TRIP

# 1. Sum all hauls (first + all others) for each of the 1,000 simulations per trip
trip_simulations_scen1_profit <- as.data.frame(welfare_output_scen1_profit) %>%
  mutate(trip_id = main_data$trip_id) %>%
  group_by(trip_id) %>%
  summarise(across(starts_with("V"), sum))

# 2. Calculate the mean and SD across the entire datasets
welfare_per_trip_scen1_profit <- trip_simulations_scen1_profit %>%
  tidyr::pivot_longer(cols = -trip_id, names_to = "sim", values_to = "welfare") %>%
  mutate(welfare_loss = -welfare) %>% 
  summarise(
    mean_welfare_loss = mean(welfare_loss, na.rm=TRUE),
    sd_welfare_loss = sd(welfare_loss, na.rm = TRUE)
  ) %>%
  mutate(mui = "profit")

# View results
head(welfare_per_trip_scen1_profit)

} else {
  
  # SKIP ANALYSIS: Create dummy row with 0s
  message("Scenario 1: No overlap found for this port. Filling with zeros.")
  
  welfare_per_trip_scen1_profit <- data.frame(
    mean_welfare_loss = 0,
    sd_welfare_loss = 0,
    mui = "profit"
  )
}

## SCENARIO 2

scenario_name <- "closure_2"

filename <- paste0(locoutput(project), scenario_name, "_closures.yaml")

scenario_2 <- readRDS("~/Documents/GitHub/FishSET/data/non-confidential/other/scen_2.rds")

yaml::write_yaml(scenario_2, filename)

closures <- yaml::read_yaml(filename)

closed_zones <- gsub("Zone_", "", closures$GRID5KM_ID[closures$scenario == scenario_name])

unique_zones <- unique(main_data$new_zoneID)
closed_zones_scen2 <- intersect(closed_zones, unique_zones)

# --- Predict redistributed probabilities under the closure ---

# First item in list is redistributed probabilities for each trip
# Second item is redistributed probabilities for each zone

redistributed_probabilities_scen2_profit <- predict_redistributed_probs(results_profit_vessel, all_covariates, closed_zones_scen2)

# Save data
redist_probs_zones_scen2_profit <- redistributed_probabilities_scen2_profit[[2]]
saveRDS(redist_probs_zones_scen2_profit, file=here::here("data", "confidential", "FishSETfolder", "NEWP", "output", "redist_probs_zones_scen2_profit_NEWP.rds"))

# --- Run welfare analysis ---
# The losses are reported as positive values in the output (thus, negative 
# values would indicate gains)

## profit

if (length(closed_zones_scen1) > 0) {
  
welfare_output_scen2_profit <- calculate_welfare_change(model_fit = results_profit_vessel,
                                                        covariate_list = all_covariates,
                                                        closed_zones = closed_zones_scen2,
                                                        first_haul_idx = 1,           # Position of profit (First Haul) Beta
                                                        other_haul_idx = 2,           # Position of profit (Subsequent) Beta
                                                        is_first_haul_vec = main_data$is_first, 
                                                        beta_samples = 1000)

# Calculate mean welfare loss per TRIP

# 1. Sum all hauls (first + all others) for each of the 1,000 simulations per trip
trip_simulations_scen2_profit <- as.data.frame(welfare_output_scen2_profit) %>%
  mutate(trip_id = main_data$trip_id) %>%
  group_by(trip_id) %>%
  summarise(across(starts_with("V"), sum))

# 2. Calculate the mean and SD across the entire datasets
welfare_per_trip_scen2_profit <- trip_simulations_scen2_profit %>%
  tidyr::pivot_longer(cols = -trip_id, names_to = "sim", values_to = "welfare") %>%
  mutate(welfare_loss = -welfare) %>% 
  summarise(
    mean_welfare_loss = mean(welfare_loss, na.rm=TRUE),
    sd_welfare_loss = sd(welfare_loss, na.rm = TRUE)
  ) %>%
  mutate(mui = "profit")

# View results
head(welfare_per_trip_scen2_profit)

} else {
  
  # SKIP ANALYSIS: Create dummy row with 0s
  message("Scenario 2: No overlap found for this port. Filling with zeros.")
  
  welfare_per_trip_scen2_profit <- data.frame(
    mean_welfare_loss = 0,
    sd_welfare_loss = 0,
    mui = "profit"
  )
}

# --- Save welfare results ---

welfare_results_profit_vessel <- bind_rows(
  mutate(welfare_per_trip_scen1_profit, scenario = 1, closed_val = length(closed_zones_scen1)),
  mutate(welfare_per_trip_scen2_profit, scenario = 2, closed_val = length(closed_zones_scen2))
) %>%
  mutate(
    IOPAC_PORT_GROUP = "NEWP",
    beta_samples     = 1000,
    unique_zones     = length(zOut$table$n),
    closed_zones     = closed_val,
    prop_closed      = closed_zones / unique_zones,
    model            = "profit"
  ) %>%
  dplyr::select(-closed_val)

# Save the combined file
saveRDS(welfare_results_profit_vessel, 
        file = here::here("data", "confidential", "FishSETfolder", "NEWP", "output", "complete_welfare_outputs_profit_NEWP.rds"))
