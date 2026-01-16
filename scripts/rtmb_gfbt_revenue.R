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
cond_logit_model <- function(response_matrix, covariate_list, start_params) {
  # Define the negative log-likelihood function dynamically
  nll_func <- function(par) {
    # Calculate utility by multiplying each covariate matrix by its corresponding beta
    # and summing the results.
    utility_components <- mapply(function(cov, beta) cov * beta,
                                 covariate_list, par, SIMPLIFY = FALSE)
    utility <- Reduce('+', utility_components)
    
    # Calculate the log-likelihood
    utility_of_chosen <- rowSums(response_matrix * utility)
    log_sum_exp_utility <- log(rowSums(exp(utility)))
    nll <- -sum(utility_of_chosen - log_sum_exp_utility)
    
    return(nll)
  }
  
  # Prepare data and parameter lists for RTMB
  data_list <- c(list(Y = response_matrix), covariate_list)
  
  # Create and fit the RTMB object
  obj <- RTMB::MakeADFun(func = nll_func, 
                         data = data_list, 
                         parameters = start_params, 
                         silent = TRUE)
  
  fit <- nlminb(start = obj$par, 
                objective = obj$fn, 
                gradient = obj$gr)
  
  sdr <- sdreport(obj)
  
  # Format and return results
  summary_sdr <- summary(sdr)
  z_scores <- summary_sdr[,"Estimate"]/summary_sdr[,"Std. Error"]
  results <- data.frame(
    Estimate = summary_sdr[,"Estimate"],
    Std_error = summary_sdr[,"Std. Error"],
    z_scores = z_scores,
    p_value = 2 * pnorm(-abs(z_scores))
  )
  
  return(list(fit = fit, sdr = sdr, results = results))
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
calculate_welfare_change <- function(model_fit, 
                                     covariate_list, 
                                     closed_zones, 
                                     cost_variable_index,
                                     is_cost_variable = FALSE,  # <--- NEW ARGUMENT HERE
                                     beta_samples = 20) {
  
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("The 'MASS' package is required for this function. Please install it.", call. = FALSE)
  }
  
  beta_means <- model_fit$results$Estimate
  beta_vcov <- model_fit$sdr$cov.fixed
  beta_draws <- MASS::mvrnorm(n = beta_samples, mu = beta_means, Sigma = beta_vcov)
  
  total_welfare_changes <- numeric(beta_samples)
  all_trip_welfare_changes <- vector("list", beta_samples)
  welfare_change_per_trip <- matrix(NA, nrow=dim(covariate_list[[1]])[1], ncol=beta_samples)
  
  for (i in 1:beta_samples) {
    betas_drawn <- beta_draws[i, ]
    
    if(is_cost_variable) {
      # CASE: Fuel Cost (Negative Coefficient)
      alpha_drawn <- -betas_drawn[cost_variable_index]
      
    } else {
      # CASE: Revenue (Positive Coefficient)
      alpha_drawn <- betas_drawn[cost_variable_index]
    }
    
    # Calculate utility
    utility_components <- mapply(function(cov, beta) cov * beta, 
                                 covariate_list, betas_drawn, SIMPLIFY = FALSE)
    utility_before <- Reduce('+', utility_components)
    
    utility_after <- utility_before
    utility_after[, closed_zones] <- -Inf 
    
    # Calculate welfare
    welfare_before <- (1 / alpha_drawn) * log(rowSums(exp(utility_before)))
    welfare_after <- (1 / alpha_drawn) * log(rowSums(exp(utility_after)))
    
    welfare_change <- welfare_after - welfare_before
    welfare_change_per_trip[,i] <- welfare_change
    all_trip_welfare_changes[[i]] <- welfare_change
    total_welfare_changes[i] <- sum(welfare_change)
  }
  
  all_trip_changes_vec <- unlist(all_trip_welfare_changes)
  
  results_df <- list(
    summary = data.frame(
      Metric = c("Per Trip", "Total Sample"),
      Mean_Welfare_Loss = c(-mean(all_trip_changes_vec), -mean(total_welfare_changes)),
      Standard_Error = c(sd(all_trip_changes_vec), sd(total_welfare_changes))
    ),
    welfare_change_per_trip
  )
  return(results_df)
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
project <- "scwa"
update_folderpath()

# DATA PREPARATION --------------------------------------------------------------------------------

# Load data from the FishSET project
main_data <- table_view(paste0(project, "MainDataTable"), project)
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

# Reshape covariate data (distance, catch/revenue, and fuel) to long format and join

# distance
distance_long <- as.data.frame(mdf$distance) %>%
  mutate(haul_id = main_data$haul_id) %>%
  pivot_longer(
    cols = -c(haul_id),
    names_to = "zones",
    values_to = "distance_from_haul")

# catch/expected revenue
exp_revenue_long <- as.data.frame(mdf$gridVaryingVariables$exp1) %>%
  mutate(haul_id = main_data$haul_id) %>%
  pivot_longer(
    cols = -c(haul_id),
    names_to = "zones",
    values_to = "expected_revenue")

# Join all data together
df_long <- df %>%
  mutate(zones = as.character(zones)) %>%
  left_join(distance_long, by = c("zones" = "zones", "obsID" = "haul_id")) %>%
  left_join(exp_revenue_long, by = c("zones" = "zones", "obsID" = "haul_id")) %>%
  rename(ZoneID = zones, haul_id = obsID) %>%
  dplyr::select(ZoneID, haul_id, selected, distance_from_haul, expected_revenue) %>%
  mutate(distance_from_haul = (as.numeric(distance_from_haul))*1.60934)

# Final data cleaning: remove zones with NA distance values
zones_to_remove <- unique(df_long[which(is.na(df_long$distance_from_haul)),]$ZoneID)

df_long <- df_long %>% 
  filter(!(ZoneID %in% zones_to_remove))

# Convert data from long to wide matrix format for the model
model_matrices <- pivot_to_wide_matrices(
  data = df_long,
  id_col = "haul_id",
  names_from_col = "ZoneID",
  values_to_spread = c("selected", "expected_revenue", "distance_from_haul")
)

Y <- model_matrices$selected
revenue <- model_matrices$expected_revenue
distance <- model_matrices$distance_from_haul

# Fuel Covariate

# fuel cost
price_per_km <- main_data$price_per_km
fuel <- price_per_km * distance

# profit
profit <- revenue - fuel

# MODEL FITTING ----------------------------------------------------------------

# Revenue and Cost

# Define inputs for conditional logit model
covariates <- list(revenue = revenue, fuel = fuel)
starting_params <- list(beta_revenue = 0, beta_fuel = 0)

# Fit the conditional logit model
results <- cond_logit_model(response_matrix = Y,
                            covariate_list = covariates,
                            start_params = starting_params)

print(results)

# Profit

# Define inputs for conditional logit model
covariates <- list(profit = profit)
starting_params <- list(beta_profit = 0)

# Fit the conditional logit model
results_profit <- cond_logit_model(response_matrix = Y,
                            covariate_list = covariates,
                            start_params = starting_params)

print(results_profit)

# POLICY SIMULATION AND WELFARE ANALYSIS ---------------------------------------

# --- Predict baseline probabilities ---
# First item in list is predicted probabilities for each trip
# Second item in list is predicted probabilities for each zone
predicted_probabilities <- predict_choice_probs(results_profit, covariates)

# Save data
predicted_probs_zones <- predicted_probabilities[[2]]
saveRDS(predicted_probs_zones, file=here::here("data", "confidential", "FishSETfolder", "eureka", "output", "predicted_probs_zones.rds"))

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
redistributed_probabilities <- predict_redistributed_probs(results_profit, covariates, closed_zones_scen1)

# Save data
redist_probs_zones <- redistributed_probabilities[[2]]
saveRDS(redist_probs_zones, file=here::here("data", "confidential", "FishSETfolder", "eureka", "output", "redist_probs_zones_scen1.rds"))

# --- Run welfare analysis ---
# The losses are reported as positive values in the output (thus, negative 
# values would indicate gains)

## Revenue
welfare_output_scen1 <- calculate_welfare_change(results_profit, 
                                                     covariates, 
                                                     closed_zones_scen1, 
                                                     cost_variable_index = 1,
                                                     is_cost_variable = FALSE,
                                                     beta_samples = 20)

# Calculate mean welfare loss per TRIP
welfare_per_haul <- welfare_output_scen1[[2]]
mean_welfare_per_haul <- rowMeans(welfare_per_haul, na.rm = TRUE)

welfare_per_trip_scen1 <- main_data %>%
  dplyr::select(trip_id, haul_id) %>%
  mutate(mean_welfare_per_haul = mean_welfare_per_haul) %>%
  group_by(trip_id) %>%
  summarise(welfare_per_trip=sum(mean_welfare_per_haul)) %>%
  mutate(welfare_loss_per_trip = -1 * welfare_per_trip) %>%
  summarise(mean_welfare_per_trip = mean(welfare_loss_per_trip), sd = sd(welfare_loss_per_trip))

# positive values indicates LOSS
welfare_per_trip_scen1

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
redistributed_probabilities <- predict_redistributed_probs(results_profit, covariates, closed_zones_scen2)

# Save data
redist_probs_zones <- redistributed_probabilities[[2]]
saveRDS(redist_probs_zones, file=here::here("data", "confidential", "FishSETfolder", "eureka", "output", "redist_probs_zones_scen2.rds"))

# --- Run welfare analysis ---
# The losses are reported as positive values in the output (thus, negative 
# values would indicate gains)

## Revenue
welfare_output_scen2 <- calculate_welfare_change(results_profit, 
                                                     covariates, 
                                                     closed_zones_scen2, 
                                                     cost_variable_index = 1,
                                                     is_cost_variable = FALSE,
                                                     beta_samples = 20)

# Calculate mean welfare loss per TRIP
welfare_per_haul <- welfare_output_scen2[[2]]
mean_welfare_per_haul <- rowMeans(welfare_per_haul, na.rm = TRUE)

welfare_per_trip_scen2 <- main_data %>%
  dplyr::select(trip_id, haul_id) %>%
  mutate(mean_welfare_per_haul = mean_welfare_per_haul) %>%
  group_by(trip_id) %>%
  summarise(welfare_per_trip=sum(mean_welfare_per_haul)) %>%
  mutate(welfare_loss_per_trip = -1 * welfare_per_trip) %>%
  summarise(mean_welfare_per_trip = mean(welfare_loss_per_trip), sd = sd(welfare_loss_per_trip))

# positive values indicates LOSS
welfare_per_trip_scen2 