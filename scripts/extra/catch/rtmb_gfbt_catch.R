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
                                     beta_samples = 20) {
  
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("The 'MASS' package is required for this function. Please install it.", call. = FALSE)
  }
  
  # Extract means (point estimates) and the variance-covariance matrix for the betas
  beta_means <- model_fit$results$Estimate
  beta_vcov <- model_fit$sdr$cov.fixed
  
  # Generate random draws of the beta parameters using the MASS package
  beta_draws <- MASS::mvrnorm(n = beta_samples, mu = beta_means, Sigma = beta_vcov)
  
  # Store the total welfare change from each simulation
  total_welfare_changes <- numeric(beta_samples)
  # Store all per-trip welfare changes for overall mean and se
  all_trip_welfare_changes <- vector("list", beta_samples)
  # Store trip loss
  welfare_change_per_trip <- matrix(NA, nrow=dim(covariate_list[[1]])[1], ncol=beta_samples)
  
  # Loop through each simulated set of betas
  for (i in 1:beta_samples) {
    # Extract the i-th row, which is one full set of drawn betas
    betas_drawn <- beta_draws[i, ]
    
    # Calculate marginal utility of income (alpha) from the drawn cost coefficient
    alpha_drawn <- -betas_drawn[cost_variable_index]
    
    # Calculate utility before and after closure using the drawn betas
    utility_components <- mapply(function(cov, beta) cov * beta,
                                 covariate_list, betas_drawn, SIMPLIFY = FALSE)
    utility_before <- Reduce('+', utility_components)
    
    utility_after <- utility_before
    utility_after[, closed_zones] <- -Inf # Because exp(-Inf) = 0
    
    # Calculate welfare using the log-sum formula
    welfare_before <- (1 / alpha_drawn) * log(rowSums(exp(utility_before)))
    welfare_after <- (1 / alpha_drawn) * log(rowSums(exp(utility_after)))
    
    # Calculate and store the change for this simulation
    welfare_change <- welfare_after - welfare_before
    welfare_change_per_trip[,i] <- welfare_change
    all_trip_welfare_changes[[i]] <- welfare_change
    total_welfare_changes[i] <- sum(welfare_change)
  }
  
  # Combine all per-trip changes into a single vector to calculate mean and se
  all_trip_changes_vec <- unlist(all_trip_welfare_changes)
  
  # Summarize the results, reporting losses as positive values
  results_df <- list(
    summary = data.frame(
      Metric = c("Per Trip", "Total Sample"),
      Mean_Welfare_Loss = c(-mean(all_trip_changes_vec), -mean(total_welfare_changes)),
      Standard_Error = c(sd(all_trip_changes_vec), sd(total_welfare_changes))
    ),
    
    welfare_change_per_trip
  )
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
project <- "eureka"

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
exp_catch_long <- as.data.frame(mdf$gridVaryingVariables$exp1) %>%
  mutate(haul_id = main_data$haul_id) %>%
  pivot_longer(
    cols = -c(haul_id),
    names_to = "zones",
    values_to = "expected_catch")

# Join all data together
df_long <- df %>%
  mutate(zones = as.character(zones)) %>%
  left_join(distance_long, by = c("zones" = "zones", "obsID" = "haul_id")) %>%
  left_join(exp_catch_long, by = c("zones" = "zones", "obsID" = "haul_id")) %>%
  rename(ZoneID = zones, haul_id = obsID) %>%
  dplyr::select(ZoneID, haul_id, selected, distance_from_haul, expected_catch) %>%
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
  values_to_spread = c("selected", "expected_catch", "distance_from_haul")
)

Y <- model_matrices$selected
catch <- model_matrices$expected_catch
distance <- model_matrices$distance_from_haul

# Add in additional covariates

# fuel cost
price_per_km <- main_data$price_per_km
fuel <- price_per_km * distance

# fishery fidelity x catch/revenue interaction

# center fishery fidelity
ff <- main_data %>%
  dplyr::select(vessel_id, fishery_fidelity) %>%
  group_by(vessel_id) %>%
  unique()

ff_mean <- mean(ff$fishery_fidelity)
fishery_fidelity <- (main_data$fishery_fidelity)-ff_mean

ff_catch <- fishery_fidelity*catch

# vessel length x catch/revenue interaction

# center vessel length
vl <- main_data %>%
  dplyr::select(vessel_id, vessel_length) %>%
  group_by(vessel_id) %>%
  unique()

length_mean <- mean(vl$vessel_length)
vessel_length <- (main_data$vessel_length)-length_mean

length_catch <- vessel_length*catch

# vessel length x fuel cost interaction

# center vessel length
vl <- main_data %>%
  dplyr::select(vessel_id, vessel_length) %>%
  group_by(vessel_id) %>%
  unique()

length_mean <- mean(vl$vessel_length)
vessel_length <- (main_data$vessel_length)-length_mean

length_cost <- vessel_length*fuel

# MODEL FITTING ----------------------------------------------------------------

# Define inputs for conditional logit model
covariates <- list(catch = catch, fuel = fuel, ff_catch = ff_catch, length_catch = length_catch, length_cost = length_cost) #
starting_params <- list(beta_catch = 0, beta_fuel = 0, beta_ff_catch = 0, beta_length_catch = 0, beta_length_cost = 0) #

# Fit the conditional logit model
results <- cond_logit_model(response_matrix = Y,
                            covariate_list = covariates,
                            start_params = starting_params)

print(results)

# SAVE MODEL OUTPUTS -----------------------------------------------------------

# --- Model fit table ---
LL <- results$fit$objective
k <- length(results$fit$par)
n <- dim(Y)[1]
# Null model to calculate PseudoR2
n_choices <- dim(Y)[2]  # Number of choices (zones)
LL_null <- -n * log(n_choices)

AIC <- round((-2 * LL) + (2 * k), 3)
AICc <- round(AIC + (((2 * k) * (k + 1)) / (n - k - 1)), 3)
BIC <- round((-2 * LL) + (k * log(n)), 3)
PseudoR2 <- round(1 - (-LL / LL_null), 3)

# Save dataframe
mod.out <- data.frame(matrix(NA, nrow = 4, ncol = 1))
mod.out[, 1] <- c(AIC, AICc, BIC, PseudoR2)
rownames(mod.out) <- c("AIC", "AICc", "BIC", "PseudoR2")
colnames(mod.out) <- mdf$mod.name

# --- Model output table ---
ModelOut <- list(
  name = mdf$mod.name,
  errorExplain = NULL, # Just filling in NULL value here when running RTMB  
  OutLogit = data.frame(
    estimate = results$results$Estimate,
    std_error = results$results$Std_error,
    t_value = results$results$Estimate/results$results$Std_error),
  optoutput = list(
    counts = results$fit$evaluations,
    convergence = results$fit$convergence,
    optim_message = NULL), # Just fill in NULL value here when running RTMB
  seoutmat2 = results$results$Std_error,
  MCM = list(
    AIC = AIC,
    AICc = AICc,
    BIC = BIC,
    PseudoR2 = PseudoR2),
  H1 = unname(results$sdr$cov.fixed),
  choice.table = data.frame(
    choice = colnames(Y)[max.col(Y)]),
  params = unname(results$fit$par),
  modTime = 10 # Just filling in random time here - this isn't used for anything  
)

# Save ModelOut table to outputs folder

MCM <- ModelOut[["MCM"]]
write.csv(MCM, 
          file = here::here("data", "confidential", "FishSETfolder", "eureka", "output", "MCM_mod4.csv"),
          row.names = FALSE)

params <- results[["results"]]
write.csv(params, 
          file=here::here("data", "confidential", "FishSETfolder", "eureka", "output", "params_mod4.csv"))


# POLICY SIMULATION AND WELFARE ANALYSIS ---------------------------------------

# --- Predict baseline probabilities ---
# First item in list is predicted probabilities for each trip
# Second item in list is predicted probabilities for each zone
predicted_probabilities <- predict_choice_probs(results, covariates)

# --- Load and process the zone closure scenario ---

scenario_name <- "closure_2"

filename <- paste0(locoutput(project), scenario_name, "_closures.yaml")

scenario_2 <- readRDS("~/Documents/GitHub/FishSET/data/non-confidential/other/scen_2.rds")

yaml::write_yaml(scenario_2, filename)

closures <- yaml::read_yaml(filename)

closed_zones <- gsub("Zone_", "", closures$GRID5KM_ID[closures$scenario == scenario_name])

unique_zones <- unique(main_data$new_zoneID)
closed_zones <- intersect(closed_zones, unique_zones)

# --- Predict redistributed probabilities under the closure ---
# First item in list is redistributed probabilities for each trip
# Second item is redistributed probabilities for each zone
redistributed_probabilities <- predict_redistributed_probs(results, covariates, closed_zones)

# --- Run welfare analysis ---
# The losses are reported as positive values in the output (thus, negative 
# values would indicate gains)
welfare_output <- calculate_welfare_change(results, 
                                           covariates, 
                                           closed_zones, 
                                           cost_variable_index = 1,
                                           beta_samples = 20)

# SAVE RESULTS -----------------------------------------------------------------

redist_probs <- redistributed_probabilities[["trip_probabilities_redist"]]
saveRDS(redist_probs, file=here::here("data", "confidential", "FishSETfolder", "eureka", "output", "redist_probs.rds"))

welfare_per_haul <- welfare_output[[2]]
saveRDS(welfare_per_haul, file=here::here("data", "confidential", "FishSETfolder", "eureka", "output", "welfare_per_haul_scen2.rds"))

# SCENARIO 1

scenario_name <- "closure_1"

filename <- paste0(locoutput(project), scenario_name, "_closures.yaml")

scenario_1 <- readRDS("~/Documents/GitHub/FishSET/data/non-confidential/other/scen_1.rds")

yaml::write_yaml(scenario_1, filename)

closures <- yaml::read_yaml(filename)

closed_zones <- gsub("Zone_", "", closures$GRID5KM_ID[closures$scenario == scenario_name])

unique_zones <- unique(main_data$new_zoneID)
closed_zones <- intersect(closed_zones, unique_zones)

# --- Predict redistributed probabilities under the closure ---
# First item in list is redistributed probabilities for each trip
# Second item is redistributed probabilities for each zone
redistributed_probabilities <- predict_redistributed_probs(results, covariates, closed_zones)

# --- Run welfare analysis ---
# The losses are reported as positive values in the output (thus, negative 
# values would indicate gains)
welfare_output <- calculate_welfare_change(results, 
                                           covariates, 
                                           closed_zones, 
                                           cost_variable_index = 1,
                                           beta_samples = 20)

# SAVE RESULTS -----------------------------------------------------------------

redist_probs <- redistributed_probabilities[["trip_probabilities_redist"]]
saveRDS(redist_probs, file=here::here("data", "confidential", "FishSETfolder", "eureka", "output", "redist_probs.rds"))

welfare_per_haul <- welfare_output[[2]]
saveRDS(welfare_per_haul, file=here::here("data", "confidential", "FishSETfolder", "eureka", "output", "welfare_per_haul_scen1.rds"))

# TESTING WELFARE VALUES -------------------------------------------------------

zones <- zOut$table

scenario_3 <- zones %>%
  filter(n > 30) %>% #30 (38), 8 (same number of zones as closure), 1 (153 closed zones)
  rename(GRID5KM_ID = new_zoneID) %>%
  dplyr::select(-n) %>%
  mutate(scenario = "closure_3")
 
scenario_name <- "closure_3"

filename <- paste0(locoutput(project), scenario_name, "_closures.yaml")

yaml::write_yaml(scenario_3, filename)

closures <- yaml::read_yaml(filename)

closed_zones <- gsub("Zone_", "", closures$GRID5KM_ID[closures$scenario == scenario_name])

unique_zones <- unique(main_data$new_zoneID)
closed_zones <- intersect(closed_zones, unique_zones)

# --- Predict redistributed probabilities under the closure ---
# First item in list is redistributed probabilities for each trip
# Second item is redistributed probabilities for each zone
redistributed_probabilities <- predict_redistributed_probs(results, covariates, closed_zones)

# --- Run welfare analysis ---
# The losses are reported as positive values in the output (thus, negative 
# values would indicate gains)
welfare_output <- calculate_welfare_change(results, 
                                           covariates, 
                                           closed_zones, 
                                           cost_variable_index = 2,
                                           beta_samples = 20)

#####

mean_rev <- main_data %>%
  summarise(mean_rev_per_haul = mean(tow_r), sd_rev = sd(tow_r))

mean_fuel <- mean(fuel)
sd_fuel <- sd(fuel)

