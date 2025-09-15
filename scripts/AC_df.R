
# DATA PREP ---------------------------------------------------------------------------------------
# Set project
project <- "scallop_rtmb"

# Load data
main_data <- table_view(paste0(project, "MainDataTable"), project)
altc_data <- unserialize_table(paste0(project,"AltMatrix"), project)

# Save unique zones based on alternative catch with specified min haul threshold
unique_zones <- unique(altc_data$greaterNZ)

# Filter main data for 
main_data <- main_data %>%
  filter(ZoneID %in% unique_zones)

# Create a long format dataframe for RTMB
# TRIPID in this example is the unique row/observation ID
df <- data.frame(zones = rep(unique_zones, length(main_data$TRIPID)),
                 obsID = rep(main_data$TRIPID, each = length(unique_zones)))

head(df)

# Combine zone and obsID to link to observations in the main_data
df$zone_obs <- paste0(df$zones, df$obsID)

selected <- paste0(main_data$ZoneID, main_data$TRIPID)

# Initialize selected zone for each observation
df$selected <- 0

# Update with the select zone for each observation
df$selected[which(df$zone_obs %in% selected)] <- 1

# double check that this is correct
sum(df$selected) == nrow(main_data)

# Remove the combo variable
df$zone_obs <- NULL

# Now get the distance matrix
mdf <- model_design_list(project)
mdf <- mdf[[1]] # Just grab the first model design
distance <- as.data.frame(mdf$distance)

# Add trip ID to distance to pivot long
distance$TRIPID <- main_data$TRIPID
distance_long <- pivot_longer(
  data = distance,
  cols = -c(TRIPID),
  names_to = "zones",
  values_to = "distance_from_port"
)

# Prep for left join
names(distance_long) <- c("obsID", "zones", "distance_from_port")
distance_long$zones <- as.numeric(distance_long$zones)

df <- left_join(
  df,
  distance_long,
  by = c("zones", "obsID")
)

head(df)

# Get expected catch
exp_catch <- as.data.frame(mdf$gridVaryingVariables$exp1)
exp_catch$TRIPID <- main_data$TRIPID

exp_catch_long <- pivot_longer(
  data = exp_catch,
  cols = -c(TRIPID),
  names_to = "zones",
  values_to = "expected_catch"
)

# Prep for left join
names(exp_catch_long) <- c("obsID", "zones", "expected_catch")
exp_catch_long$zones <- as.numeric(exp_catch_long$zones)

df <- left_join(
  df,
  exp_catch_long,
  by = c("zones", "obsID") 
)

head(df)
names(df) <- c("ZoneID", "TRIPID", "selected", "distance_from_port", "expected_catch")

# saveRDS(df, "scallop_long.rds")

df$distance_from_port <- as.numeric(df$distance_from_port)

# Final data check
# Two zones had all NAs for distance, so filter zones out of data
zones_to_remove <- unique(df[which(is.na(df$distance_from_port)),]$ZoneID)
df <- df %>% 
  filter(!(ZoneID %in% zones_to_remove))


# RTMB --------------------------------------------------------------------------------------------
library(RTMB)

# Specify negative log-likelihood function (conditional logit model)
nll_func <- function(par) {
  # 1. Get parameters
  beta_catch <- par$beta_catch
  beta_distance <- par$beta_distance
  
  # 2. Calculate systematic utility for all choices
  # These are our data matrices defined in the global environment
  utility <- beta_catch * catch + beta_distance * distance
  
  # 3. Calculate the log-likelihood
  # This calculates the utility of ONLY the chosen alternative for each fisher
  utility_of_chosen <- rowSums(Y * utility)
  
  # This is the log-sum-exp term, which is the denominator in the logit formula
  log_sum_exp_utility <- log(rowSums(exp(utility)))
  
  # The total negative log-likelihood is the sum over all fishers
  nll <- -sum(utility_of_chosen - log_sum_exp_utility)
  
  return(nll)
}

# Reformat data
head(df)
Y <- df %>%
  select(TRIPID, ZoneID, selected) %>%
  pivot_wider(id_cols = TRIPID, names_from = ZoneID, values_from = selected)
Y <- as.matrix(Y[,-1])

catch <- df %>%
  select(TRIPID, ZoneID, expected_catch) %>%
  pivot_wider(id_cols = TRIPID, names_from = ZoneID, values_from = expected_catch)
catch <- as.matrix(catch[,-1])

distance <- df %>%
  select(TRIPID, ZoneID, distance_from_port) %>%
  pivot_wider(id_cols = TRIPID, names_from = ZoneID, values_from = distance_from_port)
distance <- as.matrix(distance[,-1])

# List of data to be used by the model function
data_list <- list(
  Y = Y,
  catch = catch,
  distance = distance
)

# List of parameters with their starting values for the optimization
params <- list(
  beta_catch = 0,    # Start at 0
  beta_distance = 0  # Start at 0
)

# Compile TMB object
# Translates model to C++ representation, automatically differentiates likelihood function
obj <- RTMB::MakeADFun(
  func = nll_func, 
  data = data_list, 
  parameters = params
)

# Optimize
fit <- nlminb(
  start = obj$par,
  objective = obj$fn,
  gradient = obj$gr
)

# Get coefficients and std error
sdr <- sdreport(obj)
summary <- summary(sdr)

# Calculate p-values
estimates <- summary[,"Estimate"]
std_errors <- summary[,"Std. Error"]

# Calculate z-scores
z_scores <- estimates/std_errors

# Calculate two-tailed p-values
p_values <- 2 * pnorm(-abs(z_scores))

# Final table
results <- data.frame(
  Estimate = estimates,
  Std_error = std_errors,
  z_value = z_scores,
  p_value = p_values
)

print(results)


# ### TEST ABOVE USING RTMB ####
# head(df)
# tmp_Y <- df %>%
#   select(TRIPID, ZoneID, selected) %>%
#   pivot_wider(id_cols = TRIPID, names_from = ZoneID, values_from = selected)
# tmp_Y <- as.matrix(tmp_Y[,-1])
# 
# tmp_catch <- df %>%
#   select(TRIPID, ZoneID, expected_catch) %>%
#   pivot_wider(id_cols = TRIPID, names_from = ZoneID, values_from = expected_catch)
# tmp_catch <- as.matrix(tmp_catch[,-1])
# 
# tmp_distance <- df %>%
#   select(TRIPID, ZoneID, distance_from_port) %>%
#   pivot_wider(id_cols = TRIPID, names_from = ZoneID, values_from = distance_from_port)
# tmp_distance <- as.matrix(tmp_distance[,-1])
# 
# nll_func <- function(tmp_par) {
#   # 1. Get parameters
#   tmp_beta_catch <- tmp_par$tmp_beta_catch
#   tmp_beta_distance <- tmp_par$tmp_beta_distance
# 
#   # 2. Calculate systematic utility for all choices
#   # These are our data matrices defined in the global environment
#   tmp_utility <- tmp_beta_catch * tmp_catch + tmp_beta_distance * tmp_distance
# 
#   # 3. Calculate the log-likelihood
#   # 'Y' is our 0/1 response matrix
#   # This calculates the utility of ONLY the chosen alternative for each fisher
#   tmp_utility_of_chosen <- rowSums(tmp_Y * tmp_utility)
# 
#   # This is the log-sum-exp term, which is the denominator in the logit formula
#   tmp_log_sum_exp_utility <- log(rowSums(exp(tmp_utility)))
# 
#   # The total negative log-likelihood is the sum over all fishers
#   tmp_nll <- -sum(tmp_utility_of_chosen - tmp_log_sum_exp_utility)
# 
#   return(tmp_nll)
# }
# 
# # List of data to be used by the model function
# tmp_data_list <- list(
#   tmp_Y = tmp_Y,
#   tmp_catch = tmp_catch,
#   tmp_distance = tmp_distance
# )
# 
# # List of parameters with their starting values for the optimization
# tmp_params <- list(
#   tmp_beta_catch = 0,    # Start at 0
#   tmp_beta_distance = 0  # Start at 0
# )
# 
# # Create the RTMB object
# # This compiles the R function into C++ and links the data/parameters
# obj <- RTMB::MakeADFun(
#   func = tmp_nll_func,
#   data = tmp_data_list,
#   parameters = tmp_params
# )
# 
# # Fit the model using an optimization routine (like nlminb)
# fit <- nlminb(
#   start = obj$par,
#   objective = obj$fn,
#   gradient = obj$gr
# )
# 
# sdr <- sdreport(obj)
# 
# # View the summary of the results
# summary(sdr)


# cpp_code <- '
# // [[Rcpp::depends(RTMB)]]
# #include <rtmb.hpp>
# 
# template <class Type>
# Type conditional_logit_rcpp(Rcpp::List CORE) {
#   using namespace rtmb;
#   
#   vector<Type> beta = Rcpp::as<vector<Type> >(CORE["beta"]);
#   matrix<Type> X = Rcpp::as<matrix<Type> >(CORE["X"]);
#   vector<int> selected_int = Rcpp::as<vector<int> >(CORE["selected"]);
#   vector<int> trip_id = Rcpp::as<vector<int> >(CORE["trip_id"]);
#   vector<Type> selected = selected_int.template cast<Type>();
#   
#   vector<Type> V = X * beta;
#   vector<Type> log_sum_exp_V = log(GROUPED_SUM(exp(V), trip_id));
#   Type selected_V_sum = sum(V * selected);
#   Type nll = sum(log_sum_exp_V) - selected_V_sum;
# 
#   return nll;
# }
# 
# RTMB_FUNC(conditional_logit_rcpp)
# '
# writeLines(cpp_code, "conditional_logit.cpp")
# sourceCpp("conditional_logit.cpp")

# 
# 
# 
# 
# #### RTMB EXAMPLE ####
# library(RTMB)
# set.seed(123)
# 
# n_fishers <- 500
# n_locations <- 5
# 
# # True parameter values we want to recover
# true_beta_catch <- 1.5   # Fishers prefer high catch rates
# true_beta_distance <- -0.8 # Fishers dislike long distances
# 
# # Create location-specific attributes (these don't change per fisher)
# # Matrix where rows are fishers, columns are locations
# catch <- matrix(rlnorm(n_fishers * n_locations, meanlog = 2, sdlog = 0.5),
#                 nrow = n_fishers, ncol = n_locations, byrow = TRUE)
# distance <- matrix(runif(n_fishers * n_locations, min = 5, max = 100),
#                    nrow = n_fishers, ncol = n_locations, byrow = TRUE)
# 
# utility <- true_beta_catch * catch + true_beta_distance * distance
# 
# # Convert utility to choice probabilities using the logit formula
# # P_ij = exp(U_ij) / sum(exp(U_ik)) for k=1 to 5
# probs <- t(apply(utility, 1, function(u) exp(u) / sum(exp(u))))
# 
# # Simulate the actual choice for each fisher based on the probabilities
# # This gives us a vector of chosen location indices (1, 2, 3, 4, or 5)
# choices <- apply(probs, 1, function(p) sample(1:n_locations, size = 1, prob = p))
# 
# # Convert the choice vector to a 0/1 response matrix (dummy variable format)
# # This is the required format for the model's likelihood function
# Y <- matrix(0, nrow = n_fishers, ncol = n_locations)
# Y[cbind(1:n_fishers, choices)] <- 1
# 
# head(catch[,1:5])
# 
# nll_func <- function(par) {
#   # 1. Get parameters
#   beta_catch <- par$beta_catch
#   beta_distance <- par$beta_distance
#   
#   # 2. Calculate systematic utility for all choices
#   # These are our data matrices defined in the global environment
#   utility <- beta_catch * catch + beta_distance * distance
#   
#   # 3. Calculate the log-likelihood
#   # This calculates the utility of ONLY the chosen alternative for each fisher
#   utility_of_chosen <- rowSums(Y * utility)
#   
#   # This is the log-sum-exp term, which is the denominator in the logit formula
#   log_sum_exp_utility <- log(rowSums(exp(utility)))
#   
#   # The total negative log-likelihood is the sum over all fishers
#   nll <- -sum(utility_of_chosen - log_sum_exp_utility)
#   
#   return(nll)
# }
# 
# # List of data to be used by the model function
# data_list <- list(
#   Y = Y,
#   catch = catch,
#   distance = distance
# )
# 
# # List of parameters with their starting values for the optimization
# params <- list(
#   beta_catch = 0,    # Start at 0
#   beta_distance = 0  # Start at 0
# )
# 
# # Create the RTMB object
# # This compiles the R function into C++ and links the data/parameters
# obj <- RTMB::MakeADFun(
#   func = nll_func, 
#   data = data_list, 
#   parameters = params
# )
# 
# # Fit the model using an optimization routine (like nlminb)
# fit <- nlminb(
#   start = obj$par,
#   objective = obj$fn,
#   gradient = obj$gr
# )
# 
# sdr <- sdreport(obj)
# 
# 
# # View the summary of the results
# summary(sdr)
# 
# 

