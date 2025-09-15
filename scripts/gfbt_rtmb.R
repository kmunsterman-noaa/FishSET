
## Note: Run all FishSET function prior to this; up until discretefish_subroutine()

# DATA PREP ---------------------------------------------------------------------------------------
# Set project # This should already be set
project <- "astoria"

# Load data
main_data <- table_view(paste0(project, "MainDataTable"), project)
altc_data <- unserialize_table(paste0(project,"AltMatrix"), project)

# Save unique zones based on alternative catch with specified min haul threshold
unique_zones <- unique(altc_data$greaterNZ)

# Filter main data for 
main_data <- main_data %>%
  filter(ZoneID %in% unique_zones)

# Create a long format dataframe for RTMB
# 'haul_id' in this example is the unique row/observation ID
df <- data.frame(zones = rep(unique_zones, length(main_data$haul_id)),
                 obsID = rep(main_data$haul_id, each = length(unique_zones)))

head(df)

# Combine zone and obsID to link to observations in the main_data
df$zone_obs <- paste0(df$zones, df$obsID)

selected <- paste0(main_data$ZoneID, main_data$haul_id)

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

# Add haul ID to distance to pivot long
distance$haul_id <- main_data$haul_id
distance_long <- pivot_longer(
  data = distance,
  cols = -c(haul_id),
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
exp_catch$haul_id <- main_data$haul_id

exp_catch_long <- pivot_longer(
  data = exp_catch,
  cols = -c(haul_id),
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
names(df) <- c("ZoneID", "haul_id", "selected", "distance_from_port", "expected_catch")

# saveRDS(df, "astoria_long.rds")

df$distance_from_port <- as.numeric(df$distance_from_port)

# Final data check
# Check to see if any zones had all NAs for distance, so filter zones out of data
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
  dplyr::select(haul_id, ZoneID, selected) %>%
  pivot_wider(id_cols = haul_id, names_from = ZoneID, values_from = selected)
Y <- as.matrix(Y[,-1])

catch <- df %>%
  dplyr::select(haul_id, ZoneID, expected_catch) %>%
  pivot_wider(id_cols = haul_id, names_from = ZoneID, values_from = expected_catch)
catch <- as.matrix(catch[,-1])

distance <- df %>%
  dplyr::select(haul_id, ZoneID, distance_from_port) %>%
  pivot_wider(id_cols = haul_id, names_from = ZoneID, values_from = distance_from_port)
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

tic()
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

toc()

print(results)
