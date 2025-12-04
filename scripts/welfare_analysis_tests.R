
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
