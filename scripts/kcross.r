

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