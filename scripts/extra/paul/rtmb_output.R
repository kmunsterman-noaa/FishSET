# This is the script to add to rtmb_gfbt.R to save model outputs


# SAVE RESULTS ------------------------------------------------------------------------------------

# --- Save outputs to FishSET project ---
# Create connection to db
fishset_db <- DBI::dbConnect(RSQLite::SQLite(), locdatabase(project = project))


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

# Table name
mft_tab_nm <- paste0(project, "ModelFit") 

# Check if the table exists, if it does we want to combine results
if (table_exists(mft_tab_nm, project)) {
  mft <- model_fit(project, CV = FALSE)
  
  # If model is in the table, replace old value
  if (names(mod.out) %in% names(mft)){
    # Get the name of the column to replace
    col_to_replace <- names(mod.out)
    # Remove the old column
    mft <- mft[, !names(mft) %in% col_to_replace]
  }
  
  DBI::dbWriteTable(fishset_db, mft_tab_nm, cbind(mft, mod.out), overwrite = TRUE)
  
} else {
  # If not, create the table and save model fit
  DBI::dbWriteTable(fishset_db, mft_tab_nm, mod.out)
}

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

mot_tab_nm <- paste0(project, "ModelOut")  
mot_exists <- table_exists(mot_tab_nm, project)

if (mot_exists) {
  mot <- model_out_view(project, CV = FALSE)
  table_remove(mot_tab_nm, project)
  mot_n <- vapply(mot, function(x) x$name, character(1)) 
  if (mdf$mod.name %in% mot_n) {
    #replace existing model
    mot[[which(mot_n %in% mdf$mod.name)]] <- ModelOut
    
  } else {
    # add new model
    mot[[length(mot) + 1]] <- ModelOut
  }
  
  mot_to_save <- mot
  
} else {
  # new mot
  mot_to_save <- list(ModelOut)
}

# Save ModelOut table to db
raw_sql <- paste("INSERT INTO", mot_tab_nm, "VALUES (:data)")
DBI::dbExecute(fishset_db, 
               paste("CREATE TABLE IF NOT EXISTS", mot_tab_nm, "(data ModelOut)"))
DBI::dbExecute(fishset_db, raw_sql, 
               params = list(data = list(serialize(mot_to_save, NULL))))

# Save ModelOut table to outputs folder

MCM <- ModelOut[["MCM"]]
saveRDS(MCM, file=here::here("data", "confidential", "FishSETfolder", "eureka", "output", "MCM.rds"))

params <- ModelOut[["params"]]
saveRDS(params, file=here::here("data", "confidential", "FishSETfolder", "eureka", "output", "params.rds"))













#### KSM Updates


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
          file = here::here("data", "confidential", "FishSETfolder", "eureka", "output", "MCM_mod1.csv"),
          row.names = FALSE)

params <- results[["results"]]
write.csv(params, 
          file=here::here("data", "confidential", "FishSETfolder", "eureka", "output", "params_mod1.csv"))







