# SET UP ---------------------------------------------------------------------------------------

# Load packages

# Install the FishSET package (Note: Local Install worked best -- reinstalled 09/2025)
# install.packages("~/Documents/FishSET/FishSET-1.1.0.tar.gz", repos = NULL, type = "source")
# devtools::install_github("noaa-nwfsc/FishSET")

# any updates to fishSET functions must be run and then reload the package using this::
# devtools::load_all()

library(FishSET)
library(tidyverse)
library(maps)
library(readr)
library(here)
library(sf)
library(rnaturalearth)
library(adehabitatHR)
library(mapview)
library(purrr)
library(magrittr)
library(lwgeom)
library(tictoc)
library(future)
library(furrr)

# Set working directory
setwd("~/Documents/GitHub/FishSET/data/confidential/FishSETFolder")

# Load spatial data
# (1) 5x5km grid cell
spat <- "~/Documents/GitHub/FishSET/data/non-confidential/shape_files/5km_grid/master_5km_grid_tmer.shp"

# (2) Port coordinates
ports <- "~/Documents/GitHub/FishSET/data/non-confidential/other/port_coords.csv"

# (2) WEA Spatial Scenarios
scen_1 <- "~/Documents/GitHub/FishSET/data/non-confidential/shape_files/weas/scen_1_grid/scen_1_grid.shp"
scen_2 <- "~/Documents/GitHub/FishSET/data/non-confidential/shape_files/weas/scen_2_grid/scen_2_grid.shp"

# eureka -----------------------------------------------------------------------------------------

# Set project variables
project <- "eureka"

# LOAD DATA ---------------------------------------------------------------------------------------

# Load main data
eureka_data <- "~/Documents/GitHub/FishSET/data/confidential/rds/iopac_port/EUREKA.rds"
EUREKA <- readRDS("~/Documents/GitHub/FishSET/data/confidential/rds/iopac_port/EUREKA.rds")
update_folderpath()
load_maindata(eureka_data, project = "eureka", over_write = TRUE)
eurekaMainDataTable <- table_view("eurekaMainDataTable", 
                                    project = "eureka")

# Load spatial data
load_spatial(spat, name = "5x5", project = "eureka")
eureka5x5SpatTable <- table_view("eureka5x5SpatTable",
                                    project = "eureka")

# Load port data
load_port(ports, port_name = "port_code", project = "eureka")
eurekaPortTable <- table_view("eurekaPortTable",
                                project = "eureka")

# Load closure data
load_spatial(scen_1, name = "closure_1", project = "eureka")
load_spatial(scen_2, name = "closure_2", project = "eureka")

# DATA PREP ---------------------------------------------------------------------------------------
# Scale catch data to tens
eurekaMainDataTable <- create_var_num(dat = eurekaMainDataTable, 
                                        project = project, 
                                        x = "tow_lb",
                                        y = 1000, 
                                        method = 'division', 
                                        name = 'tow_lbs_thousands')
# I think I can delete this?
# Assign zone ID for primary data
#eurekaMainDataTable <- assignment_column(dat = eurekaMainDataTable, 
                                           #project = project, 
                                           #spat = eureka5x5SpatTable,
                                           #lon.dat = "centro_lon", 
                                           #lat.dat = "centro_lat", 
                                           #cat = "GRID5KM_ID", 
                                           #name = "new_ZoneID")

# Plot zone summary
zone_summary(dat = eurekaMainDataTable, 
                    spat = eureka5x5SpatTable, 
                    project = project, 
                    zone.dat = "ZoneID", 
                    zone.spat = "GRID5KM_ID", 
                    output = "plot")

# Create centroid
create_centroid(spat = eureka5x5SpatTable, 
                project = project, 
                spatID = "GRID5KM_ID", 
                type = "zonal centroid", 
                output = "centroid table")

# Create starting location
eurekaMainDataTable <- create_startingloc(dat = eurekaMainDataTable, 
                                            project = project, 
                                            spat = eureka5x5SpatTable, 
                                            port = eurekaPortTable, 
                                            port_name = "Port_Name",
                                            port_lon = "Port_Long", 
                                            port_lat = "Port_Lat", 
                                            trip_id = "trip_id", 
                                            haul_order = "haul_counter", 
                                            starting_port = "depart_port", 
                                            zoneID = "ZoneID", 
                                            spatID = "GRID5KM_ID", 
                                            name = "start_loc")

# DATA QA/QC --------------------------------------------------------------------------------------

# Check NAs 
eurekaMainDataTable <- na_filter(eurekaMainDataTable, 
                                   project = project, 
                                   x = c("tow_lbs_thousands"),
                                   remove = TRUE,
                                   over_write = TRUE)

# Alternative choice
create_alternative_choice(dat = eurekaMainDataTable, 
                          project = project, 
                          occasion = "zonal centroid", 
                          occasion_var = "start_loc",
                          alt_var = "zonal centroid", 
                          min.haul = 1, 
                          zoneID = "ZoneID", 
                          zone.cent.name = "eurekaZoneCentroid")

z_ind <- which(alt_choice_list(project)$dataZoneTrue == 1)

zOut <- zone_summary(dat = eurekaMainDataTable[z_ind, ], 
                     spat = eureka5x5SpatTable, 
                     project = project, 
                     zone.dat = "ZoneID",
                     zone.spat = "GRID5KM_ID", 
                     output = "tab_plot")

zOut$table
zOut$plot

# Create expectations
create_expectations(dat = eurekaMainDataTable, 
                    project = project, 
                    catch = "tow_lbs_thousands",
                    temp.var = "date_time", 
                    temp.window = 7, 
                    temp.lag = 1, 
                    year.lag = 0,
                    temporal = 'daily', 
                    empty.catch = NA, 
                    empty.expectation = 1e-14,
                    default.exp = FALSE, 
                    replace.output = TRUE)


# Check model data
check_model_data(eurekaMainDataTable, 
                 project = project, 
                 uniqueID = "haul_id", 
                 latlon = c("centro_lon","centro_lat"))

# MODEL DESIGN ------------------------------------------------------------------------------------

# Make the model design for conditional logit model
make_model_design(project = project, 
                  catchID = "tow_lbs_thousands",
                  likelihood = "logit_c", 
                  initparams = c(0,0),
                  startloc = "start_loc",
                  mod.name = "logit_c_mod1",
                  vars1 = "price_per_km",
                  expectcatchmodels = list('individual'))

# Design closure
FishSET::zone_closure(project = project,
             spatdat = eureka5x5SpatTable,
             cat = "GRID5KM_ID")

# STOP HERE AND RUN RTMB

# RUN MODEL ---------------------------------------------------------------------------------------
tic()
discretefish_subroutine(project = project, 
                        explorestarts = TRUE,
                        run = "logit_c_mod1")
toc()

# MODEL OUTPUTS -----------------------------------------------------------------------------------
model_params(project = project, 
             output = "print")

# CROSS VALIDATION --------------------------------------------------------------------------------
cross_validation(
  project,
  mod.name="logit_c_mod1",
  zone.dat = "new_ZoneID",
  groups = "Observations",
  k = 5,
  time_var = NULL,
  use.scalers = FALSE,
  scaler.func = NULL
)

# astoria -----------------------------------------------------------------------------------------

# Set project variables
project <- "astoria"

# LOAD DATA ---------------------------------------------------------------------------------------

# Load main data
astoria_data <- "~/Documents/GitHub/FishSET/data/confidential/rds/iopac_port/ASTORIA.rds"
load_maindata(astoria_data, project = "astoria", over_write = TRUE)
astoriaMainDataTable <- table_view("astoriaMainDataTable", 
                                  project = "astoria")

# Load spatial data
load_spatial(spat, name = "5x5", project = "astoria")
astoria5x5SpatTable <- table_view("astoria5x5SpatTable",
                                 project = "astoria")

# Load port data
load_port(ports, port_name = "port_code", project = "astoria")
astoriaPortTable <- table_view("astoriaPortTable",
                              project = "astoria")

# DATA PREP ---------------------------------------------------------------------------------------
# Scale catch data to tens
astoriaMainDataTable <- create_var_num(dat = astoriaMainDataTable, 
                                      project = project, 
                                      x = "tow_lb",
                                      y = 1000, 
                                      method = 'division', 
                                      name = 'tow_lbs_thousands')

# Assign zone ID for primary data
astoriaMainDataTable <- assignment_column(dat = astoriaMainDataTable, 
                                         project = project, 
                                         spat = astoria5x5SpatTable,
                                         lon.dat = "centro_lon", 
                                         lat.dat = "centro_lat", 
                                         cat = "GRID5KM_ID", 
                                         name = "new_ZoneID")

# Plot zone summary
zone_summary(dat = astoriaMainDataTable, 
             spat = astoria5x5SpatTable, 
             project = project, 
             zone.dat = "new_ZoneID", 
             zone.spat = "GRID5KM_ID", 
             output = "plot")

# Create centroid
create_centroid(spat = astoria5x5SpatTable, 
                project = project, 
                spatID = "GRID5KM_ID", 
                type = "zonal centroid", 
                output = "centroid table")

# Create starting location
astoriaMainDataTable <- create_startingloc(dat = astoriaMainDataTable, 
                                          project = project, 
                                          spat = astoria5x5SpatTable, 
                                          port = astoriaPortTable, 
                                          port_name = "Port_Name",
                                          port_lon = "Port_Long", 
                                          port_lat = "Port_Lat", 
                                          trip_id = "trip_id", 
                                          haul_order = "haul_counter", 
                                          starting_port = "depart_port", 
                                          zoneID = "new_ZoneID", 
                                          spatID = "GRID5KM_ID", 
                                          name = "start_loc")

# DATA QA/QC --------------------------------------------------------------------------------------

# Check NAs 
astoriaMainDataTable <- na_filter(astoriaMainDataTable, 
                                 project = project, 
                                 x = "tow_lb",
                                 remove = TRUE,
                                 over_write = TRUE)

# Alternative choice
create_alternative_choice(dat = astoriaMainDataTable, 
                          project = project, 
                          occasion = "zonal centroid", 
                          occasion_var = "start_loc",
                          alt_var = "zonal centroid", 
                          min.haul = 1, 
                          zoneID = "new_ZoneID", 
                          zone.cent.name = "astoriaZoneCentroid")

z_ind <- which(alt_choice_list(project)$dataZoneTrue == 1)

zOut <- zone_summary(dat = astoriaMainDataTable[z_ind, ], 
                     spat = astoria5x5SpatTable, 
                     project = project, 
                     zone.dat = "new_ZoneID",
                     zone.spat = "GRID5KM_ID", 
                     output = "tab_plot")

zOut$table
zOut$plot

# Create expected catch matrices

create_expectations(dat = astoriaMainDataTable, 
                    project = project, 
                    catch = "tow_lbs_thousands",
                    temp.var = "date_time", 
                    temp.window = 7, 
                    temp.lag = 1, 
                    year.lag = 0,
                    temporal = 'daily', 
                    empty.catch = NA, 
                    empty.expectation = 1e-14,
                    default.exp = FALSE, 
                    replace.output = TRUE)

# Check model data
check_model_data(astoriaMainDataTable, 
                 project = project, 
                 uniqueID = "haul_id", 
                 latlon = c("centro_lon","centro_lat"))

# MODEL DESIGN ------------------------------------------------------------------------------------

# Make the model design for conditional logit model
make_model_design(project = project, 
                  catchID = "tow_lbs_thousands",
                  likelihood = "logit_c", 
                  initparams = c(0,0),
                  startloc = "start_loc",
                  mod.name = "logit_c_mod1", 
                  expectcatchmodels = list('individual'))

# RUN MODEL ---------------------------------------------------------------------------------------

discretefish_subroutine(project = project, 
                        explorestarts = TRUE,
                        run = "logit_c_mod3")

# MODEL OUTPUTS -----------------------------------------------------------------------------------
model_params(project = project, 
             output = "print")

# morro bay ----------------------------------------------------------------------------------------

# Set project variables
project <- "morro"

# LOAD DATA ---------------------------------------------------------------------------------------

# Load main data
morro_data <- "~/Documents/GitHub/FishSET/data/confidential/rds/iopac_port/MORRO.rds"
load_maindata(morro_data, project = "morro", over_write = TRUE)
morroMainDataTable <- table_view("morroMainDataTable", 
                                   project = "morro")

# Load spatial data
load_spatial(spat, name = "5x5", project = "morro")
morro5x5SpatTable <- table_view("morro5x5SpatTable",
                                  project = "morro")

# Load port data
load_port(ports, port_name = "port_code", project = "morro")
morroPortTable <- table_view("morroPortTable",
                               project = "morro")

# DATA PREP ---------------------------------------------------------------------------------------
# Scale catch data to tens
morroMainDataTable <- create_var_num(dat = morroMainDataTable, 
                                       project = project, 
                                       x = "tow_lb",
                                       y = 1000, 
                                       method = 'division', 
                                       name = 'tow_lbs_thousands')

# Assign zone ID for primary data
morroMainDataTable <- assignment_column(dat = morroMainDataTable, 
                                          project = project, 
                                          spat = morro5x5SpatTable,
                                          lon.dat = "centro_lon", 
                                          lat.dat = "centro_lat", 
                                          cat = "GRID5KM_ID", 
                                          name = "new_ZoneID")

# Plot zone summary
zone_summary(dat = morroMainDataTable, 
             spat = morro5x5SpatTable, 
             project = project, 
             zone.dat = "new_ZoneID", 
             zone.spat = "GRID5KM_ID", 
             output = "plot")

# Create centroid
create_centroid(spat = morro5x5SpatTable, 
                project = project, 
                spatID = "GRID5KM_ID", 
                type = "zonal centroid", 
                output = "centroid table")

# Create starting location
morroMainDataTable <- create_startingloc(dat = morroMainDataTable, 
                                           project = project, 
                                           spat = morro5x5SpatTable, 
                                           port = morroPortTable, 
                                           port_name = "Port_Name",
                                           port_lon = "Port_Long", 
                                           port_lat = "Port_Lat", 
                                           trip_id = "trip_id", 
                                           haul_order = "haul_counter", 
                                           starting_port = "depart_port", 
                                           zoneID = "new_ZoneID", 
                                           spatID = "GRID5KM_ID", 
                                           name = "start_loc")

# DATA QA/QC --------------------------------------------------------------------------------------

# Check NAs 
morroMainDataTable <- na_filter(morroMainDataTable, 
                                  project = project, 
                                  x = "tow_lb",
                                  remove = TRUE,
                                  over_write = TRUE)

# Alternative choice
create_alternative_choice(dat = morroMainDataTable, 
                          project = project, 
                          occasion = "zonal centroid", 
                          occasion_var = "start_loc",
                          alt_var = "zonal centroid", 
                          min.haul = 1, 
                          zoneID = "new_ZoneID", 
                          zone.cent.name = "morroZoneCentroid",
                          spatname = morro5x5SpatTable,
                          spatID = "GRID5KM_ID")

z_ind <- which(alt_choice_list(project)$dataZoneTrue == 1)

zOut <- zone_summary(dat = morroMainDataTable[z_ind, ], 
                     spat = morro5x5SpatTable, 
                     project = project, 
                     zone.dat = "new_ZoneID",
                     zone.spat = "GRID5KM_ID", 
                     output = "tab_plot")

zOut$table
zOut$plot

# Create expected catch matrices
# Run in parallel across multiple cores

# Create a vector of the 'temp.window' values to iterate over
temp.window_values <- c(7, 14, 30, 90, 180)

# Configure parallel plan - with one core free for OS to run smoothly
workers_to_use <- parallel::detectCores() - 1
plan(multicore, workers = workers_to_use)

# Run the function in parallel
message(paste("Starting parallel processing across", workers_to_use, "cores..."))

tic()
future_walk(temp.window_values, ~ {
  
  # makes sure  FishSET package is loaded on each parallel worker
  library(FishSET)
  
  message(paste("Processing temp.window =", .x))
  
  create_expectations(
    dat = morroMainDataTable,
    project = project,
    catch = "tow_lbs_thousands",
    temp.var = "date_time",
    temp.window = .x,
    temp.lag = 1,
    year.lag = 0,
    temporal = 'daily',
    empty.catch = NA,
    empty.expectation = 1e-14,
    default.exp = FALSE,
    replace.output = FALSE
  )
  
  message(paste("Finished temp.window =", .x))
  
})
toc()


##############################################

# Create expectations
create_expectations(dat = morroMainDataTable, 
                    project = project, 
                    catch = "tow_lbs_thousands",
                    temp.var = "date_time", 
                    temp.window = 14, 
                    temp.lag = 1, 
                    year.lag = 0,
                    temporal = 'daily', 
                    empty.catch = NA, 
                    empty.expectation = 1e-14,
                    default.exp = FALSE, 
                    replace.output = TRUE)


# Check model data
check_model_data(morroMainDataTable, 
                 project = project, 
                 uniqueID = "haul_id", 
                 latlon = c("centro_lon","centro_lat"))

# MODEL DESIGN ------------------------------------------------------------------------------------

# Make the model design for conditional logit model
make_model_design(project = project, 
                  catchID = "tow_lbs_thousands",
                  likelihood = "logit_c", 
                  initparams = c(0,0),
                  startloc = "start_loc",
                  mod.name = "logit_c_mod1", 
                  expectcatchmodels = list('individual')) 

# RUN MODEL ---------------------------------------------------------------------------------------
discretefish_subroutine(project = project, 
                        explorestarts = TRUE,
                        run = "logit_c_mod1")

# MODEL OUTPUTS -----------------------------------------------------------------------------------
model_params(project = project, 
             output = "print")


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

clean_dates <- morroMainDataTable %>%
mutate(cleaned_dates = sub(" .*", "\\1", date_time))

# crescent city ----------------------------------------------------------------------------------------

# Set project variables
project <- "crescent"

# LOAD DATA ---------------------------------------------------------------------------------------

# Load main data
crescent_data <- "~/Documents/GitHub/FishSET/data/confidential/rds/iopac_port/CRESCENT CITY.rds"
load_maindata(crescent_data, project = "crescent", over_write = TRUE)
crescentMainDataTable <- table_view("crescentMainDataTable", 
                                 project = "crescent")

# Load spatial data
load_spatial(spat, name = "5x5", project = "crescent")
crescent5x5SpatTable <- table_view("crescent5x5SpatTable",
                                project = "crescent")

# Load port data
load_port(ports, port_name = "port_code", project = "crescent")
crescentPortTable <- table_view("crescentPortTable",
                             project = "crescent")

# DATA PREP ---------------------------------------------------------------------------------------
# Scale catch data to tens
crescentMainDataTable <- create_var_num(dat = crescentMainDataTable, 
                                     project = project, 
                                     x = "tow_lb",
                                     y = 1000, 
                                     method = 'division', 
                                     name = 'tow_lbs_thousands')

# Assign zone ID for primary data
crescentMainDataTable <- assignment_column(dat = crescentMainDataTable, 
                                        project = project, 
                                        spat = crescent5x5SpatTable,
                                        lon.dat = "centro_lon", 
                                        lat.dat = "centro_lat", 
                                        cat = "GRID5KM_ID", 
                                        name = "new_ZoneID")

# Plot zone summary
zone_summary(dat = crescentMainDataTable, 
             spat = crescent5x5SpatTable, 
             project = project, 
             zone.dat = "new_ZoneID", 
             zone.spat = "GRID5KM_ID", 
             output = "plot")

# Create centroid
create_centroid(spat = crescent5x5SpatTable, 
                project = project, 
                spatID = "GRID5KM_ID", 
                type = "zonal centroid", 
                output = "centroid table")

# Create starting location
crescentMainDataTable <- create_startingloc(dat = crescentMainDataTable, 
                                         project = project, 
                                         spat = crescent5x5SpatTable, 
                                         port = crescentPortTable, 
                                         port_name = "Port_Name",
                                         port_lon = "Port_Long", 
                                         port_lat = "Port_Lat", 
                                         trip_id = "trip_id", 
                                         haul_order = "haul_counter", 
                                         starting_port = "depart_port", 
                                         zoneID = "new_ZoneID", 
                                         spatID = "GRID5KM_ID", 
                                         name = "start_loc")

# DATA QA/QC --------------------------------------------------------------------------------------

# Check NAs 
crescentMainDataTable <- na_filter(crescentMainDataTable, 
                                project = project, 
                                x = "tow_lb",
                                remove = TRUE,
                                over_write = TRUE)

# Alternative choice
create_alternative_choice(dat = crescentMainDataTable, 
                          project = project, 
                          occasion = "zonal centroid", 
                          occasion_var = "start_loc",
                          alt_var = "zonal centroid", 
                          min.haul = 1, 
                          zoneID = "new_ZoneID", 
                          zone.cent.name = "crescentZoneCentroid")

z_ind <- which(alt_choice_list(project)$dataZoneTrue == 1)

zOut <- zone_summary(dat = crescentMainDataTable[z_ind, ], 
                     spat = crescent5x5SpatTable, 
                     project = project, 
                     zone.dat = "new_ZoneID",
                     zone.spat = "GRID5KM_ID", 
                     output = "tab_plot")

zOut$table
zOut$plot

# Create expected catch matrices
# Run in parallel across multiple cores

# Create a vector of the 'temp.window' values to iterate over
temp.window_values <- c(7, 14, 30, 90, 180)

# Configure parallel plan - with one core free for OS to run smoothly
workers_to_use <- parallel::detectCores() - 1
plan(multicore, workers = workers_to_use)

# Run the function in parallel
message(paste("Starting parallel processing across", workers_to_use, "cores..."))

tic()
future_walk(temp.window_values, ~ {
  
  # makes sure  FishSET package is loaded on each parallel worker
  library(FishSET)
  
  message(paste("Processing temp.window =", .x))
  
  create_expectations(
    dat = crescentMainDataTable,
    project = project,
    catch = "tow_lbs_thousands",
    temp.var = "date_time",
    temp.window = .x,
    temp.lag = 1,
    year.lag = 0,
    temporal = 'daily',
    empty.catch = NA,
    empty.expectation = 1e-14,
    default.exp = FALSE,
    replace.output = FALSE
  )
  
  message(paste("Finished temp.window =", .x))
  
})
toc()


# Check model data
check_model_data(crescentMainDataTable, 
                 project = project, 
                 uniqueID = "haul_id", 
                 latlon = c("centro_lon","centro_lat"))

# MODEL DESIGN ------------------------------------------------------------------------------------

# Make the model design for conditional logit model
make_model_design(project = project, 
                  catchID = "tow_lbs_thousands",
                  likelihood = "logit_c", 
                  initparams = c(0,0,0,0,0,0),
                  startloc = "start_loc",
                  mod.name = "logit_c_mod1", 
                  expectcatchmodels = list('individual')) 

# RUN MODEL ---------------------------------------------------------------------------------------
tic()
discretefish_subroutine(project = project, 
                        explorestarts = TRUE,
                        run = "logit_c_mod1")
toc()

# MODEL OUTPUTS -----------------------------------------------------------------------------------
model_params(project = project, 
             output = "print")

# brookings ----------------------------------------------------------------------------------------

# Set project variables
project <- "brookings"

# LOAD DATA ---------------------------------------------------------------------------------------

# Load main data
brookings_data <- "~/Documents/GitHub/FishSET/data/confidential/rds/iopac_port/BROOKINGS.rds"
load_maindata(brookings_data, project = "brookings", over_write = TRUE)
brookingsMainDataTable <- table_view("brookingsMainDataTable", 
                                    project = "brookings")

# Load spatial data
load_spatial(spat, name = "5x5", project = "brookings")
brookings5x5SpatTable <- table_view("brookings5x5SpatTable",
                                   project = "brookings")

# Load port data
load_port(ports, port_name = "port_code", project = "brookings")
brookingsPortTable <- table_view("brookingsPortTable",
                                project = "brookings")

# DATA PREP ---------------------------------------------------------------------------------------
# Scale catch data to tens
brookingsMainDataTable <- create_var_num(dat = brookingsMainDataTable, 
                                        project = project, 
                                        x = "tow_lb",
                                        y = 1000, 
                                        method = 'division', 
                                        name = 'tow_lbs_thousands')

# Assign zone ID for primary data
brookingsMainDataTable <- assignment_column(dat = brookingsMainDataTable, 
                                           project = project, 
                                           spat = brookings5x5SpatTable,
                                           lon.dat = "centro_lon", 
                                           lat.dat = "centro_lat", 
                                           cat = "GRID5KM_ID", 
                                           name = "new_ZoneID")

# Plot zone summary
zone_summary(dat = brookingsMainDataTable, 
             spat = brookings5x5SpatTable, 
             project = project, 
             zone.dat = "new_ZoneID", 
             zone.spat = "GRID5KM_ID", 
             output = "plot")

# Create centroid
create_centroid(spat = brookings5x5SpatTable, 
                project = project, 
                spatID = "GRID5KM_ID", 
                type = "zonal centroid", 
                output = "centroid table")

# Create starting location
brookingsMainDataTable <- create_startingloc(dat = brookingsMainDataTable, 
                                            project = project, 
                                            spat = brookings5x5SpatTable, 
                                            port = brookingsPortTable, 
                                            port_name = "Port_Name",
                                            port_lon = "Port_Long", 
                                            port_lat = "Port_Lat", 
                                            trip_id = "trip_id", 
                                            haul_order = "haul_counter", 
                                            starting_port = "depart_port", 
                                            zoneID = "new_ZoneID", 
                                            spatID = "GRID5KM_ID", 
                                            name = "start_loc")

# DATA QA/QC --------------------------------------------------------------------------------------

# Check NAs 
brookingsMainDataTable <- na_filter(brookingsMainDataTable, 
                                   project = project, 
                                   x = "tow_lb",
                                   remove = TRUE,
                                   over_write = TRUE)

# Alternative choice
create_alternative_choice(dat = brookingsMainDataTable, 
                          project = project, 
                          occasion = "zonal centroid", 
                          occasion_var = "start_loc",
                          alt_var = "zonal centroid", 
                          min.haul = 1, 
                          zoneID = "new_ZoneID", 
                          zone.cent.name = "brookingsZoneCentroid")

z_ind <- which(alt_choice_list(project)$dataZoneTrue == 1)

zOut <- zone_summary(dat = brookingsMainDataTable[z_ind, ], 
                     spat = brookings5x5SpatTable, 
                     project = project, 
                     zone.dat = "new_ZoneID",
                     zone.spat = "GRID5KM_ID", 
                     output = "tab_plot")

zOut$table
zOut$plot

# Create expected catch matrices
# Run in parallel across multiple cores

# Create a vector of the 'temp.window' values to iterate over
temp.window_values <- c(7, 14, 30, 90, 180)

# Configure parallel plan - with one core free for OS to run smoothly
workers_to_use <- parallel::detectCores() - 1
plan(multicore, workers = workers_to_use)

# Run the function in parallel
message(paste("Starting parallel processing across", workers_to_use, "cores..."))

tic()
future_walk(temp.window_values, ~ {
  
  # makes sure  FishSET package is loaded on each parallel worker
  library(FishSET)
  
  message(paste("Processing temp.window =", .x))
  
  create_expectations(
    dat = brookingsMainDataTable,
    project = project,
    catch = "tow_lbs_thousands",
    temp.var = "date_time",
    temp.window = .x,
    temp.lag = 1,
    year.lag = 0,
    temporal = 'daily',
    empty.catch = NA,
    empty.expectation = 1e-14,
    default.exp = FALSE,
    replace.output = FALSE
  )
  
  message(paste("Finished temp.window =", .x))
  
})
toc()


# Check model data
check_model_data(brookingsMainDataTable, 
                 project = project, 
                 uniqueID = "haul_id", 
                 latlon = c("centro_lon","centro_lat"))

# MODEL DESIGN ------------------------------------------------------------------------------------

# Make the model design for conditional logit model
make_model_design(project = project, 
                  catchID = "tow_lbs_thousands",
                  likelihood = "logit_c", 
                  initparams = c(0,0,0,0,0,0),
                  startloc = "start_loc",
                  mod.name = "logit_c_mod1", 
                  expectcatchmodels = list('individual')) 

# RUN MODEL ---------------------------------------------------------------------------------------
tic()
discretefish_subroutine(project = project, 
                        explorestarts = TRUE,
                        run = "logit_c_mod1")
toc()

# MODEL OUTPUTS -----------------------------------------------------------------------------------
model_params(project = project, 
             output = "print")


# fort bragg ----------------------------------------------------------------------------------------

# Set project variables
project <- "fort"

# LOAD DATA ---------------------------------------------------------------------------------------

# Load main data
fort_data <- "~/Documents/GitHub/FishSET/data/confidential/rds/iopac_port/FORT BRAGG.rds"
load_maindata(fort_data, project = "fort", over_write = TRUE)
fortMainDataTable <- table_view("fortMainDataTable", 
                                     project = "fort")

# Load spatial data
load_spatial(spat, name = "5x5", project = "fort")
fort5x5SpatTable <- table_view("fort5x5SpatTable",
                                    project = "fort")

# Load port data
load_port(ports, port_name = "port_code", project = "fort")
fortPortTable <- table_view("fortPortTable",
                                 project = "fort")

# DATA PREP ---------------------------------------------------------------------------------------
# Scale catch data to tens
fortMainDataTable <- create_var_num(dat = fortMainDataTable, 
                                         project = project, 
                                         x = "tow_lb",
                                         y = 1000, 
                                         method = 'division', 
                                         name = 'tow_lbs_thousands')

# Assign zone ID for primary data
fortMainDataTable <- assignment_column(dat = fortMainDataTable, 
                                            project = project, 
                                            spat = fort5x5SpatTable,
                                            lon.dat = "centro_lon", 
                                            lat.dat = "centro_lat", 
                                            cat = "GRID5KM_ID", 
                                            name = "new_ZoneID")

# Plot zone summary
zone_summary(dat = fortMainDataTable, 
             spat = fort5x5SpatTable, 
             project = project, 
             zone.dat = "new_ZoneID", 
             zone.spat = "GRID5KM_ID", 
             output = "plot")

# Create centroid
create_centroid(spat = fort5x5SpatTable, 
                project = project, 
                spatID = "GRID5KM_ID", 
                type = "zonal centroid", 
                output = "centroid table")

# Create starting location
fortMainDataTable <- create_startingloc(dat = fortMainDataTable, 
                                             project = project, 
                                             spat = fort5x5SpatTable, 
                                             port = fortPortTable, 
                                             port_name = "Port_Name",
                                             port_lon = "Port_Long", 
                                             port_lat = "Port_Lat", 
                                             trip_id = "trip_id", 
                                             haul_order = "haul_counter", 
                                             starting_port = "depart_port", 
                                             zoneID = "new_ZoneID", 
                                             spatID = "GRID5KM_ID", 
                                             name = "start_loc")

# DATA QA/QC --------------------------------------------------------------------------------------

# Check NAs 
fortMainDataTable <- na_filter(fortMainDataTable, 
                                    project = project, 
                                    x = "tow_lb",
                                    remove = TRUE,
                                    over_write = TRUE)

# Alternative choice
create_alternative_choice(dat = fortMainDataTable, 
                          project = project, 
                          occasion = "zonal centroid", 
                          occasion_var = "start_loc",
                          alt_var = "zonal centroid", 
                          min.haul = 1, 
                          zoneID = "new_ZoneID", 
                          zone.cent.name = "fortZoneCentroid")

z_ind <- which(alt_choice_list(project)$dataZoneTrue == 1)

zOut <- zone_summary(dat = fortMainDataTable[z_ind, ], 
                     spat = fort5x5SpatTable, 
                     project = project, 
                     zone.dat = "new_ZoneID",
                     zone.spat = "GRID5KM_ID", 
                     output = "tab_plot")

zOut$table
zOut$plot

# Create expected catch matrices
# Run in parallel across multiple cores

# Create a vector of the 'temp.window' values to iterate over
temp.window_values <- c(7, 14, 30, 90, 180)

# Configure parallel plan - with one core free for OS to run smoothly
workers_to_use <- parallel::detectCores() - 1
plan(multicore, workers = workers_to_use)

# Run the function in parallel
message(paste("Starting parallel processing across", workers_to_use, "cores..."))

tic()
future_walk(temp.window_values, ~ {
  
  # makes sure  FishSET package is loaded on each parallel worker
  library(FishSET)
  
  message(paste("Processing temp.window =", .x))
  
  create_expectations(
    dat = fortMainDataTable,
    project = project,
    catch = "tow_lbs_thousands",
    temp.var = "date_time",
    temp.window = .x,
    temp.lag = 1,
    year.lag = 0,
    temporal = 'daily',
    empty.catch = NA,
    empty.expectation = 1e-14,
    default.exp = FALSE,
    replace.output = FALSE
  )
  
  message(paste("Finished temp.window =", .x))
  
})
toc()


# Check model data
check_model_data(fortMainDataTable, 
                 project = project, 
                 uniqueID = "haul_id", 
                 latlon = c("centro_lon","centro_lat"))

# MODEL DESIGN ------------------------------------------------------------------------------------

# Make the model design for conditional logit model
make_model_design(project = project, 
                  catchID = "tow_lbs_thousands",
                  likelihood = "logit_c", 
                  initparams = c(0,0,0,0,0,0),
                  startloc = "start_loc",
                  mod.name = "logit_c_mod1", 
                  expectcatchmodels = list('individual')) 

# RUN MODEL ---------------------------------------------------------------------------------------
tic()
discretefish_subroutine(project = project, 
                        explorestarts = TRUE,
                        run = "logit_c_mod1")
toc()

# MODEL OUTPUTS -----------------------------------------------------------------------------------
model_params(project = project, 
             output = "print")
