# SET UP ---------------------------------------------------------------------------------------

# Load packages

# Install the FishSET package (Note: Local Install worked best -- reinstalled 09/2025)
# install.packages("~/Documents/FishSET/FishSET-1.1.0.tar.gz", repos = NULL, type = "source")
# devtools::install_github("noaa-nwfsc/FishSET")

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

# Set working directory
setwd("~/Documents/GitHub/FishSET/data/confidential/FishSETFolder")

# Load spatial data
# (1) 5x5km grid cell

spat <- "~/Documents/GitHub/FishSET/data/non-confidential/shape_files/5km_grid/master_5km_grid_tmer.shp"

# (2) Port coordinates
ports <- "~/Documents/GitHub/FishSET/data/non-confidential/other/port_coords.csv"

# (2) WEA Spatial Scenarios
scen_1 <- "~/Documents/GitHub/FishSET/data/non-confidential/shape_files/weas/OSW_Scen_1.rds"
scen_2 <- "~/Documents/GitHub/FishSET/data/non-confidential/shape_files/weas/OSW_Scen_2.rds"

# eureka -----------------------------------------------------------------------------------------

# Set project variables
project <- "eureka"

# LOAD DATA ---------------------------------------------------------------------------------------

# Load main data
eureka_data <- "~/Documents/GitHub/FishSET/data/confidential/rds/iopac_port/EUREKA.rds"
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

# DATA PREP ---------------------------------------------------------------------------------------
# Scale catch data to tens
eurekaMainDataTable <- create_var_num(dat = eurekaMainDataTable, 
                                        project = project, 
                                        x = "tow_lb",
                                        y = 1000, 
                                        method = 'division', 
                                        name = 'tow_lbs_thousands')

# Assign zone ID for primary data
eurekaMainDataTable <- assignment_column(dat = eurekaMainDataTable, 
                                           project = project, 
                                           spat = eureka5x5SpatTable,
                                           lon.dat = "centro_lon", 
                                           lat.dat = "centro_lat", 
                                           cat = "GRID5KM_ID", 
                                           name = "new_ZoneID")

# Plot zone summary
zone_summary(dat = eurekaMainDataTable, 
                    spat = eureka5x5SpatTable, 
                    project = project, 
                    zone.dat = "new_ZoneID", 
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
                                            zoneID = "new_ZoneID", 
                                            spatID = "GRID5KM_ID", 
                                            name = "start_loc")

# DATA QA/QC --------------------------------------------------------------------------------------

# Check NAs 
eurekaMainDataTable <- na_filter(eurekaMainDataTable, 
                                   project = project, 
                                   x = "tow_lb",
                                   remove = TRUE,
                                   over_write = TRUE)

# Alternative choice
create_alternative_choice(dat = eurekaMainDataTable, 
                          project = project, 
                          occasion = "zonal centroid", 
                          occasion_var = "start_loc",
                          alt_var = "zonal centroid", 
                          min.haul = 1, 
                          zoneID = "new_ZoneID", 
                          zone.cent.name = "eurekaZoneCentroid")

z_ind <- which(alt_choice_list(project)$dataZoneTrue == 1)

zOut <- zone_summary(dat = eurekaMainDataTable[z_ind, ], 
                     spat = eureka5x5SpatTable, 
                     project = project, 
                     zone.dat = "new_ZoneID",
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

# Create expectations
create_expectations(dat = eurekaMainDataTable, 
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

# Create expectations
create_expectations(dat = eurekaMainDataTable, 
                    project = project, 
                    catch = "tow_lbs_thousands",
                    temp.var = "date_time", 
                    temp.window = 30, 
                    temp.lag = 1, 
                    year.lag = 0,
                    temporal = 'daily', 
                    empty.catch = NA, 
                    empty.expectation = 1e-14,
                    default.exp = FALSE, 
                    replace.output = TRUE)

# Create expectations
create_expectations(dat = eurekaMainDataTable, 
                    project = project, 
                    catch = "tow_lbs_thousands",
                    temp.var = "date_time", 
                    temp.window = 180, 
                    temp.lag = 1, 
                    year.lag = 0,
                    temporal = 'daily', 
                    empty.catch = NA, 
                    empty.expectation = 1e-14,
                    default.exp = FALSE, 
                    replace.output = TRUE)

# Create expectations
create_expectations(dat = eurekaMainDataTable, 
                    project = project, 
                    catch = "tow_lbs_thousands",
                    temp.var = "date_time", 
                    temp.window = 365, 
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
                  expectcatchmodels = list('individual'))

# RUN MODEL ---------------------------------------------------------------------------------------
discretefish_subroutine(project = project, 
                        explorestarts = TRUE,
                        run = "logit_c_mod1")

# MODEL OUTPUTS -----------------------------------------------------------------------------------
model_params(project = project, 
             output = "print")


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
                          min.haul = 15, 
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

# Create expectations
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
                  mod.name = "logit_c_mod3", 
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
                          zone.cent.name = "morroZoneCentroid")

z_ind <- which(alt_choice_list(project)$dataZoneTrue == 1)

zOut <- zone_summary(dat = morroMainDataTable[z_ind, ], 
                     spat = morro5x5SpatTable, 
                     project = project, 
                     zone.dat = "new_ZoneID",
                     zone.spat = "GRID5KM_ID", 
                     output = "tab_plot")

zOut$table
zOut$plot

# Create expectations
create_expectations(dat = morroMainDataTable, 
                    project = project, 
                    catch = "tow_lbs_thousands",
                    temp.var = "date_time", 
                    temp.window = 180, 
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
                  mod.name = "logit_c_mod3", 
                  expectcatchmodels = list('individual')) 

# RUN MODEL ---------------------------------------------------------------------------------------
discretefish_subroutine(project = project, 
                        explorestarts = TRUE,
                        run = "logit_c_mod3")

# MODEL OUTPUTS -----------------------------------------------------------------------------------
model_params(project = project, 
             output = "print")



