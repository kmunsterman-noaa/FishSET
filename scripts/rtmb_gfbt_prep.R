# ==============================================================================
# File: rtmb_gfbt_prep.R
# Desc: FishSET model preparation for RTMB script
# Author: Katrina Sky Munsterman
# ==============================================================================

# LOAD LIBRARIES ---------------------------------------------------------------

devtools::load_all()
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

# LOAD UNIVERSAL DATA ----------------------------------------------------------

# Set working directory
setwd("~/Documents/GitHub/FishSET/data/confidential/FishSETFolder")

#   5x5 km grid cells
#   VMS-pipleline::Blake Feist
spat <- "~/Documents/GitHub/FishSET/data/non-confidential/shape_files/5km_grid/master_5km_grid_tmer.shp"

#   Port Coordinates
#   PacFIN
ports <- "~/Documents/GitHub/FishSET/data/non-confidential/other/port_coords.csv"

# ==============================================================================
# astoria -----------------------------------------------------------------------
# ==============================================================================

# Set project variables

project <- "astoria"

# LOAD DATA --------------------------------------------------------------------

#   Main Data
#   fisher-behavior-displacement::fishet_prep.R

astoria_data <- "~/Documents/GitHub/FishSET/data/confidential/rds/iopac_port/ASTORIA.rds"
update_folderpath()
load_maindata(astoria_data, project = "astoria", over_write = TRUE)
astoriaMainDataTable <- table_view("astoriaMainDataTable", 
                                  project = "astoria")

#   Spatial Data
#   5x5 km grid

load_spatial(spat, name = "5x5", project = "astoria")
astoria5x5SpatTable <- table_view("astoria5x5SpatTable",
                                 project = "astoria")

#   Port Data
#   Port coordinates

load_port(ports, port_name = "port_code", project = "astoria")
astoriaPortTable <- table_view("astoriaPortTable",
                              project = "astoria")

# DATA PREP --------------------------------------------------------------------

# Scale revenue data to tens

astoriaMainDataTable <- create_var_num(dat = astoriaMainDataTable, 
                                      project = project, 
                                      x = "tow_r",
                                      y = 1000, 
                                      method = 'division', 
                                      name = 'tow_r_thousands')

# Assign zone ID for primary data

astoriaMainDataTable <- assignment_column(dat = astoriaMainDataTable, 
                                          project = project, 
                                          spat = astoria5x5SpatTable,
                                          lon.dat = "centro_lon", 
                                          lat.dat = "centro_lat", 
                                          cat = "GRID5KM_ID", 
                                          name = "astoria_ZoneID")

# Plot zone summary

zone_summary(dat = astoriaMainDataTable, 
             spat = astoria5x5SpatTable, 
             project = project, 
             zone.dat = "astoria_ZoneID", 
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
                                           zoneID = "astoria_ZoneID", 
                                           spatID = "GRID5KM_ID", 
                                           name = "start_loc")

# DATA QA/QC -------------------------------------------------------------------

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
                          zoneID = "astoria_ZoneID", 
                          zone.cent.name = "astoriaZoneCentroid")

z_ind <- which(alt_choice_list(project)$dataZoneTrue == 1)

zOut <- zone_summary(dat = astoriaMainDataTable[z_ind, ], 
                     spat = astoria5x5SpatTable, 
                     project = project, 
                     zone.dat = "astoria_ZoneID",
                     zone.spat = "GRID5KM_ID", 
                     output = "tab_plot")

zOut$table
zOut$plot

# Create expected catch matrices

create_expectations(dat = astoriaMainDataTable, 
                    project = project,
                    name = "exp1",
                    catch = "tow_r_thousands",
                    temp_var = "date_time", 
                    temp_window = 10, 
                    day_lag = 1, 
                    year_lag = 0,
                    temporal = 'daily', 
                    empty_catch = NA, 
                    empty_expectation = 1e-14,
                    dummy_exp = FALSE)

# Check model data

check_model_data(astoriaMainDataTable, 
                 project = project, 
                 uniqueID = "haul_id", 
                 latlon = c("centro_lon","centro_lat"))

# MODEL DESIGN -----------------------------------------------------------------

# Make the model design for conditional logit model

make_model_design(project = project, 
                  catchID = "tow_lbs_thousands",
                  likelihood = "logit_c", 
                  initparams = c(0,0),
                  startloc = "start_loc",
                  mod.name = "logit_c_mod1", 
                  expectcatchmodels = list('individual'))

