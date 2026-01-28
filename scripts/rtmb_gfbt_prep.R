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
# XXXX -----------------------------------------------------------------------
# ==============================================================================

# Set project variables

project <- "XXXX"

# LOAD DATA --------------------------------------------------------------------

#   Main Data
#   fisher-behavior-displacement::fishet_prep.R

XXXX_data <- "~/Documents/GitHub/FishSET/data/confidential/rds/iopac_return/XXXX.rds"
update_folderpath()
load_maindata(XXXX_data, project = "XXXX", over_write = TRUE)
XXXXMainDataTable <- table_view("XXXXMainDataTable", 
                                project = "XXXX")

#   Spatial Data
#   5x5 km grid

load_spatial(spat, name = "5x5", project = "XXXX")
XXXX5x5SpatTable <- table_view("XXXX5x5SpatTable",
                               project = "XXXX")

#   Port Data
#   Port coordinates

load_port(ports, port_name = "port_code", project = "XXXX")
XXXXPortTable <- table_view("XXXXPortTable",
                            project = "XXXX")

# DATA PREP --------------------------------------------------------------------

# Assign zone ID for primary data

XXXXMainDataTable <- assignment_column(dat = XXXXMainDataTable, 
                                       project = project, 
                                       spat = XXXX5x5SpatTable,
                                       lon.dat = "centro_lon", 
                                       lat.dat = "centro_lat", 
                                       cat = "GRID5KM_ID", 
                                       name = "new_zoneID")

# Plot zone summary

zone_summary(dat = XXXXMainDataTable, 
             spat = XXXX5x5SpatTable, 
             project = project, 
             zone.dat = "new_zoneID", 
             zone.spat = "GRID5KM_ID", 
             output = "plot")

# Create centroid

create_centroid(spat = XXXX5x5SpatTable, 
                project = project, 
                spatID = "GRID5KM_ID", 
                type = "zonal centroid", 
                output = "centroid table")

# Starting location

XXXXMainDataTable <- change_class(dat = XXXXMainDataTable, 
                                  project = project, 
                                  x = "startingloc", 
                                  new_class = 'character', 
                                  save = TRUE)

# DATA QA/QC -------------------------------------------------------------------

# Check NAs 
XXXXMainDataTable <- nan_filter(XXXXMainDataTable, 
                                project = project, 
                                remove = TRUE,
                                over_write = TRUE)

XXXXMainDataTable <- na_filter(XXXXMainDataTable, 
                               project = project,
                               x = "tow_r",
                               remove = TRUE,
                               over_write = TRUE)

XXXXMainDataTable <- na_filter(XXXXMainDataTable, 
                               project = project,
                               x = "tow_lb",
                               remove = TRUE,
                               over_write = TRUE)

XXXXMainDataTable <- na_filter(XXXXMainDataTable, 
                               project = project,
                               x = "port_zoneID",
                               remove = TRUE,
                               over_write = TRUE)

XXXXMainDataTable <- na_filter(XXXXMainDataTable, 
                                   project = project,
                                   x = "vessel_id",
                                   remove = TRUE,
                                   over_write = TRUE)

# Alternative choice
create_alternative_choice(dat = XXXXMainDataTable, 
                          project = project, 
                          occasion = "zonal centroid", 
                          occasion_var = "startingloc",
                          alt_var = "zonal centroid", 
                          min.haul = 1, 
                          zoneID = "new_zoneID", 
                          zone.cent.name = "XXXXZoneCentroid")

z_ind <- which(alt_choice_list(project)$dataZoneTrue == 1)

zOut <- zone_summary(dat = XXXXMainDataTable[z_ind, ], 
                     spat = XXXX5x5SpatTable, 
                     project = project, 
                     zone.dat = "new_zoneID",
                     zone.spat = "GRID5KM_ID", 
                     output = "tab_plot")

zOut$table
zOut$plot

# Create expected catch matrices

create_expectations(dat = XXXXMainDataTable, 
                    project = project,
                    name = "exp1",
                    catch = "tow_r",
                    temp_var = "date_time", 
                    temp_window = 10, 
                    day_lag = 1, 
                    year_lag = 0,
                    temporal = 'daily', 
                    empty_catch = NA, 
                    empty_expectation = 1e-14,
                    dummy_exp = FALSE)

# Check model data

check_model_data(XXXXMainDataTable, 
                 project = project, 
                 uniqueID = "haul_id", 
                 latlon = c("centro_lon","centro_lat"))

# MODEL DESIGN -----------------------------------------------------------------

# Make the model design for conditional logit model

make_model_design(project = project, 
                  catchID = "tow_r",
                  likelihood = "logit_c", 
                  initparams = c(0,0),
                  startloc = "startingloc",
                  mod.name = "logit_c_mod1", 
                  expectcatchmodels = list('individual'))
