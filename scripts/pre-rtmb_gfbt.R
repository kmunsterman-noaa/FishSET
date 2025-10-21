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

# eureka -----------------------------------------------------------------------------------------

# Set project variables
project <- "eureka"

# LOAD DATA ---------------------------------------------------------------------------------------

# Load main data
eureka_data <- "~/Documents/GitHub/FishSET/data/confidential/rds/iopac_port/EUREKA.rds"
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
                                  x = "tow_lbs",
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

# Create expected catch matrices

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
                  expectcatchmodels = list('individual'))

zone_closure(project = "eureka",
             spatdat = eureka5x5SpatTable,
             cat = "GRID5KM_ID")
