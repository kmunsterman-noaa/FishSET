# Groundfish bottom trawl

# EUREKA ##########################################################################################
# Set project variables
project <- "gfbt_erk"

## LOAD DATA --------------------------------------------------------------------------------------
# Load main data
erk_data <- "C:\\Users\\Paul.Carvalho\\Documents\\fishset_projects\\groundfish\\GFBTdata\\EUREKA.rds"
load_maindata(erk_data, project = "gfbt_erk", over_write = TRUE) # noaa laptop
gfbt_erkMainDataTable <- table_view(paste0(project, "MainDataTable"), 
                                    project = project)

# Load spatial data
erk_spat <- "C:\\Users\\Paul.Carvalho\\Documents\\fishset_projects\\groundfish\\master_5km_grid_tmer.shp"
load_spatial(erk_spat, name = "spat", project = project) # noaa laptop
gfbt_erkspatSpatTable <- table_view(paste0(project, "spatSpatTable"), 
                                    project = project)

# Load port data
gfbt_ports <- "C:\\Users\\Paul.Carvalho\\Documents\\fishset_projects\\groundfish\\port_coords.csv"
load_port(gfbt_ports, port_name = "port_code", project = project)
gfbt_erkPortTable <- table_view(paste0(project, "PortTable"),
                                project = project)

## DATA PREP --------------------------------------------------------------------------------------
# Scale catch data to tens
gfbt_erkMainDataTable <- create_var_num(dat = gfbt_erkMainDataTable, 
                                        project = project, 
                                        x = "tow_lb",
                                        y = 1000, 
                                        method = 'division', 
                                        name = 'tow_lbs_thousands')

# Assign zone ID for primary data
gfbt_erkMainDataTable <- assignment_column(dat = gfbt_erkMainDataTable, 
                                           project = project, 
                                           spat = gfbt_erkspatSpatTable,
                                           lon.dat = "centro_lon", 
                                           lat.dat = "centro_lat", 
                                           cat = "GRID5KM_ID", 
                                           name = "new_ZoneID")

# Plot zone summary
tmp <- zone_summary(dat = gfbt_erkMainDataTable, 
                    spat = gfbt_erkspatSpatTable, 
                    project = project, 
                    zone.dat = "new_ZoneID", 
                    zone.spat = "GRID5KM_ID", 
                    output = "plot")

# Create centroid
create_centroid(spat = gfbt_erkspatSpatTable, 
                project = project, 
                spatID = "GRID5KM_ID", 
                type = "zonal centroid", 
                output = "centroid table")

# Create starting location
gfbt_erkMainDataTable <- create_startingloc(dat = gfbt_erkMainDataTable, 
                                            project = project, 
                                            spat = gfbt_erkspatSpatTable, 
                                            port = gfbt_erkPortTable, 
                                            port_name = "Port_Name",
                                            port_lon = "Port_Long", 
                                            port_lat = "Port_Lat", 
                                            trip_id = "trip_id", 
                                            haul_order = "haul_counter", 
                                            starting_port = "depart_port", 
                                            zoneID = "new_ZoneID", 
                                            spatID = "GRID5KM_ID", 
                                            name = "start_loc1")

## DATA QA/QC -------------------------------------------------------------------------------------
# Check NAs 
gfbt_erkMainDataTable <- na_filter(gfbt_erkMainDataTable, 
                                   project = project, 
                                   x = "tow_lb",
                                   remove = TRUE,
                                   over_write = TRUE)

# Alternative choice
create_alternative_choice(dat = gfbt_erkMainDataTable, 
                          project = project, 
                          occasion = "zonal centroid", 
                          occasion_var = "start_loc1",
                          alt_var = "zonal centroid", 
                          min.haul = 5, 
                          zoneID = "new_ZoneID", 
                          zone.cent.name = "gfbt_erkZoneCentroid")

z_ind <- which(alt_choice_list(project)$dataZoneTrue == 1)

zOut <- zone_summary(dat = gfbt_erkMainDataTable[z_ind, ], 
                     spat = gfbt_erkspatSpatTable, 
                     project = project, 
                     zone.dat = "new_ZoneID",
                     zone.spat = "GRID5KM_ID", 
                     output = "tab_plot")

zOut$table
zOut$plot

# Create expectations
create_expectations(dat = gfbt_erkMainDataTable, 
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

# check model data
check_model_data(gfbt_erkMainDataTable, 
                 project = project, 
                 uniqueID = "haul_id", 
                 latlon = c("centro_lon","centro_lat"))

## MODEL DESIGN -----------------------------------------------------------------------------------
# Make the model design for conditional logit model
make_model_design(project = project, 
                  catchID = "tow_lbs_thousands",
                  likelihood = "logit_c", 
                  initparams = c(0,0),
                  startloc = "start_loc1",
                  mod.name = "logit_c_mod1", 
                  expectcatchmodels = list('individual'))

## RUN MODEL --------------------------------------------------------------------------------------
discretefish_subroutine(project = project, 
                        explorestarts = TRUE)

## MODEL OUTPUTS ----------------------------------------------------------------------------------
model_params(project = project, 
             output = "print")

## MODEL PREDICTIONS ------------------------------------------------------------------------------
# 1. plot grid
# 2. click on areas to 'close'
# 3. 'add closure'
# 4. 'save closure'
# 5. close the shiny app
zone_closure(project = project,
             spatdat = gfbt_erkspatSpatTable,
             cat = "GRID5KM_ID")

policy_data <- run_policy(project = project,
                          mod.name = "logit_c_mod1",
                          policy.name = "closure_1",
                          betadraws = 500,
                          marg_util_income = "V1",
                          income_cost = TRUE,
                          zone.dat = "new_ZoneID")

pred_prob_outputs(project = project,
                  mod.name = "logit_c_mod1",
                  policy.name = "closure_1",
                  output_option = "model_fig")

prob_map <- predict_map(project = project,
                        mod.name = "logit_c_mod1",
                        # Setting policy.name as the mod.name gives plot without closures
                        policy.name = "logit_c_mod1", 
                        spat = gfbt_erkspatSpatTable,
                        zone.spat = "GRID5KM_ID",
                        plot_type = "static") # Can also create 'dynamic' interactive maps

closure_prob_map <- predict_map(project = project,
                                mod.name = "logit_c_mod1",
                                policy.name = "closure_1", 
                                spat = gfbt_erkspatSpatTable,
                                zone.spat = "GRID5KM_ID",
                                plot_type = "static")



# ASTORIA #########################################################################################
# Set project variables
project <- "gfbt_ast"

# LOAD DATA ---------------------------------------------------------------------------------------
# Load main data
ast_data <- "C:\\Users\\Paul.Carvalho\\Documents\\fishset_projects\\groundfish\\GFBTdata\\ASTORIA.rds"
load_maindata(ast_data, project = "gfbt_ast", over_write = TRUE) # noaa laptop
gfbt_astMainDataTable <- table_view(paste0(project, "MainDataTable"), 
                                    project = project)

# Load spatial data
ast_spat <- "C:\\Users\\Paul.Carvalho\\Documents\\fishset_projects\\groundfish\\master_5km_grid_tmer.shp"
load_spatial(ast_spat, name = "spat", project = project) # noaa laptop
gfbt_astspatSpatTable <- table_view(paste0(project, "spatSpatTable"), 
                                    project = project)

# Load port data
gfbt_ports <- "C:\\Users\\Paul.Carvalho\\Documents\\fishset_projects\\groundfish\\port_coords.csv"
load_port(gfbt_ports, port_name = "port_code", project = project)
gfbt_astPortTable <- table_view(paste0(project, "PortTable"),
                                project = project)

# DATA PREP ---------------------------------------------------------------------------------------
# Scale catch data to tens
gfbt_astMainDataTable <- create_var_num(dat = gfbt_astMainDataTable, 
                                        project = project, 
                                        x = "tow_lb",
                                        y = 1000, 
                                        method = 'division', 
                                        name = 'tow_lbs_thousands')

# Assign zone ID for primary data
gfbt_astMainDataTable <- assignment_column(dat = gfbt_astMainDataTable, 
                                           project = project, 
                                           spat = gfbt_astspatSpatTable,
                                           lon.dat = "centro_lon", 
                                           lat.dat = "centro_lat", 
                                           cat = "GRID5KM_ID", 
                                           name = "new_ZoneID")

# Plot zone summary
tmp <- zone_summary(dat = gfbt_astMainDataTable, 
                    spat = gfbt_astspatSpatTable, 
                    project = project, 
                    zone.dat = "new_ZoneID", 
                    zone.spat = "GRID5KM_ID", 
                    output = "plot")

# Create centroid
create_centroid(spat = gfbt_astspatSpatTable, 
                project = project, 
                spatID = "GRID5KM_ID", 
                type = "zonal centroid", 
                output = "centroid table")

# Create starting location
gfbt_astMainDataTable <- create_startingloc(dat = gfbt_astMainDataTable, 
                                            project = project, 
                                            spat = gfbt_astspatSpatTable, 
                                            port = gfbt_astPortTable, 
                                            port_name = "Port_Name",
                                            port_lon = "Port_Long", 
                                            port_lat = "Port_Lat", 
                                            trip_id = "trip_id", 
                                            haul_order = "haul_counter", 
                                            starting_port = "depart_port", 
                                            zoneID = "new_ZoneID", 
                                            spatID = "GRID5KM_ID", 
                                            name = "start_loc1")

# DATA QA/QC --------------------------------------------------------------------------------------
# Check NAs 
gfbt_astMainDataTable <- na_filter(gfbt_astMainDataTable, 
                                   project = project, 
                                   x = c("tow_lb","tow_lbs_thousands"),
                                   remove = TRUE,
                                   over_write = TRUE)

# Alternative choice
create_alternative_choice(dat = gfbt_astMainDataTable, 
                          project = project, 
                          occasion = "zonal centroid", 
                          occasion_var = "start_loc1",
                          alt_var = "zonal centroid", 
                          min.haul = 5, 
                          zoneID = "new_ZoneID", 
                          zone.cent.name = "gfbt_astZoneCentroid")

z_ind <- which(alt_choice_list(project)$dataZoneTrue == 1)

zOut <- zone_summary(dat = gfbt_astMainDataTable[z_ind, ], 
                     spat = gfbt_astspatSpatTable, 
                     project = project, 
                     zone.dat = "new_ZoneID",
                     zone.spat = "GRID5KM_ID", 
                     output = "tab_plot")

zOut$table
zOut$plot

# Create expectations
create_expectations(dat = gfbt_astMainDataTable, 
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

# check model data
check_model_data(gfbt_astMainDataTable, 
                 project = project, 
                 uniqueID = "haul_id", 
                 latlon = c("centro_lon","centro_lat"))

# MODEL DESIGN ------------------------------------------------------------------------------------
# Make the model design for conditional logit model
make_model_design(project = project, 
                  catchID = "tow_lbs_thousands",
                  likelihood = "logit_c", 
                  initparams = c(0,0),
                  startloc = "start_loc1",
                  mod.name = "logit_c_mod1", 
                  expectcatchmodels = list('individual'))

# RUN MODEL ---------------------------------------------------------------------------------------
tic()
discretefish_subroutine(project = project, 
                        explorestarts = TRUE)
toc()

# MODEL OUTPUTS -----------------------------------------------------------------------------------
model_params(project = project, 
             output = "print")











