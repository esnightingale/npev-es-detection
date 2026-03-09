################################################################################
# Run full analysis workflow for configured country and time period
#
# 1. Edit R/config.R: set country and tp
# 2. Ensure raw POLIS data in Dropbox: .../polis/csv/<country>/
# 3. For new countries: add WorldPop rasters to data/<country>/WorldPop/
# 4. For watersheds: run R/setup/setup_watersheds.Rmd (requires Novel-T data)
# 5. Source and run this script
################################################################################

library(here)
source(here("R/config.R"))

message("Running analysis for ", country, " (", tp[1], " to ", tp[2], ")")

# Step 1: Clean raw POLIS data
message("\n--- Step 1: Cleaning POLIS data ---")
source(here("R/setup/clean_polis.R"))

# Step 2: Set up analysis data (clears env; config is re-sourced inside)
message("\n--- Step 2: Analysis data setup ---")
source(here("R/setup/analysis_data_setup.R"))

# Re-load config (analysis_data_setup clears environment)
source(here("R/config.R"))

# Step 3: Prevalence model setup
message("\n--- Step 3: Prevalence model setup ---")
source(here("R/prevalence/prev_model_setup.R"))

# Step 4: Fit prevalence model
message("\n--- Step 4: Fitting prevalence model ---")
source(here("R/prevalence/prev_model_fit.R"))

# Step 5: Extract prevalence predictions
message("\n--- Step 5: Prevalence predictions ---")
source(here("R/prevalence/prev_model_predict.R"))

# Step 6: ES detection model - requires es_npev_xy from es_data_setup.Rmd
message("\n--- Step 6: ES detection model ---")
message("Run R/detection/es_data_setup.Rmd to create es_npev_xy.rds, then:")
source(here("R/detection/es_model_fit.R"))

message("\nDone. Outputs in output/prevalence/", country, " and output/detection/", country)
