################################################################################
# Run analysis workflow - setup, prevalence, and/or detection
#
# Usage:
#   source("run_analysis.R")                    # run all components
#   source("run_analysis.R", local = new.env()) # with run_* = TRUE/FALSE set below
#
#   Rscript run_analysis.R                      # run all
#   Rscript run_analysis.R setup                # setup only
#   Rscript run_analysis.R prevalence           # prevalence only
#   Rscript run_analysis.R detection            # detection only
#   Rscript run_analysis.R setup prevalence     # setup + prevalence
#
# Prerequisites: Edit R/config.R (country, tp), raw POLIS data, WorldPop rasters
################################################################################

# Parse command-line args (when run via Rscript)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  run_setup      <- TRUE
  run_prevalence <- TRUE
  run_detection  <- TRUE
} else {
  run_setup      <- "setup"     %in% args
  run_prevalence <- "prevalence" %in% args
  run_detection  <- "detection" %in% args
}

# Alternatively, set these directly when sourcing interactively:
# run_setup <- TRUE
# run_prevalence <- FALSE
# run_detection <- TRUE

library(here)
source(here("R/config.R"))

message("========================================")
message("NPEV analysis workflow: ", country)
message("Time period: ", tp[1], " to ", tp[2])
message("Components: ", 
        if (run_setup) "setup " else "",
        if (run_prevalence) "prevalence " else "",
        if (run_detection) "detection" else "")
message("========================================\n")

################################################################################
# 1. SETUP
# -----------------------------------------------------------------------------

if (run_setup) {
  message("--- Step 1: Setup ---")
  message("  Cleaning POLIS data")
  source(here("R/setup/clean_polis.R"))

  message("  Analysis data setup")
  source(here("R/setup/analysis_data_setup.R"))

  source(here("R/config.R"))  # re-load (analysis_data_setup clears env)
  # Re-set component flags (analysis_data_setup clears workspace)
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) {
    run_setup <- TRUE; run_prevalence <- TRUE; run_detection <- TRUE
  } else {
    run_setup <- "setup" %in% args
    run_prevalence <- "prevalence" %in% args
    run_detection <- "detection" %in% args
  }
}

################################################################################
# 2. PREVALENCE
# -----------------------------------------------------------------------------

if (run_prevalence) {
  message("\n--- Step 2: Prevalence ---")
  if (!run_setup) source(here("R/config.R"))

  message("  Model setup")
  source(here("R/prevalence/prev_model_setup.R"))

  message("  Fitting model")
  source(here("R/prevalence/prev_model_fit.R"))

  message("  Extracting predictions")
  source(here("R/prevalence/prev_model_predict.R"))
}

################################################################################
# 3. DETECTION
# -----------------------------------------------------------------------------

if (run_detection) {
  message("\n--- Step 3: Detection ---")
  if (!run_setup && !run_prevalence) source(here("R/config.R"))

  message("  ES data setup (link samples to district prevalence, creates es_ev_xy.rds)")
  source(here("R/setup/es_data_setup.R"))

  message("  Fitting ES detection model")
  source(here("R/detection/es_model_fit.R"))
}

################################################################################

message("\n========================================")
message("Done.")
message("Outputs: output/prevalence/", country, ", output/detection/", country)
message("========================================")
