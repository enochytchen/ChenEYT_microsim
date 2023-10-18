## Filename: 00_master
## Purpose: List the order of the syntaxes

#' Set the wd as where this R file is
if (require("rstudioapi") && isAvailable()) {
  wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(wd)
}

##############################################################
##============================================================
## Install/read required packages
##============================================================
##############################################################
source("00_packages.R") 

##############################################################
##============================================================
## Read data
##============================================================
##############################################################
source("01a_fischer2016.R")   # CLL-8 trial: 8 years of follow-up 
source("01a_williams2017.R")  # Extrapolation from Williams et al.
source("01b_popmort.R")       # Population mortality file
source("02_msmcancer.R")      # Claim 4 years of survival data into a multistate structure 

##############################################################
##============================================================
## Relative survival framework
##============================================================
##############################################################
source("03_rel_survmod.R")      # Fit survival models
source("04_rel_survextrap.R")   # Evaluate survival extrapolation
source("04b_table_survextrap")  # Tabulate the observed and the extrapolated survival
source("04c_figure_survextrap") # Plot the observed and the extrapolated survival
source("05_rel_microsim.R")     # Microsimulation model for cost-effectiveness analysis
source("06_rel_onewaysen.R")    # One-way sensitivity analysis
source("07_rel_PSA.R")          # Probabilistic sensitivity analysis 
source("08_rel_CEplane_CEAC.R") # Plot cost-effectiveness plane and cost-effectiveness acceptability curve
