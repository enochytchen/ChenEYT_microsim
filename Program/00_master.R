## Filename: 00_master
## Purpose: List the order of the syntaxes

## Set the wd as where this R file is
if (require("rstudioapi") && isAvailable()) {
  wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
} else wd <- "~/src/R/ChenEYT_microsim/Program/"
setwd(wd)

##############################################################
##============================================================
## Install required packages
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
## Survival models
##============================================================
##############################################################
source("03a_survmod_fpm.R")      # Fit flexible parametric models using rstpm2/flexsurv packages
source("03a_survmod_williams.R") # Fit standard parametric models used in Williams et al. 2017
source("03b_compare_hazard.R")   # Compare the extrapolated hazard functions

##############################################################
##============================================================
## Compare survival extrapolations
##============================================================
##############################################################
source("04a_extrap_microsim_fpm_ac.R")       # Microsimulation model using flexible parametric models within an all-cause survival framework
source("04a_extrap_microsim_fpm_rel.R")      # Microsimulation model using flexible parametric models within a relative survival framework
source("04a_extrap_microsim_fpm_rel(hesim)") # Microsimulation model using flexible parametric models within a relative survival framework, using hesim package
source("04a_extrap_microsim_williams_ac.R")  # Microsimulation model using standard parametric models, used in Williams et al. 2017, within an all-cause survival framework
source("04a_extrap_semiMarkov_fpm_ac.R")     # Semi-Markov model using flexible parametric models within an all-cause survival framework
source("04b_table_survextrap.R")             # Tabulate the observed and the extrapolated survival
source("04c_figure_survextrap.R")            # Plot the observed and the extrapolated survival

##############################################################
##============================================================
## Full health economics model using microsimulation
## flexible parametric models within a relative survival framework
##============================================================
##############################################################
source("05_microsim_fpm_rel.R")        # Microsimulation model for cost-effectiveness analysis
source("06_onewaysen.R")               # One-way sensitivity analysis
source("07_PSA.R")                     # Probabilistic sensitivity analysis 
source("08_CEplane_CEAC.R")            # Plot cost-effectiveness plane and cost-effectiveness acceptability curve

################################################################
# Copyright 2023 Chen EYT. All Rights Reserved.
# A microsimulation model incorporating relative survival extrapolation and 
# multiple timescales for health technology assessment
# 
# Licensed to the Apache Software Foundation (ASF) under one or more
# contributor license agreements.  See the NOTICE file distributed with
# this work for additional information regarding copyright ownership.
# The ASF licenses this file to You under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance with
# the License.  You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
