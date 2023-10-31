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

################################################################
setwd(original_wd)  # Reset to the original working directory
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