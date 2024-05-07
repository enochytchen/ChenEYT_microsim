## Filename: 00_master
## Purpose: List the order of the syntaxes
## Reference: Williams C, Lewsey JD, Briggs AH, Mackay DF. 
##            Cost-effectiveness Analysis in R Using a Multi-state Modeling Survival Analysis Framework: A Tutorial. 
##            Med Decis Making. 2017 May;37(4):340–52.
##            Cost-effectiveness Analysis in R Using a Multi-state Modeling Survival Analysis Framework: A Tutorial © 2017 by Williams C. et al is licensed under CC BY 3.0. 

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
source("04a1_Markov_fpm_ac_microsim.R")           # Markov (clock-forward) model using flexible parametric models within an all-cause survival framework (microsimulation package)
source("04a1_Markov_fpm_ac_mstate.R")             # Markov (clock-forward) model using flexible parametric models within an all-cause survival framework (mstate package)
source("04a2_Markov_fpm_rel_microsim.R")          # Markov (clock-forward) model using flexible parametric models within a relative survival framework (microsimulation package)
source("04a3_Markov_williams_ac_microsim.R")      # Markov (clock-forward) model using standard parametric models within an all-cause survival framework (microsimulation package)
source("04a5_semiMarkov_fpm_ac_microsim.R")       # semi-Markov (clock-reset) model using flexible parametric models within an all-cause survival framework (microsimulation package)
source("04a6_semiMarkov_fpm_rel_hesim.R")         # Proposed method: semi-Markov (clock-reset) model using flexible parametric models within a relative survival framework (hesim package)
source("04a6_semiMarkov_fpm_rel_microsim.R")      # Proposed method: semi-Markov (clock-reset) model using flexible parametric models within a relative survival framework (microsimulation package)
source("04a7_semiMarkov_williams_ac_microsim.R")  # semi-Markov (clock-reset) model using standard parametric models within an all-cause survival framework (microsimulation package)
source("04a7_semiMarkov_williams_ac_mstate.R")    # Williams et al. 2017: semi-Markov (clock-reset) model using standard parametric models within an all-cause survival framework (mstate package)
source("04b_table_survextrap.R")                  # Tabulate the observed and the extrapolated survival
source("04c_figure_survextrap.R")                 # Plot the observed and the extrapolated survival

##############################################################
##============================================================
## Full health economics model using microsimulation
## flexible parametric models within a relative survival framework
##============================================================
##############################################################
source("05_fullmod_semiMarkov_fpm_rel_microsim.R")        # Full model (including costs) of the proposed method (04a6_semiMarkov_fpm_rel_microsim.R)
source("06_fullmod_onewaysen.R")                          # One-way sensitivity analysis using the full model
source("07_fullmod_PSA.R")                                # Probabilistic sensitivity analysis using the full model
source("08_CEplane_CEAC.R")                               # Plot cost-effectiveness plane and cost-effectiveness acceptability curve

################################################################
# Copyright 2024 Chen EYT. All Rights Reserved.
# A Multistate Model Incorporating Relative Survival Extrapolation and 
# Mixed Time Scales for Health Technology Assessment
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
