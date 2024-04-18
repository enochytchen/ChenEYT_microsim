## Filename: 04a_extrap_semiMarkov_fpm
## Purpose: Run the semi-Markov model using flexible parametric models (flexsurvreg)
## Notes: Must have installed the lastest Rtools and R4.2.2+

##############################################################
##============================================================
## Read the packages
##============================================================
##############################################################
library(mstate)

##############################################################
##============================================================
## Read the models
##============================================================
##############################################################
## This file also reads in popmort file from 01b_popmort.R
source("03b_survmod_fpm.R")

##############################################################
##============================================================
## Construct a multistate model from a set of parametric survival models
##============================================================
##############################################################
## Define transition matrix
tmat <- rbind("PFS"=c(NA, 1, 2), "Prog"=c(NA, NA, 3), "Death"=c(NA, NA, NA))

## Put all the fitted models into one object
crfs <- fmsm(m1_flexfpm_ac, m2_flexfpm_ac, m3_flexfpm_ac, trans = tmat)

## Define a multistate model
ms_RFC<- msfit.flexsurvreg(crfs, t = seq(0, 15, by = 0.1), trans = tmat, newdata=data.frame(treat=1))
ms_FC<- msfit.flexsurvreg(crfs, t = seq(0, 15, by = 0.1), trans = tmat, newdata=data.frame(treat=0))

## Estimate probability at times
fpm1fpm2fpm3sim_RFC<- mssample(ms_RFC$Haz, trans = tmat, clock = "reset", M = 5000,
                               tvec =  seq(0, 15, by = 0.1))
fpm1fpm2fpm3sim_FC<- mssample(ms_RFC$Haz, trans = tmat, clock = "reset", M = 5000,
                              tvec =  seq(0, 15, by = 0.1))

##############################################################
##============================================================
## Organize the results for comparison
##============================================================
##############################################################
## Create a variable for undiscounted survival
## that is, probability of being at PFS (pstate1) + Progression (pstate2)
fpm1fpm2fpm3sim_RFC$surv <- fpm1fpm2fpm3sim_RFC$pstate1 + fpm1fpm2fpm3sim_RFC$pstate2
fpm1fpm2fpm3sim_FC$surv <- fpm1fpm2fpm3sim_FC$pstate1 + fpm1fpm2fpm3sim_FC$pstate2

## Create a variable to indicate treatments
fpm1fpm2fpm3sim_RFC$treat <- 1
fpm1fpm2fpm3sim_FC$treat <- 0

## Combine data for saving 
semiMarkov_fpm <- bind_rows(fpm1fpm2fpm3sim_RFC, fpm1fpm2fpm3sim_FC)

## Plot survival curve to check
ggplot(semiMarkov_fpm, aes(x = time, y = surv, color=as.factor(treat))) +
  geom_line()

write.table(semiMarkov_fpm, file = "../Data/04a_extrap_semiMarkov_fpm.txt", sep = " ",
            row.names = FALSE, col.names = TRUE,  quote = FALSE)
##############################################################
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