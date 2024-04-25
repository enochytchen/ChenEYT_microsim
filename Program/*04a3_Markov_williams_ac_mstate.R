## Filename: 04a4_Markov_williams_ac_mstate
## Purpose: Run the Markov model using standard parametric models (flexsurvreg)
## Notes: Must have installed the lastest Rtools and R4.2.2+

##############################################################
##============================================================
## Read the packages
##============================================================
##############################################################
library(mstate)
library(flexsurv)
library(tidyverse)

##############################################################
##============================================================
## Read the models
##============================================================
##############################################################
## This file also reads in popmort file from 01b_popmort.R
source("03a_survmod_williams.R")

##############################################################
##============================================================
## Construct a multistate model from a set of parametric survival models
##============================================================
##############################################################
## Define transition matrix
tmat <- rbind("PFS"=c(NA, 1, 2), "Prog"=c(NA, NA, 3), "Death"=c(NA, NA, NA))

## Put all the fitted models into one object
crfs <- fmsm(m1_gom, m2_gam, m3_Markov_gom, trans = tmat)

## Define a multistate model
t = c(seq(0, 4, 1/12), seq(4+1/144, 15, 1/144), seq(15+1/144, 20, 1/144))

ms_RFC <- msfit.flexsurvreg(crfs, t = t, trans = tmat, newdata=data.frame(treat=1))
ms_FC  <- msfit.flexsurvreg(crfs, t = t, trans = tmat, newdata=data.frame(treat=0))

## Estimate probability at times
gom1gam2gom3sim_RFC<- mssample(ms_RFC$Haz, trans = tmat, clock = "forward", M = 1000,
                               tvec =  t)
gom1gam2gom3sim_FC<- mssample(ms_FC$Haz, trans = tmat, clock = "forward", M = 1000,
                              tvec =  t)

##############################################################
##============================================================
## Organize the results for comparison
##============================================================
##############################################################
## Create a variable for undiscounted survival
## that is, probability of being at PFS (pstate1) + Progression (pstate2)
gom1gam2gom3sim_RFC$surv <- gom1gam2gom3sim_RFC$pstate1 + gom1gam2gom3sim_RFC$pstate2
gom1gam2gom3sim_FC$surv <- gom1gam2gom3sim_FC$pstate1 + gom1gam2gom3sim_FC$pstate2

## Create a variable to indicate treatments
gom1gam2gom3sim_RFC$treat <- 1
gom1gam2gom3sim_FC$treat <- 0

## Combine data for saving 
Markov_williams <- bind_rows(gom1gam2gom3sim_RFC, gom1gam2gom3sim_FC)

## Plot survival curve to check
ggplot(Markov_williams, aes(x = time, y = surv, color=as.factor(treat))) +
  geom_line()

write.table(Markov_williams, file = "../Data/04a4_Markov_williams_ac_mstate.txt", sep = " ",
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