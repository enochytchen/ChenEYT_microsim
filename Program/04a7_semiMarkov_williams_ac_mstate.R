## Filename: 04a7_semiMarkov_williams_ac_mstate
## Purpose: Run semi-Markov (clock-reset) model using mstate package with
##          standard parametric models within an all-cause survival framework
##          Models (m1: gompertz; m2: g-gamma; m3: gompertz) same as the base-case
##          analysis as Williams et al. 2017. 
##          This analysis is essentially the same one as the base case analysis
##          in Williams et al. 2017.
## Reference: Williams C, Lewsey JD, Briggs AH, Mackay DF. 
##            Cost-effectiveness Analysis in R Using a Multi-state Modeling Survival Analysis Framework: A Tutorial. 
##            Med Decis Making. 2017 May;37(4):340–52.
##            Cost-effectiveness Analysis in R Using a Multi-state Modeling Survival Analysis Framework: A Tutorial © 2017 by Williams C. et al is licensed under CC BY 3.0. 

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
crfs <- fmsm(m1_gom, m2_gam, m3_semiMarkov_gom, trans = tmat)

## Define a multistate model
t = c(seq(0,4,1/12), seq(4+1/288, 19, 1/288))

ms_RFC <- msfit.flexsurvreg(crfs, t = t, trans = tmat, newdata=data.frame(treat=1))
ms_FC  <- msfit.flexsurvreg(crfs, t = t, trans = tmat, newdata=data.frame(treat=0))

## Estimate probability at times
gom1gam2gom3sim_RFC<- mssample(ms_RFC$Haz, trans = tmat, clock = "reset", M = 1e4,
                               tvec =  t)
gom1gam2gom3sim_FC<- mssample(ms_FC$Haz, trans = tmat, clock = "reset", M = 1e4,
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
semiMarkov_williams <- bind_rows(gom1gam2gom3sim_RFC, gom1gam2gom3sim_FC)

## Plot survival curve to check
ggplot(semiMarkov_williams, aes(x = time, y = surv, color=as.factor(treat))) +
  geom_line()

## Save results
## Prevent from changing the results. We put # here.
# write.table(semiMarkov_williams, file = "../Data/04a7_semiMarkov_williams_ac_mstate.txt", sep = " ",
#             row.names = FALSE, col.names = TRUE,  quote = FALSE)
##############################################################
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