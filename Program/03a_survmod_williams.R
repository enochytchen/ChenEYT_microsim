## Filename: 03a_williams_survmod
## Purpose: Refit Williams et al.'s models for hazard comparisons
## Reference: Williams C, Lewsey JD, Briggs AH, Mackay DF. 
##            Cost-effectiveness Analysis in R Using a Multi-state Modeling Survival Analysis Framework: A Tutorial. 
##            Med Decis Making. 2017 May;37(4):340–52.
##            Cost-effectiveness Analysis in R Using a Multi-state Modeling Survival Analysis Framework: A Tutorial © 2017 by Williams C. et al is licensed under CC BY 3.0. 

##############################################################
##============================================================
## Read packages
##============================================================
##############################################################
library(ggplot2)
library(ggpubr)
library(flexsurv)

##############################################################
##============================================================
## Read the survival data
##============================================================
##############################################################
## This file also reads in popmort file from 01b_popmort.R
source("02_msmcancer.R")

##############################################################
##============================================================
## Make datasets for plotting later
##============================================================
##############################################################
plotdataRFC_empty = expand.grid(treat = 1, Tstop=seq(0, 15, length.out = 180))
plotdataFC_empty = expand.grid(treat = 0, Tstop=seq(0, 15, length.out = 180))

##############################################################
##============================================================
## Transition 1 PFS to prog: Gompertz (Williams et al.)
##============================================================
##############################################################
## Make empty datasets for predictions later
plotdataRFC_m1_w <- plotdataRFC_empty
plotdataFC_m1_w <- plotdataFC_empty

## Fit gompertz model
m1_gom <- flexsurvreg(Surv(Tstart, Tstop, status) ~ treat,
                      dist="gompertz",
                      data=msmcancer1)

## Predict hazards
## RFC
haz_m1_gom <- predict(m1_gom, newdata = plotdataRFC_m1_w, type = "haz", times = plotdataRFC_m1_w$Tstop)[[1]][[1]]$".pred_hazard"
plotdataRFC_m1_w <- cbind(plotdataRFC_m1_w, haz_m1_gom)

## FC
haz_m1_gom <- predict(m1_gom, newdata = plotdataFC_m1_w, type = "haz", times = plotdataFC_m1_w$Tstop)[[1]][[1]]$".pred_hazard"
plotdataFC_m1_w <- cbind(plotdataFC_m1_w, haz_m1_gom)

## Plot the hazards
plot(plotdataRFC_m1_w$Tstop, plotdataRFC_m1_w$haz_m1_gom, col = "blue", type = "l", lty="solid",
     xlim=c(0,15), ylim=c(0,150), xlab = "Time since study (years)", ylab = "Hazard rate (events/person-year)")
lines(plotdataFC_m1_w$Tstop, plotdataFC_m1_w$haz_m1_gom, col = "red", type = "l", lty="solid")
legend(0.2, 150, c("Gompertz (RFC)", "Gompertz (FC)"), bty="n", 
       lty=c("solid", "solid"), 
       col = c("blue", "red"), cex=1)
title(main="Progression-free -> progression", cex.main=1)

#################################################################
##===============================================================
## Transition 2 PFS to death: generalized gamma (Williams et al.)
##===============================================================
#################################################################
## Check the events
with(msmcancer2, table(status==1))
## Maximum follow-up
max(msmcancer2$time)

## Make empty datasets for predictions later
plotdataRFC_m2_w <- plotdataRFC_empty
plotdataFC_m2_w <- plotdataFC_empty

## Fit generalized gamma model
m2_gam <- flexsurvreg(Surv(Tstart, Tstop, status) ~ treat,
                      dist="gengamma",
                      data=msmcancer2)

## Predict hazards
## RFC
haz_m2_gam <- predict(m2_gam, newdata = plotdataRFC_m2_w, type = "haz", times = plotdataRFC_m2_w$Tstop)[[1]][[1]]$".pred_hazard"
plotdataRFC_m2_w <- cbind(plotdataRFC_m2_w, haz_m2_gam)

## FC
haz_m2_gam <- predict(m2_gam, newdata = plotdataFC_m2_w, type = "haz", times = plotdataFC_m2_w$Tstop)[[1]][[1]]$".pred_hazard"
plotdataFC_m2_w <- cbind(plotdataFC_m2_w, haz_m2_gam)

## Plot the hazards
plot(plotdataRFC_m2_w$Tstop, plotdataRFC_m2_w$haz_m2_gam, col = "blue", type = "l", lty="solid",
     xlim=c(0,15), ylim=c(0,0.1), xlab = "Time since study (years)", ylab = "Hazard rate (events/person-year)")
lines(plotdataFC_m2_w$Tstop, plotdataFC_m2_w$haz_m2_gam, col = "red", type = "l", lty="solid")
legend(0.1, 0.1, c("Ggamma (RFC)", "Ggamma (FC)"), bty="n", 
       lty=c("solid", "solid"), 
       col = c("blue", "red"), cex=1)
title(main="Progression-free -> death", cex.main=1)

################################################################
##==============================================================
## Transition 3 Progression to death: Gompertz (Williams et al.)
##==============================================================
################################################################
## Make a dataset for plotting later
## Caution: the time variable name here is "time"
##          prepare Markov-renewal
plotdataRFC_empty = expand.grid(treat = 1, time=seq(0, 15, length.out = 180))
plotdataFC_empty = expand.grid(treat = 0, time=seq(0, 15, length.out = 180))

## Make empty datasets for predictions later
plotdataRFC_m3_w <- plotdataRFC_empty
plotdataFC_m3_w <- plotdataFC_empty

## Fit Gompertz model
m3_gom <- flexsurvreg(Surv(time, status)~ treat,
                         dist="gompertz",
                         data=msmcancer3)

## Predict hazards
## RFC
haz_m3_gom <- predict(m3_gom, newdata = plotdataRFC_m3_w, type = "haz", times = plotdataRFC_m3_w$time)[[1]][[1]]$".pred_hazard"
plotdataRFC_m3_w <- cbind(plotdataRFC_m3_w, haz_m3_gom)

## FC
haz_m3_gom <- predict(m3_gom, newdata = plotdataFC_m3_w, type = "haz", times = plotdataFC_m3_w$time)[[1]][[1]]$".pred_hazard"
plotdataFC_m3_w <- cbind(plotdataFC_m3_w, haz_m3_gom)

## Plot the hazards
plot(plotdataRFC_m3_w$time, plotdataRFC_m3_w$haz_m3_gom, col = "blue", type = "l", lty="solid",
     xlim=c(0,15), ylim=c(0.0,0.3), xlab = "Time since study (years)", ylab = "Hazard rate (events/person-year)")
lines(plotdataFC_m3_w$time, plotdataFC_m3_w$haz_m3_gom, col = "red", type = "l", lty="solid")
legend(5, 0.05, c("Gompertz (RFC)", "Gompertz (FC)"), bty="n", 
       lty=c("solid", "solid"), 
       col = c("blue", "red"), cex=1)
title(main="Progression -> death", cex.main=1)

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
