## Filename: 03_rel_survmod
## Purpose: Fit parametric models for multistate-structed data
##          Relative survival framework

## Set the wd as where this R file is
if (require("rstudioapi") && isAvailable()) {
  original_wd <- getwd()  # Store the original working directory
  wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(wd)
}

## Set seed for coherency
set.seed(12345)

##############################################################
##============================================================
## Read packages
##============================================================
##############################################################
library(rstpm2)
library(ggplot2)
library(ggpubr)

##############################################################
##============================================================
## Read the survival data
##============================================================
##############################################################
source("02_msmcancer.R")

##############################################################
##============================================================
## Transition 1 PFS to prog
##============================================================
##############################################################
## Make a dataset for plotting later
plotdataRFC_empty = expand.grid(treat = 1, Tstop=seq(0, 15, length.out = 180))
plotdataRFC <- plotdataRFC_empty
plotdataFC_empty = expand.grid(treat = 0, Tstop=seq(0, 15, length.out = 180))
plotdataFC <- plotdataFC_empty

## Check the events
with(msmcancer1, table(status==1))
## Maximum follow-up
max(msmcancer1$prog_ty)

## Build models
models <- lapply(c(2, 3, 4), function(df) {
  mod <- stpm2(Surv(Tstart, Tstop, status == 1) ~ treat, data = msmcancer1, df = df)
  mod_tvc <- lapply(c(2, 3), function(tvc) {
    stpm2(Surv(Tstart, Tstop, status == 1) ~ treat, data = msmcancer1, tvc = list(treat = tvc), df = df)
  })
  c(mod, mod_tvc)
})

aic_values <- unlist(lapply(models, function(model) {
  unlist(lapply(model, AIC))
}))

## List the AIC of all fitted mods
## The lowest is the first one
print(aic_values)

## Select model for transition 1 based on the lowest AIC
m1 <- models[[1]][[1]]
AIC(m1)
summary(m1)

## Predict 15-year horizon based on the fitted model
plotdataRFC$cumhaz <- predict(m1, plotdataRFC, type = "cumhaz")
plotdataFC$cumhaz <- predict(m1, plotdataFC, type = "cumhaz")

## Plot
## 4 years
## Cumulative hazard functions
plot(survfit(Surv(time, status) ~ 1, data = msmcancer1RFC)$time, 
     -log(survfit(Surv(time, status) ~ 1, data = msmcancer1RFC)$surv),
     type = "l", lty = "solid",  
     xlab = "Time", ylab = "Cumulative hazard",
     xaxs = "i", yaxs = "i",
     las = 1, cex.lab = 1, cex.axis = 1, col = "navy", 
     xlim=c(0, 4), ylim=c(0,2))
lines(survfit(Surv(time,status)~ 1, data=msmcancer1FC)$time, 
      -log(survfit(Surv(time,status)~ 1, data=msmcancer1FC)$surv), col="firebrick4", lty="solid")
lines(plotdataRFC$Tstop, plotdataRFC$cumhaz, col = "blue", lty="longdash")
lines(plotdataFC$Tstop, plotdataFC$cumhaz, col = "red", lty="longdash")  
legend(0.5, 6, c("Observed--RFC", "Observed--FC", "FPM--RFC", "FPM--FC"), bty="n", 
       lty=c("solid","solid","longdash","longdash"), 
       col = c("navy", "firebrick4", "blue", "red"),cex=1)
title(main="Progression-free -> progression", cex.main=1)

## 15 years
## Cumulative hazard functions
plot(survfit(Surv(time, status) ~ 1, data = msmcancer1RFC)$time, 
     -log(survfit(Surv(time, status) ~ 1, data = msmcancer1RFC)$surv),
     type = "l", lty = "solid",  
     xlab = "Time", ylab = "Cumulative hazard",
     xaxs = "i", yaxs = "i",
     las = 1, cex.lab = 1, cex.axis = 1, col = "navy", 
     xlim=c(0, 15), ylim=c(0,10))
lines(survfit(Surv(time,status)~ 1, data=msmcancer1FC)$time, 
      -log(survfit(Surv(time,status)~ 1, data=msmcancer1FC)$surv), col="firebrick4", lty="solid")
lines(plotdataRFC$Tstop, plotdataRFC$cumhaz, col = "blue", lty="longdash")
lines(plotdataFC$Tstop, plotdataFC$cumhaz, col = "red", lty="longdash")  
legend(0.5, 10, c("Observed--RFC", "Observed--FC", "FPM--RFC", "FPM--FC"), bty="n", 
       lty=c("solid","solid","longdash","longdash"), 
       col = c("navy", "firebrick4", "blue", "red"),cex=1)
title(main="Progression-free -> progression", cex.main=1)

##############################################################
##============================================================
## Transition 2 PFS to death
##============================================================
##############################################################
## Make a dataset for plotting later
plotdataRFC = plotdataRFC_empty
plotdataFC = plotdataFC_empty

## Check the events
with(msmcancer2, table(status==1))
## Maximum follow-up
max(msmcancer2$time)

## Build models
models <- lapply(c(2, 3, 4), function(df) {
  mod <- stpm2(Surv(Tstart, Tstop, status == 1) ~ treat + bhazard(rate), 
               data = msmcancer2, df = df)
  mod_tvc <- lapply(c(2, 3), function(tvc) {
    stpm2(Surv(Tstart, Tstop, status == 1) ~ treat  + bhazard(rate), 
          data = msmcancer2, tvc = list(treat = tvc), df = df)
  })
  c(mod, mod_tvc)
})

aic_values <- unlist(lapply(models, function(model) {
  unlist(lapply(model, AIC))
}))

## List the AIC of all fitted mods
## The lowest is the first one
print(aic_values)

## Select model for transition 1 based on the lowest AIC
m2 <- models[[1]][[1]]
AIC(m2)
summary(m2)

## Predict 15-year horizon based on the fitted model
## Excess hazard functions
plotdataRFC$haz <- predict(m2, plotdataRFC, type = "haz")
plotdataFC$haz <- predict(m2, plotdataFC, type = "haz")
## Cumulative excess hazard functions
plotdataRFC$cumhaz <- predict(m2, plotdataRFC, type = "cumhaz")
plotdataFC$cumhaz <- predict(m2, plotdataFC, type = "cumhaz")

## Plot
## 4 years
## Excess hazard functions
plot(plotdataRFC$Tstop, plotdataRFC$haz, col = "blue", type = "l", lty="longdash",
     xlim=c(0,4), ylim=c(0,0.2))
lines(plotdataFC$Tstop, plotdataFC$haz, col = "red", type = "l", lty="longdash")  
legend(0.2, 0.2, c("FPRSM--RFC", "FPRSM--FC"), bty="n", 
       lty=c("longdash","longdash"), 
       col = c("blue", "red"), cex=1)
title(main="Progression-free -> death", cex.main=1)

## Cumulative excess hazard functions
plot(plotdataRFC$Tstop, plotdataRFC$cumhaz, col = "blue", type = "l", lty="longdash",
     xlim=c(0,4), ylim=c(0,0.2))
lines(plotdataFC$Tstop, plotdataFC$cumhaz, col = "red", type = "l", lty="longdash")  
legend(0.2, 0.2, c("FPRSM--RFC", "FPRSM--FC"), bty="n", 
       lty=c("longdash","longdash"), 
       col = c("blue", "red"), cex=1)
title(main="Progression-free -> death", cex.main=1)

## 15 years
## Excess hazard functions
plot(plotdataRFC$Tstop, plotdataRFC$haz, col = "blue", type = "l", lty="longdash",
     xlim=c(0,15), ylim=c(0,0.2))
lines(plotdataFC$Tstop, plotdataFC$haz, col = "red", type = "l", lty="longdash")  
legend(0.2, 0.2, c("FPRSM--RFC", "FPRSM--FC"), bty="n", 
       lty=c("longdash","longdash"), 
       col = c("blue", "red"), cex=1)
title(main="Progression-free -> death", cex.main=1)

## Cumulative excess hazard functions
plot(plotdataRFC$Tstop, plotdataRFC$cumhaz, col = "blue", type = "l", lty="longdash",
     xlim=c(0,15), ylim=c(0,0.5))
lines(plotdataFC$Tstop, plotdataFC$cumhaz, col = "red", type = "l", lty="longdash")  
legend(0.2, 0.5, c("FPRSM--RFC", "FPRSM--FC"), bty="n", 
       lty=c("longdash","longdash"), 
       col = c("blue", "red"), cex=1)
title(main="Progression-free -> death", cex.main=1)

##############################################################
##============================================================
## Transition 3 Progression to death
##============================================================
##############################################################
## Check the events
with(msmcancer3, table(status==1))
## Maximum follow-up
max(msmcancer3$time)
## Make a dataset for plotting later
## Caution: the time variable name here is "time"
##          prepare Markov-renewal
plotdataRFC = expand.grid(treat = 1, time=seq(0, 15, length.out = 180))
plotdataFC = expand.grid(treat = 0, time=seq(0, 15, length.out = 180))

## Build models
models <- lapply(c(2, 3, 4), function(df) {
  mod <- stpm2(Surv(time, status == 1) ~ treat + bhazard(rate), data = msmcancer3,
               df = df)
  mod_tvc <- lapply(c(2, 3), function(tvc) {
    stpm2(Surv(time, status == 1) ~ treat + bhazard(rate), data = msmcancer3,
          tvc = list(treat = tvc), df = df)
  })
  c(mod, mod_tvc)
})

aic_values <- unlist(lapply(models, function(model) {
  unlist(lapply(model, AIC))
}))

## List the AIC of all fitted mods
## The lowest is the first one
print(aic_values)

## Select model for transition 1 based on the lowest AIC
m3 <- models[[2]][[1]]
AIC(m3)
summary(m3)

## Predict 15-year horizon based on the fitted model
## Excess hazard functions
plotdataRFC$haz <- predict(m3, plotdataRFC, type = "haz")
plotdataFC$haz <- predict(m3, plotdataFC, type = "haz")
## Cumulative excess hazard functions
plotdataRFC$cumhaz <- predict(m3, plotdataRFC, type = "cumhaz")
plotdataFC$cumhaz <- predict(m3, plotdataFC, type = "cumhaz")

## Plot
## 4 years
## Excess hazard functions
plot(plotdataRFC$time, plotdataRFC$haz, col = "blue", type = "l", lty="longdash",
     xlim=c(0,4), ylim=c(0,1))
lines(plotdataFC$time, plotdataFC$haz, col = "red", type = "l", lty="longdash")  
legend(0.2, 1, c("FPRSM--RFC", "FPRSM--FC"), bty="n", 
       lty=c("longdash","longdash"), 
       col = c("blue", "red"), cex=1)
title(main="Progression -> death", cex.main=1)

## Cumulative excess hazard functions
plot(plotdataRFC$time, plotdataRFC$cumhaz, col = "blue", type = "l", lty="longdash",
     xlim=c(0,4), ylim=c(0,2))
lines(plotdataFC$time, plotdataFC$cumhaz, col = "red", type = "l", lty="longdash")  
legend(0.2, 2, c("FPRSM--RFC", "FPRSM--FC"), bty="n", 
       lty=c("longdash","longdash"), 
       col = c("blue", "red"), cex=1)
title(main="Progression -> death", cex.main=1)

## 15 years
## Excess hazard functions
plot(plotdataRFC$time, plotdataRFC$haz, col = "blue", type = "l", lty="longdash",
     xlim=c(0,15), ylim=c(0,2))
lines(plotdataFC$time, plotdataFC$haz, col = "red", type = "l", lty="longdash")  
legend(0.2, 2, c("FPRSM--RFC", "FPRSM--FC"), bty="n", 
       lty=c("longdash","longdash"), 
       col = c("blue", "red"), cex=1)
title(main="Progression -> death", cex.main=1)

## Cumulative excess hazard functions
plot(plotdataRFC$time, plotdataRFC$cumhaz, col = "blue", type = "l", lty="longdash",
     xlim=c(0,15), ylim=c(0,5))
lines(plotdataFC$time, plotdataFC$cumhaz, col = "red", type = "l", lty="longdash")  
legend(0.2, 5, c("FPRSM--RFC", "FPRSM--FC"), bty="n", 
       lty=c("longdash","longdash"), 
       col = c("blue", "red"), cex=1)
title(main="Progression -> death", cex.main=1)

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