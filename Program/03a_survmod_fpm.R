## Filename: 03a_survmod_fpm
## Purpose: Fit flexible parametric parametric models in 
##          1) an all-cause survival framework
##          2) a relative survival framework
##          using rstpm2 and flexsurv packages

## Set seed for coherency
set.seed(12345)

##############################################################
##============================================================
## Read packages
##============================================================
##############################################################
library(rstpm2)
library(flexsurv)
library(ggplot2)
library(ggpubr)

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
plotdataRFC_empty = expand.grid(treat = 1, Tstop=seq(0, 30, length.out = 180))
plotdataFC_empty = expand.grid(treat = 0, Tstop=seq(0, 30, length.out = 180))

## Generate a variable for merging with popmort file later
plotdataRFC_empty$t_floor <- floor(plotdataRFC_empty$Tstop)
plotdataFC_empty$t_floor <- floor(plotdataFC_empty$Tstop)

##############################################################
##============================================================
## Transition 1 PFS to prog
##============================================================
##############################################################
## Make empty datasets for predictions later 
plotdataRFC_m1_fpm <- plotdataRFC_empty
plotdataFC_m1_fpm <- plotdataFC_empty

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
m1_fpm <- models[[1]][[1]]
AIC(m1_fpm)
summary(m1_fpm)

## The model above is equivalent to this one using flexsurv::flexsurvspline below
m1_flexfpm_ac <- flexsurvspline(Surv(Tstart, Tstop, status == 1) ~ treat, 
                                data=msmcancer1, scale="hazard", k=1)
m1_flexfpm_ac

## Predict hazards
plotdataRFC_m1_fpm$haz_m1_fpm <- predict(m1_fpm, plotdataRFC_m1_fpm, type = "haz")
plotdataFC_m1_fpm$haz_m1_fpm <- predict(m1_fpm, plotdataFC_m1_fpm, type = "haz")

## Plot the hazards
plot(plotdataRFC_m1_fpm$Tstop, plotdataRFC_m1_fpm$haz, col = "blue", type = "l", lty="solid",
     xlim=c(0,30), ylim=c(0,1), xlab = "Time since study (years)", ylab = "Hazard rate (events/person-year)")
lines(plotdataFC_m1_fpm$Tstop, plotdataFC_m1_fpm$haz, col = "red", type = "l", lty="solid")
legend(0.2, 1, c("FPM (RFC)", "FPM (FC)"), bty="n", 
       lty=c("solid", "solid"), 
       col = c("blue", "red"), cex=1)
title(main="Progression-free -> progression", cex.main=1)

##############################################################
##============================================================
## Transition 2 PFS to death
##============================================================
##############################################################
## Check the events
with(msmcancer2, table(status==1))
## Maximum follow-up
max(msmcancer2$time)

##############################################
##=======All-cause survival framework=======##
##############################################
## Make a dataset for plotting later
plotdataRFC_m2_fpm_ac <- plotdataRFC_empty
plotdataFC_m2_fpm_ac <- plotdataFC_empty

## Build models
models <- lapply(c(2, 3, 4), function(df) {
  mod <- stpm2(Surv(Tstart, Tstop, status == 1) ~ treat, 
               data = msmcancer2, df = df)
  mod_tvc <- lapply(c(2, 3), function(tvc) {
    stpm2(Surv(Tstart, Tstop, status == 1) ~ treat, 
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
m2_fpm_ac <- models[[1]][[1]]
AIC(m2_fpm_ac)
summary(m2_fpm_ac)

## The model above is equivalent to this one using flexsurv::flexsurvspline below
m2_flexfpm_ac <- flexsurvspline(Surv(Tstart, Tstop, status == 1) ~ treat, data=msmcancer2, scale="hazard", k=1)
m2_flexfpm_ac

## Predict cumulative hazards and hazards
plotdataRFC_m2_fpm_ac$haz <- predict(m2_fpm_ac, plotdataRFC_m2_fpm_ac, type = "haz")
plotdataFC_m2_fpm_ac$haz <- predict(m2_fpm_ac, plotdataFC_m2_fpm_ac, type = "haz")

##############################################
##=======Relative survival framework========##
##############################################
## Make a dataset for plotting later
plotdataRFC_m2_fpm_rel <- plotdataRFC_empty
plotdataFC_m2_fpm_rel <- plotdataFC_empty

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
m2_fpm_rel <- models[[1]][[1]]
AIC(m2_fpm_rel)
summary(m2_fpm_rel)

## The model above is equivalent to this one using flexsurv::flexsurvspline below
m2_flexfpm_rel <- flexsurvspline(Surv(Tstart, Tstop, status == 1) ~ treat, 
                                bhazard = rate, data=msmcancer2, scale="hazard", k=1)
m2_flexfpm_rel

## Predict excess hazards
plotdataRFC_m2_fpm_rel$exhaz <- predict(m2_fpm_rel, plotdataRFC_m2_fpm_rel, type = "haz")
plotdataFC_m2_fpm_rel$exhaz <- predict(m2_fpm_rel, plotdataFC_m2_fpm_rel, type = "haz")

## To obtain all-cause hazard, it is required to estimate the expected hazard.
## Use the interrelationship (all-cause hazard = excess hazard + expected hazard)
## Expected hazard function of age 61 year 2005
popmort2 <- transform(popmort, cohort=att_year-att_age) %>%
            filter(cohort == 2005 - 61) %>%
            mutate(t_floor = att_age - 61) %>%
            filter(t_floor >= 0)

## rate here is the expected hazard
plot(popmort2$t_floor, popmort2$rate, col = "blue", type = "l", lty="longdash",
     xlim=c(0,100), ylim=c(0,0.30))

## Merge the expected hazard with the original data
plotdataRFC_m2_fpm_rel <- merge(plotdataRFC_m2_fpm_rel, popmort2, by=c("t_floor"))
plotdataRFC_m2_fpm_rel <- plotdataRFC_m2_fpm_rel %>%
                          mutate(achaz = exhaz + rate)

plotdataFC_m2_fpm_rel <- merge(plotdataFC_m2_fpm_rel, popmort2, by=c("t_floor"))
plotdataFC_m2_fpm_rel <- plotdataFC_m2_fpm_rel %>%
                         mutate(achaz = exhaz + rate)

## Plot all the all-cause hazards predicted in either an all-cause survival framework or
## a relative survival framework 
plot(plotdataRFC_m2_fpm_ac$Tstop, plotdataRFC_m2_fpm_ac$haz, col = "blue", type = "l", lty="solid",
     xlim=c(0,30), ylim=c(0,0.1), xlab = "Time since study (years)", ylab = "Hazard rate (events/person-year)")
lines(plotdataFC_m2_fpm_ac$Tstop, plotdataFC_m2_fpm_ac$haz, col = "red", type = "l", lty="solid")
lines(plotdataRFC_m2_fpm_rel$Tstop, plotdataRFC_m2_fpm_rel$achaz, col = "blue", type = "l", lty="longdash")  
lines(plotdataFC_m2_fpm_rel$Tstop, plotdataFC_m2_fpm_rel$achaz, col = "red", type = "l", lty="longdash")  
legend(0.2, 0.1, c("FPM all-cause (RFC)", "FPM all-cause (FC)", "FPM relative (RFC)", "FPM relative (FC)"), bty="n", 
       lty=c("solid", "solid", "longdash","longdash"), 
       col = c("blue", "red", "blue", "red"), cex=1)
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

##############################################
##=======All-cause survival framework=======##
##############################################
## Make a dataset for plotting later
## Caution: semi-Markov--time variable is "time"
plotdataRFC_m3_fpm_ac = expand.grid(treat = 1, time=seq(0, 30, length.out = 180))
plotdataRFC_m3_fpm_ac$Tstop <- plotdataRFC_m3_fpm_ac$time

plotdataFC_m3_fpm_ac = expand.grid(treat = 0, time=seq(0, 30, length.out = 180))
plotdataFC_m3_fpm_ac$Tstop <- plotdataFC_m3_fpm_ac$time

## Build models
models <- lapply(c(2, 3, 4), function(df) {
  mod <- stpm2(Surv(time, status == 1) ~ treat, data = msmcancer3,
               df = df)
  mod_tvc <- lapply(c(2, 3), function(tvc) {
    stpm2(Surv(time, status == 1) ~ treat, data = msmcancer3,
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
## Semi-Markov
m3_semiMarkov_fpm_ac <- models[[2]][[1]]
AIC(m3_semiMarkov_fpm_ac)
summary(m3_semiMarkov_fpm_ac)

## The model above is equivalent to this one using flexsurv::flexsurvspline below
## Semi-Markov
m3_semiMarkov_flexfpm_ac <-flexsurvspline(Surv(time, status == 1) ~ treat, 
                                          data=msmcancer3, scale="hazard", k=2)
m3_semiMarkov_flexfpm_ac

## Predict hazards
plotdataRFC_m3_fpm_ac$haz_semiMarkov <- predict(m3_semiMarkov_fpm_ac, plotdataRFC_m3_fpm_ac, type = "haz")
plotdataRFC_m3_fpm_ac <- plotdataRFC_m3_fpm_ac[plotdataRFC_m3_fpm_ac$time>0,]

plotdataFC_m3_fpm_ac$haz_semiMarkov <- predict(m3_semiMarkov_fpm_ac, plotdataFC_m3_fpm_ac, type = "haz")
plotdataFC_m3_fpm_ac <- plotdataFC_m3_fpm_ac[plotdataFC_m3_fpm_ac$time>0,]

##############################################
##=======Relative survival framework========##
##############################################
## Make a dataset for plotting later
## Caution: semi-Markov--time variable is "time"
plotdataRFC_m3_fpm_rel = expand.grid(treat = 1, time=seq(0, 30, length.out = 180))
plotdataRFC_m3_fpm_rel$Tstop <- plotdataRFC_m3_fpm_rel$time
plotdataRFC_m3_fpm_rel$t_floor <- floor(plotdataRFC_m3_fpm_rel$time)

plotdataFC_m3_fpm_rel = expand.grid(treat = 0, time=seq(0, 30, length.out = 180))
plotdataFC_m3_fpm_rel$Tstop <- plotdataFC_m3_fpm_rel$time
plotdataFC_m3_fpm_rel$t_floor <- floor(plotdataFC_m3_fpm_rel$time)

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
## Semi-Markov
m3_semiMarkov_fpm_rel <- models[[2]][[1]]
AIC(m3_semiMarkov_fpm_rel)
summary(m3_semiMarkov_fpm_rel)

## The model above is equivalent to this one using flexsurv::flexsurvspline below
## Semi-Markov
m3_semiMarkov_flexfpm_rel <-flexsurvspline(Surv(time, status == 1) ~ treat, bhazard = rate,
                                           data=msmcancer3, scale="hazard", k=2)
m3_semiMarkov_flexfpm_rel

## Predict excess hazards
plotdataRFC_m3_fpm_rel$exhaz_semiMarkov <- predict(m3_semiMarkov_fpm_rel, plotdataRFC_m3_fpm_rel, type = "haz")
plotdataFC_m3_fpm_rel$exhaz_semiMarkov <- predict(m3_semiMarkov_fpm_rel, plotdataFC_m3_fpm_rel, type = "haz")

## To obtain all-cause hazard, it is required to estimate the expected hazard.
## Use the interrelationship (all-cause hazard = excess hazard + expected hazard)
## Expected hazard function of age 61 year 2005
## Merge the expected hazard with the original data
plotdataRFC_m3_fpm_rel <- merge(plotdataRFC_m3_fpm_rel, popmort2, by=c("t_floor"))
plotdataRFC_m3_fpm_rel <- plotdataRFC_m3_fpm_rel %>%
                          mutate(achaz_semiMarkov = exhaz_semiMarkov + rate) 
plotdataRFC_m3_fpm_rel <- plotdataRFC_m3_fpm_rel[plotdataRFC_m3_fpm_rel$time>0,]

plotdataFC_m3_fpm_rel <- merge(plotdataFC_m3_fpm_rel, popmort2, by=c("t_floor"))
plotdataFC_m3_fpm_rel <- plotdataFC_m3_fpm_rel %>%
                         mutate(achaz_semiMarkov = exhaz_semiMarkov + rate)
plotdataFC_m3_fpm_rel <- plotdataFC_m3_fpm_rel[plotdataFC_m3_fpm_rel$time>0,]

## Plot all the all-cause hazards predicted in either an all-cause survival framework or
## a relative survival framework 
plot(plotdataRFC_m3_fpm_ac$time, plotdataRFC_m3_fpm_ac$haz_semiMarkov, col = "blue", type = "l", lty="solid",
     xlim=c(0,30), ylim=c(0,1), xlab = "Time since progression (years)", ylab = "Hazard rate (events/person-year)")
lines(plotdataFC_m3_fpm_ac$time, plotdataFC_m3_fpm_ac$haz_semiMarkov, col = "red", type = "l", lty="solid")
lines(plotdataRFC_m3_fpm_rel$time, plotdataRFC_m3_fpm_rel$achaz_semiMarkov, col = "blue", type = "l", lty="longdash")  
lines(plotdataFC_m3_fpm_rel$time, plotdataFC_m3_fpm_rel$achaz_semiMarkov, col = "red", type = "l", lty="longdash")  
legend(0.2, 1, c("FPM all-cause (RFC)", "FPM all-cause (FC)", "FPM relative (RFC)", "FPM relative (FC)"), bty="n", 
       lty=c("solid", "solid", "longdash","longdash"), 
       col = c("blue", "red", "blue", "red"), cex=1)
title(main="Progression -> death", cex.main=1)

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