## Filename: 04b_table_survextrap
## Purpose: Caluclate area under curve for selected time points
##          Observed: Hallek 2010, Fischer 2016
##          Extrapolated: Semi-Markov with stand parametric models within an all-cause survival framework (Williams 2017)  
##                        Semi-Markov with flexible parametric models within an all-cause survival framework
##                        Semi-Markov with flexible parametric models within a relative survival framework 

## Reference: Williams C, Lewsey JD, Briggs AH, Mackay DF. 
##            Cost-effectiveness Analysis in R Using a Multi-state Modeling Survival Analysis Framework: A Tutorial. 
##            Med Decis Making. 2017 May;37(4):340–52.
##            Cost-effectiveness Analysis in R Using a Multi-state Modeling Survival Analysis Framework: A Tutorial © 2017 by Williams C. et al is licensed under CC BY 3.0. 

##############################################################
##============================================================
## Install/Read packages
##============================================================
##############################################################
library(tidyverse)
library(ggplot2)
library(survival)
library(survRM2)
library(pracma)
library(flextable)

##############################################################
##============================================================
## Read data
##============================================================
###############################################################
#################
## Hallek 2010 ##
#################
hallek2010 <- read.table("../Data/data.txt", header = TRUE)
colnames(hallek2010)
hallek2010$progdeath_ty <- hallek2010$progdeath_t/12
hallek2010$death_ty <- hallek2010$death_t/12

mfit1 <- survfit(Surv(death_ty, death==1) ~ treat, data = hallek2010)
summary(mfit1)

mfit1 <- survfit(Surv(progdeath_ty, progdeath==1) ~ treat, data = hallek2010)
summary(mfit1)

##################
## Fischer 2016 ##
##################
fischer2016 <- read.table("../Data/01a_fischer2016.txt", header = TRUE)
fischer2016$time_y <- fischer2016$time/12
mfit2 <- survfit(Surv(time_y, event==1) ~ treat, data = fischer2016)

fischer_FC_4yr_surv <- summary(mfit2, times=4)$surv[1]
fischer_RFC_4yr_surv <- summary(mfit2, times=4)$surv[2]

fischer_FC_8yr_surv <- summary(mfit2, times=8)$surv[1]
fischer_RFC_8yr_surv <- summary(mfit2, times=8)$surv[2]

######################
## Fischer 2016 PFS ##
######################
fischer2016_pfs <- read.table("../Data/01a_fischer2016_pfs.txt", header = TRUE)
fischer2016_pfs$time_y <- fischer2016_pfs$time/12

# PFS
mfit2_pfs <- survfit(Surv(time_y, event==1) ~ treat, data = fischer2016_pfs)

fischer_FC_4yr_pfs <- summary(mfit2_pfs, times=4)$surv[1]
fischer_RFC_4yr_pfs <- summary(mfit2_pfs, times=4)$surv[2]

fischer_FC_8yr_pfs <- summary(mfit2_pfs, times=8)$surv[1]
fischer_RFC_8yr_pfs <- summary(mfit2_pfs, times=8)$surv[2]

#######################
## Fischer 2016 Prog ##
#######################
fischer_FC_4yr_prog <-  fischer_FC_4yr_surv - fischer_FC_4yr_pfs
fischer_RFC_4yr_prog <-  fischer_RFC_4yr_surv - fischer_RFC_4yr_pfs

fischer_FC_8yr_prog <-  fischer_FC_8yr_surv - fischer_FC_8yr_pfs
fischer_RFC_8yr_prog <-  fischer_RFC_8yr_surv - fischer_RFC_8yr_pfs

#######################
#######################
##    Semi-Markov    ##
#######################
#######################

###############################################################################################
## Semi-Markov with standard parametric models, all-cause survival framework (Williams 2017) ##
###############################################################################################
## Use the data for extrapolation until lifetime instead
williams2017 <- read.table("../Data/04a3_semiMarkov_williams_ac_mstate.txt", header = TRUE)
colnames(williams2017)

williams_FC_4yr_surv <- subset(williams2017, treat == 0 & time == 4)$surv
williams_FC_4yr_surv_PFS <- subset(williams2017, treat == 0 & time == 4)$pstate1
williams_FC_4yr_surv_Prog <- subset(williams2017, treat == 0 & time == 4)$pstate2

williams_RFC_4yr_surv <- subset(williams2017, treat == 1 & time == 4)$surv
williams_RFC_4yr_surv_PFS <- subset(williams2017, treat == 1 & time == 4)$pstate1
williams_RFC_4yr_surv_Prog <- subset(williams2017, treat == 1 & time == 4)$pstate2

williams_FC_8yr_surv <- subset(williams2017, treat == 0 & time == 8)$surv
williams_FC_8yr_surv_PFS <- subset(williams2017, treat == 0 & time == 8)$pstate1
williams_FC_8yr_surv_Prog <- subset(williams2017, treat == 0 & time == 8)$pstate2

williams_RFC_8yr_surv <- subset(williams2017, treat == 1 & time == 8)$surv
williams_RFC_8yr_surv_PFS <- subset(williams2017, treat == 1 & time == 8)$pstate1
williams_RFC_8yr_surv_Prog <- subset(williams2017, treat == 1 & time == 8)$pstate2

williams_FC_15yr_surv <- subset(williams2017, treat == 0 & time == 15)$surv
williams_FC_15yr_surv_PFS <- subset(williams2017, treat == 0 & time == 15)$pstate1
williams_FC_15yr_surv_Prog <- subset(williams2017, treat == 0 & time == 15)$pstate2

williams_RFC_15yr_surv <- subset(williams2017, treat == 1 & time == 15)$surv
williams_RFC_15yr_surv_PFS <- subset(williams2017, treat == 1 & time == 15)$pstate1
williams_RFC_15yr_surv_Prog <- subset(williams2017, treat == 1 & time == 15)$pstate2

williams_FC_lifetime_surv <- 0
williams_FC_lifetime_surv_PFS <- 0
williams_FC_lifetime_surv_Prog <- 0

williams_RFC_lifetime_surv <- 0
williams_RFC_lifetime_surv_PFS <- 0
williams_RFC_lifetime_surv_Prog <- 0

###############################################################################
## Semi-Markov with standard parametric models, all-cause survival framework ##
###############################################################################
## This is to show the results are identical to Williams 2017 .
## Read the results
semiMarkov_williams_ac_FC <- readRDS("../Data/04a3_semiMarkov_williams_ac_microsim_FC.rds")
semiMarkov_williams_ac_RFC <- readRDS("../Data/04a3_semiMarkov_williams_ac_microsim_RFC.rds")

## Obtain survival at times 4, 8, 15, and 50 years/lifetime
prev_FC <- as.data.frame(semiMarkov_williams_ac_FC$prev)
prev_FC$prob <- prev_FC$number / semiMarkov_williams_ac_FC$n
prev_FC$time <- round(prev_FC$time, digits = 1)

prev_RFC <- as.data.frame(semiMarkov_williams_ac_RFC$prev)
prev_RFC$prob <- prev_RFC$number / semiMarkov_williams_ac_RFC$n
prev_RFC$time <- round(prev_RFC$time, digits = 1)

semiMarkov_williams_ac_FC_4yr_surv_PFS <- prev_FC$prob[prev_FC$time == 4 & prev_FC$state == "PFS"]
semiMarkov_williams_ac_FC_4yr_surv_Prog <- prev_FC$prob[prev_FC$time == 4 & prev_FC$state == "Prog"]
semiMarkov_williams_ac_FC_4yr_surv <- semiMarkov_williams_ac_FC_4yr_surv_PFS + semiMarkov_williams_ac_FC_4yr_surv_Prog

semiMarkov_williams_ac_RFC_4yr_surv_PFS <- prev_RFC$prob[prev_RFC$time == 4 & prev_RFC$state == "PFS"]
semiMarkov_williams_ac_RFC_4yr_surv_Prog <- prev_RFC$prob[prev_RFC$time == 4 & prev_RFC$state == "Prog"]
semiMarkov_williams_ac_RFC_4yr_surv <- semiMarkov_williams_ac_RFC_4yr_surv_PFS + semiMarkov_williams_ac_RFC_4yr_surv_Prog

semiMarkov_williams_ac_FC_8yr_surv_PFS <- prev_FC$prob[prev_FC$time == 8 & prev_FC$state == "PFS"]
semiMarkov_williams_ac_FC_8yr_surv_Prog <- prev_FC$prob[prev_FC$time == 8 & prev_FC$state == "Prog"]
semiMarkov_williams_ac_FC_8yr_surv <- semiMarkov_williams_ac_FC_8yr_surv_PFS + semiMarkov_williams_ac_FC_8yr_surv_Prog

semiMarkov_williams_ac_RFC_8yr_surv_PFS <- prev_RFC$prob[prev_RFC$time == 8 & prev_RFC$state == "PFS"]
semiMarkov_williams_ac_RFC_8yr_surv_Prog <- prev_RFC$prob[prev_RFC$time == 8 & prev_RFC$state == "Prog"]
semiMarkov_williams_ac_RFC_8yr_surv <- semiMarkov_williams_ac_RFC_8yr_surv_PFS + semiMarkov_williams_ac_RFC_8yr_surv_Prog

semiMarkov_williams_ac_FC_15yr_surv_PFS <- prev_FC$prob[prev_FC$time == 15 & prev_FC$state == "PFS"]
semiMarkov_williams_ac_FC_15yr_surv_PFS <- 0
semiMarkov_williams_ac_FC_15yr_surv_Prog <- prev_FC$prob[prev_FC$time == 15 & prev_FC$state == "Prog"]
semiMarkov_williams_ac_FC_15yr_surv <- semiMarkov_williams_ac_FC_15yr_surv_PFS + semiMarkov_williams_ac_FC_15yr_surv_Prog

semiMarkov_williams_ac_RFC_15yr_surv_PFS <- prev_RFC$prob[prev_RFC$time == 15 & prev_RFC$state == "PFS"]
semiMarkov_williams_ac_RFC_15yr_surv_PFS <- 0
semiMarkov_williams_ac_RFC_15yr_surv_Prog <- prev_RFC$prob[prev_RFC$time == 15 & prev_RFC$state == "Prog"]
semiMarkov_williams_ac_RFC_15yr_surv <- semiMarkov_williams_ac_RFC_15yr_surv_PFS + semiMarkov_williams_ac_RFC_15yr_surv_Prog

semiMarkov_williams_ac_FC_lifetime_surv <- 0
semiMarkov_williams_ac_FC_lifetime_surv_PFS <- 0
semiMarkov_williams_ac_FC_lifetime_surv_Prog <- 0

semiMarkov_williams_ac_RFC_lifetime_surv <- 0
semiMarkov_williams_ac_RFC_lifetime_surv_PFS <- 0
semiMarkov_williams_ac_RFC_lifetime_surv_Prog <- 0

###############################################################################
## Semi-Markov with flexible parametric models, all-cause survival framework ##
###############################################################################
## Read the results
semiMarkov_fpm_ac_FC <- readRDS("../Data/04a1_semiMarkov_fpm_ac_microsim_FC.rds")
semiMarkov_fpm_ac_RFC <- readRDS("../Data/04a1_semiMarkov_fpm_ac_microsim_RFC.rds")

## Obtain survival at times 4, 8, 15, and 50 years/lifetime
prev_FC <- as.data.frame(semiMarkov_fpm_ac_FC$prev)
prev_FC$prob <- prev_FC$number / semiMarkov_fpm_ac_FC$n
prev_FC$time <- round(prev_FC$time, digits = 1)

prev_RFC <- as.data.frame(semiMarkov_fpm_ac_RFC$prev)
prev_RFC$prob <- prev_RFC$number / semiMarkov_fpm_ac_RFC$n
prev_RFC$time <- round(prev_RFC$time, digits = 1)

semiMarkov_fpm_ac_FC_4yr_surv_PFS <- prev_FC$prob[prev_FC$time == 4 & prev_FC$state == "PFS"]
semiMarkov_fpm_ac_FC_4yr_surv_Prog <- prev_FC$prob[prev_FC$time == 4 & prev_FC$state == "Prog"]
semiMarkov_fpm_ac_FC_4yr_surv <- semiMarkov_fpm_ac_FC_4yr_surv_PFS + semiMarkov_fpm_ac_FC_4yr_surv_Prog

semiMarkov_fpm_ac_RFC_4yr_surv_PFS <- prev_RFC$prob[prev_RFC$time == 4 & prev_RFC$state == "PFS"]
semiMarkov_fpm_ac_RFC_4yr_surv_Prog <- prev_RFC$prob[prev_RFC$time == 4 & prev_RFC$state == "Prog"]
semiMarkov_fpm_ac_RFC_4yr_surv <- semiMarkov_fpm_ac_RFC_4yr_surv_PFS + semiMarkov_fpm_ac_RFC_4yr_surv_Prog

semiMarkov_fpm_ac_FC_8yr_surv_PFS <- prev_FC$prob[prev_FC$time == 8 & prev_FC$state == "PFS"]
semiMarkov_fpm_ac_FC_8yr_surv_Prog <- prev_FC$prob[prev_FC$time == 8 & prev_FC$state == "Prog"]
semiMarkov_fpm_ac_FC_8yr_surv <- semiMarkov_fpm_ac_FC_8yr_surv_PFS + semiMarkov_fpm_ac_FC_8yr_surv_Prog

semiMarkov_fpm_ac_RFC_8yr_surv_PFS <- prev_RFC$prob[prev_RFC$time == 8 & prev_RFC$state == "PFS"]
semiMarkov_fpm_ac_RFC_8yr_surv_Prog <- prev_RFC$prob[prev_RFC$time == 8 & prev_RFC$state == "Prog"]
semiMarkov_fpm_ac_RFC_8yr_surv <- semiMarkov_fpm_ac_RFC_8yr_surv_PFS + semiMarkov_fpm_ac_RFC_8yr_surv_Prog

semiMarkov_fpm_ac_FC_15yr_surv_PFS <- prev_FC$prob[prev_FC$time == 15 & prev_FC$state == "PFS"]
semiMarkov_fpm_ac_FC_15yr_surv_PFS <- 0
semiMarkov_fpm_ac_FC_15yr_surv_Prog <- prev_FC$prob[prev_FC$time == 15 & prev_FC$state == "Prog"]
semiMarkov_fpm_ac_FC_15yr_surv <- semiMarkov_fpm_ac_FC_15yr_surv_PFS + semiMarkov_fpm_ac_FC_15yr_surv_Prog

semiMarkov_fpm_ac_RFC_15yr_surv_PFS <- prev_RFC$prob[prev_RFC$time == 15 & prev_RFC$state == "PFS"]
semiMarkov_fpm_ac_RFC_15yr_surv_Prog <- prev_RFC$prob[prev_RFC$time == 15 & prev_RFC$state == "Prog"]
semiMarkov_fpm_ac_RFC_15yr_surv <- semiMarkov_fpm_ac_RFC_15yr_surv_PFS + semiMarkov_fpm_ac_RFC_15yr_surv_Prog

semiMarkov_fpm_ac_FC_50yr_surv_PFS <- prev_FC$prob[prev_FC$time == 50 & prev_FC$state == "PFS"]
semiMarkov_fpm_ac_FC_50yr_surv_PFS <- 0
semiMarkov_fpm_ac_FC_50yr_surv_Prog <- prev_FC$prob[prev_FC$time == 50 & prev_FC$state == "Prog"]
semiMarkov_fpm_ac_FC_50yr_surv <- semiMarkov_fpm_ac_FC_50yr_surv_PFS + semiMarkov_fpm_ac_FC_50yr_surv_Prog

semiMarkov_fpm_ac_RFC_50yr_surv_PFS <- prev_RFC$prob[prev_RFC$time == 50 & prev_RFC$state == "PFS"]
semiMarkov_fpm_ac_RFC_50yr_surv_PFS <- 0
semiMarkov_fpm_ac_RFC_50yr_surv_Prog <- prev_RFC$prob[prev_RFC$time == 50 & prev_RFC$state == "Prog"]
semiMarkov_fpm_ac_RFC_50yr_surv <- semiMarkov_fpm_ac_RFC_50yr_surv_PFS + semiMarkov_fpm_ac_RFC_50yr_surv_Prog

semiMarkov_fpm_ac_FC_lifetime_surv <- semiMarkov_fpm_ac_FC_50yr_surv
semiMarkov_fpm_ac_FC_lifetime_surv_PFS <- semiMarkov_fpm_ac_FC_50yr_surv_PFS
semiMarkov_fpm_ac_FC_lifetime_surv_Prog <- semiMarkov_fpm_ac_FC_50yr_surv_Prog

semiMarkov_fpm_ac_RFC_lifetime_surv <- semiMarkov_fpm_ac_RFC_50yr_surv
semiMarkov_fpm_ac_RFC_lifetime_surv_PFS <- semiMarkov_fpm_ac_RFC_50yr_surv_PFS
semiMarkov_fpm_ac_RFC_lifetime_surv_Prog <- semiMarkov_fpm_ac_RFC_50yr_surv_Prog


###################################################################################
## Semi-Markov with flexible parametric models, relative survival framework ##
###################################################################################
## Read the results
semiMarkov_fpm_rel_FC <- readRDS("../Data/04a2_semiMarkov_fpm_rel_microsim_FC.rds")
semiMarkov_fpm_rel_RFC <- readRDS("../Data/04a2_semiMarkov_fpm_rel_microsim_RFC.rds")

## Obtain survival at times 4, 8, 15, and 50 years/lifetime
prev_FC <- as.data.frame(semiMarkov_fpm_rel_FC$prev)
prev_FC$prob <- prev_FC$number / semiMarkov_fpm_rel_FC$n
prev_FC$time <- round(prev_FC$time, digits = 1)

prev_RFC <- as.data.frame(semiMarkov_fpm_rel_RFC$prev)
prev_RFC$prob <- prev_RFC$number / semiMarkov_fpm_rel_RFC$n
prev_RFC$time <- round(prev_RFC$time, digits = 1)

semiMarkov_fpm_rel_FC_4yr_surv_PFS <- prev_FC$prob[prev_FC$time == 4 & prev_FC$state == "PFS"]
semiMarkov_fpm_rel_FC_4yr_surv_Prog <- prev_FC$prob[prev_FC$time == 4 & prev_FC$state == "Prog"]
semiMarkov_fpm_rel_FC_4yr_surv <- semiMarkov_fpm_rel_FC_4yr_surv_PFS + semiMarkov_fpm_rel_FC_4yr_surv_Prog

semiMarkov_fpm_rel_RFC_4yr_surv_PFS <- prev_RFC$prob[prev_RFC$time == 4 & prev_RFC$state == "PFS"]
semiMarkov_fpm_rel_RFC_4yr_surv_Prog <- prev_RFC$prob[prev_RFC$time == 4 & prev_RFC$state == "Prog"]
semiMarkov_fpm_rel_RFC_4yr_surv <- semiMarkov_fpm_rel_RFC_4yr_surv_PFS + semiMarkov_fpm_rel_RFC_4yr_surv_Prog

semiMarkov_fpm_rel_FC_8yr_surv_PFS <- prev_FC$prob[prev_FC$time == 8 & prev_FC$state == "PFS"]
semiMarkov_fpm_rel_FC_8yr_surv_Prog <- prev_FC$prob[prev_FC$time == 8 & prev_FC$state == "Prog"]
semiMarkov_fpm_rel_FC_8yr_surv <- semiMarkov_fpm_rel_FC_8yr_surv_PFS + semiMarkov_fpm_rel_FC_8yr_surv_Prog

semiMarkov_fpm_rel_RFC_8yr_surv_PFS <- prev_RFC$prob[prev_RFC$time == 8 & prev_RFC$state == "PFS"]
semiMarkov_fpm_rel_RFC_8yr_surv_Prog <- prev_RFC$prob[prev_RFC$time == 8 & prev_RFC$state == "Prog"]
semiMarkov_fpm_rel_RFC_8yr_surv <- semiMarkov_fpm_rel_RFC_8yr_surv_PFS + semiMarkov_fpm_rel_RFC_8yr_surv_Prog

semiMarkov_fpm_rel_FC_15yr_surv_PFS <- prev_FC$prob[prev_FC$time == 15 & prev_FC$state == "PFS"]
semiMarkov_fpm_rel_FC_15yr_surv_Prog <- prev_FC$prob[prev_FC$time == 15 & prev_FC$state == "Prog"]
semiMarkov_fpm_rel_FC_15yr_surv <- semiMarkov_fpm_rel_FC_15yr_surv_PFS + semiMarkov_fpm_rel_FC_15yr_surv_Prog

semiMarkov_fpm_rel_RFC_15yr_surv_PFS <- prev_RFC$prob[prev_RFC$time == 15 & prev_RFC$state == "PFS"]
semiMarkov_fpm_rel_RFC_15yr_surv_Prog <- prev_RFC$prob[prev_RFC$time == 15 & prev_RFC$state == "Prog"]
semiMarkov_fpm_rel_RFC_15yr_surv <- semiMarkov_fpm_rel_RFC_15yr_surv_PFS + semiMarkov_fpm_rel_RFC_15yr_surv_Prog

semiMarkov_fpm_rel_FC_lifetime_surv <- 0
semiMarkov_fpm_rel_FC_lifetime_surv_PFS <- 0
semiMarkov_fpm_rel_FC_lifetime_surv_Prog <- 0

semiMarkov_fpm_rel_RFC_lifetime_surv <- 0
semiMarkov_fpm_rel_RFC_lifetime_surv_PFS <- 0
semiMarkov_fpm_rel_RFC_lifetime_surv_Prog <- 0

##############################################################
##============================================================
## Calculate restricted mean survival time (RMST)
##============================================================
##############################################################
var_list <- list()

##################
## Fischer 2016 ##
##################
## 4 years OS (PFS + Progression) 
rmst_mfit2 <- rmst2(fischer2016$time_y, fischer2016$event, arm = fischer2016$treat, tau = 4)
var_list[["fischer_RFC_4yr"]] <- rmst_mfit2$RMST.arm1$rmst[[1]]
print(paste("Fischer 2016,", "RFC,", "AUC for time <= 4 years:", var_list[["fischer_RFC_4yr"]]))
var_list[["fischer_FC_4yr"]] <- rmst_mfit2$RMST.arm0$rmst[[1]]
print(paste("Fischer 2016,", "FC,", "AUC for time <= 4 years:", var_list[["fischer_FC_4yr"]]))

## 4 years PFS
rmst_mfit2_pfs <- rmst2(fischer2016_pfs$time_y, fischer2016_pfs$event, arm = fischer2016_pfs$treat, tau = 4)
var_list[["fischer_RFC_4yr_pfs"]] <- rmst_mfit2_pfs$RMST.arm1$rmst[[1]]
print(paste("Fischer 2016,", "RFC,", "AUC for PFS for time <= 4 years:", var_list[["fischer_RFC_4yr_pfs"]]))
var_list[["fischer_FC_4yr_pfs"]] <- rmst_mfit2_pfs$RMST.arm0$rmst[[1]]
print(paste("Fischer 2016,", "FC,", "AUC for PFS for time <= 4 years:", var_list[["fischer_FC_4yr_pfs"]]))

## 4 years Progression
## We calculate time by subtracting PFS from OS
var_list[["fischer_RFC_4yr_prog"]] <- var_list[["fischer_RFC_4yr"]] - var_list[["fischer_RFC_4yr_pfs"]]
var_list[["fischer_FC_4yr_prog"]] <- var_list[["fischer_FC_4yr"]] - var_list[["fischer_FC_4yr_pfs"]]

## 8 years OS (PFS + Progression)
rmst_mfit3 <- rmst2(fischer2016$time_y, fischer2016$event, arm = fischer2016$treat, tau = 8)
var_list[["fischer_RFC_8yr"]] <- rmst_mfit3$RMST.arm1$rmst[[1]]
print(paste("Fischer 2016,", "RFC,", "AUC for time <= 3.94 years:", var_list[["fischer_RFC_8yr"]]))
var_list[["fischer_FC_8yr"]] <- rmst_mfit3$RMST.arm0$rmst[[1]]
print(paste("Fischer 2016,", "FC,", "AUC for time <= 3.94 years:", var_list[["fischer_FC_8yr"]]))

## 8 years PFS
rmst_mfit3_pfs <- rmst2(fischer2016_pfs$time_y, fischer2016_pfs$event, arm = fischer2016_pfs$treat, tau = 8)
var_list[["fischer_RFC_8yr_pfs"]] <- rmst_mfit3_pfs$RMST.arm1$rmst[[1]]
print(paste("Fischer 2016,", "RFC,", "AUC for PFS for time <= 8 years:", var_list[["fischer_RFC_8yr_pfs"]]))
var_list[["fischer_FC_8yr_pfs"]] <- rmst_mfit3_pfs$RMST.arm0$rmst[[1]]
print(paste("Fischer 2016,", "FC,", "AUC for PFS for time <= 8 years:", var_list[["fischer_FC_8yr_pfs"]]))

## 8 years Progression
## We calculate time by subtracting PFS from OS
var_list[["fischer_RFC_8yr_prog"]] <- var_list[["fischer_RFC_8yr"]] - var_list[["fischer_RFC_8yr_pfs"]]
var_list[["fischer_FC_8yr_prog"]] <- var_list[["fischer_FC_8yr"]] - var_list[["fischer_FC_8yr_pfs"]]

###############################################################################################
## Semi-Markov with standard parametric models, all-cause survival framework (Williams 2017) ##
###############################################################################################
## Define the time intervals
time_intervals <- c(4, 8, 15, 19)

## Ensure data is sorted by time
williams2017 <- williams2017[order(c(williams2017$treat,williams2017$time)), ]

## Loop over the time intervals
for (i in time_intervals) {
  ## subset data to time <= i
  williams2017_subset <- williams2017[williams2017$time <= i, ]
  
  ## further split into FC and RFC subsets based on 'treat' column
  williams2017_subset_FC <- williams2017_subset[grepl(0, williams2017_subset$treat), ]
  williams2017_subset_RFC <- williams2017_subset[grepl(1, williams2017_subset$treat), ]
  
  ## calculate area under the curve for 'surv' with time <= i for RFC 
  value <- trapz(williams2017_subset_RFC$time, williams2017_subset_RFC$surv)
  var <- paste("williams_RFC_",i,"yr")
  var_list[[var]] <- value
  print(paste("Williams 2017", "RFC", "Area under curve for time <= ", i, "and treat = RFC: ", value))
  
  ## calculate area under the curve for 'surv' with time <= i for FC
  value <- trapz(williams2017_subset_FC$time, williams2017_subset_FC$surv)
  var <- paste("williams_FC_",i,"yr")
  var_list[[var]] <- value
  print(paste("Williams 2017", "FC", "Area under curve for time <= ", i, "and treat = FC: ", value))
}

# PFS
for (i in time_intervals) {
  ## subset data to time <= i
  williams2017_subset <- williams2017[williams2017$time <= i, ]
  
  ## further split into FC and RFC subsets based on 'treat' column
  williams2017_subset_FC <- williams2017_subset[grepl(0, williams2017_subset$treat), ]
  williams2017_subset_RFC <- williams2017_subset[grepl(1, williams2017_subset$treat), ]
  
  ## calculate area under the curve for 'surv' with time <= i for RFC 
  value <- trapz(williams2017_subset_RFC$time, williams2017_subset_RFC$pstate1)
  var <- paste("williams_RFC_",i,"yr","_pfs")
  var_list[[var]] <- value
  print(paste("Williams 2017", "RFC", "Area under curve (PFS) for time <= ", i, "and treat = RFC: ", value))
  
  ## calculate area under the curve for 'surv' with time <= i for FC
  value <- trapz(williams2017_subset_FC$time, williams2017_subset_FC$pstate1)
  var <- paste("williams_FC_",i,"yr","_pfs")
  var_list[[var]] <- value
  print(paste("Williams 2017", "FC", "Area under curve  (PFS) for time <= ", i, "and treat = FC: ", value))
}

# Progression
for (i in time_intervals) {
  ## subset data to time <= i
  williams2017_subset <- williams2017[williams2017$time <= i, ]
  
  ## further split into FC and RFC subsets based on 'treat' column
  williams2017_subset_FC <- williams2017_subset[grepl(0, williams2017_subset$treat), ]
  williams2017_subset_RFC <- williams2017_subset[grepl(1, williams2017_subset$treat), ]
  
  ## calculate area under the curve for 'surv' with time <= i for RFC 
  value <- trapz(williams2017_subset_RFC$time, williams2017_subset_RFC$pstate2)
  var <- paste("williams_RFC_",i,"yr","_prog")
  var_list[[var]] <- value
  print(paste("Williams 2017", "RFC", "Area under curve (Progression) for time <= ", i, "and treat = RFC: ", value))
  
  ## calculate area under the curve for 'surv' with time <= i for FC
  value <- trapz(williams2017_subset_FC$time, williams2017_subset_FC$pstate2)
  var <- paste("williams_FC_",i,"yr","_prog")
  var_list[[var]] <- value
  print(paste("Williams 2017", "FC", "Area under curve  (Progression) for time <= ", i, "and treat = FC: ", value))
}

###############################################################################
## Semi-Markov with standard parametric models, all-cause survival framework ##
###############################################################################
## Define the time intervals
time_intervals <- c(4, 8, 15, 50)

## loop over the time intervals
## OS
for (i in time_intervals) {
  report <- lapply(list(semiMarkov_williams_ac_FC, semiMarkov_williams_ac_RFC), function(sim){
    # create a subset of the data where 'time' is less than or equal to i
    sim$pt_subset <- subset(sim$pt, time <= i)
    
    # calculate LY using the subset and print it
    LY <- (sum(sim$pt_subset$pt) - sum(subset(sim$pt_subset, state %in% c("AcD"))$pt)) / sim$n
  })
  
  ## FC, treat = 0
  value <- report[[1]]
  var <- paste("semiMarkov_williams_ac_FC",i,"yr")
  var_list[[var]] <- value
  print(paste("Semi-Markov: SPMs, ASF", "FC", "Area under curve for time <=", i, ": ", value))
  
  ## RFC, treat = 1
  value <- report[[2]]
  var <- paste("semiMarkov_williams_ac_RFC",i,"yr")
  var_list[[var]] <- value
  print(paste("Semi-Markov: SPMs, ASF", "RFC", "Area under curve for time <=", i, ": ", value))
}

## PFS 
## loop over the time intervals
for (i in time_intervals) {
  report <- lapply(list(semiMarkov_williams_ac_FC, semiMarkov_williams_ac_RFC), function(sim){
    # create a subset of the data where 'time' is less than or equal to i
    sim$pt_subset <- subset(sim$pt, time <= i)
    
    # calculate LY using the subset and print it
    LY_pfs <- (sum(subset(sim$pt_subset, state %in% c("PFS"))$pt)) / sim$n
  })
  
  ## FC, treat = 0
  value <- report[[1]]
  var <- paste("semiMarkov_williams_ac_FC",i,"yr","_pfs")
  var_list[[var]] <- value
  print(paste("Semi-Markov: SPMs, ASF", "FC", "Area under curve (PFS) for time <=", i, ": ", value))
  
  ## RFC, treat = 1
  value <- report[[2]]
  var <- paste("semiMarkov_williams_ac_RFC",i,"yr","_pfs")
  var_list[[var]] <- value
  print(paste("Semi-Markov: SPMs, ASF", "RFC", "Area under curve (PFS) for time <=", i, ": ", value))
}

## Progression
## loop over the time intervals
for (i in time_intervals) {
  report <- lapply(list(semiMarkov_williams_ac_FC, semiMarkov_williams_ac_RFC), function(sim){
    # create a subset of the data where 'time' is less than or equal to i
    sim$pt_subset <- subset(sim$pt, time <= i)
    
    # calculate LY using the subset and print it
    LY_prog <- (sum(subset(sim$pt_subset, state %in% c("Prog"))$pt)) / sim$n
  })
  
  ## FC, treat = 0
  value <- report[[1]]
  var <- paste("semiMarkov_williams_ac_FC",i,"yr","_prog")
  var_list[[var]] <- value
  print(paste("Semi-Markov: SPMs, ASF", "FC", "Area under curve (Progression) for time <=", i, ": ", value))
  
  ## RFC, treat = 1
  value <- report[[2]]
  var <- paste("semiMarkov_williams_ac_RFC",i,"yr","_prog")
  var_list[[var]] <- value
  print(paste("Semi-Markov: SPMs, ASF", "RFC", "Area under curve (Progression) for time <=", i, ": ", value))
}

###############################################################################
## Semi-Markov with flexible parametric models, all-cause survival framework ##
###############################################################################
## Define the time intervals
time_intervals <- c(4, 8, 15, 50)

## loop over the time intervals
## OS
for (i in time_intervals) {
  report <- lapply(list(semiMarkov_fpm_ac_FC, semiMarkov_fpm_ac_RFC), function(sim){
    # create a subset of the data where 'time' is less than or equal to i
    sim$pt_subset <- subset(sim$pt, time <= i)
    
    # calculate LY using the subset and print it
    LY <- (sum(sim$pt_subset$pt) - sum(subset(sim$pt_subset, state %in% c("AcD"))$pt)) / sim$n
  })
  
  ## FC, treat = 0
  value <- report[[1]]
  var <- paste("semiMarkov_fpm_ac_FC",i,"yr")
  var_list[[var]] <- value
  print(paste("Semi-Markov: FPMs, ASF", "FC", "Area under curve for time <=", i, ": ", value))
  
  ## RFC, treat = 1
  value <- report[[2]]
  var <- paste("semiMarkov_fpm_ac_RFC",i,"yr")
  var_list[[var]] <- value
  print(paste("Semi-Markov: FPMs, ASF", "RFC", "Area under curve for time <=", i, ": ", value))
}


## PFS
for (i in time_intervals) {
  report <- lapply(list(semiMarkov_fpm_ac_FC, semiMarkov_fpm_ac_RFC), function(sim){
    # create a subset of the data where 'time' is less than or equal to i
    sim$pt_subset <- subset(sim$pt, time <= i)
    
    # calculate LY using the subset and print it
    LY <- (sum(subset(sim$pt_subset, state %in% c("PFS"))$pt)) / sim$n
  })
  
  ## FC, treat = 0
  value <- report[[1]]
  var <- paste("semiMarkov_fpm_ac_FC",i,"yr","_pfs")
  var_list[[var]] <- value
  print(paste("Semi-Markov: FPMs, ASF", "FC", "Area under curve (PFS) for time <=", i, ": ", value))
  
  ## RFC, treat = 1
  value <- report[[2]]
  var <- paste("semiMarkov_fpm_ac_RFC",i,"yr","_pfs")
  var_list[[var]] <- value
  print(paste("Semi-Markov: FPMs, ASF", "RFC", "Area under curve (PFS) for time <=", i, ": ", value))
}

## Progression
for (i in time_intervals) {
  report <- lapply(list(semiMarkov_fpm_ac_FC, semiMarkov_fpm_ac_RFC), function(sim){
    # create a subset of the data where 'time' is less than or equal to i
    sim$pt_subset <- subset(sim$pt, time <= i)
    
    # calculate LY using the subset and print it
    LY <- (sum(subset(sim$pt_subset, state %in% c("Prog"))$pt)) / sim$n
  })
  
  ## FC, treat = 0
  value <- report[[1]]
  var <- paste("semiMarkov_fpm_ac_FC",i,"yr","_prog")
  var_list[[var]] <- value
  print(paste("Semi-Markov: FPMs, ASF", "FC", "Area under curve (Progression) for time <=", i, ": ", value))
  
  ## RFC, treat = 1
  value <- report[[2]]
  var <- paste("semiMarkov_fpm_ac_RFC",i,"yr","_prog")
  var_list[[var]] <- value
  print(paste("Semi-Markov: FPMs, ASF", "RFC", "Area under curve (Progression) for time <=", i, ": ", value))
}

###############################################################################
## Semi-Markov with flexible parametric models, relative survival framework ##
###############################################################################
## Define the time intervals
time_intervals <- c(4, 8, 15, 50)

## OS
for (i in time_intervals) {
  report <- lapply(list(semiMarkov_fpm_rel_FC, semiMarkov_fpm_rel_RFC), function(sim){
    # create a subset of the data where 'time' is less than or equal to i
    sim$pt_subset <- subset(sim$pt, time <= i)
    
    # calculate LY using the subset and print it
    LY <- (sum(sim$pt_subset$pt) - sum(subset(sim$pt_subset, state %in% c("ExpD", "ExcD"))$pt)) / sim$n
  })
  
  ## FC, treat = 0
  value <- report[[1]]
  var <- paste("semiMarkov_fpm_rel_FC",i,"yr")
  var_list[[var]] <- value
  print(paste("Semi-Markov: FPMs, RSF", "FC", "Area under curve for time <=", i, ": ", value))
  
  ## RFC, treat = 1
  value <- report[[2]]
  var <- paste("semiMarkov_fpm_rel_RFC",i,"yr")
  var_list[[var]] <- value
  print(paste("Semi-Markov: FPMs, RSF", "RFC", "Area under curve for time <=", i, ": ", value))
}

## PFS
for (i in time_intervals) {
  report <- lapply(list(semiMarkov_fpm_rel_FC, semiMarkov_fpm_rel_RFC), function(sim){
    # create a subset of the data where 'time' is less than or equal to i
    sim$pt_subset <- subset(sim$pt, time <= i)
    
    # calculate LY using the subset and print it
    LY_pfs <- (sum(subset(sim$pt_subset, state %in% c("PFS"))$pt)) / sim$n
  })
  
  ## FC, treat = 0
  value <- report[[1]]
  var <- paste("semiMarkov_fpm_rel_FC",i,"yr","_pfs")
  var_list[[var]] <- value
  print(paste("Semi-Markov: FPMs, RSF", "FC", "Area under curve (PFS) for time <=", i, ": ", value))
  
  ## RFC, treat = 1
  value <- report[[2]]
  var <- paste("semiMarkov_fpm_rel_RFC",i,"yr","_pfs")
  var_list[[var]] <- value
  print(paste("Semi-Markov: FPMs, RSF", "RFC", "Area under curve (PFS) for time <=", i, ": ", value))
}

## Progression
for (i in time_intervals) {
  report <- lapply(list(semiMarkov_fpm_rel_FC, semiMarkov_fpm_rel_RFC), function(sim){
    # create a subset of the data where 'time' is less than or equal to i
    sim$pt_subset <- subset(sim$pt, time <= i)
    
    # calculate LY using the subset and print it
    LY_prog <- (sum(subset(sim$pt_subset, state %in% c("Prog"))$pt)) / sim$n
  })
  
  ## FC, treat = 0
  value <- report[[1]]
  var <- paste("semiMarkov_fpm_rel_FC",i,"yr","_prog")
  var_list[[var]] <- value
  print(paste("Semi-Markov: FPMs, RSF", "FC", "Area under curve (Progression) for time <=", i, ": ", value))
  
  ## RFC, treat = 1
  value <- report[[2]]
  var <- paste("semiMarkov_fpm_rel_RFC",i,"yr","_prog")
  var_list[[var]] <- value
  print(paste("Semi-Markov: FPMs, RSF", "RFC", "Area under curve (Progression) for time <=", i, ": ", value))
}

###################################
## Store the results in one list ##
###################################
variables <- c("fischer_FC_4yr_surv", "fischer_RFC_4yr_surv", "fischer_FC_8yr_surv", "fischer_RFC_8yr_surv",
               "williams_FC_4yr_surv", "williams_RFC_4yr_surv", "williams_FC_8yr_surv", "williams_RFC_8yr_surv", "williams_FC_15yr_surv", "williams_RFC_15yr_surv", "williams_FC_lifetime_surv", "williams_RFC_lifetime_surv",
               "semiMarkov_williams_ac_FC_4yr_surv", "semiMarkov_williams_ac_RFC_4yr_surv", "semiMarkov_williams_ac_FC_8yr_surv", "semiMarkov_williams_ac_RFC_8yr_surv", 
                  "semiMarkov_williams_ac_FC_15yr_surv", "semiMarkov_williams_ac_RFC_15yr_surv", "semiMarkov_williams_ac_FC_lifetime_surv", "semiMarkov_williams_ac_RFC_lifetime_surv", 
               "semiMarkov_fpm_ac_FC_4yr_surv", "semiMarkov_fpm_ac_RFC_4yr_surv", "semiMarkov_fpm_ac_FC_8yr_surv", "semiMarkov_fpm_ac_RFC_8yr_surv", 
                      "semiMarkov_fpm_ac_FC_15yr_surv", "semiMarkov_fpm_ac_RFC_15yr_surv", "semiMarkov_fpm_ac_FC_lifetime_surv", "semiMarkov_fpm_ac_RFC_lifetime_surv",
               "semiMarkov_fpm_rel_FC_4yr_surv", "semiMarkov_fpm_rel_RFC_4yr_surv", "semiMarkov_fpm_rel_FC_8yr_surv", "semiMarkov_fpm_rel_RFC_8yr_surv", 
                      "semiMarkov_fpm_rel_FC_15yr_surv", "semiMarkov_fpm_rel_RFC_15yr_surv", "semiMarkov_fpm_rel_FC_lifetime_surv", "semiMarkov_fpm_rel_RFC_lifetime_surv",
               
               "fischer_FC_4yr_pfs", "fischer_RFC_4yr_pfs", "fischer_FC_8yr_pfs", "fischer_RFC_8yr_pfs",
               "williams_FC_4yr_surv_PFS", "williams_RFC_4yr_surv_PFS", "williams_FC_8yr_surv_PFS", "williams_RFC_8yr_surv_PFS", "williams_FC_15yr_surv_PFS", "williams_RFC_15yr_surv_PFS", "williams_FC_lifetime_surv_PFS", "williams_RFC_lifetime_surv_PFS",
               "semiMarkov_williams_ac_FC_4yr_surv_PFS", "semiMarkov_williams_ac_RFC_4yr_surv_PFS", "semiMarkov_williams_ac_FC_8yr_surv_PFS", "semiMarkov_williams_ac_RFC_8yr_surv_PFS", 
               "semiMarkov_williams_ac_FC_15yr_surv_PFS", "semiMarkov_williams_ac_RFC_15yr_surv_PFS", "semiMarkov_williams_ac_FC_lifetime_surv_PFS", "semiMarkov_williams_ac_RFC_lifetime_surv_PFS", 
               "semiMarkov_fpm_ac_FC_4yr_surv_PFS", "semiMarkov_fpm_ac_RFC_4yr_surv_PFS", "semiMarkov_fpm_ac_FC_8yr_surv_PFS", "semiMarkov_fpm_ac_RFC_8yr_surv_PFS", 
               "semiMarkov_fpm_ac_FC_15yr_surv_PFS", "semiMarkov_fpm_ac_RFC_15yr_surv_PFS", "semiMarkov_fpm_ac_FC_lifetime_surv_PFS", "semiMarkov_fpm_ac_RFC_lifetime_surv_PFS",
               "semiMarkov_fpm_rel_FC_4yr_surv_PFS", "semiMarkov_fpm_rel_RFC_4yr_surv_PFS", "semiMarkov_fpm_rel_FC_8yr_surv_PFS", "semiMarkov_fpm_rel_RFC_8yr_surv_PFS", 
               "semiMarkov_fpm_rel_FC_15yr_surv_PFS", "semiMarkov_fpm_rel_RFC_15yr_surv_PFS", "semiMarkov_fpm_rel_FC_lifetime_surv_PFS", "semiMarkov_fpm_rel_RFC_lifetime_surv_PFS",
               
               "fischer_FC_4yr_prog", "fischer_RFC_4yr_prog", "fischer_FC_8yr_prog", "fischer_RFC_8yr_prog",
               "williams_FC_4yr_surv_Prog", "williams_RFC_4yr_surv_Prog", "williams_FC_8yr_surv_Prog", "williams_RFC_8yr_surv_Prog", "williams_FC_15yr_surv_Prog", "williams_RFC_15yr_surv_Prog", "williams_FC_lifetime_surv_Prog", "williams_RFC_lifetime_surv_Prog",
               "semiMarkov_williams_ac_FC_4yr_surv_Prog", "semiMarkov_williams_ac_RFC_4yr_surv_Prog", "semiMarkov_williams_ac_FC_8yr_surv_Prog", "semiMarkov_williams_ac_RFC_8yr_surv_Prog", 
               "semiMarkov_williams_ac_FC_15yr_surv_Prog", "semiMarkov_williams_ac_RFC_15yr_surv_Prog", "semiMarkov_williams_ac_FC_lifetime_surv_Prog", "semiMarkov_williams_ac_RFC_lifetime_surv_Prog", 
               "semiMarkov_fpm_ac_FC_4yr_surv_Prog", "semiMarkov_fpm_ac_RFC_4yr_surv_Prog", "semiMarkov_fpm_ac_FC_8yr_surv_Prog", "semiMarkov_fpm_ac_RFC_8yr_surv_Prog", 
               "semiMarkov_fpm_ac_FC_15yr_surv_Prog", "semiMarkov_fpm_ac_RFC_15yr_surv_Prog", "semiMarkov_fpm_ac_FC_lifetime_surv_Prog", "semiMarkov_fpm_ac_RFC_lifetime_surv_Prog",
               "semiMarkov_fpm_rel_FC_4yr_surv_Prog", "semiMarkov_fpm_rel_RFC_4yr_surv_Prog", "semiMarkov_fpm_rel_FC_8yr_surv_Prog", "semiMarkov_fpm_rel_RFC_8yr_surv_Prog", 
               "semiMarkov_fpm_rel_FC_15yr_surv_Prog", "semiMarkov_fpm_rel_RFC_15yr_surv_Prog", "semiMarkov_fpm_rel_FC_lifetime_surv_Prog", "semiMarkov_fpm_rel_RFC_lifetime_surv_Prog"
               )
               
surv_list <- mget(variables)

surv_list <- lapply(surv_list, function(x) {
  formatted <- format(round(x, digits = 2), nsmall = 2)
  return(formatted)
})

##############################################################
##============================================================
## Tabulate the results
##============================================================
##############################################################
## Round to 2 decimals
var_list <- lapply(var_list, function(x) {
  formatted <- format(round(x, digits = 2), nsmall = 2)
  return(formatted)
})

var_matrix <- matrix(
  c("FC","OS", "Observed", var_list$fischer_FC_4yr, var_list$fischer_FC_8yr, NA, NA, surv_list$fischer_FC_4yr_surv, surv_list$fischer_FC_8yr_surv, NA, NA,
    NA,NA, "SPMs, ASF (Williams et al.)", var_list$`williams_FC_ 4 yr`, var_list$`williams_FC_ 8 yr`, var_list$`williams_FC_ 15 yr`,  var_list$`williams_FC_ 19 yr`, surv_list$williams_FC_4yr_surv, surv_list$williams_FC_8yr_surv, surv_list$williams_FC_15yr_surv,  surv_list$williams_FC_lifetime_surv,
    NA,NA, "FPMs, ASF", var_list$`semiMarkov_fpm_ac_FC 4 yr`, var_list$`semiMarkov_fpm_ac_FC 8 yr`, var_list$`semiMarkov_fpm_ac_FC 15 yr`, var_list$`semiMarkov_fpm_ac_FC 50 yr`, surv_list$semiMarkov_fpm_ac_FC_4yr_surv, surv_list$semiMarkov_fpm_ac_FC_8yr_surv, surv_list$semiMarkov_fpm_ac_FC_15yr_surv, surv_list$semiMarkov_fpm_ac_FC_lifetime_surv,
    NA,NA, "FPMs, RSF (Proposed method)", var_list$`semiMarkov_fpm_rel_FC 4 yr`, var_list$`semiMarkov_fpm_rel_FC 8 yr`, var_list$`semiMarkov_fpm_rel_FC 15 yr`, var_list$`semiMarkov_fpm_rel_FC 50 yr`, surv_list$semiMarkov_fpm_rel_FC_4yr_surv, surv_list$semiMarkov_fpm_rel_FC_8yr_surv, surv_list$semiMarkov_fpm_rel_FC_15yr_surv, surv_list$semiMarkov_fpm_rel_FC_lifetime_surv,
    NA, "Progression-free", "Observed", var_list$fischer_FC_4yr_pfs, var_list$fischer_FC_8yr_pfs, NA, NA, surv_list$fischer_FC_4yr_pfs, surv_list$fischer_FC_8yr_pfs, NA, NA,
    NA,NA, "SPMs, ASF", var_list$`williams_FC_ 4 yr _pfs`, var_list$`williams_FC_ 8 yr _pfs`, var_list$`williams_FC_ 15 yr _pfs`,  var_list$`williams_FC_ 19 yr _pfs`, surv_list$williams_FC_4yr_surv_PFS, surv_list$williams_FC_8yr_surv_PFS, surv_list$williams_FC_15yr_surv_PFS,  surv_list$williams_FC_lifetime_surv_PFS,
    NA,NA, "FPMs, ASF", var_list$`semiMarkov_fpm_ac_FC 4 yr _pfs`, var_list$`semiMarkov_fpm_ac_FC 8 yr _pfs`, var_list$`semiMarkov_fpm_ac_FC 15 yr _pfs`, var_list$`semiMarkov_fpm_ac_FC 50 yr _pfs`, surv_list$semiMarkov_fpm_ac_FC_4yr_surv_PFS, surv_list$semiMarkov_fpm_ac_FC_8yr_surv_PFS, surv_list$semiMarkov_fpm_ac_FC_15yr_surv_PFS, surv_list$semiMarkov_fpm_ac_FC_lifetime_surv_PFS,
    NA,NA, "FPMs, RSF", var_list$`semiMarkov_fpm_rel_FC 4 yr _pfs`, var_list$`semiMarkov_fpm_rel_FC 8 yr _pfs`, var_list$`semiMarkov_fpm_rel_FC 15 yr _pfs`, var_list$`semiMarkov_fpm_rel_FC 50 yr _pfs`, surv_list$semiMarkov_fpm_rel_FC_4yr_surv_PFS, surv_list$semiMarkov_fpm_rel_FC_8yr_surv_PFS, surv_list$semiMarkov_fpm_rel_FC_15yr_surv_PFS, surv_list$semiMarkov_fpm_rel_FC_lifetime_surv_PFS,
    NA,"Progression", "Observed", var_list$fischer_FC_4yr_prog, var_list$fischer_FC_8yr_prog, NA, NA,  surv_list$fischer_FC_4yr_prog,  surv_list$fischer_FC_4yr_prog, NA, NA,
    NA,NA, "SPMs, ASF", var_list$`williams_FC_ 4 yr _prog`, var_list$`williams_FC_ 8 yr _prog`, var_list$`williams_FC_ 15 yr _prog`,  var_list$`williams_FC_ 19 yr _prog`, surv_list$williams_FC_4yr_surv_Prog, surv_list$williams_FC_8yr_surv_Prog, surv_list$williams_FC_15yr_surv_Prog, surv_list$williams_FC_lifetime_surv_Prog,
    NA,NA, "FPMs, ASF", var_list$`semiMarkov_fpm_ac_FC 4 yr _prog`, var_list$`semiMarkov_fpm_ac_FC 8 yr _prog`, var_list$`semiMarkov_fpm_ac_FC 15 yr _prog`, var_list$`semiMarkov_fpm_ac_FC 50 yr _prog`,surv_list$semiMarkov_fpm_ac_FC_4yr_surv_Prog, surv_list$semiMarkov_fpm_ac_FC_8yr_surv_Prog, surv_list$semiMarkov_fpm_ac_FC_15yr_surv_Prog, surv_list$semiMarkov_fpm_ac_FC_lifetime_surv_Prog,
    NA,NA, "FPMs, RSF", var_list$`semiMarkov_fpm_rel_FC 4 yr _prog`, var_list$`semiMarkov_fpm_rel_FC 8 yr _prog`, var_list$`semiMarkov_fpm_rel_FC 15 yr _prog`, var_list$`semiMarkov_fpm_rel_FC 50 yr _prog`, surv_list$semiMarkov_fpm_rel_FC_4yr_surv_Prog, surv_list$semiMarkov_fpm_rel_FC_8yr_surv_Prog, surv_list$semiMarkov_fpm_rel_FC_15yr_surv_Prog, surv_list$semiMarkov_fpm_rel_FC_lifetime_surv_Prog,
    "RFC","OS", "Observed", var_list$fischer_RFC_4yr, var_list$fischer_RFC_8yr, NA, NA, surv_list$fischer_RFC_4yr_surv, surv_list$fischer_RFC_8yr_surv, NA, NA,
    NA,NA, "SPMs, ASF", var_list$`williams_RFC_ 4 yr`, var_list$`williams_RFC_ 8 yr`, var_list$`williams_RFC_ 15 yr`,  var_list$`williams_RFC_ 19 yr`, surv_list$williams_RFC_4yr_surv, surv_list$williams_RFC_8yr_surv, surv_list$williams_RFC_15yr_surv,  surv_list$williams_RFC_lifetime_surv,
    NA,NA, "FPMs, ASF", var_list$`semiMarkov_fpm_ac_RFC 4 yr`, var_list$`semiMarkov_fpm_ac_RFC 8 yr`, var_list$`semiMarkov_fpm_ac_RFC 15 yr`, var_list$`semiMarkov_fpm_ac_RFC 50 yr`, surv_list$semiMarkov_fpm_ac_RFC_4yr_surv, surv_list$semiMarkov_fpm_ac_RFC_8yr_surv, surv_list$semiMarkov_fpm_ac_RFC_15yr_surv, surv_list$semiMarkov_fpm_ac_RFC_lifetime_surv,
    NA,NA, "FPMs, RSF", var_list$`semiMarkov_fpm_rel_RFC 4 yr`, var_list$`semiMarkov_fpm_rel_RFC 8 yr`, var_list$`semiMarkov_fpm_rel_RFC 15 yr`, var_list$`semiMarkov_fpm_rel_RFC 50 yr`, surv_list$semiMarkov_fpm_rel_RFC_4yr_surv, surv_list$semiMarkov_fpm_rel_RFC_8yr_surv, surv_list$semiMarkov_fpm_rel_RFC_15yr_surv, surv_list$semiMarkov_fpm_rel_RFC_lifetime_surv,
    NA,"Progression-free", "Observed", var_list$fischer_RFC_4yr_pfs, var_list$fischer_RFC_8yr_pfs, NA, NA, surv_list$fischer_RFC_4yr_pfs, surv_list$fischer_RFC_8yr_pfs, NA, NA,
    NA,NA, "SPMs, ASF", var_list$`williams_RFC_ 4 yr _pfs`, var_list$`williams_RFC_ 8 yr _pfs`, var_list$`williams_RFC_ 15 yr _pfs`,  var_list$`williams_RFC_ 19 yr _pfs`, surv_list$williams_RFC_4yr_surv_PFS, surv_list$williams_RFC_8yr_surv_PFS, surv_list$williams_RFC_15yr_surv_PFS,  surv_list$williams_RFC_lifetime_surv_PFS,
    NA,NA, "FPMs, ASF", var_list$`semiMarkov_fpm_ac_RFC 4 yr _pfs`, var_list$`semiMarkov_fpm_ac_RFC 8 yr _pfs`, var_list$`semiMarkov_fpm_ac_RFC 15 yr _pfs`, var_list$`semiMarkov_fpm_ac_RFC 50 yr _pfs`, surv_list$semiMarkov_fpm_ac_RFC_4yr_surv_PFS, surv_list$semiMarkov_fpm_ac_RFC_8yr_surv_PFS, surv_list$semiMarkov_fpm_ac_RFC_15yr_surv_PFS, surv_list$semiMarkov_fpm_ac_RFC_lifetime_surv_PFS,
    NA,NA, "FPMs, RSF", var_list$`semiMarkov_fpm_rel_RFC 4 yr _pfs`, var_list$`semiMarkov_fpm_rel_RFC 8 yr _pfs`, var_list$`semiMarkov_fpm_rel_RFC 15 yr _pfs`, var_list$`semiMarkov_fpm_rel_RFC 50 yr _pfs`, surv_list$semiMarkov_fpm_rel_RFC_4yr_surv_PFS, surv_list$semiMarkov_fpm_rel_RFC_8yr_surv_PFS, surv_list$semiMarkov_fpm_rel_RFC_15yr_surv_PFS, surv_list$semiMarkov_fpm_rel_RFC_lifetime_surv_PFS,
    NA,"Progression", "Observed", var_list$fischer_RFC_4yr_prog, var_list$fischer_RFC_8yr_prog, NA, NA, surv_list$fischer_RFC_4yr_prog,  surv_list$fischer_RFC_4yr_prog, NA, NA,
    NA,NA, "SPMs, ASF", var_list$`williams_RFC_ 4 yr _prog`, var_list$`williams_RFC_ 8 yr _prog`, var_list$`williams_RFC_ 15 yr _prog`,  var_list$`williams_RFC_ 19 yr _prog`, surv_list$williams_RFC_4yr_surv_Prog, surv_list$williams_RFC_8yr_surv_Prog, surv_list$williams_RFC_15yr_surv_Prog, surv_list$williams_RFC_lifetime_surv_Prog,
    NA,NA, "FPMs, ASF", var_list$`semiMarkov_fpm_ac_RFC 4 yr _prog`, var_list$`semiMarkov_fpm_ac_RFC 8 yr _prog`, var_list$`semiMarkov_fpm_ac_RFC 15 yr _prog`, var_list$`semiMarkov_fpm_ac_RFC 50 yr _prog`,surv_list$semiMarkov_fpm_ac_RFC_4yr_surv_Prog, surv_list$semiMarkov_fpm_ac_RFC_8yr_surv_Prog, surv_list$semiMarkov_fpm_ac_RFC_15yr_surv_Prog, surv_list$semiMarkov_fpm_ac_RFC_lifetime_surv_Prog,
    NA,NA, "FPMs, RSF", var_list$`semiMarkov_fpm_rel_RFC 4 yr _prog`, var_list$`semiMarkov_fpm_rel_RFC 8 yr _prog`, var_list$`semiMarkov_fpm_rel_RFC 15 yr _prog`, var_list$`semiMarkov_fpm_rel_RFC 50 yr _prog`, surv_list$semiMarkov_fpm_rel_RFC_4yr_surv_Prog, surv_list$semiMarkov_fpm_rel_RFC_8yr_surv_Prog, surv_list$semiMarkov_fpm_rel_RFC_15yr_surv_Prog, surv_list$semiMarkov_fpm_rel_RFC_lifetime_surv_Prog)
    ,nrow = 24, byrow = TRUE)


## Make data frame
df <- as.data.frame(var_matrix)

## print df
print(df)

## Export as docx
## Make a flextable and export it as docx 
ft <- flextable(df)  
ft <- autofit(ft)
ft <- set_header_labels(ft, V1 = "Treatment", V2 = "State", V3 = "Source", 
                            V4 = "4", V5 = "8", V6 = "15", V7 = "50*",
                            V8 = "4", V9 = "8", V10 = "15", V11 = "50*")
ft <- add_header_row(ft, values = c("","Lengths of stay in years at times (years)", "Survival proportions at times (years)"),
                     colwidths = c(3, 4, 4))

ft <- align(ft, align = c("right"), i = 1, j = NULL, part = "header")
ft <- align(ft, align = c("center"), i = 2, j = NULL, part = "header")
ft <- align(ft, align = c("left"), i = 1:24, j = 1:3, part = "body")
ft <- align(ft, align = c("right"), i = 1:24, j = 4:11, part = "body")

ft <- bold(ft, i = c(1,5,9,13,17,21), j = 1:11, bold = TRUE, part = "body")

ft <- add_footer_lines(ft,"LE, life expectancy; FC, fludarabine and cyclophosphamide; RFC, rituximab, fludarabine, and cyclophosphamide; SPMs, standard parametric models; FPMs, flexible parametric models; ASF, all-cause survival framework; RSF, relative survival framework; OS, overall survival.")
ft <- autofit(ft)
ft

## Save results
## Prevent from changing the results. We put # here.
# save_as_docx(ft, path="../Output/04b_table_survextrap.docx")
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


