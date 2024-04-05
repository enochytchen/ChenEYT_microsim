## Filename: 04b_table_survextrap
## Purpose: Caluclate area under curve for selected time points
##          Observed: Hallek 2010, Fischer 2016
##          Extrapolated: Williams 2017, Microsimulation
## Reference: Williams C, Lewsey JD, Briggs AH, Mackay DF. 
##            Cost-effectiveness Analysis in R Using a Multi-state Modeling Survival Analysis Framework: A Tutorial. 
##            Med Decis Making. 2017 May;37(4):340–52.
##            Cost-effectiveness Analysis in R Using a Multi-state Modeling Survival Analysis Framework: A Tutorial © 2017 by Williams C. et al is licensed under CC BY 3.0. 

## Set seed for coherency
set.seed(12345)

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
##############################################################
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

###################
## Williams 2017 ##
###################
williams2017 <- read.table("../Data/01a_williams2017.txt", header = TRUE)
colnames(williams2017)

williams_FC_4yr_surv <- subset(williams2017, treat == 0 & time == 4)$surv
williams_RFC_4yr_surv <- subset(williams2017, treat == 1 & time == 4)$surv

williams_FC_8yr_surv <- subset(williams2017, treat == 0 & time == 8)$surv
williams_RFC_8yr_surv <- subset(williams2017, treat == 1 & time == 8)$surv

williams_FC_15yr_surv <- subset(williams2017, treat == 0 & time == 15)$surv
williams_RFC_15yr_surv <- subset(williams2017, treat == 1 & time == 15)$surv

#####################
## Microsimulation ##
#####################
## Run the microsimulation model to obtain the results
## source("04a_rel_survextrap.R")
## Read the results
results_FC <- readRDS("../Data/04a_rel_survextrap_results_FC.rds")
results_RFC <- readRDS("../Data/04a_rel_survextrap_results_RFC.rds")

## Obtain survival at times 4, 8, and 15 years
prev_FC <- as.data.frame(results_FC$prev)
prev_FC$prob <- prev_FC$number / results_FC$n
prev_FC$time <- round(prev_FC$time, digits = 1)

prev_RFC <- as.data.frame(results_RFC$prev)
prev_RFC$prob <- prev_RFC$number / results_RFC$n
prev_RFC$time <- round(prev_RFC$time, digits = 1)

microsim_FC_4yr_surv_PFS <- prev_FC$prob[prev_FC$time == 4 & prev_FC$state == "PFS"]
microsim_FC_4yr_surv_Prog <- prev_FC$prob[prev_FC$time == 4 & prev_FC$state == "Prog"]
microsim_FC_4yr_surv <- microsim_FC_4yr_surv_PFS + microsim_FC_4yr_surv_Prog

microsim_RFC_4yr_surv_PFS <- prev_RFC$prob[prev_RFC$time == 4 & prev_RFC$state == "PFS"]
microsim_RFC_4yr_surv_Prog <- prev_RFC$prob[prev_RFC$time == 4 & prev_RFC$state == "Prog"]
microsim_RFC_4yr_surv <- microsim_RFC_4yr_surv_PFS + microsim_RFC_4yr_surv_Prog

microsim_FC_8yr_surv_PFS <- prev_FC$prob[prev_FC$time == 8 & prev_FC$state == "PFS"]
microsim_FC_8yr_surv_Prog <- prev_FC$prob[prev_FC$time == 8 & prev_FC$state == "Prog"]
microsim_FC_8yr_surv <- microsim_FC_8yr_surv_PFS + microsim_FC_8yr_surv_Prog

microsim_RFC_8yr_surv_PFS <- prev_RFC$prob[prev_RFC$time == 8 & prev_RFC$state == "PFS"]
microsim_RFC_8yr_surv_Prog <- prev_RFC$prob[prev_RFC$time == 8 & prev_RFC$state == "Prog"]
microsim_RFC_8yr_surv <- microsim_RFC_8yr_surv_PFS + microsim_RFC_8yr_surv_Prog

microsim_FC_15yr_surv_PFS <- prev_FC$prob[prev_FC$time == 15 & prev_FC$state == "PFS"]
microsim_FC_15yr_surv_Prog <- prev_FC$prob[prev_FC$time == 15 & prev_FC$state == "Prog"]
microsim_FC_15yr_surv <- microsim_FC_15yr_surv_PFS + microsim_FC_15yr_surv_Prog

microsim_RFC_15yr_surv_PFS <- prev_RFC$prob[prev_RFC$time == 15 & prev_RFC$state == "PFS"]
microsim_RFC_15yr_surv_Prog <- prev_RFC$prob[prev_RFC$time == 15 & prev_RFC$state == "Prog"]
microsim_RFC_15yr_surv <- microsim_RFC_15yr_surv_PFS + microsim_RFC_15yr_surv_Prog

## Store the results in one list
variables <- c("fischer_FC_4yr_surv", "fischer_RFC_4yr_surv", "fischer_FC_8yr_surv", "fischer_RFC_8yr_surv",
               "williams_FC_4yr_surv", "williams_RFC_4yr_surv", "williams_FC_8yr_surv", "williams_RFC_8yr_surv", "williams_FC_15yr_surv", "williams_RFC_15yr_surv",
               "microsim_FC_4yr_surv", "microsim_RFC_4yr_surv", "microsim_FC_8yr_surv", "microsim_RFC_8yr_surv", "microsim_FC_15yr_surv", "microsim_RFC_15yr_surv")

surv_list <- mget(variables)

surv_list <- lapply(surv_list, function(x) {
  formatted <- format(round(x, digits = 2), nsmall = 2)
  return(formatted)
})

##############################################################
##============================================================
## Calculate restricted mean survival time (RMST)
##============================================================
##############################################################
var_list <- list()

##################
## Fischer 2016 ##
##################
## 4 years (Hallek 2010's max follow-up was 4 years)
rmst_mfit2 <- rmst2(fischer2016$time_y, fischer2016$event, arm = fischer2016$treat, tau = 4)
var_list[["fischer_RFC_4yr"]] <- rmst_mfit2$RMST.arm1$rmst[[1]]
print(paste("Fischer 2016,", "RFC,", "AUC for time <= 4 years:", var_list[["fischer_RFC_4yr"]]))
var_list[["fischer_FC_4yr"]] <- rmst_mfit2$RMST.arm0$rmst[[1]]
print(paste("Fischer 2016,", "FC,", "AUC for time <= 4 years:", var_list[["fischer_FC_4yr"]]))

## 8 years
rmst_mfit3 <- rmst2(fischer2016$time_y, fischer2016$event, arm = fischer2016$treat, tau = 8)
var_list[["fischer_RFC_8yr"]] <- rmst_mfit3$RMST.arm1$rmst[[1]]
print(paste("Fischer 2016,", "RFC,", "AUC for time <= 3.94 years:", var_list[["fischer_RFC_8yr"]]))
var_list[["fischer_FC_8yr"]] <- rmst_mfit3$RMST.arm0$rmst[[1]]
print(paste("Fischer 2016,", "FC,", "AUC for time <= 3.94 years:", var_list[["fischer_FC_8yr"]]))

###################
## Williams 2017 ##
###################
## Define the time intervals
time_intervals <- c(4, 8, 15)

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
  print(paste("Williams 2017", "FC", "Area under curve for time <= ", i, "and treat = RFC: ", value))
  
  ## calculate area under the curve for 'surv' with time <= i for FC
  value <- trapz(williams2017_subset_FC$time, williams2017_subset_FC$surv)
  var <- paste("williams_FC_",i,"yr")
  var_list[[var]] <- value
  print(paste("Williams 2017", "FC", "Area under curve for time <= ", i, "and treat = FC: ", value))
}

#####################
## Microsimulation ##
#####################
## Define the time intervals
time_intervals <- c(4, 8, 15)

## loop over the time intervals
for (i in time_intervals) {
    report <- lapply(list(results_FC, results_RFC), function(sim){
      # create a subset of the data where 'time' is less than or equal to i
      sim$pt_subset <- subset(sim$pt, time <= i)
      
      # calculate LY using the subset and print it
      LY <- (sum(sim$pt_subset$pt) - sum(subset(sim$pt_subset, state %in% c("ExpD", "ExcD"))$pt)) / sim$n
})

## FC, treat = 0
value <- report[[1]]
var <- paste("microsim_FC_",i,"yr")
var_list[[var]] <- value
print(paste("Microsimulation", "FC", "Area under curve for time <=", i, ": ", value))

## RFC, treat = 1
value <- report[[2]]
var <- paste("microsim_RFC_",i,"yr")
var_list[[var]] <- value
print(paste("Microsimulation", "RFC", "Area under curve for time <=", i, ": ", value))
}

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
  c("FC", "Observed", var_list$fischer_FC_4yr, var_list$fischer_FC_8yr, NA, surv_list$fischer_FC_4yr_surv, surv_list$fischer_FC_8yr_surv, NA, 
    NA, "Semi-Markov", var_list$`williams_FC_ 4 yr`, var_list$`williams_FC_ 8 yr`, var_list$`williams_FC_ 15 yr`, surv_list$williams_FC_4yr_surv, surv_list$williams_FC_8yr_surv, surv_list$williams_FC_15yr_surv,
    NA, "Microsimulation", var_list$`microsim_FC_ 4 yr`, var_list$`microsim_FC_ 8 yr`, var_list$`microsim_FC_ 15 yr`, surv_list$microsim_FC_4yr_surv, surv_list$microsim_FC_8yr_surv, surv_list$microsim_FC_15yr_surv,
    "RFC", "Observed", var_list$fischer_RFC_4yr, var_list$fischer_RFC_8yr, NA, surv_list$fischer_RFC_4yr_surv, surv_list$fischer_RFC_8yr_surv, NA, 
    NA, "Semi-Markov", var_list$`williams_RFC_ 4 yr`, var_list$`williams_RFC_ 8 yr`, var_list$`williams_RFC_ 15 yr`, surv_list$williams_RFC_4yr_surv, surv_list$williams_RFC_8yr_surv, surv_list$williams_RFC_15yr_surv,
    NA, "Microsimulation", var_list$`microsim_RFC_ 4 yr`, var_list$`microsim_RFC_ 8 yr`, var_list$`microsim_RFC_ 15 yr`, surv_list$microsim_RFC_4yr_surv, surv_list$microsim_RFC_8yr_surv, surv_list$microsim_RFC_15yr_surv),
    nrow = 6, byrow = TRUE)


## Make data frame
df <- as.data.frame(var_matrix)

## print df
print(df)

## Export as docx
## Make a flextable and export it as docx 
ft <- flextable(df)  
ft <- autofit(ft)
ft <- set_header_labels(ft, V1 = "Treatment", V2 = "Source", 
                            V3 = "4", V4 = "8", V5 = "15",
                            V6 = "4", V7 = "8", V8 = "15")
ft <- add_header_row(ft, values = c("","RMST at times (years)", "Survival proportions at times (years)"),
                     colwidths = c(1,4, 3))

ft <- align(ft, align = c("right"), i = 1, j = NULL, part = "header")
ft <- align(ft, align = c("center"), i = 2, j = NULL, part = "header")

ft <- add_footer_lines(ft,"RMST, restricted mean survival time; FC, fludarabine and cyclophosphamide; RFC, rituximab, fludarabine, and cyclophosphamide.")
ft

save_as_docx(ft, path="../Output/04b_table_survextrap.docx")
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


