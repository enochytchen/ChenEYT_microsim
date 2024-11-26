## Filename: 01a_fischer2016_PFS
## Purpose: Organise digitised data from Fischer et al., 2016
##          8 years of follow-up from CLL-8 trial
## Reference: Fischer K, Bahlo J, Fink AM, Goede V, Herling CD, Cramer P, et al. 
##            Long-term remissions after FCR chemoimmunotherapy in previously untreated patients with CLL: updated results of the CLL8 trial.
##            Blood. 2016 Jan 14;127(2):208–15.

## Note: Start from making digitised data
##       the difference between PFS and OS is actually Prog :))


##############################################################
##============================================================
## Read packages
##============================================================
##############################################################
library(survHE)
library(survival)
library(dplyr)

##############################################################
##============================================================
## Read and manage digitised data
##============================================================
##############################################################
## Digitised data from Fischer et al., Figure 1B
surv_RFC <- read.csv("../Data/Fischer2016_Fig1B_RFC.csv", header=FALSE)
surv_FC <- read.csv("../Data/Fischer2016_Fig1B_FC.csv", header=FALSE)

##############################################################
##============================================================
## Manage data
##============================================================
##############################################################
## Rename column names
colnames(surv_RFC) <- c("time", "surv")
colnames(surv_FC) <- c("time", "surv")

## Sort the data by time
surv_RFC <- surv_RFC[order(surv_RFC$time),]
surv_FC <- surv_FC[order(surv_FC$time),]

## Add row id 
surv_RFC$id <- seq_along(surv_RFC[,1])
surv_FC$id <- seq_along(surv_FC[,1])

## Order variables
surv_RFC <- surv_RFC[, c("id", "time", "surv")]
surv_FC <- surv_FC[, c("id", "time", "surv")]

## Check whether the survival is monotonically decreasing
## Print TRUE if the variable is monotonically decreasing
is_decreasing <- all(diff(surv_RFC$surv) <= 0)
print(is_decreasing)  
is_decreasing <- all(diff(surv_FC$surv) <= 0)
print(is_decreasing)  

## Save the data
## Prevent from changing the results. We put # here.
# write.table(surv_RFC, "../Data/01a_fischer2016_surv_RFC.txt", row.names = FALSE)
# write.table(surv_FC, "../Data/01a_fischer2016_surv_FC.txt", row.names = FALSE)

##############################################################
##============================================================
## Numbers at risk over time
##============================================================
##############################################################
#### RFC
## Totally 9 time points
i <- c(1:9)

## Number at risk at times
trisk <- c(seq(0, 96, by=12))
surv_RFC$time
lower <- c(1,11,22,31,40,49,60,74,79)
upper <- lower - 1
upper <- append(upper[-1],79)
print(upper)

nrisk <- c(408, 384, 363, 342, 318, 290, 134, 41, 2)

nrisk_RFC <- data.frame(i, trisk, lower, upper, nrisk)  

#### FC
## Totally 9 time points
i <- c(1:9)

## Number at risk at times
trisk <- c(seq(0,96,by=12))
surv_FC$time
lower <- c(1,10,20,29,39,49,59,70,76)
upper <- lower - 1
upper <- append(upper[-1],76)
print(upper)

nrisk<- c(409, 360, 232, 297, 262, 220, 100, 33, 1)
nrisk_FC <- data.frame(i, trisk, lower, upper, nrisk)  

## Save the data 
## Prevent from changing the results. We put # here.
# write.table(nrisk_RFC, "../Data/01a_fischer2016_nrisk_RFC.txt", row.names = FALSE)
# write.table(nrisk_FC, "../Data/01a_fischer2016_nrisk_FC.txt", row.names = FALSE)

##############################################################
##============================================================
## Generate digitised data
##============================================================
##############################################################
filenames <- c("RFC", "FC")
for (filename in filenames) {
    surv <- paste("../Data/01a_fischer2016_surv_", filename, ".txt", sep="")
    nrisk <- paste("../Data/01a_fischer2016_nrisk_", filename, ".txt", sep="")
    km <- paste("../Data/01a_fischer2016_km_", filename, ".txt", sep="")
    ipd <- paste("../Data/01a_fischer2016_ipd_", filename,".txt", sep="")
  
  digitise(surv_inp = surv, nrisk = nrisk, km_output = km, ipd_output = ipd)
  
  make.ipd(ipd, ctr = 1, var.labs = c("time", "event"))
}

##############################################################
##============================================================
## Combine data into one dataset
##============================================================
##############################################################
## Read data
ipd_RFC <- read.table("../Data/01a_fischer2016_ipd_RFC.txt", header = TRUE)
ipd_RFC$arm <- 1
ipd_FC <- read.table("../Data/01a_fischer2016_ipd_FC.txt", header = TRUE)
ipd_FC$arm <- 0

## Combine 2 datasets into one data
fischer2016 <- rbind(ipd_RFC, ipd_FC)
fischer2016 <- rename(fischer2016, treat = arm)

## Plot K-M curves
mfit <- survfit(Surv(time, event==1) ~ treat, data = fischer2016)

## Print Kaplan-Meier table
summary(mfit)  

## Plot Kaplan-Meier curve
plot(mfit,                                                     
     ylab="S(t)",
     xlab="Time (months)",
     main = "Kaplan-Meier estimates")

## Save the data
## Prevent from changing the results. We put # here.
# write.table(fischer2016, file = "../Data/01a_fischer2016.txt", sep = " ",
#             row.names = FALSE, col.names = TRUE,  quote = FALSE)
################################################################
# Chen EYT, Dickman PW, Clements MS. 
# A Multistate Model Incorporating Relative Survival Extrapolation and 
# Mixed Time Scales for Health Technology Assessment. 
# PharmacoEconomics. 2024. (doi:https://doi.org/10.1007/s40273-024-01457-w).
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
