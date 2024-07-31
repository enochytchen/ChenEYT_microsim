## Filename: ../Data/01a_fischer2016_PFS
## Purpose: Organise progression-free survival digitised data from Fischer et al., 2016
##          8 years of follow-up from CLL-8 trial
## Reference: Fischer K, Bahlo J, Fink AM, Goede V, Herling CD, Cramer P, et al. 
##            Long-term remissions after FCR chemoimmunotherapy in previously untreated patients with CLL: updated results of the CLL8 trial.
##            Blood. 2016 Jan 14;127(2):208–15.

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
## Digitised data from Fischer et al., Figure 1A
pfs_RFC <- read.csv("../Data/Fischer2016_Fig1A_RFC.csv", header=FALSE)
pfs_FC <- read.csv("../Data/Fischer2016_Fig1A_FC.csv", header=FALSE)

# Add a row for calculation later
pfs_FC <- rbind(pfs_FC, data.frame(V1 = 96, V2 = tail(pfs_FC$V2, 1)))

##############################################################
##============================================================
## Manage data
##============================================================
##############################################################
## Rename column names
colnames(pfs_RFC) <- c("time", "surv")
colnames(pfs_FC) <- c("time", "surv")

## Sort the data by time
pfs_RFC <- pfs_RFC[order(pfs_RFC$time),]
pfs_FC <- pfs_FC[order(pfs_FC$time),]

## Add row id 
pfs_RFC$id <- seq_along(pfs_RFC[,1])
pfs_FC$id <- seq_along(pfs_FC[,1])

## Order variables
pfs_RFC <- pfs_RFC[, c("id", "time", "surv")]
pfs_FC <- pfs_FC[, c("id", "time", "surv")]

## Check whether the survival is monotonically decreasing
## Print TRUE if the variable is monotonically decreasing
is_decreasing <- all(diff(pfs_RFC$surv) <= 0)
print(is_decreasing)  
is_decreasing <- all(diff(pfs_FC$surv) <= 0)
print(is_decreasing)  

## Save the data
## Prevent from changing the results. We put # here.
write.table(pfs_RFC, "../Data/01a_fischer2016_pfs_RFC.txt", row.names = FALSE)
write.table(pfs_FC, "../Data/01a_fischer2016_pfs_FC.txt", row.names = FALSE)

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
pfs_RFC$time
lower <- c(1,11,21,31,41,51,60,71,72)
upper <- lower - 1
upper <- append(upper[-1],72)
print(upper)

nrisk <- c(408, 358, 310, 261, 222, 178, 76, 18, 1)

nrisk_RFC <- data.frame(i, trisk, lower, upper, nrisk)  
nrisk_RFC
pfs_RFC$time

#### FC
## Totally 9 time points
i <- c(1:9)

## Number at risk at times
trisk <- c(seq(0,96,by=12))
pfs_FC$time
lower <- c(1,13,27,40,49,58,67,70,72)
upper <- lower - 1
upper <- append(upper[-1],72)
print(upper)

nrisk<- c(409, 232, 236, 167, 119, 86, 39, 13, 0)
nrisk_FC <- data.frame(i, trisk, lower, upper, nrisk)  
nrisk_FC
pfs_FC$time

## Save the data 
## Prevent from changing the results. We put # here.
write.table(nrisk_RFC, "../Data/01a_fischer2016_pfs_nrisk_RFC.txt", row.names = FALSE)
write.table(nrisk_FC, "../Data/01a_fischer2016_pfs_nrisk_FC.txt", row.names = FALSE)

##############################################################
##============================================================
## Generate digitised data
##============================================================
##############################################################
filenames <- c("RFC", "FC")
for (filename in filenames) {
    surv <- paste("../Data/01a_fischer2016_pfs_", filename, ".txt", sep="")
    nrisk <- paste("../Data/01a_fischer2016_pfs_nrisk_", filename, ".txt", sep="")
    km <- paste("../Data/01a_fischer2016_pfs_km_", filename, ".txt", sep="")
    ipd <- paste("../Data/01a_fischer2016_pfs_ipd_", filename,".txt", sep="")
  
  digitise(surv_inp = surv, nrisk = nrisk, km_output = km, ipd_output = ipd)
  
  make.ipd(ipd, ctr = 1, var.labs = c("time", "event"))
}

##############################################################
##============================================================
## Combine data into one dataset
##============================================================
##############################################################
## Read data
ipd_RFC <- read.table("../Data/01a_fischer2016_pfs_ipd_RFC.txt", header = TRUE)
ipd_RFC$arm <- 1
ipd_FC <- read.table("../Data/01a_fischer2016_pfs_ipd_FC.txt", header = TRUE)
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
write.table(fischer2016, file = "../Data/01a_fischer2016_pfs.txt", sep = " ",
            row.names = FALSE, col.names = TRUE,  quote = FALSE)
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
