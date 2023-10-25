## Filename:  01a_williams2017
## Purpose:   read Williams et al. 2017's results for saving the extrapolation results
## Notes:     The codes were modified from Williams et al. 2017's supplementary materials
## Reference: Williams C, Lewsey JD, Briggs AH, Mackay DF. 
##            Cost-effectiveness Analysis in R Using a Multi-state Modeling Survival Analysis Framework: A Tutorial. 
##            Med Decis Making. 2017 May;37(4):340–52.
##            Cost-effectiveness Analysis in R Using a Multi-state Modeling Survival Analysis Framework: A Tutorial © 2017 by Williams C. et al is licensed under CC BY 3.0. 

## Set the wd as where this R file is
if (require("rstudioapi") && isAvailable()) {
  original_wd <- getwd()  # Store the original working directory
  wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(wd)
}
setwd("../Data")
getwd()

##############################################################
##============================================================
## Install/Read packages
##============================================================
##############################################################
library(tidyverse)
library(ggplot2)

##############################################################
##============================================================
## Read data
##============================================================
##############################################################
## Load the base case analysis extrapolation results provided by Williams et al.
load("basecase15gom1gam2gom3.RData")

## Select the one chosen by Wiiliams et al. 2017 (Appendix1)
## Progression -> death: Gompertz
## Progression-free -> death: Generalized gamma
## Progression–free -> progression: Gompertz
## The basecase analysis was gom1gam2gom3
basecase_FC <-  gom1gam2gom3simFC
basecase_RFC <- gom1gam2gom3simRFC

## Create a variable for undiscounted survival
## that is, probability of being at PFS (pstate1) + Progression (pstate2)
basecase_FC$surv <- basecase_FC$pstate1 + basecase_FC$pstate2
basecase_RFC$surv <- basecase_RFC$pstate1 + basecase_RFC$pstate2

## Create a variable to indicate treatments
basecase_FC$treat <- 0
basecase_RFC$treat <- 1

## Combine data for 
basecase <- bind_rows(basecase_FC, basecase_RFC)

## Plot survival curve
ggplot(basecase, aes(x = time, y = surv, color=as.factor(treat))) +
       geom_line()

##############################################################
##============================================================
## Save extrapolation results
##============================================================
##############################################################
write.table(basecase, file = "01a_williams2017.txt", sep = " ",
            row.names = FALSE, col.names = TRUE,  quote = FALSE)

################################################################
setwd(original_wd)  # Reset to the original working directory
################################################################
## A microsimulation model incorporating relative survival extrapolation and multiple timescales for health technology assessment © 2023 by Chen EYT, Dickman PW, Clements MS is licensed under CC BY 4.0
