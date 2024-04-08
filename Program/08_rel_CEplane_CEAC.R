## Filename: 08_rel_CEplane_CEAC
## Purpose: Plot cost-effectiveness plane and cost-effectiveness acceptability curve
##          for probabilistic sensitivity analysis
##          Relative survival framework
## Caution: Must have installed the lastest Rtools and R4.2.2+
## Reference: Williams C, Lewsey JD, Briggs AH, Mackay DF. 
##            Cost-effectiveness Analysis in R Using a Multi-state Modeling Survival Analysis Framework: A Tutorial. 
##            Med Decis Making. 2017 May;37(4):340–52.
##            Cost-effectiveness Analysis in R Using a Multi-state Modeling Survival Analysis Framework: A Tutorial © 2017 by Williams C. et al is licensed under CC BY 3.0. 

##############################################################
##============================================================
## Install/Read packages
##============================================================
##############################################################
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(mc2d) # for beta pert distribution
library(pracma)

##############################################################
##============================================================
## Read PSA results 
##============================================================
##############################################################
## Read data
load("../Data/07_rel_PSA.RData")

##############################################################
##============================================================
## Plot of cost-effectiveness plane 
##============================================================
##############################################################
plot <- ggplot() +
        geom_point(data = PSA_results, aes(x = disQALY_inc, y = disCost_inc)) +
        scale_x_continuous(breaks = seq(-2, 2, 1), limits = c(-2, 2), expand = c(0, 0)) +
        scale_y_continuous(breaks = seq(0, 20000, 5000), expand = c(0, 0), limits = c(0, 20000)) +
        labs(x = "Incremental QALY", y = "Incremental Cost (£)") +
        geom_hline(yintercept = 0, linetype = "dashed") +
        geom_vline(xintercept = 0, linetype = "dashed") + 
        geom_abline(intercept = 0, slope = 30000, linetype = "dashed", color = "indianred1") +
        annotate(geom="text", x = 1, y = 18000, label = "CET = £30,000", color = "indianred1", size = 4)+
        theme_bw() +  # Set a white background
        theme(panel.grid = element_blank(), 
              panel.background = element_rect(fill = "white"))
plot

ggsave("../Output/08c_rel_CEplane_microsim.png", plot, width = 10, height = 7, units = "in", dpi = 300)

##############################################################
##============================================================
## Plot of cost-effectiveness plane 
##============================================================
##############################################################
## Create Threshold Vector
## Initialize an empty character vector for T_Names
T_Names <- character(0)

## Specify the maximum threshold (100000 in this case) and the increment (100)
max_threshold <- 100000
increment <- 100

## Create T_Names using a loop
for (i in seq(100, max_threshold, by = increment)) {
  T_Names <- c(T_Names, paste0("T", i))
}

Threshold<-c(seq(100, 100000, by=100))

## Create a data frame
id<-seq(1,1000)
df<-matrix(NA, nrow=1000, ncol=length(T_Names), dimnames=list(NULL, WTP=T_Names))

## For loop for Incremental Net benefits Calculation
for(i in 1:1000){
  df[,i]<-Threshold[i]*PSA_results$disQALY_inc-PSA_results$disCost_inc  
}

## Create a vector to store the values
Prob_CE <- numeric(length(T_Names))

## Calculate and store the values
for (i in 1:length(T_Names)) {
  Prob_CE[i] <- sum(df[, i] >= 0) / 1000
}

## CEAC Dataframe
CEAC_Data<-data.frame(Threshold,Prob_CE)
CEAC_Data

## Plot CEAC
plot <- ggplot(CEAC_Data)+
        geom_line(aes(x=Threshold, y=Prob_CE),color="blue",linewidth=1.0)+
        scale_x_continuous(breaks = seq(0, 100000, 20000), limits = c(0, 100000), 
                           expand = c(0.05, 0.05), labels = scales::comma) +
        ylim(0,1) +
        xlab("Cost-effectiveness threshold  (£)") +
        ylab("Probability of RFC being cost-effective") +
        # ggtitle("The Cost Effectiveness Acceptability Curve") +
        # theme(plot.title = element_text(hjust = 0.5, size = 14))
        theme_bw() +  # Set a white background
        theme(panel.grid = element_blank(), 
              panel.background = element_rect(fill = "white"))
plot 

ggsave("../Output/08c_rel_CEAC_microsim.png", plot, width = 10, height = 7, units = "in", dpi = 300)

##############################################################
##============================================================
## Williams et al. 2017 PSA results
##============================================================
##############################################################
## Run the codes provided by Williams et al.
load("../Data/PSAdata15gom1gam2gom3.RData")

## Preload user-defined functions by Williams et al.
source("williams_function PSAQALY.R")
source("williams_function PSAmeanLY.R")
source("williams_function PSAprob.R")

      ##########################################################################################
      ##========================================================================================
      ## PSA: REMOVING THE DRAWS WITH COMPUTATIONAL DIFFICULTIES 
      ##========================================================================================
      ##########################################################################################
      
      ### identify the error ids for either RFC or FC
      
      errid_<-sort(unique(c(PSASMSARFCwitherrid[[1]],PSASMSAFCwitherrid[[1]])))
      length(errid_) #102 : no. of errors
      temp<-seq(1,1000,1)
      temp2<-temp[-errid_] # position of the 898 usuable draws
      
      temp<-seq(1,1000,1)
      temp<-temp[-PSASMSARFCwitherrid[[1]]] # currently 961 draws for rfc 
      
      ##### remove the appropriate draws from rfc  to obtain the 898 usuable draws 
      temp3<-rep(NA,length(errid_) )
      for (i in 1:length(errid_))
        if(length(which(temp==errid_[i]))>0 ) 
          temp3[i]<-which(temp==errid_[i])
      
      temp3<- temp3[!is.na(temp3)]
      PSASMSARFC<-PSASMSARFCwitherrid[[2]][-temp3]
      
      ##### obtain the usuable draws for fc 
      temp<-seq(1,1000,1)
      temp<-temp[-PSASMSAFCwitherrid[[1]]] 
      
      temp3<-rep(NA,length(errid_) )
      for (i in 1:length(errid_))
        if(length(which(temp==errid_[i]))>0 ) 
          temp3[i]<-which(temp==errid_[i])
      
      temp3<- temp3[!is.na(temp3)]
      PSASMSAFC<-PSASMSAFCwitherrid[[2]][-temp3]
      
      ##########################################################################################
      ##========================================================================================
      ## PSA: MEAN LIFE YEARS 
      ##========================================================================================
      ##########################################################################################
      
      ### RFC treatment arm 
      ### mean life years in initial progression-free state 
      (PSARFCmeanLYPFSdis<-PSAmeanLY(object=PSASMSARFC, state=1,instate=TRUE,discounted=TRUE, rate=0.035))
      
      ### mean life years in progression
      (PSARFCmeanLYprogdis<-PSAmeanLY(object=PSASMSARFC,state=2,instate=TRUE,discounted=TRUE, rate=0.035))
      
      ### FC treatment arm 
      ### mean life years in initial progression-free state 
      (PSAFCmeanLYPFSdis<-PSAmeanLY(object=PSASMSAFC, state=1,instate=TRUE,discounted=TRUE, rate=0.035))
      
      ### mean life years in progression
      (PSAFCmeanLYprogdis<-PSAmeanLY(object=PSASMSAFC, state=2,instate=TRUE,discounted=TRUE, rate=0.035))
      
      
      
      ##########################################################################################
      ##========================================================================================
      ## PSA: MEAN QALYs
      ##========================================================================================
      ##########################################################################################
      
      set.seed(12345)
      utility1=rbeta(1000,800,200)
      set.seed(12345)
      utility2=rbeta(1000,600,400)
      
      ### RFC treatment arm 
      ### mean QALYs in initial progression-free state 
      PSARFCQALYPFSdis<-PSAQALY(object=PSASMSARFC, utility=utility1[1:898], state=1,discounted=TRUE, dis1yronwards=TRUE, rate=0.035)
      mean(PSARFCQALYPFSdis)
      
      ### mean QALYs in progression
      PSARFCQALYprogdis<-PSAQALY(object=PSASMSARFC, utility=utility2[1:898],state=2,discounted=TRUE, rate=0.035)
      mean(PSARFCQALYprogdis)
      
      ### FC treatment arm 
      ### mean QALYs in initial progression-free state 
      PSAFCQALYPFSdis<-PSAQALY(object=PSASMSAFC, utility=utility1[1:898],state=1,discounted=TRUE, rate=0.035)
      mean(PSAFCQALYPFSdis)
      
      ### mean QALYs in progression
      PSAFCQALYprogdis<-PSAQALY(object=PSASMSAFC, utility=utility2[1:898],state=2,discounted=TRUE, rate=0.035)
      mean(PSAFCQALYprogdis)
      
      
      ##########################################################################################
      ##========================================================================================
      ## PSA: INCORPORATING COSTS AND MEASURING COST-EFFECTIVENESS 
      ##========================================================================================
      ##########################################################################################
      
      ### costs probabilistic
      
      nruns<-1000     
      set.seed(12345)   
      cost_PFS_support<- rpert(nruns,14,28,42)
      
      temp1<-84+5179/((PSARFCmeanLYprogdis+PSAFCmeanLYprogdis)/2*12)
      cost_prog_support_dis_<- rpert(nruns,0.5*temp1,temp1,1.5*temp1)
      
      cost_admin_oral <- rpert(nruns, 174, 280, 482)
      cost_admin_complex <- rpert(nruns, 210, 430, 795)
      cost_bmt<-rpert(nruns, 34318.25, 47565.05, 54646.47)
      cost_transfusion<-rpert(nruns, 173.84, 289.73, 405.62)
      cost_blood_unit<-rpert(nruns, 96.67, 161.11, 225.26)
      
      cost_prog_RFC_dis<-cost_prog_support_dis_*12*PSARFCmeanLYprogdis
      cost_prog_FC_dis<-cost_prog_support_dis_*12*PSAFCmeanLYprogdis
      
      cost_PFS_RFC_dis<-10113+1222+2776+1109+21+1109+
        cost_PFS_support*12*PSARFCmeanLYPFSdis+
        cost_bmt*5/402+(cost_transfusion*318+cost_blood_unit*1025)/402
      
      cost_PFS_FC_dis<-2790+1115+22+1115+cost_PFS_support*12*PSAFCmeanLYPFSdis+
        cost_bmt*3/396+(cost_transfusion*269+cost_blood_unit*762)/396
      
      total_cost_RFC_dis<-cost_prog_RFC_dis+cost_PFS_RFC_dis
      total_cost_FC_dis<-cost_prog_FC_dis+cost_PFS_FC_dis
      
      mean_cost_prog_RFC_dis<-mean(cost_prog_RFC_dis)
      mean_cost_prog_FC_dis<-mean(cost_prog_FC_dis)
      
      mean_cost_PFS_RFC_dis<-mean(cost_PFS_RFC_dis)
      mean_cost_PFS_FC_dis<-mean(cost_PFS_FC_dis)
      
      mean_total_cost_RFC_dis<-mean(total_cost_RFC_dis)
      mean_total_cost_FC_dis<-mean(total_cost_FC_dis)
      
      mean_cost_prog_RFC_dis
      mean_cost_prog_FC_dis
      
      mean_cost_PFS_RFC_dis
      mean_cost_PFS_FC_dis
      
      mean_total_cost_RFC_dis
      mean_total_cost_FC_dis
      
      
      ##########################################################################################
      ##========================================================================================
      ## PSA: COST-EFFECTIVENESS PLANE 
      ##========================================================================================
      ##########################################################################################
      
      ### incremental discounted qaly     
      
      incQALY<-PSARFCQALYPFSdis+PSARFCQALYprogdis - 
        PSAFCQALYPFSdis-PSAFCQALYprogdis
      
      ### incremental costs
      incCost<-total_cost_RFC_dis - total_cost_FC_dis    
      
      incCost<-incCost[1:898]
      
      ## Make dataset
      williams2017_PSA <- data.frame(incQALY, incCost)
      
      ##########################################################################################
      ##========================================================================================
      ## PSA: COST-EFFECTIVENESS ACCEPTABILITY CURVE 
      ##========================================================================================
      ##########################################################################################
      # Create Threshold Vector
      # Initialize an empty character vector for T_Names
      T_Names <- character(0)
      
      # Specify the maximum threshold (100000 in this case) and the increment (100)
      max_threshold <- 100000
      increment <- 100
      
      # Create T_Names using a loop
      for (i in seq(100, max_threshold, by = increment)) {
        T_Names <- c(T_Names, paste0("T", i))
      }
      
      Threshold<-c(seq(100, 100000, by=100))
      
      #Create Dataframe
      id<-seq(1,1000)
      df<-matrix(NA, nrow=898, ncol=length(T_Names), dimnames=list(NULL, WTP=T_Names))
      
      #For loop for Incremental Net benefits Calculation
      for(i in 1:1000){
        df[,i]<-Threshold[i]*williams2017_PSA$incQALY-williams2017_PSA$incCost  
      }
      
      ## Create a vector to store the values
      Prob_CE_williams <- numeric(length(T_Names))
      
      # Calculate and store the values
      for (i in 1:length(T_Names)) {
        Prob_CE_williams[i] <- sum(df[, i] >= 0) / 1000
      }
      
      #CEAC Dataframe
      williams2017_PSA_CEAC <-data.frame(Threshold, Prob_CE_williams)

##############################################################
##============================================================
## Comparison: microsim vs. Williams et al. 2017
##============================================================
##############################################################
## Plot of cost-effectiveness plane for both results
plot1 <- ggplot() +
           geom_point(data = PSA_results, aes(x = disQALY_inc, y = disCost_inc, color="Microsimulation"), alpha = 0.8, size = 2) +
           geom_point(data = williams2017_PSA, aes(x = incQALY, y = incCost, color="Semi-Markov"), alpha = 0.8, size = 2) +
           scale_x_continuous(breaks = seq(-2, 2, 1), limits = c(-2, 2), expand = c(0, 0)) +
           scale_y_continuous(breaks = seq(0, 20000, 5000), expand = c(0, 0), limits = c(0, 20000)) +
           scale_color_manual(values = c("Semi-Markov" = "orange", "Microsimulation" = "cornflowerblue")) +
           labs(x = "Incremental QALY", y = "Incremental Cost (£)") +
           labs(color="") +
           geom_hline(yintercept = 0, linetype = "dashed") +
           geom_vline(xintercept = 0, linetype = "dashed") + 
           geom_abline(intercept = 0, slope = 30000, linetype = "dashed", color = "indianred1") +
           annotate(geom="text", x = 1, y = 18000, label = "CET = £30,000", color = "indianred1", size = 4)+
           theme_bw() +  # Set a white background
           theme(panel.grid = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  axis.text = element_text(size = 12),
                  axis.title.x = element_text(size = 12),
                  axis.title.y = element_text(size = 12),
                  legend.text = element_text(size = 12),
                  legend.key.size = unit(0.4, 'cm'),
                  legend.position="bottom", 
                  legend.box="vertical",
                  legend.key = element_rect(colour = NA, fill = NA)) +
        ggtitle("Cost-effectiveness plane")
      
plot1
ggsave("../Output/08c_rel_CEplane_comparison.png", plot1, width = 10, height = 7, units = "in", dpi = 300)
      
### Plot of cost-effectiveness acceptability curves 
plot2 <- ggplot()+
              geom_line(data = CEAC_Data, aes(x = Threshold, y = Prob_CE, color = "Microsimulation"), linewidth = 1.0) +
              geom_line(data = williams2017_PSA_CEAC, aes(x = Threshold, y = Prob_CE_williams, color = "Semi-Markov"), linewidth = 1.0) +
              scale_x_continuous(breaks = seq(0, 100000, 20000), limits = c(0, 100000), 
                                 expand = c(0.05, 0.05), labels = scales::comma) +
              scale_color_manual(values = c("Microsimulation" = "cornflowerblue", "Semi-Markov" = "orange")) +
              ylim(0, 1) +
              xlab("Cost-effectiveness threshold (£)") +
              ylab("Probability of RFC being cost-effective") +
              labs(color = "") +
              guides(
                color = guide_legend(
                  keywidth = 2,  # Adjust the width of the legend key
                  keyheight = 1,  # Adjust the height of the legend key
                  override.aes = list(
                    linetype = c("solid", "solid")
                  )
                )
              ) +
              theme_bw() + 
              theme(panel.grid = element_blank(),
                    plot.title = element_text(hjust = 0.5),
                    axis.text = element_text(size = 12),
                    axis.title.x = element_text(size = 12),
                    axis.title.y = element_text(size = 12),
                    legend.text = element_text(size = 12),
                    legend.key.size = unit(0.4, 'cm'),
                    legend.position="bottom", 
                    legend.box="vertical",
                    legend.key = element_rect(colour = NA, fill = NA)) +
              ggtitle("Cost-effectiveness acceptability curve")

plot2 
      
ggsave("../Output/08c_rel_CEAC_comparison.png", plot2, width = 10, height = 7, units = "in", dpi = 300)

## Plot two figures together
plot <- ggarrange(plot1, plot2, ncol = 2,
                  common.legend = FALSE,
                  align = "h") 
plot
ggsave("../Output/08c_rel_CEplane_CEAC.png", plot, width = 10, height = 7, units = "in", dpi = 300)
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

