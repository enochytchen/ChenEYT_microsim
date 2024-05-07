## Filename: 04c_figure_survextrap
## Purpose: Plot the observed and the extrapolated survival functions
##          Observed: Hallek 2010, Fischer 2016
##          Extrapolated: Semi-Markov with stand parametric models within an all-cause survival framework (Williams et al.)  
##                        Semi-Markov with flexible parametric models within an all-cause survival framework
##                        Semi-Markov with flexible parametric models within a relative survival framework 
##                        Markov with stand parametric models within an all-cause survival framework 
##                        Markov with flexible parametric models within an all-cause survival framework
##                        Markov with flexible parametric models within a relative survival framework 
## Reference: Williams C, Lewsey JD, Briggs AH, Mackay DF. 
##            Cost-effectiveness Analysis in R Using a Multi-state Modeling Survival Analysis Framework: A Tutorial. 
##            Med Decis Making. 2017 May;37(4):340–52.
##            Cost-effectiveness Analysis in R Using a Multi-state Modeling Survival Analysis Framework: A Tutorial © 2017 by Williams C. et al is licensed under CC BY 3.0. 

##############################################################
##============================================================
## Read packages
##============================================================
##############################################################
library(survival)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggh4x)

##############################################################
##============================================================
## Read Hallek2010 and Fischer2016 as comparison
##============================================================
##############################################################
#################
## Hallek 2010 ##
#################
data<-read.table("../Data/data.txt",sep="\t",skip=1) 
names(data)=c("patid","treat","prog","death", "prog_t","death_t",
              "progdeath", "progdeath_t")
data$death_ty <- data$death_t/12
data <- data %>%
        mutate(death = ifelse(death_t>=47, 0 , death))

## Fit to plot K-M curve
mfit1 <- survfit(Surv(death_ty, death==1) ~ treat, data = subset(data, death_t <= 47))
plot(mfit1, conf.int = FALSE, mark.time = FALSE, lty = 1, col = c("red","blue"),
     xlim = c(0, 8), ylim=c(0,1))

## Make a data frame
mfit1_sum <- summary(mfit1)
mfit1_data <- data.frame(time = mfit1_sum$time,
                         surv = mfit1_sum$surv,
                         strata = mfit1_sum$strata)

## Subset the survival data for treat == 0 and 1
mfit1_data_FC <- mfit1_data[grepl("treat=0", mfit1_data$strata), ]
mfit1_data_RFC <- mfit1_data[grepl("treat=1", mfit1_data$strata), ]

## Add a time point for censoring, for the purpose of plotting
mfit1_data_FC[nrow(mfit1_data_FC)+1,] <- mfit1_data_FC[nrow(mfit1_data_FC),]
mfit1_data_FC[nrow(mfit1_data_FC), "time"] <- 4  
mfit1_data_RFC[nrow(mfit1_data_RFC)+1,] <- mfit1_data_RFC[nrow(mfit1_data_RFC),]
mfit1_data_RFC[nrow(mfit1_data_RFC), "time"] <- 4  

##################
## Fischer 2016 ##
##################
fischer2016 <- read.table("../Data/01a_fischer2016.txt", header = TRUE)
fischer2016$time_y <- fischer2016$time/12

## Fit to plot K-M curve
mfit2 <- survfit(Surv(time_y, event==1) ~ treat, data = fischer2016)
plot(mfit2, conf.int = FALSE, mark.time = FALSE, lty = 1, col = c("red","blue"),
     xlim = c(0, 8), ylim=c(0,1))

## Make a data frame
mfit2_sum <- summary(mfit2)
mfit2_data <- data.frame(time = mfit2_sum$time,
                         surv = mfit2_sum$surv,
                         strata = mfit2_sum$strata)

## Subset the survival data for treat == 0 and 1
mfit2_data_FC <- mfit2_data[grepl("treat=0", mfit2_data$strata), ]
mfit2_data_RFC <- mfit2_data[grepl("treat=1", mfit2_data$strata), ]

## Add a time point for censoring, for the purpose of plotting
mfit2_data_FC[nrow(mfit2_data_FC)+1,] <- mfit2_data_FC[nrow(mfit2_data_FC),]
mfit2_data_FC[nrow(mfit2_data_FC), "time"] <- 8  
mfit2_data_RFC[nrow(mfit2_data_RFC)+1,] <- mfit2_data_RFC[nrow(mfit2_data_RFC),]
mfit2_data_RFC[nrow(mfit2_data_RFC), "time"] <- 8  

###############################################################################################
## Semi-Markov with standard parametric models, all-cause survival framework (Williams et al.) ##
###############################################################################################
semiMarkov_williams_ac_FC <- readRDS("../Data/04a7_semiMarkov_williams_ac_microsim_FC.rds")
semiMarkov_williams_ac_RFC <- readRDS("../Data/04a7_semiMarkov_williams_ac_microsim_RFC.rds")

###############################################################################
## Semi-Markov with flexible parametric models, all-cause survival framework ##
###############################################################################
semiMarkov_fpm_ac_FC <- readRDS("../Data/04a5_semiMarkov_fpm_ac_microsim_FC.rds")
semiMarkov_fpm_ac_RFC <- readRDS("../Data/04a5_semiMarkov_fpm_ac_microsim_RFC.rds")

###################################################################################
## Semi-Markov with flexible parametric models, relative survival framework ##
###################################################################################
semiMarkov_fpm_rel_FC <- readRDS("../Data/04a6_semiMarkov_fpm_rel_microsim_FC.rds")
semiMarkov_fpm_rel_RFC <- readRDS("../Data/04a6_semiMarkov_fpm_rel_microsim_RFC.rds")

##########################################################################
## Markov with standard parametric models, all-cause survival framework ##
##########################################################################
Markov_williams_ac_FC <- readRDS("../Data/04a3_Markov_williams_ac_microsim_FC.rds")
Markov_williams_ac_RFC <- readRDS("../Data/04a3_Markov_williams_ac_microsim_RFC.rds")

##########################################################################
## Markov with flexible parametric models, all-cause survival framework ##
##########################################################################
Markov_fpm_ac_FC <- readRDS("../Data/04a1_Markov_fpm_ac_microsim_FC.rds")
Markov_fpm_ac_RFC <- readRDS("../Data/04a1_Markov_fpm_ac_microsim_RFC.rds")

#########################################################################
## Markov with flexible parametric models, relative survival framework ##
#########################################################################
Markov_fpm_rel_FC <- readRDS("../Data/04a2_Markov_fpm_rel_microsim_FC.rds")
Markov_fpm_rel_RFC <- readRDS("../Data/04a2_Markov_fpm_rel_microsim_RFC.rds")

##############################################################
##============================================================
## Calculate overall survival of all the microsim models
##============================================================
##############################################################
# List of data frames
names <- c("semiMarkov_williams_ac_RFC", "semiMarkov_williams_ac_FC", 
           "semiMarkov_fpm_ac_RFC", "semiMarkov_fpm_ac_FC", 
           "semiMarkov_fpm_rel_RFC", "semiMarkov_fpm_rel_FC",
           "Markov_williams_ac_RFC", "Markov_williams_ac_FC", 
           "Markov_fpm_ac_RFC", "Markov_fpm_ac_FC", 
           "Markov_fpm_rel_RFC", "Markov_fpm_rel_FC")

# Loop over each data frame
for (name in names) {
  # Claim data frame
  prev_microsim <- as.data.frame(get(paste0(name))$prev)
  
  # Calculate the probabilities
  prev_microsim$prob <- prev_microsim$number / get(paste0(name))$n
  
  # Calculate overall survival
  OS_microsim <- prev_microsim %>%
                 filter(state %in% c("PFS", "Prog")) %>%  
                 group_by(time) %>%                       
                 summarise(prob = sum(prob))
  
  # Optionally, you can assign the result to a new variable or do something else with it
  assign(paste0("OS_", name), OS_microsim)
}

##############################################################
##============================================================
## Plot the observed and the extrapolated results
##============================================================
##############################################################
# Plot
plot_RFC <- ggplot() +
            # Semi-Markov: SPMs, ASF (Williams2017)
            geom_line(data = OS_semiMarkov_williams_ac_RFC, aes(x = time, y = prob*100, color = "l3", linetype = "l3"), inherit.aes = FALSE,
                      linewidth = 1.2, alpha = 1)  +
            # Semi-Markov: FPMs, ASF
            geom_line(data = OS_semiMarkov_fpm_ac_RFC, aes(x = time, y = prob*100, color = "l4", linetype = "l4"), inherit.aes = FALSE,
                      linewidth = 1.2, alpha = 1)  +
            # Semi-Markov: FPMs, RSF (Proposed method)
            geom_line(data = OS_semiMarkov_fpm_rel_RFC, aes(x = time, y = prob*100, color = "l5", linetype = "l5"), inherit.aes = FALSE,
                      linewidth = 1.2, alpha = 0.8) +
            # Markov: SPMs, ASF
            geom_line(data = OS_Markov_williams_ac_RFC, aes(x = time, y = prob*100, color = "l6", linetype = "l6"), inherit.aes = FALSE,
                      linewidth = 1.2, alpha = 0.8) +
            # Markov: FPMs, ASF
            geom_line(data = OS_Markov_fpm_ac_RFC, aes(x = time, y = prob*100, color = "l7", linetype = "l7"), inherit.aes = FALSE,
                      linewidth = 1.2, alpha = 0.8) +
            # Markov: FPMs, RSF
            geom_line(data = OS_Markov_fpm_rel_RFC, aes(x = time, y = prob*100, color = "l8", linetype = "l8"), inherit.aes = FALSE,
                      linewidth = 1.2, alpha = 0.8) +
            # Observed 8 years  (Fischer2016)
            geom_step(data = mfit2_data_RFC, aes(x = time, y = surv*100, color = "l2", linetype = "l2"), inherit.aes = FALSE,
            direction = "hv", linewidth = 1.2, alpha = 0.8)  +
            # Observed 4 years (Hallek2010)
            geom_step(data = mfit1_data_RFC, aes(x = time, y = surv*100, color = "l1", linetype = "l1"), inherit.aes = FALSE,
            direction = "hv", linewidth = 1.2, alpha = 0.8) +
            labs(x = "Time (years)", y = "Survival probability (%)", fill = "") +
            scale_x_continuous(breaks = seq(0, 50, by = 5), limits = c(0, 50),
                               labels = seq(0, 50, by = 5), minor_breaks = seq(0, 50, by = 1)) +
            scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0, 100),
                               labels = seq(0, 100, by = 20), minor_breaks = seq(0, 100, by = 5)) +
            scale_color_manual(values = c("l1" = "black", "l2" = "darkolivegreen3",
                                          "l3" = "royalblue2", "l4" = "mediumorchid3", 
                                          "l5" = "gold", "l6" = "darkorange1",
                                          "l7" = "firebrick2", "l8" = "cyan"), 
                               labels = c("Observed 4-year OS", "Observed 8-year OS", 
                                          "Semi-Markov: SPMs, ASF (Williams et al.)", "Semi-Markov: FPMs, ASF", "Semi-Markov: FPMs, RSF (Proposed method)", 
                                          "Markov: SPMs, ASF", "Markov: FPMs, ASF", "Markov: FPMs, RSF"), 
                               name = "") +
            scale_linetype_manual(values = c("l1" = "solid", "l2" = "solid", "l3" = "dashed", "l4" = "twodash", "l5" = "dotdash", "l6" = "longdash", "l7" = "dashed", "l8" = "twodash"), 
                                  labels = c("Observed 4-year OS", "Observed 8-year OS", "Semi-Markov: SPMs, ASF (Williams et al.)", 
                                             "Semi-Markov: FPMs, ASF", "Semi-Markov: FPMs, RSF (Proposed method)", 
                                             "Markov: SPMs, ASF", "Markov: FPMs, ASF", "Markov: FPMs, RSF"), 
                                  name = "") +
            guides(x = "axis_minor", y = "axis_minor",
                   color = guide_legend(nrow = 4, override.aes = list(linewidth = 1.2)), linetype = guide_legend(nrow = 4), keyheight = unit(2, "lines")) +
            theme(legend.position = "bottom",
                  legend.title = element_text(size = 24),
                  legend.key = element_rect(colour = NA, fill = NA),
                  legend.key.width = unit(1.5, "cm"),
                  panel.grid = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(color = "black"),
                  plot.title = element_text(hjust = 0.5, size = 24),
                  axis.text = element_text(size = 12),
                  axis.title.x = element_text(size = 12),
                  axis.title.y = element_text(size = 12),
                  legend.text = element_text(size = 12)) +
            ggtitle("RFC")
plot_RFC  

plot_FC <- ggplot() +
  # Semi-Markov: SPMs, ASF (Williams2017)
  geom_line(data = OS_semiMarkov_williams_ac_FC, aes(x = time, y = prob*100, color = "l3", linetype = "l3"), inherit.aes = FALSE,
            linewidth = 1.2, alpha = 1)  +
  # Semi-Markov: FPMs, ASF
  geom_line(data = OS_semiMarkov_fpm_ac_FC, aes(x = time, y = prob*100, color = "l4", linetype = "l4"), inherit.aes = FALSE,
            linewidth = 1.2, alpha = 1)  +
  # Semi-Markov: FPMs, RSF (Proposed method)
  geom_line(data = OS_semiMarkov_fpm_rel_FC, aes(x = time, y = prob*100, color = "l5", linetype = "l5"), inherit.aes = FALSE,
            linewidth = 1.2, alpha = 0.8) +
  # Markov: SPMs, ASF
  geom_line(data = OS_Markov_williams_ac_FC, aes(x = time, y = prob*100, color = "l6", linetype = "l6"), inherit.aes = FALSE,
            linewidth = 1.2, alpha = 0.8) +
  # Markov: FPMs, ASF
  geom_line(data = OS_Markov_fpm_ac_FC, aes(x = time, y = prob*100, color = "l7", linetype = "l7"), inherit.aes = FALSE,
            linewidth = 1.2, alpha = 0.8) +
  # Markov: FPMs, RSF
  geom_line(data = OS_Markov_fpm_rel_FC, aes(x = time, y = prob*100, color = "l8", linetype = "l8"), inherit.aes = FALSE,
            linewidth = 1.2, alpha = 0.8) +
  # Observed 8 years  (Fischer2016)
  geom_step(data = mfit2_data_FC, aes(x = time, y = surv*100, color = "l2", linetype = "l2"), inherit.aes = FALSE,
            direction = "hv", linewidth = 1.2, alpha = 0.8)  +
  # Observed 4 years (Hallek2010)
  geom_step(data = mfit1_data_FC, aes(x = time, y = surv*100, color = "l1", linetype = "l1"), inherit.aes = FALSE,
            direction = "hv", linewidth = 1.2, alpha = 0.8) +
  labs(x = "Time (years)", y = "Survival probability (%)", fill = "") +
  scale_x_continuous(breaks = seq(0, 50, by = 5), limits = c(0, 50),
                     labels = seq(0, 50, by = 5), minor_breaks = seq(0, 50, by = 1)) +
  scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0, 100),
                     labels = seq(0, 100, by = 20), minor_breaks = seq(0, 100, by = 5)) +
  scale_color_manual(values = c("l1" = "black", "l2" = "darkolivegreen3",
                                "l3" = "royalblue2", "l4" = "mediumorchid3", 
                                "l5" = "gold", "l6" = "darkorange1",
                                "l7" = "firebrick2", "l8" = "cyan"), 
                     labels = c("Observed 4-year OS", "Observed 8-year OS", 
                                "Semi-Markov: SPMs, ASF (Williams et al.)", "Semi-Markov: FPMs, ASF", "Semi-Markov: FPMs, RSF (Proposed method)", 
                                "Markov: SPMs, ASF", "Markov: FPMs, ASF", "Markov: FPMs, RSF"), 
                     name = "") +
  scale_linetype_manual(values = c("l1" = "solid", "l2" = "solid", "l3" = "dashed", "l4" = "twodash", "l5" = "dotdash", "l6" = "longdash", "l7" = "dashed", "l8" = "twodash"), 
                        labels = c("Observed 4-year OS", "Observed 8-year OS", "Semi-Markov: SPMs, ASF (Williams et al.)", 
                                   "Semi-Markov: FPMs, ASF", "Semi-Markov: FPMs, RSF (Proposed method)", 
                                   "Markov: SPMs, ASF", "Markov: FPMs, ASF", "Markov: FPMs, RSF"), 
                        name = "") +
  guides(x = "axis_minor", y = "axis_minor",
         color = guide_legend(nrow = 4, override.aes = list(linewidth = 1.2)), linetype = guide_legend(nrow = 4), keyheight = unit(2, "lines")) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 24),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.width = unit(1.5, "cm"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 24),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  ggtitle("FC")
plot_FC  

##############################################################
##============================================================
## Combine RFC and FC plots into one figure and export
##============================================================
##############################################################
## Plot two figures together
plot <- ggarrange(plot_RFC, plot_FC, ncol = 2,
                  common.legend = TRUE, legend="bottom",
                  align = "h") 
plot

## Save results
## Prevent from changing the results. We put # here.
# ggsave("../Output/04c_figure_survextrap.png", plot, width = 10, height = 7, units = "in", dpi = 300)
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

