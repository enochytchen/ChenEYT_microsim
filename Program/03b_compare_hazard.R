## Filename: 03b_compare_hazard
## Purpose: Compare the hazard functions of (progression–free -> progression; progression -> death; progression-free -> death)
##          1) Williams et al.'s hazard models (Gompertz; generalized gamma; Gompertz) in an all-cause survival framework
##          2) Flexible parametric models in an all-cause survival framework
##          3) Flexible parametric models in an relative survival framework
##
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
library(biostat3)
library(ggplot2)

##############################################################
##============================================================
## Read the survival data
##============================================================
##############################################################
source("03a_survmod_fpm.R")
source("03a_survmod_williams.R")

##############################################################
##============================================================
## Transition 1 PFS to prog
##============================================================
##############################################################
## Plot
## Estimate hazard functions of 4 years
haz_1RFC <- muhaz2(Surv(time, status) ~ 1, data = msmcancer1RFC)
haz_1FC <- muhaz2(Surv(time, status) ~ 1, data = msmcancer1FC)

plotdataRFC_haz <- data.frame(time = haz_1RFC$est.grid, haz= haz_1RFC$haz.est)
plotdataFC_haz <- data.frame(time = haz_1FC$est.grid, haz= haz_1FC$haz.est)

## Plot using ggplot
## RFC
trans1_RFC <- ggplot() + 
              geom_line(data=plotdataRFC_haz, aes(x=time, y=haz, color = "Observed"), linetype = "solid", size = 1.5) +
              geom_line(data=plotdataRFC_m1_fpm, aes(x=Tstop, y=haz_m1_fpm, color = "FPM (ASF)"), linetype = "dashed", size = 1.5) + 
              geom_line(data=plotdataRFC_m1_w, aes(x=Tstop, y=haz_m1_gom, color = "Gompertz"), linetype = "dashed", size = 1.5) +
              geom_line(data=plotdataRFC_m1_w, aes(x=0, y=0, color = "G-gamma"), linetype = "dashed", size = 1.5) +   # for legend
              geom_line(data=plotdataRFC_m1_w, aes(x=0, y=0, color = "FPM (RSF)"), linetype = "dashed", size = 1.5) + # for legend
              scale_x_continuous(breaks = seq(0, 15, by = 1), limits = c(0, 15),
                                 labels = seq(0, 15, by = 1)) +
              scale_y_continuous(breaks = seq(0, 2, by = 0.5), limits = c(0, 2), 
                                 labels = seq(0, 2, by = 0.5)) +
              labs(x = "", y = "Hazard",
                   subtitle = "Progression-free -> progression", 
                   title = "RFC") +
              theme_minimal() + 
              scale_color_manual(values = c("Observed" = "navy", "FPM (ASF)" = "blue", "FPM (RSF)" = "cyan", 
                                            "Gompertz" = "darkolivegreen3", "G-gamma" = "darkorange"),
                                 breaks = c("Observed", "FPM (ASF)", "Gompertz", "FPM (RSF)", "G-gamma"),
                                 name = NULL) +
              theme(legend.position = "bottom",
                    legend.key = element_rect(colour = NA, fill = NA),
                    panel.grid = element_blank(),
                    panel.background = element_blank(),
                    axis.line = element_line(color = "black"),
                    plot.title = element_text(hjust = 0.5, size = 24),
                    plot.subtitle = element_text(hjust = 0.5, size = 18),
                    axis.text = element_text(size = 12),
                    axis.title.x = element_text(size = 18),
                    axis.title.y = element_text(size = 18),
                    legend.text = element_text(size = 18)) + 
              guides(color = guide_legend(nrow = 1))
trans1_RFC

## FC
trans1_FC <- ggplot() + 
              geom_line(data=plotdataFC_haz, aes(x=time, y=haz, color = "Observed"), linetype = "solid", size = 1.5) +
              geom_line(data=plotdataFC_m1_fpm, aes(x=Tstop, y=haz_m1_fpm, color = "FPM (ASF)"), linetype = "dashed", size = 1.5) + 
              geom_line(data=plotdataFC_m1_w, aes(x=Tstop, y=haz_m1_gom, color = "Gompertz"), linetype = "dashed", size = 1.5) +
              scale_x_continuous(breaks = seq(0, 15, by = 1), limits = c(0, 15),
                                 labels = seq(0, 15, by = 1)) +
              scale_y_continuous(breaks = seq(0, 2, by = 0.5), limits = c(0, 2), 
                                 labels = seq(0, 2, by = 0.5)) +
              labs(x = "", y = "Hazard",
                   subtitle = "Progression-free -> progression",
                   title = "FC") +
              theme_minimal() +
              scale_color_manual(values = c("Observed" = "navy", "FPM (ASF)" = "blue", 
                                            "Gompertz" = "darkolivegreen3"), 
                                 name = NULL) +
              theme(legend.position = "bottom",
                    legend.key = element_rect(colour = NA, fill = NA),
                    panel.grid = element_blank(),
                    panel.background = element_blank(),
                    axis.line = element_line(color = "black"),
                    plot.title = element_text(hjust = 0.5, size = 24),
                    plot.subtitle = element_text(hjust = 0.5, size = 18),
                    axis.text = element_text(size = 12),
                    axis.title.x = element_text(size = 18),
                    axis.title.y = element_text(size = 18),
                    legend.text = element_text(size = 18))
trans1_FC

##############################################################
##============================================================
## Transition 2 PFS to death
##============================================================
##############################################################
## Plot
## Estimate hazard functions of 4 years
haz_2RFC <- muhaz2(Surv(time, status) ~ 1, data = msmcancer2RFC)
haz_2FC <- muhaz2(Surv(time, status) ~ 1, data = msmcancer2FC)

plotdataRFC_haz <- data.frame(time = haz_2RFC$est.grid, haz= haz_2RFC$haz.est)
plotdataFC_haz <- data.frame(time = haz_2FC$est.grid, haz= haz_2FC$haz.est)

## Plot using ggplot
## RFC
trans2_RFC <- ggplot() + 
              geom_line(data=plotdataRFC_haz, aes(x=time, y=haz, color = "Observed"), linetype = "solid", size = 1.5) +
              geom_line(data=plotdataRFC_m2_fpm_ac, aes(x=Tstop, y=haz, color = "FPM (ASF)"), linetype = "dashed", size = 1.5) +
              geom_line(data=plotdataRFC_m2_fpm_rel, aes(x=Tstop, y=achaz, color = "FPM (RSF)"), linetype = "dashed", size = 1.5) + 
              geom_line(data=plotdataRFC_m2_w, aes(x=Tstop, y=haz_m2_gam, color = "G-gamma"), linetype = "dashed", size = 1.5) + 
              scale_x_continuous(breaks = seq(0, 15, by = 1), limits = c(0, 15),
                                 labels = seq(0, 15, by = 1)) +
              scale_y_continuous(breaks = seq(0, 0.1, by = 0.02), limits = c(0, 0.1), 
                                 labels = seq(0, 0.1, by = 0.02)) +
              labs(x = "", y = "Hazard",
                   subtitle = "Progression-free -> death") +
              theme_minimal() + 
              scale_color_manual(values = c("Observed" = "navy", "FPM (ASF)" = "blue", "FPM (RSF)" = "cyan", "G-gamma" = "darkorange"), 
                                 name = NULL) +
              theme(legend.position = "bottom",
                    legend.key = element_rect(colour = NA, fill = NA),
                    panel.grid = element_blank(),
                    panel.background = element_blank(),
                    axis.line = element_line(color = "black"),
                    plot.title = element_text(hjust = 0.5, size = 24),
                    plot.subtitle = element_text(hjust = 0.5, size = 18),
                    axis.text = element_text(size = 12),
                    axis.title.x = element_text(size = 18),
                    axis.title.y = element_text(size = 18),
                    legend.text = element_text(size = 18))
trans2_RFC

## FC
trans2_FC <-  ggplot() + 
              geom_line(data=plotdataFC_haz, aes(x=time, y=haz, color = "Observed"), linetype = "solid", size = 1.5) +
              geom_line(data=plotdataFC_m2_fpm_ac, aes(x=Tstop, y=haz, color = "FPM (ASF)"), linetype = "dashed", size = 1.5) +
              geom_line(data=plotdataFC_m2_fpm_rel, aes(x=Tstop, y=achaz, color = "FPM (RSF)"), linetype = "dashed", size = 1.5) + 
              geom_line(data=plotdataFC_m2_w, aes(x=Tstop, y=haz_m2_gam, color = "G-gamma"), linetype = "dashed", size = 1.5) + 
              scale_x_continuous(breaks = seq(0, 15, by = 1), limits = c(0, 15),
                                 labels = seq(0, 15, by = 1)) +
              scale_y_continuous(breaks = seq(0, 0.1, by = 0.02), limits = c(0, 0.1), 
                                 labels = seq(0, 0.1, by = 0.02)) +
              labs(x = "", y = "Hazard",
                   subtitle = "Progression-free -> death") +
              theme_minimal() + 
              scale_color_manual(values = c("Observed" = "navy", "FPM (ASF)" = "blue", "FPM (RSF)" = "cyan", "G-gamma" = "darkorange"), 
                                 name = NULL) +
              theme(legend.position = "bottom",
                    legend.key = element_rect(colour = NA, fill = NA),
                    panel.grid = element_blank(),
                    panel.background = element_blank(),
                    axis.line = element_line(color = "black"),
                    plot.title = element_text(hjust = 0.5, size = 24),
                    plot.subtitle = element_text(hjust = 0.5, size = 18),
                    axis.text = element_text(size = 12),
                    axis.title.x = element_text(size = 18),
                    axis.title.y = element_text(size = 18),
                    legend.text = element_text(size = 18))
trans2_FC

##############################################################
##============================================================
## Transition 3 Progression to death
##============================================================
##############################################################
## Plot
## Estimate hazard functions of 4 years
haz_3RFC <- muhaz2(Surv(time, status) ~ 1, data = msmcancer3RFC)
haz_3FC <- muhaz2(Surv(time, status) ~ 1, data = msmcancer3FC)

plotdataRFC_haz <- data.frame(time = haz_3RFC$est.grid, haz= haz_3RFC$haz.est)
plotdataFC_haz <- data.frame(time = haz_3FC$est.grid, haz= haz_3FC$haz.est)

## Plot using ggplot
## RFC
trans3_RFC <- ggplot() + 
              geom_line(data=plotdataRFC_haz, aes(x=time, y=haz, color = "Observed"), linetype = "solid", size = 1.5) +
              geom_line(data=plotdataRFC_m3_fpm_ac, aes(x=time, y=haz, color = "FPM (ASF)"), linetype = "dashed", size = 1.5) +
              geom_line(data=plotdataRFC_m3_fpm_rel, aes(x=time, y=achaz, color = "FPM (RSF)"), linetype = "dashed", size = 1.5) + 
              geom_line(data=plotdataRFC_m3_w, aes(x=time, y=haz_m3_gom, color = "Gompertz"), linetype = "dashed", size = 1.5) + 
              scale_x_continuous(breaks = seq(0, 15, by = 1), limits = c(0, 15),
                                 labels = seq(0, 15, by = 1)) +
              scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1), 
                                 labels = seq(0, 1, by = 0.2)) +
              labs(x = "Time (years)", y = "Hazard",
                   subtitle = "Progression -> death") +
              theme_minimal() + 
              scale_color_manual(values = c("Observed" = "navy", "FPM (ASF)" = "blue", "FPM (RSF)" = "cyan", "Gompertz" = "darkolivegreen3"), 
                                 name = NULL) +
              theme(legend.position = "bottom",
                    legend.key = element_rect(colour = NA, fill = NA),
                    panel.grid = element_blank(),
                    panel.background = element_blank(),
                    axis.line = element_line(color = "black"),
                    plot.title = element_text(hjust = 0.5, size = 24),
                    plot.subtitle = element_text(hjust = 0.5, size = 18),
                    axis.text = element_text(size = 12),
                    axis.title.x = element_text(size = 18),
                    axis.title.y = element_text(size = 18),
                    legend.text = element_text(size = 18))
trans3_RFC

## FC
trans3_FC <-  ggplot() + 
              geom_line(data=plotdataFC_haz, aes(x=time, y=haz, color = "Observed"), linetype = "solid", size = 1.5) +
              geom_line(data=plotdataFC_m3_fpm_ac, aes(x=time, y=haz, color = "FPM (ASF)"), linetype = "dashed", size = 1.5) +
              geom_line(data=plotdataFC_m3_fpm_rel, aes(x=time, y=achaz, color = "FPM (RSF)"), linetype = "dashed", size = 1.5) + 
              geom_line(data=plotdataFC_m3_w, aes(x=time, y=haz_m3_gom, color = "Gompertz"), linetype = "dashed", size = 1.5) + 
              scale_x_continuous(breaks = seq(0, 15, by = 1), limits = c(0, 15),
                                 labels = seq(0, 15, by = 1)) +
              scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1), 
                                 labels = seq(0, 1, by = 0.2)) +
              labs(x = "Time (years)", y = "Hazard",
                   subtitle = "Progression -> death") +
              theme_minimal() + 
              scale_color_manual(values = c("Observed" = "navy", "FPM (ASF)" = "blue", "FPM (RSF)" = "cyan", "Gompertz" = "darkolivegreen3"), 
                                 name = NULL) +
              theme(legend.position = "bottom",
                    legend.key = element_rect(colour = NA, fill = NA),
                    panel.grid = element_blank(),
                    panel.background = element_blank(),
                    axis.line = element_line(color = "black"),
                    plot.title = element_text(hjust = 0.5, size = 24),
                    plot.subtitle = element_text(hjust = 0.5, size = 18),
                    axis.text = element_text(size = 12),
                    axis.title.x = element_text(size = 18),
                    axis.title.y = element_text(size = 18),
                    legend.text = element_text(size = 18))
trans3_FC

################################################################
## Combine all plots together
################################################################
plot <- ggarrange(trans1_RFC, trans1_FC, trans2_RFC, trans2_FC, trans3_RFC, trans3_FC, 
                  ncol = 2, nrow = 3, 
                  common.legend = TRUE, legend="bottom")
plot

ggsave(filename = "../Output/03b_compare_hazard.pdf", plot, width = 10, height = 16, dpi = 1000)

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