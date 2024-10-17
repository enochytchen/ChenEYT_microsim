## Filename: 10_test
## Purpose: Plot cost-effectiveness plane for probabilistic sensitivity analysis
##          All-cause survival framework vs. Relative survival framework
## Caution: Must have installed the lastest Rtools and R4.2.2+
## Reference: Williams C, Lewsey JD, Briggs AH, Mackay DF. 
##            Med Decis Making. 2017 May;37(4):340–52.
##            Cost-effectiveness Analysis in R Using a Multi-state Modeling Survival Analysis Framework: A Tutorial © 2017 by Williams C. et al is licensed under CC BY 3.0. 

##############################################################
##============================================================
## Read packages
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
## Results from the relative survival framework
load("../Output/07_fullmod_fpm_rel_PSA.RData")
PSA_fpm_rel <- output_df
## Results from the all-cause survival framework
load("../Output/07_fullmod_fpm_ac_PSA.RData")
PSA_fpm_ac <- output_df
## Results from the all-cause survival framework, but using weibull distributions
load("../Output/07_fullmod_fpm_ac2_PSA.RData")
PSA_fpm_ac2 <- output_df
## Results from the all-cause survival framework, but using weibull distributions with time-dependent effect
load("../Output/07_fullmod_fpm_ac3_PSA.RData")
PSA_fpm_ac3 <- output_df

##############################################################
##============================================================
## Plot of cost-effectiveness plane to compare the spread
##============================================================
##############################################################
## Compare using the same models but within an all-cause or a relative survival framework
## Answer: different survival frameworks do not make a difference
plot <- ggplot() +
  geom_point(data = PSA_fpm_rel, aes(x = disQALY_inc, y = disCost_inc, color="FPMs, RSF"), alpha = 0.8, size = 2) +
  geom_point(data = PSA_fpm_ac, aes(x = disQALY_inc, y = disCost_inc, color="FPMs, ASF"), alpha = 0.8, size = 4) +
  scale_color_manual(values = c("FPMs, ASF" = "red2", "FPMs, RSF" = "gold")) +
  scale_x_continuous(breaks = seq(-2, 2, 1), limits = c(-2, 2), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 20000, 5000), expand = c(0, 0), limits = c(0, 20000)) +
  labs(x = "Incremental QALY", y = "Incremental Cost (£)") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  geom_abline(intercept = 0, slope = 30000, linetype = "dashed", color = "indianred1") +
  annotate(geom="text", x = 1, y = 18000, label = "CET = £30,000", color = "indianred1", size = 4) +
  theme_bw() +  # Set a white background
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"))
plot

## Compare using FPMs with time-dependent effect vs. df == 1 (Weibull)
## Both are within an all-cause survival framework
## Answer: Models with time-dependent effect have more spread out results than a simple modek
plot <- ggplot() +
  geom_point(data = PSA_fpm_ac, aes(x = disQALY_inc, y = disCost_inc, color="FPMs, ASF"), alpha = 0.8, size = 4) +
  geom_point(data = PSA_fpm_ac2, aes(x = disQALY_inc, y = disCost_inc, color="FPMs (df=1), ASF"), alpha = 0.8, size = 4) +
  scale_color_manual(values = c("FPMs, ASF" = "red2", "FPMs (df=1), ASF" = "blue")) +
  scale_x_continuous(breaks = seq(-2, 2, 1), limits = c(-2, 2), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 20000, 5000), expand = c(0, 0), limits = c(0, 20000)) +
  labs(x = "Incremental QALY", y = "Incremental Cost (£)") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  geom_abline(intercept = 0, slope = 30000, linetype = "dashed", color = "indianred1") +
  annotate(geom="text", x = 1, y = 18000, label = "CET = £30,000", color = "indianred1", size = 4) +
  theme_bw() +  # Set a white background
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"))
plot

## Compare using FPMs with df == 1 (Weibull) within an all-cause survival framework
## with or without time-dependent effects
## Answer: Models with time-dependent effect have more spread out results than a simple modek
plot <- ggplot() +
  geom_point(data = PSA_fpm_ac3, aes(x = disQALY_inc, y = disCost_inc, color="FPMs (df=1, dftvc=3), ASF"), alpha = 0.8, size = 4) +
  geom_point(data = PSA_fpm_ac2, aes(x = disQALY_inc, y = disCost_inc, color="FPMs (df=1), ASF"), alpha = 0.8, size = 4) +
  scale_color_manual(values = c("FPMs (df=1), ASF" = "red2", "FPMs (df=1, dftvc=3), ASF" = "blue")) +
  scale_x_continuous(breaks = seq(-2, 2, 1), limits = c(-2, 2), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 20000, 5000), expand = c(0, 0), limits = c(0, 20000)) +
  labs(x = "Incremental QALY", y = "Incremental Cost (£)") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  geom_abline(intercept = 0, slope = 30000, linetype = "dashed", color = "indianred1") +
  annotate(geom="text", x = 1, y = 18000, label = "CET = £30,000", color = "indianred1", size = 4) +
  theme_bw() +  # Set a white background
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"))
plot


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

