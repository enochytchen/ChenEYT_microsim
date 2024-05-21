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
library(ggh4x)

##############################################################
##============================================================
## Read the survival data
##============================================================
##############################################################
source("03a_survmod_fpm.R")
source("03a_survmod_williams.R")

##############################################################
##============================================================
## Transition 1 progression-free to prog
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
              geom_line(data=plotdataRFC_m1_fpm, aes(x=Tstop, y=haz_m1_fpm,  color = "lt2", linetype = "lt2"), linewidth = 1.5, alpha = 0.8) + 
              geom_line(data=plotdataRFC_m1_w, aes(x=Tstop, y=haz_m1_gom,  color = "lt4", linetype = "lt4"), linewidth = 1.5, alpha = 0.8) +
              geom_line(data=plotdataRFC_m1_w, aes(x=0, y=0,  color = "lt3", linetype = "lt3"), linewidth = 1.5, alpha = 0.8) +   # for legend
              geom_line(data=plotdataRFC_m1_w, aes(x=0, y=0,  color = "lt5", linetype = "lt5"), linewidth = 1.5, alpha = 0.8) + # for legend
              geom_line(data=plotdataRFC_haz, aes(x=time, y=haz, color = "lt1", linetype = "lt1"), linewidth = 1.5, alpha = 0.8) +
              geom_vline(xintercept = 4, linetype = "dotted", color = "darkgray") +  
              geom_vline(xintercept = 8, linetype = "dotted", color = "darkgray") +  
              geom_vline(xintercept = 15, linetype = "dotted", color = "darkgray") +
              scale_x_continuous(breaks = seq(0, 30, by = 5), limits = c(0, 30),
                                 labels = seq(0, 30, by = 5), minor_breaks = seq(0, 30, by = 1)) +
              scale_y_continuous(breaks = seq(0, 2, by = 0.5), limits = c(0, 2), 
                                 labels = seq(0, 2, by = 0.5)) +
              labs(x = "Time since study entry  (years)", y = "Hazard", title = "RFC") +
              scale_color_manual(values = c("lt1" = "black", "lt2" = "blue", "lt3" = "cyan", "lt4" = "darkolivegreen3", "lt5" = "darkorange"),
                                 labels = c("Observed", "FPM (ASF)", "FPM (RSF)", "Gompertz", "Generalized gamma"),
                                 name = "") +
              scale_linetype_manual(values = c("lt1" = "solid", "lt2" = "dashed", "lt3" = "twodash", "lt4" = "dotdash", "lt5" = "longdash"), 
                                    labels = c("Observed", "FPM (ASF)", "FPM (RSF)", "Gompertz", "Generalized gamma"),
                                    name = "") +
              guides(x = "axis_minor", y = "axis_minor",
                     color = guide_legend(nrow = 1, override.aes = list(linewidth = 1.2)), 
                     linetype = guide_legend(nrow = 1), keyheight = unit(2, "lines")) +
              theme(legend.position = "bottom",
                    legend.key = element_rect(colour = NA, fill = NA),
                    legend.key.width = unit(1.5, "cm"),
                    panel.grid = element_blank(),
                    panel.background = element_blank(),
                    axis.line = element_line(color = "black"),
                    plot.title = element_text(hjust = 0.5, size = 18),
                    plot.subtitle = element_text(hjust = 0.5, size = 18),
                    axis.text = element_text(size = 12),
                    axis.title.x = element_text(size = 18),
                    axis.title.y = element_text(size = 18),
                    legend.text = element_text(size = 18),
                    legend.box.spacing = unit(1, "pt")) 
trans1_RFC

## FC
trans1_FC <- ggplot() + 
              geom_line(data=plotdataFC_m1_fpm, aes(x=Tstop, y=haz_m1_fpm,  color = "lt2", linetype = "lt2"), linewidth = 1.5, alpha = 0.8) + 
              geom_line(data=plotdataFC_m1_w, aes(x=Tstop, y=haz_m1_gom,  color = "lt4", linetype = "lt4"), linewidth = 1.5, alpha = 0.8) +
              geom_line(data=plotdataFC_m1_w, aes(x=0, y=0,  color = "lt3", linetype = "lt3"), linewidth = 1.5, alpha = 0.8) +   # for legend
              geom_line(data=plotdataFC_m1_w, aes(x=0, y=0,  color = "lt5", linetype = "lt5"), linewidth = 1.5, alpha = 0.8) + # for legend
              geom_line(data=plotdataFC_haz, aes(x=time, y=haz, color = "lt1", linetype = "lt1"), linewidth = 1.5, alpha = 0.8) +
              geom_vline(xintercept = 4, linetype = "dotted", color = "darkgray") +  
              geom_vline(xintercept = 8, linetype = "dotted", color = "darkgray") +  
              geom_vline(xintercept = 15, linetype = "dotted", color = "darkgray") +
              scale_x_continuous(breaks = seq(0, 30, by = 5), limits = c(0, 30),
                                 labels = seq(0, 30, by = 5), minor_breaks = seq(0, 30, by = 1)) +
              scale_y_continuous(breaks = seq(0, 2, by = 0.5), limits = c(0, 2), 
                                 labels = seq(0, 2, by = 0.5)) +
              labs(x = "Time since study entry  (years)", y = "Hazard",
                   title = "FC") +
              scale_color_manual(values = c("lt1" = "black", "lt2" = "blue", "lt3" = "cyan", "lt4" = "darkolivegreen3", "lt5" = "darkorange"),
                                 labels = c("Observed", "FPM (ASF)", "FPM (RSF)", "Gompertz", "Generalized gamma"),
                                 name = "") +
              scale_linetype_manual(values = c("lt1" = "solid", "lt2" = "dashed", "lt3" = "twodash", "lt4" = "dotdash", "lt5" = "longdash"), 
                                    labels = c("Observed", "FPM (ASF)", "FPM (RSF)", "Gompertz", "Generalized gamma"),
                                    name = "") +
              guides(x = "axis_minor", y = "axis_minor",
                     color = guide_legend(nrow = 1, override.aes = list(linewidth = 1.2)), 
                     linetype = guide_legend(nrow = 1), keyheight = unit(2, "lines")) +
              theme(legend.position = "bottom",
                    legend.key = element_rect(colour = NA, fill = NA),
                    legend.key.width = unit(1.5, "cm"),
                    panel.grid = element_blank(),
                    panel.background = element_blank(),
                    axis.line = element_line(color = "black"),
                    plot.title = element_text(hjust = 0.5, size = 18),
                    plot.subtitle = element_text(hjust = 0.5, size = 18),
                    axis.text = element_text(size = 12),
                    axis.title.x = element_text(size = 18),
                    axis.title.y = element_text(size = 18),
                    legend.text = element_text(size = 18)) 
trans1_FC

plot1 <- ggarrange(trans1_RFC, trans1_FC, ncol = 2, nrow = 1, common.legend = TRUE, legend="none")
plot1 <- annotate_figure(plot1, top = text_grob("(A) Transition 1: Progression-free -> progression", 
                                      color = "black", face = "bold", size = 24))
plot1

##############################################################
##============================================================
## Transition 2 progression-free to death
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
              geom_line(data=plotdataRFC_m2_fpm_ac, aes(x=Tstop, y=haz, color = "lt2", linetype = "lt2"), linewidth = 1.5, alpha = 0.8) +
              geom_line(data=plotdataRFC_m2_fpm_rel, aes(x=Tstop, y=achaz, color = "lt3", linetype = "lt3"), linewidth = 1.5, alpha = 0.8) + 
              geom_line(data=plotdataRFC_m2_w, aes(x=Tstop, y=haz_m2_gam, color = "lt5", linetype = "lt5"), linewidth = 1.5, alpha = 0.8) + 
              geom_line(data=plotdataRFC_haz, aes(x=time, y=haz, color = "lt1", linetype = "lt1"), linewidth = 1.5, alpha = 0.8) +
              geom_vline(xintercept = 4, linetype = "dotted", color = "darkgray") +  
              geom_vline(xintercept = 8, linetype = "dotted", color = "darkgray") +  
              geom_vline(xintercept = 15, linetype = "dotted", color = "darkgray") +
              scale_x_continuous(breaks = seq(0, 30, by = 5), limits = c(0, 30),
                                 labels = seq(0, 30, by = 5), minor_breaks = seq(0, 30, by = 1)) +
              scale_y_continuous(breaks = seq(0, 0.1, by = 0.02), limits = c(0, 0.1), 
                     labels = seq(0, 0.1, by = 0.02)) +
              labs(x = "Time since study entry  (years)", y = "Hazard", title = "RFC") +
              scale_color_manual(values = c("lt1" = "black", "lt2" = "blue", "lt3" = "cyan", "lt5" = "darkorange"),
                                 labels = c("Observed", "FPM (ASF)", "FPM (RSF)", "Generalized gamma"),
                                 name = "") +
              scale_linetype_manual(values = c("lt1" = "solid", "lt2" = "dashed", "lt3" = "twodash", "lt5" = "longdash"), 
                                    labels = c("Observed", "FPM (ASF)", "FPM (RSF)", "Generalized gamma"),
                                    name = "") +
              guides(x = "axis_minor", y = "axis_minor",
                     color = guide_legend(nrow = 1, override.aes = list(linewidth = 1.2)), 
                     linetype = guide_legend(nrow = 1), keyheight = unit(2, "lines")) +
              theme(legend.position = "bottom",
                    legend.key = element_rect(colour = NA, fill = NA),
                    legend.key.width = unit(1.5, "cm"),
                    panel.grid = element_blank(),
                    panel.background = element_blank(),
                    axis.line = element_line(color = "black"),
                    plot.title = element_text(hjust = 0.5, size = 18),
                    plot.subtitle = element_text(hjust = 0.5, size = 18),
                    axis.text = element_text(size = 12),
                    axis.title.x = element_text(size = 18),
                    axis.title.y = element_text(size = 18),
                    legend.text = element_text(size = 18)) 
trans2_RFC

## FC
trans2_FC <- ggplot() + 
              geom_line(data=plotdataFC_m2_fpm_ac, aes(x=Tstop, y=haz, color = "lt2", linetype = "lt2"), linewidth = 1.5, alpha = 0.8) +
              geom_line(data=plotdataFC_m2_fpm_rel, aes(x=Tstop, y=achaz, color = "lt3", linetype = "lt3"), linewidth = 1.5, alpha = 0.8) + 
              geom_line(data=plotdataFC_m2_w, aes(x=Tstop, y=haz_m2_gam, color = "lt5", linetype = "lt5"), linewidth = 1.5, alpha = 0.8) + 
              geom_line(data=plotdataFC_haz, aes(x=time, y=haz, color = "lt1", linetype = "lt1"), linewidth = 1.5, alpha = 0.8) +                geom_vline(xintercept = 4, linetype = "dotted", color = "darkgray") +  
              geom_vline(xintercept = 8, linetype = "dotted", color = "darkgray") +  
              geom_vline(xintercept = 15, linetype = "dotted", color = "darkgray") +
              scale_x_continuous(breaks = seq(0, 30, by = 5), limits = c(0, 30),
                                 labels = seq(0, 30, by = 5), minor_breaks = seq(0, 30, by = 1)) +
              scale_y_continuous(breaks = seq(0, 0.1, by = 0.02), limits = c(0, 0.1), 
                                 labels = seq(0, 0.1, by = 0.02)) +
              labs(x = "Time since study entry  (years)", y = "Hazard", title = "FC") +
              scale_color_manual(values = c("lt1" = "black", "lt2" = "blue", "lt3" = "cyan", "lt5" = "darkorange"),
                                 labels = c("Observed", "FPM (ASF)", "FPM (RSF)", "Generalized gamma"),
                                 name = "") +
              scale_linetype_manual(values = c("lt1" = "solid", "lt2" = "dashed", "lt3" = "twodash", "lt5" = "longdash"), 
                                    labels = c("Observed", "FPM (ASF)", "FPM (RSF)", "Generalized gamma"),
                                    name = "") +
              guides(x = "axis_minor", y = "axis_minor",
                     color = guide_legend(nrow = 1, override.aes = list(linewidth = 1.2)), 
                     linetype = guide_legend(nrow = 1), keyheight = unit(2, "lines")) +
              theme(legend.position = "bottom",
                    legend.key = element_rect(colour = NA, fill = NA),
                    legend.key.width = unit(1.5, "cm"),
                    panel.grid = element_blank(),
                    panel.background = element_blank(),
                    axis.line = element_line(color = "black"),
                    plot.title = element_text(hjust = 0.5, size = 18),
                    plot.subtitle = element_text(hjust = 0.5, size = 18),
                    axis.text = element_text(size = 12),
                    axis.title.x = element_text(size = 18),
                    axis.title.y = element_text(size = 18),
                    legend.text = element_text(size = 18)) 
trans2_FC

plot2 <- ggarrange(trans2_RFC, trans2_FC, ncol = 2, nrow = 1, legend="none")
plot2 <- annotate_figure(plot2, top = text_grob("(B) Transition 2: Progression-free -> death", 
                                                color = "black", face = "bold", size = 24))
plot2

##############################################################
##============================================================
## Transition 3 Progression to death (semi-Markov)
##============================================================
##############################################################
## Plot
## Estimate hazard functions of 4 years
haz_3RFC <- muhaz2(Surv(time, status) ~ 1, data = msmcancer3RFC)
haz_3FC <- muhaz2(Surv(time, status) ~ 1, data = msmcancer3FC)

plotdataRFC_haz <- data.frame(time = haz_3RFC$est.grid, haz_semiMarkov= haz_3RFC$haz.est)
plotdataFC_haz<- data.frame(time = haz_3FC$est.grid, haz_semiMarkov= haz_3FC$haz.est)

## Set the legend order
breaks = names(colors)[c(4,2,1,3)]
## Plot using ggplot
## RFC
trans3_RFC <- ggplot() + 
              geom_line(data=plotdataRFC_m3_fpm_ac, aes(x=time, y=haz_semiMarkov, color = "lt2", linetype = "lt2"), linewidth = 1.5, alpha = 0.8) +
              geom_line(data=plotdataRFC_m3_fpm_rel, aes(x=time, y=achaz_semiMarkov, color = "lt3", linetype = "lt3"), linewidth = 1.5, alpha = 0.8) + 
              geom_line(data=plotdataRFC_m3_w, aes(x=time, y=haz_m3_gom_semiMarkov, color = "lt4", linetype = "lt4"), linewidth = 1.5, alpha = 0.8) +
              geom_line(data=plotdataRFC_m3_w, aes(x=0, y=0,  color = "lt5", linetype = "lt5"), linewidth = 1.5, alpha = 0.8) + # for legend
              geom_line(data=plotdataRFC_haz, aes(x=time, y=haz_semiMarkov, color = "lt1", linetype = "lt1"), linewidth = 1.5, alpha = 0.8) +                geom_vline(xintercept = 4, linetype = "dotted", color = "darkgray") +  
              geom_vline(xintercept = 8, linetype = "dotted", color = "darkgray") +  
              geom_vline(xintercept = 15, linetype = "dotted", color = "darkgray") +
              scale_x_continuous(breaks = seq(0, 30, by = 5), limits = c(0, 30),
                                 labels = seq(0, 30, by = 5), minor_breaks = seq(0, 30, by = 1)) +
              scale_y_continuous(breaks = seq(0, 0.5, by = 0.1), limits = c(0, 0.5), 
                                 labels = seq(0, 0.5, by = 0.1)) +
              labs(x = "Time since progression (years)", y = "Hazard", title = "RFC") +
              scale_color_manual(values = c("lt1" = "black", "lt2" = "blue", "lt3" = "cyan", "lt4" = "darkolivegreen3", "lt5" = "darkorange")[c(1,4,2,5,3)],
                                 labels = c("Observed", "FPM (ASF)", "FPM (RSF)", "Gompertz", "Generalized gamma"),
                                 name = "") +
              scale_linetype_manual(values = c("lt1" = "solid", "lt2" = "dashed", "lt3" = "twodash", "lt4" = "dotdash", "lt5" = "longdash"), 
                                    labels = c("Observed", "FPM (ASF)", "FPM (RSF)", "Gompertz", "Generalized gamma"),
                                    name = "") +
              guides(x = "axis_minor", y = "axis_minor",
                     color = guide_legend(nrow = 1, override.aes = list(linewidth = 1.2)), 
                     linetype = guide_legend(nrow = 1), keyheight = unit(2, "lines")) +
              theme(legend.position = "bottom",
                    legend.key = element_rect(colour = NA, fill = NA),
                    legend.key.width = unit(1.5, "cm"),
                    panel.grid = element_blank(),
                    panel.background = element_blank(),
                    axis.line = element_line(color = "black"),
                    plot.title = element_text(hjust = 0.5, size = 18),
                    plot.subtitle = element_text(hjust = 0.5, size = 18),
                    axis.text = element_text(size = 12),
                    axis.title.x = element_text(size = 18),
                    axis.title.y = element_text(size = 18),
                    legend.text = element_text(size = 18))
trans3_RFC

## FC
trans3_FC <- ggplot() + 
              geom_line(data=plotdataFC_m3_fpm_ac, aes(x=time, y=haz_semiMarkov, color = "lt2", linetype = "lt2"), linewidth = 1.5, alpha = 0.8) +
              geom_line(data=plotdataFC_m3_fpm_rel, aes(x=time, y=achaz_semiMarkov, color = "lt3", linetype = "lt3"), linewidth = 1.5, alpha = 0.8) + 
              geom_line(data=plotdataFC_m3_w, aes(x=time, y=haz_m3_gom_semiMarkov, color = "lt4", linetype = "lt4"), linewidth = 1.5, alpha = 0.8) + 
              geom_line(data=plotdataFC_m3_w, aes(x=0, y=0,  color = "lt5", linetype = "lt5"), linewidth = 1.5, alpha = 0.8) + # for legend
              geom_line(data=plotdataFC_haz, aes(x=time, y=haz_semiMarkov, color = "lt1", linetype = "lt1"), linewidth = 1.5, alpha = 0.8) +
              geom_vline(xintercept = 4, linetype = "dotted", color = "darkgray") +  
              geom_vline(xintercept = 8, linetype = "dotted", color = "darkgray") +  
              geom_vline(xintercept = 15, linetype = "dotted", color = "darkgray") +
              scale_x_continuous(breaks = seq(0, 30, by = 5), limits = c(0, 30),
                                 labels = seq(0, 30, by = 5), minor_breaks = seq(0, 30, by = 1)) +
              scale_y_continuous(breaks = seq(0, 0.5, by = 0.1), limits = c(0, 0.5), 
                                 labels = seq(0, 0.5, by = 0.1)) +
              labs(x = "Time since progression (years)", y = "Hazard", title = "FC") +
              scale_color_manual(values = c("lt1" = "black", "lt2" = "blue", "lt3" = "cyan", "lt4" = "darkolivegreen3", "lt5" = "darkorange"),
                                 labels = c("Observed", "FPM (ASF)", "FPM (RSF)", "Gompertz", "Generalized gamma"),
                                 name = "") +
              scale_linetype_manual(values = c("lt1" = "solid", "lt2" = "dashed", "lt3" = "twodash", "lt4" = "dotdash", "lt5" = "longdash"), 
                                    labels = c("Observed", "FPM (ASF)", "FPM (RSF)", "Gompertz", "Generalized gamma"),
                                    name = "") +
              guides(x = "axis_minor", y = "axis_minor",
                     color = guide_legend(nrow = 1, override.aes = list(linewidth = 1.2)), 
                     linetype = guide_legend(nrow = 1), keyheight = unit(2, "lines")) +
              theme(legend.position = "bottom",
                    legend.key = element_rect(colour = NA, fill = NA),
                    legend.key.width = unit(1.5, "cm"),
                    panel.grid = element_blank(),
                    panel.background = element_blank(),
                    axis.line = element_line(color = "black"),
                    plot.title = element_text(hjust = 0.5, size = 18),
                    plot.subtitle = element_text(hjust = 0.5, size = 18),
                    axis.text = element_text(size = 12),
                    axis.title.x = element_text(size = 18),
                    axis.title.y = element_text(size = 18),
                    legend.text = element_text(size = 18))
trans3_FC

plot3 <- ggarrange(trans3_RFC, trans3_FC, ncol = 2, nrow = 1, common.legend=TRUE, legend="none")
plot3 <- annotate_figure(plot3, top = text_grob("(C) Transition 3: Progression -> death (semi-Markov)", 
                                                color = "black", face = "bold", size = 24))
plot3


##############################################################
##============================================================
## Transition 3 Progression to death (Markov)
##============================================================
##############################################################
## Plot
## Estimate hazard functions of 4 years
haz_3RFC <- bshazard::bshazard(Surv(Tstart, Tstop, status) ~ 1, data = msmcancer3RFC)
haz_3FC <- bshazard::bshazard(Surv(Tstart, Tstop, status)~ 1, data = msmcancer3FC)

plotdataRFC_haz <- data.frame(time = haz_3RFC$time, haz_Markov= haz_3RFC$hazard)
plotdataFC_haz<- data.frame(time = haz_3FC$time, haz_Markov= haz_3FC$hazard)

## Set the legend order
breaks = names(colors)[c(4,2,1,3)]
## Plot using ggplot
## RFC
trans3_RFC <- ggplot() + 
  geom_line(data=plotdataRFC_m3_fpm_ac, aes(x=time, y=haz_Markov, color = "lt2", linetype = "lt2"), linewidth = 1.5, alpha = 0.8) +
  geom_line(data=plotdataRFC_m3_fpm_rel, aes(x=time, y=achaz_Markov, color = "lt3", linetype = "lt3"), linewidth = 1.5, alpha = 0.8) + 
  geom_line(data=plotdataRFC_m3_w, aes(x=time, y=haz_m3_gom_Markov, color = "lt4", linetype = "lt4"), linewidth = 1.5, alpha = 0.8) +
  geom_line(data=plotdataRFC_m3_w, aes(x=0, y=0,  color = "lt5", linetype = "lt5"), linewidth = 1.5, alpha = 0.8) + # for legend
  geom_line(data=plotdataRFC_haz, aes(x=time, y=haz_Markov, color = "lt1", linetype = "lt1"), linewidth = 1.5, alpha = 0.8) +                geom_vline(xintercept = 4, linetype = "dotted", color = "darkgray") +  
  geom_vline(xintercept = 8, linetype = "dotted", color = "darkgray") +  
  geom_vline(xintercept = 15, linetype = "dotted", color = "darkgray") +
  scale_x_continuous(breaks = seq(0, 30, by = 5), limits = c(0, 30),
                     labels = seq(0, 30, by = 5), minor_breaks = seq(0, 30, by = 1)) +
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.1), limits = c(0, 0.5), 
                     labels = seq(0, 0.5, by = 0.1)) +
  labs(x = "Time since study entry  (years)", y = "Hazard", title = "RFC") +
  scale_color_manual(values = c("lt1" = "black", "lt2" = "blue", "lt3" = "cyan", "lt4" = "darkolivegreen3", "lt5" = "darkorange")[c(1,4,2,5,3)],
                     labels = c("Observed", "FPM (ASF)", "FPM (RSF)", "Gompertz", "Generalized gamma"),
                     name = "") +
  scale_linetype_manual(values = c("lt1" = "solid", "lt2" = "dashed", "lt3" = "twodash", "lt4" = "dotdash", "lt5" = "longdash"), 
                        labels = c("Observed", "FPM (ASF)", "FPM (RSF)", "Gompertz", "Generalized gamma"),
                        name = "") +
  guides(x = "axis_minor", y = "axis_minor",
         color = guide_legend(nrow = 1, override.aes = list(linewidth = 1.2)), 
         linetype = guide_legend(nrow = 1), keyheight = unit(2, "lines")) +
  theme(legend.position = "bottom",
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.width = unit(1.5, "cm"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 18))
trans3_RFC

## FC
trans3_FC <- ggplot() + 
  geom_line(data=plotdataFC_m3_fpm_ac, aes(x=time, y=haz_Markov, color = "lt2", linetype = "lt2"), linewidth = 1.5, alpha = 0.8) +
  geom_line(data=plotdataFC_m3_fpm_rel, aes(x=time, y=achaz_Markov, color = "lt3", linetype = "lt3"), linewidth = 1.5, alpha = 0.8) + 
  geom_line(data=plotdataFC_m3_w, aes(x=time, y=haz_m3_gom_Markov, color = "lt4", linetype = "lt4"), linewidth = 1.5, alpha = 0.8) + 
  geom_line(data=plotdataFC_m3_w, aes(x=0, y=0,  color = "lt5", linetype = "lt5"), linewidth = 1.5, alpha = 0.8) + # for legend
  geom_line(data=plotdataFC_haz, aes(x=time, y=haz_Markov, color = "lt1", linetype = "lt1"), linewidth = 1.5, alpha = 0.8) +
  geom_vline(xintercept = 4, linetype = "dotted", color = "darkgray") +  
  geom_vline(xintercept = 8, linetype = "dotted", color = "darkgray") +  
  geom_vline(xintercept = 15, linetype = "dotted", color = "darkgray") +
  scale_x_continuous(breaks = seq(0, 30, by = 5), limits = c(0, 30),
                     labels = seq(0, 30, by = 5), minor_breaks = seq(0, 30, by = 1)) +
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.1), limits = c(0, 0.5), 
                     labels = seq(0, 0.5, by = 0.1)) +
  labs(x = "Time since study entry  (years)", y = "Hazard", title = "FC") +
  scale_color_manual(values = c("lt1" = "black", "lt2" = "blue", "lt3" = "cyan", "lt4" = "darkolivegreen3", "lt5" = "darkorange"),
                     labels = c("Observed", "FPM (ASF)", "FPM (RSF)", "Gompertz", "Generalized gamma"),
                     name = "") +
  scale_linetype_manual(values = c("lt1" = "solid", "lt2" = "dashed", "lt3" = "twodash", "lt4" = "dotdash", "lt5" = "longdash"), 
                        labels = c("Observed", "FPM (ASF)", "FPM (RSF)", "Gompertz", "Generalized gamma"),
                        name = "") +
  guides(x = "axis_minor", y = "axis_minor",
         color = guide_legend(nrow = 1, override.aes = list(linewidth = 1.2)), 
         linetype = guide_legend(nrow = 1), keyheight = unit(2, "lines")) +
  theme(legend.position = "bottom",
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.width = unit(1.5, "cm"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 18))
trans3_FC

plot4 <- ggarrange(trans3_RFC, trans3_FC, ncol = 2, nrow = 1, common.legend=TRUE, legend="bottom")
plot4 <- annotate_figure(plot4, top = text_grob("(D) Transition 3: Progression -> death (Markov)", 
                                                color = "black", face = "bold", size = 24))
plot4

################################################################
## Combine all plots together
################################################################
plot <- ggarrange(plot1, plot2, plot3, plot4,
                  ncol = 1, nrow = 4,
                  widths = c(0.05, 0.05, 0.05, 0.05),
                  common.legend = TRUE, legend="bottom")
plot

## Prevent from changing the results. We put # here.
# ggsave(filename = "../Output/03b_compare_hazard.png", plot, bg = "white",
#        width = 12, height = 20, dpi = 300)
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