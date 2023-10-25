## Filename: 04c_figure_survextrap
## Purpose: Plot the observed and the extrapolated survival functions
##          Observed: Observed--4 years (Hallek 2010), Observed--8 years (Fischer 2016)
##          Extrapolated: Semi-Markov (Williams 2017), Microsimulation
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

##############################################################
##============================================================
## Install/Read packages
##============================================================
##############################################################
library(survival)
library(tidyverse)
library(ggplot2)
library(ggpubr)

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

###################
## Williams 2017 ##
###################
williams2017 <- read.table("../Data/01a_williams2017.txt", header = TRUE)
## Subset the survival data for treat == 0 and 1
williams_data_FC <- williams2017[grepl(0, williams2017$treat), ]
williams_data_RFC <- williams2017[grepl(1, williams2017$treat), ]

#####################
## Microsimulation ##
#####################
## Run the microsimulation model to obtain the results
## source("04a_rel_survextrap.R")
## Read the results
results_FC <- readRDS("../Data/04a_rel_survextrap_results_FC.rds")
results_RFC <- readRDS("../Data/04a_rel_survextrap_results_RFC.rds")

##############################################################
##============================================================
## Plot the observed and the extrapolated results
##============================================================
##############################################################
calculateAndPlot_prev <- function(sims, mfit1_data, mfit2_data, williams_data) {
  
  # Claim data frame
  prevalence <- as.data.frame(sims$prev)
  
  # Calculate the probabilities
  prevalence$prob <- prevalence$number / sims$n
  
  prevalence <- prevalence %>%
                mutate(state = str_replace(state, "Prog", "Progression"),
                       state = str_replace(state, "ExcD", "Excess death"),
                       state = str_replace(state, "ExpD", "Expected death"))
  
  # Define the desired order of the "state" variable
  order <- c("Expected death", "Excess death", "Progression", "PFS")
  prevalence$state <- factor(prevalence$state, levels = order)
  
  # Prepare legend through scale_custom
  scale_custom <- list(
    scale_fill_manual(values = c("Expected death" = "gray75", 
                                 "Excess death" = "gray50",  
                                 "Progression" = "indianred1", 
                                 "PFS" = "cornflowerblue")),
    scale_color_manual(values = c("Observed--4 years" = "black", 
                                  "Observed--8 years" = "chartreuse",
                                  "Semi-Markov" = "gold1")))
  
  # Plot
  plot <- ggplot(data = prevalence, aes(x = time, y = prob, fill = state)) +
          # Microsimulation
          geom_bar(stat = "identity", width = 0.2) +
          labs(x = "Time (years)", y = "Probability", fill = "") +
          scale_x_continuous(breaks = seq(0, 15, 2), limits = c(0, 15), expand = c(0, 0)) +
          scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0, 0), limits = c(0, 1)) +
          theme(panel.grid = element_blank(),
                axis.text = element_text(size = 12),
                axis.title.x = element_text(size = 12),
                axis.title.y = element_text(size = 12),
                legend.text = element_text(size = 18),
                legend.key.size = unit(0.4, 'cm'),
                legend.position="bottom", 
                legend.box="vertical",
                legend.key = element_rect(fill = "transparent")) + 
          
          # Add mfit1 (Hallek2010)
          geom_step(data = mfit1_data, aes(x = time, y = surv, color = "Observed--4 years"), inherit.aes = FALSE,
                    linetype = "solid", direction = "hv", linewidth = 1.2) +
          # Add mfit2 (Fischer2016)
          geom_step(data = mfit2_data, aes(x = time, y = surv, color = "Observed--8 years" ), inherit.aes = FALSE,
                    linetype = "solid", direction = "hv", linewidth = 1.2)  +
          # Add williams (Williams2017)
          geom_line(data = williams_data, aes(x = time, y = surv, color = "Semi-Markov"), inherit.aes = FALSE,
                    linetype = "solid", linewidth = 1.2)  +
          # Legend
          scale_custom +
          guides(fill = guide_legend(title = ""), 
                 color = guide_legend(title = ""))
  
  return(plot)
    
  }

prevplotlist <- lapply(list(list(results_FC, mfit1_data_FC, mfit2_data_FC, williams_data_FC), 
                            list(results_RFC, mfit1_data_RFC, mfit2_data_RFC,  williams_data_RFC)), 
                       function(x) {
                         calculateAndPlot_prev(x[[1]], x[[2]], x[[3]],  x[[4]])
                       })

##############################################################
##============================================================
## Combine RFC and FC plots into one figure and export
##============================================================
##############################################################
## FC, treat == 0
prev_FC <- prevplotlist[[1]]
prev_FC <- prev_FC + labs(title = "FC") + 
           theme(plot.title = element_text(hjust = 0.5, size = 18))


## RFC, treat == 1
prev_RFC <- prevplotlist[[2]]
prev_RFC <- prev_RFC + labs(title = "RFC") +
            theme(plot.title = element_text(hjust = 0.5, size = 18))


## Plot two figures together
plot <- ggarrange(prev_FC, prev_RFC, ncol = 2,
                  common.legend = TRUE, legend="bottom",
                  align = "h") 
plot

## Assign the output route
getwd()
ggsave("../Output/04c_figure_survextrap.png", plot, width = 10, height = 7, units = "in", dpi = 300)
################################################################
setwd(original_wd)  # Reset to the original working directory
################################################################
## A microsimulation model incorporating relative survival extrapolation and multiple timescales for health technology assessment © 2023 by Chen EYT, Dickman PW, Clements MS is licensed under CC BY 4.0


