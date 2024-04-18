## Filename: rel_stackedplot
## Purpose: Plot the stacked plots of transition probabilities and QALYs each state

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
# install.packages("survival")
# install.packages("tidyverse")
# install.packages("ggplot2")
# install.packages("ggpubr")
library(survival)
library(tidyverse)
library(ggplot2)
library(ggpubr)
##############################################################
##============================================================
## Read Hallek2010 and Fischer2016 as comparison
##============================================================
##############################################################
setwd("../Data")

## Hallek 2010
data<-read.table("data.txt",sep="\t",skip=1) 
names(data)=c("patid","treat","prog","death", "prog_t","death_t",
              "progdeath", "progdeath_t")
data$death_ty <- data$death_t/12
data <- data %>%
        mutate(death = ifelse(death_t>=47, 0 , death))

## Fit to plot K-M curve
mfit1 <- survfit(Surv(death_ty, death==1) ~ treat, data = subset(data, death_t <= 47))
plot(mfit1, conf.int = FALSE, mark.time = FALSE, lty = 1, col = c("red","blue"),
     xlim = c(0, 8), ylim=c(0,1))

# Make a data frame
mfit1_sum <- summary(mfit1)
mfit1_data <- data.frame(time = mfit1_sum$time,
                         surv = mfit1_sum$surv,
                         strata = mfit1_sum$strata)

# Subset the survival data for treat == 0 and 1
mfit1_data_FC <- mfit1_data[grepl("treat=0", mfit1_data$strata), ]
mfit1_data_RFC <- mfit1_data[grepl("treat=1", mfit1_data$strata), ]

# Add a time point for censoring, for the purpose of plotting
mfit1_data_FC[nrow(mfit1_data_FC)+1,] <- mfit1_data_FC[nrow(mfit1_data_FC),]
mfit1_data_FC[nrow(mfit1_data_FC), "time"] <- 4  
mfit1_data_RFC[nrow(mfit1_data_RFC)+1,] <- mfit1_data_RFC[nrow(mfit1_data_RFC),]
mfit1_data_RFC[nrow(mfit1_data_RFC), "time"] <- 4  


## Fischer2016 data
fischer2016 <- read.table("fischer2016.txt", header = TRUE)
fischer2016$time_y <- fischer2016$time/12

## Fit to plot K-M curve
mfit2 <- survfit(Surv(time_y, event==1) ~ treat, data = fischer2016)
plot(mfit2, conf.int = FALSE, mark.time = FALSE, lty = 1, col = c("red","blue"),
     xlim = c(0, 8), ylim=c(0,1))

# Make a data frame
mfit2_sum <- summary(mfit2)
mfit2_data <- data.frame(time = mfit2_sum$time,
                         surv = mfit2_sum$surv,
                         strata = mfit2_sum$strata)

# Subset the survival data for treat == 0 and 1
mfit2_data_FC <- mfit2_data[grepl("treat=0", mfit2_data$strata), ]
mfit2_data_RFC <- mfit2_data[grepl("treat=1", mfit2_data$strata), ]

# Add a time point for censoring, for the purpose of plotting
mfit2_data_FC[nrow(mfit2_data_FC)+1,] <- mfit2_data_FC[nrow(mfit2_data_FC),]
mfit2_data_FC[nrow(mfit2_data_FC), "time"] <- 8  
mfit2_data_RFC[nrow(mfit2_data_RFC)+1,] <- mfit2_data_RFC[nrow(mfit2_data_RFC),]
mfit2_data_RFC[nrow(mfit2_data_RFC), "time"] <- 8  

## Williams 2017 data
williams2017 <- read.table("williams2017.txt", header = TRUE)
# Subset the survival data for treat == 0 and 1
williams_data_FC <- williams2017[grepl(0, williams2017$treat), ]
williams_data_RFC <- williams2017[grepl(1, williams2017$treat), ]

##############################################################
##============================================================
## Plot transition probabilities
##============================================================
##############################################################
calculateAndPlot_prev_0 <- function(sims, mfit1_data, mfit2_data, williams_data) {
  
  # Prepare legend through scale_custom
  scale_custom <- list(
                  scale_color_manual(values = c("Hallek 2010" = "black", 
                                                "Williams 2017" = "gold1")))
  
   # Plot the stacked probabilities
  ggplot() +
    # Add mfit1 (Hallek2010)
    geom_step(data = mfit1_data, aes(x = time, y = surv, colour = "Hallek 2010"), inherit.aes = FALSE,
              linetype = "solid", direction = "hv", linewidth = 1) +
    # Add williams (Williams2017)
    geom_step(data = williams_data, aes(x = time, y = surv, colour = "Williams 2017" ), inherit.aes = FALSE,
              linetype = "solid", direction = "hv", linewidth = 1)  +
    labs(x = "Time (years)", y = "(A) Probability", fill = "") +
    scale_x_continuous(breaks = seq(0, 15, 2), limits = c(0, 15), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0, 0), limits = c(0, 1)) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          legend.text = element_text(size = 18),
          legend.key.size = unit(0.4, 'cm'),
          legend.position="bottom", 
          legend.box="vertical") + 

    # Legend
    scale_custom +
    guides(fill = guide_legend(title = ""), 
           color = guide_legend(title = ""))
  }

calculateAndPlot_prev_1 <- function(sims, mfit1_data, mfit2_data, williams_data) {
  
  # Claim data frame
  prevalence <- as.data.frame(sims$prev)
  
  # Calculate the stacked probabilities
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
    scale_color_manual(values = c("Hallek 2010" = "black", 
                                  "Fischer 2016" = "chartreuse",
                                  "Williams 2017" = "gold1")))
  
  # Plot the stacked probabilities
  ggplot(data = prevalence, aes(x = time, y = prob, fill = state)) +
    geom_bar(stat = "identity", width = 0.1) +
    labs(x = "Time (years)", y = "(A) Probability", fill = "") +
    scale_x_continuous(breaks = seq(0, 15, 2), limits = c(0, 15), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0, 0), limits = c(0, 1)) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          legend.text = element_text(size = 18),
          legend.key.size = unit(0.4, 'cm'),
          legend.position="bottom", 
          legend.box="vertical") + 
    
    # Add mfit1 (Hallek2010)
    geom_step(data = mfit1_data, aes(x = time, y = surv, colour = "Hallek 2010"), inherit.aes = FALSE,
              linetype = "solid", direction = "hv", linewidth = 1) +
    # Add williams (Williams2017)
    geom_step(data = williams_data, aes(x = time, y = surv, colour = "Williams 2017" ), inherit.aes = FALSE,
              linetype = "solid", direction = "hv", linewidth = 1)  +
    # Legend
    scale_custom +
    guides(fill = guide_legend(title = ""), 
           color = guide_legend(title = ""))
}

calculateAndPlot_prev_2 <- function(sims, mfit1_data, mfit2_data, williams_data) {
  
  # Claim data frame
  prevalence <- as.data.frame(sims$prev)
  
  # Calculate the stacked probabilities
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
    scale_color_manual(values = c("Hallek 2010" = "black", 
                                  "Fischer 2016" = "chartreuse",
                                  "Williams 2017" = "gold1")))
  
  # Plot the stacked probabilities
  ggplot(data = prevalence, aes(x = time, y = prob, fill = state)) +
    geom_bar(stat = "identity", width = 0.1) +
    labs(x = "Time (years)", y = "(A) Probability", fill = "") +
    scale_x_continuous(breaks = seq(0, 15, 2), limits = c(0, 15), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0, 0), limits = c(0, 1)) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          legend.text = element_text(size = 18),
          legend.key.size = unit(0.4, 'cm'),
          legend.position="bottom", 
          legend.box="vertical") + 
    
    # Add mfit1 (Hallek2010)
    geom_step(data = mfit1_data, aes(x = time, y = surv, colour = "Hallek 2010"), inherit.aes = FALSE,
              linetype = "solid", direction = "hv", linewidth = 1) +
    # Add mfit2 (Fischer2016)
    geom_step(data = mfit2_data, aes(x = time, y = surv, colour = "Fischer 2016" ), inherit.aes = FALSE,
              linetype = "solid", direction = "hv", linewidth = 1)  +
    # Add williams (Williams2017)
    geom_step(data = williams_data, aes(x = time, y = surv, colour = "Williams 2017" ), inherit.aes = FALSE,
              linetype = "solid", direction = "hv", linewidth = 1)  +
    # Legend
    scale_custom +
    guides(fill = guide_legend(title = ""), 
           color = guide_legend(title = ""))
}

prevplotlist_0 <- lapply(list(list(sims_FC, mfit1_data_FC, mfit2_data_FC, williams_data_FC), 
                            list(sims_RFC, mfit1_data_RFC, mfit2_data_RFC,  williams_data_RFC)), 
                       function(x) {
                         calculateAndPlot_prev_0(x[[1]], x[[2]], x[[3]],  x[[4]])
                       })

prevplotlist_1 <- lapply(list(list(sims_FC, mfit1_data_FC, mfit2_data_FC, williams_data_FC), 
                              list(sims_RFC, mfit1_data_RFC, mfit2_data_RFC,  williams_data_RFC)), 
                         function(x) {
                           calculateAndPlot_prev_1(x[[1]], x[[2]], x[[3]],  x[[4]])
                         })

prevplotlist_2 <- lapply(list(list(sims_FC, mfit1_data_FC, mfit2_data_FC, williams_data_FC), 
                              list(sims_RFC, mfit1_data_RFC, mfit2_data_RFC,  williams_data_RFC)), 
                         function(x) {
                           calculateAndPlot_prev_2(x[[1]], x[[2]], x[[3]],  x[[4]])
                         })
##############################################################
##============================================================
## Combine prev and util into one figure and export
##============================================================
##############################################################
## Assign the output route
getwd()
setwd("../Output")

## FC, treat == 0
prev_FC_0 <- prevplotlist_0[[1]]
prev_FC_0
prev_FC_1 <- prevplotlist_1[[1]]
prev_FC_1
prev_FC_2 <- prevplotlist_2[[1]]
prev_FC_2

ggsave("stack_FC_0.png", prev_FC_0, width = 7, height = 7, units = "in", dpi = 1000)
ggsave("stack_FC_1.png", prev_FC_1, width = 7, height = 7, units = "in", dpi = 1000)
ggsave("stack_FC_2.png", prev_FC_2, width = 7, height = 7, units = "in", dpi = 1000)

## RFC, treat == 1
prev_RFC_0 <- prevplotlist_0[[2]]
prev_RFC_0
prev_RFC_1 <- prevplotlist_1[[2]]
prev_RFC_1
prev_RFC_2 <- prevplotlist_2[[2]]
prev_RFC_2

ggsave("stack_RFC_0.png", prev_RFC_0, width = 7, height = 7, units = "in", dpi = 1000)
ggsave("stack_RFC_1.png", prev_RFC_1, width = 7, height = 7, units = "in", dpi = 1000)
ggsave("stack_RFC_2.png", prev_RFC_2, width = 7, height = 7, units = "in", dpi = 1000)

################################################################
setwd(original_wd)  # Reset to the original working directory

