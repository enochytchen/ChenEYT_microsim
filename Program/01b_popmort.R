## Filename: 01b_popmort
## Purpose: Calculate the average death rate, by age, sex, and calendar year
##          from For 11 countries (Hallek, Lancet, 2010. NCT00281918 study)
## Notes:   Data from http://www.mortality.org/
##          Death rates (period 1x1)

##############################################################
##============================================================
## Install/Read packages
##============================================================
##############################################################
library(tidyverse)
library(ggplot2)

##############################################################
##============================================================
## Read the data
##============================================================
##############################################################
## Totally 11 countries
file_names <- c("Australia.txt", "Austria.txt", "Belgium.txt", "Czech.txt", "Denmark.txt",
                "France.txt", "Germany.txt", "Israel.txt", "Italy.txt", "NewZealand.txt", "Spain.txt")
country_names <- c("Australia", "Austria", "Belgium", "Czech", "Denmark", "France", "Germany",
                   "Israel", "Italy", "NewZealand", "Spain")

## Create an empty list to store the data frames
data_list <- list()

## Loop over the files
for (i in seq_along(file_names)) {
  file_path <- paste("../Data/popmort/", file_names[i], sep = "")
  
  # Read the data from each file
  data <- read.table(file_path, skip = 3)
  
  # Create the "country" variable
  data$country <- country_names[i]
  
  # Add the data frame to the list
  data_list[[i]] <- data
}

## Combine all data frames into a single data frame
data <- do.call(rbind.data.frame, data_list)

##############################################################
##============================================================
## Organise the data
##============================================================
##############################################################
## Rename the columns
colnames(data) <- c("year", "age", "rate2", "rate1", "rate3", "country")
## Make all the missing as NA
data[data == "."] <- NA
## Make all age 110+ as 110 for making them numberic later
data <- data %>% 
        mutate(age=ifelse(age=="110+", 110, age))
## Remove rows with NA values
data <- data %>% 
        drop_na()
data <- dplyr::mutate_at(data, vars(1:5), as.numeric) 
data <- subset(data, year >= 1950 & age <= 100)    # Keep year >= 2000
data <- data[, c("year", "age", "rate1", "rate2")] # Be aware of male and female
                                                   # Eventually, sex 1 m, 2 f.
data <- tidyr::pivot_longer(data, cols = c("rate1", "rate2"), names_to = "sex", values_to = "rate")
data$sex <- ifelse(data$sex=="rate1", 1, 2)

##############################################################
##============================================================
## Take average by country
##============================================================
##############################################################
## Average across rates
data <- data %>%
        group_by(sex, year, age) %>%
        summarize(rate = mean(rate))

## The last year of enrollment was 2006 (Hallek et al. 2010)
data <- subset(data, year <= 2006)  # Keep year <= 2006

##############################################################
##============================================================
## Extrapolation to the future
##============================================================
##############################################################
## Duplicates for age 101-130 from age 100
data <- rbind(data,
              merge(data.frame(age=101:130),
                    transform(subset(data, age==100),age=NULL)))

## Duplicates for 2007-2120 from 2006
data <- rbind(data,
              merge(data.frame(year=2007:2120),
              transform(subset(data, year==2006), year=NULL)))

## Duplicates for 1940-1949 from 1950
data <- rbind(data,
              merge(data.frame(year=1940:1949),
                    transform(subset(data, year==1950),year=NULL)))

data <- arrange(data, year, sex, age)

##############################################################
##============================================================
## Take average across sexes
##============================================================
##############################################################
## Take the average age across sexes (M to F was 74% to 26% (Hallek et al. 2010) )
data <- data %>% 
        dplyr::group_by(age, year) %>%
        dplyr::summarize(rate = sum(ifelse(sex == 1, rate*0.74, rate*0.26)) / n(), .groups = "drop")

## Generate a sex variable which indicates the average rate across sexes
data$sex <- 3

##############################################################
##============================================================
## Plot to check the data
##============================================================
##############################################################
## 60-year-old across years
plot_data <- data %>%
             filter(age == 60, sex == 3)
ggplot(plot_data, aes(x = year, y = rate)) + 
       geom_line() 

## 2050 and 2006 years across ages
## They are supposed to overlap because Year 2007 to 2120 are duplicates from 2006
plot_data1 <- data %>%
             filter(year == 2050, sex == 3)
plot_data2 <- data %>%
             filter(year == 2006, sex == 3)
ggplot() + 
       geom_line(data=plot_data1, aes(x = age, y = rate)) +
       geom_line(data=plot_data2, aes(x = age, y = rate), color="red", alpha=0.5) 
  
## Compress and label data
## Prevent from changing the results. We put # here.
# saveRDS(data, "../Data/01b_popmort.rds")
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
