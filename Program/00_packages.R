## Filename: 00_packages
## Purpose: Install all the required packages to run the codes in this project

install.packages("survHE")
install.packages("survival")
install.packages("dplyr")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("survRM2")
install.packages("pracma")
install.packages("Rcpp")
install.packages("officer")
install.packages("flextable")
install.packages("remotes")
install.packages("mc2d")
install.packages("mvtnorm")
install.packages("biostat3")

## Some packages need to be pulled from GitHub
library(remotes)
remotes::install_github("mclements/rstpm2", ref="develop")
remotes::install_github("mclements/microsimulation")
remotes::install_github("chjackson/flexsurv") # to fix a bug in relative survival modelling

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
