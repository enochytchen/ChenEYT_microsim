## A Multistate Model Incorporating Relative Survival Extrapolation and Mixed Time Scales for Health Technology Assessment

This repository contains the R code for the multistate model incorporating relative survival extrapolation and mixed time scales for health technology assessment. The model is described in the paper by Chen et al. (2024) [1]. 

## File Descriptions

| **Filename**                              | **Purpose**                                                                                           |
|-------------------------------------------|-------------------------------------------------------------------------------------------------------|
|**Install required packages**| |
| `00_packages.R`                           | Installs and loads required R packages                                                              |
|**Read data**| |
| `01a_fischer2016.R`                       | Organises progression-free survival from CLL-8 trial data with 8 years of follow-up                                                |
| `01a_fischer2016_PFS.R`                       | Processes CLL-8 trial data with 8 years of follow-up                                                |
| `01a_williams2017.R`                      | Prepares data based on Williams et al. (2017)                                                       |
| `02_msmcancer.R`                          | Converts survival data into a multistate structure
|**Survival models**| |
| `03a_survmod_fpm.R`                       | Fits flexible parametric survival models                                                            |
| `03a_survmod_williams.R`                  | Fits standard parametric survival models as per Williams et al. (2017)                              |
| `03b_compare_hazard.R`                    | Compares extrapolated hazard functions from different survival models                               |
|**Compare survival extrapolations**| |
| `04a1_semiMarkov_fpm_ac_microsim.R`       | Implements semi-Markov (clock-reset) models using flexible parametric models within an all-cause survival framework (using the `microsimulation` package)|
| `04a2_semiMarkov_fpm_rel_hesim.R`         | Proposed method: semi-Markov models with flexible parametric models in a relative survival framework (using the `hesim` package) |
| `04a2_semiMarkov_fpm_rel_microsim.R`      | Proposed method: semi-Markov models using flexible parametric models in a relative survival framework (using the `microsimulation` package) |
| `04a3_semiMarkov_williams_ac_microsim.R`  | Semi-Markov models with standard parametric models in an all-cause survival framework (using the `microsimulation` package) |
| `04a3_semiMarkov_williams_ac_mstate.R`    | Semi-Markov models with standard parametric models (using the `mstate` package)                         |
| `04b_table_survextrap.R`                  | Generates tables of observed vs. extrapolated survival                                              |
| `04c_figure_survextrap.R`                 | Plots observed vs. extrapolated survival                                                            |
|**Full multistate models for cost-effectiveness analyses**| |
| `05_fullmod_semiMarkov_fpm_rel_microsim.R`| Builds a full health economic model using semi-Markov models (flexible parametric, relative survival) |
| `06_fullmod_onewaysen.R`                  | Performs one-way sensitivity analysis on the full model                                              |
| `07_fullmod_fpm_ac_PSA.R`                        | Runs microsimulation model within an all-cause survival framework for probabilistic sensitivity analysis                                                   |
| `07_fullmod_fpm_ac2_PSA.R`                        | Runs microsimulation model within an all-cause survival framework for probabilistic sensitivity analysis, using very simple models (flexible parametric models with df == 1, weibull for all transitions)                                                   |
| `07_fullmod_fpm_ac3_PSA.R`                        | Runs microsimulation model within an all-cause survival framework for probabilistic sensitivity analysis, using very simple models(flexible parametric models with df == 1, weibull for all transitions, dftvc == 2)                                               |
| `07_fullmod_fpm_rel_PSA.R`                        | Runs microsimulation model within a relative survival framework for probabilistic sensitivity analysis |
|**Compare results using the cost-effectiveness plane and cost-effectiveness acceptability curve**| |
| `08_CEplane_CEAC.R`                       | Plots cost-effectiveness planes and cost-effectiveness acceptability curves                          |
| `10_test.R`                       | Plots cost-effectiveness plane for probabilistic sensitivity analysis: all-cause survival framework vs. relative survival framework|
|**Others**| |
|rel_stackedplot_presentation.R| Stacked plot for presentation|
|williams_function PSAmeanLY.R| Function defined by Williams et al.(2017)[2] to estimate mean life years|
|williams_function PSAprob.R|Function defined by Williams et al.(2017)[2] to estimate transition probabilities |
|williams_function PSAmeanQALY.R| Function defined by Williams et al.(2017)[2] to estimate mean QALYs|


## References
1. Chen EYT, Dickman PW, Clements MS. A Multistate Model Incorporating Relative Survival Extrapolation and Mixed Time Scales for Health Technology Assessment. *PharmacoEconomics*. 2024. (doi:https://doi.org/10.1007/s40273-024-01457-w). 
2. Williams C, Lewsey JD, Briggs AH, Mackay DF. Cost-effectiveness Analysis in R Using a Multi-state Modelling Survival Analysis Framework: A Tutorial. *Med Desis Making*. 2017; 37(4): 340-52.

## License
Copyright © 2024. All rights reserved.<br>
![License](https://img.shields.io/badge/license-Apache%202.0-blue)
![R Version](https://img.shields.io/badge/R-%3E%3D4.4.1-blue)
<br>
Licensed to the Apache Software Foundation (ASF) under one or more
contributor license agreements.  See the NOTICE file distributed with
this work for additional information regarding copyright ownership.
The ASF licenses this file to You under the Apache License, Version 2.0
(the "License"); you may not use this file except in compliance with
the License.  You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.