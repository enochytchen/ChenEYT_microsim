## Filename: 07_fullmod_fpm_ac_PSA
## Purpose: Run microsimulation model for probabilistic sensitivity analysis
##          All-cause survival framework
## Caution: Must have installed the lastest Rtools and R4.2.2+

## Set seed for coherency
set.seed(12345)
##############################################################
##============================================================
## Read packages
##============================================================
##############################################################
library(tidyverse)
library(rstpm2)
library(microsimulation)
library(Rcpp)
library(mc2d) # for beta pert distribution
library(mvtnorm) # for bootstrapping parametric models 
library(parallel)

##############################################################
##============================================================
## Run the other programs
##============================================================
##############################################################
source("02_msmcancer.R")
source("03a_survmod_fpm.R")

##############################################################
##============================================================
## Determine life years in progression state for RFC and FC
##============================================================
##############################################################
## The calculation of the costs of second-line and subsequent therapy in the progression state
## involves life years spent in progression state RFC and FC
## We first run the model just in order to determine LY_Prog_RFC and LY_Prog_FC
sourceCpp(code="
  // BH is required for older versions of microsimulation
  // RcppArmadillo is required for newer versions of microsimulation
  // [[Rcpp::depends(BH)]]
  // [[Rcpp::depends(RcppArmadillo)]]
  // [[Rcpp::depends(microsimulation)]]
  #include <RcppArmadillo.h>
  #include <microsimulation.h>
  using namespace arma;
  using namespace ssim;
  enum state_t {PFS, Prog, AcD};
  enum event_t {toPFS, toProg, toAcD, toEOS};
  typedef ssim::SummaryReport<short,short> Report;
  /**
      Utility: Run a set of simulations for a single process
  */
  void runSimulations(ssim::cProcess* process, int n) {
    for (int i = 0; i < n; i++) {
      ssim::Sim::create_process(process);
      ssim::Sim::run_simulation();
      ssim::Sim::clear();
    }
  }
  /**
      Define a class for the process
  */
  class IDM : public ssim::cProcess
  {
  public:
    int id;
    state_t state;
    Rcpp::List param;
    Report *report;
    ssim::gsm m1; 
    ssim::gsm m2;
    ssim::gsm m3;
    int treat; // treat variable indicating 1=RFC or 0=FC
    IDM(Rcpp::List param, Report *report) : id(-1), param(param), report(report) {
      // read in from param...
      m1 = ssim::gsm(as<List>(param(\"m1\")));
      m2 = ssim::gsm(as<List>(param(\"m2\")));
      m3 = ssim::gsm(as<List>(param(\"m3\")));
      treat = param(\"treat\");
      }
    void pfs();
    void init(); // to be specified
    void handleMessage(const ssim::cMessage* msg); 
    void cancelEvents(); 
  };
  void IDM::pfs() {
    state = PFS;
    report -> setUtility(0.8);
    scheduleAt(m1.randU(R::runif(0.0,1.0), now()),toProg);
    scheduleAt(m2.randU(R::runif(0.0,1.0), now()),toAcD);
  }
  /**
      Initialise a simulation run for an individual
  */
  void IDM::init() {
    id++;
    pfs();
  }
  /**
      Handle receiving self-messages
  */
  void IDM::handleMessage(const ssim::cMessage* msg) {
    if (param(\"debug\")) Rprintf(\"id: %i, state: %i, kind: %i, previous: %f, now: %f\\n\",
                       id, state, msg->kind, this->previousEventTime, ssim::now());
    report->add(state, msg->kind, this->previousEventTime, ssim::now(), id);
    cancel_events();
    scheduleAt(15, toEOS); // End of study--Time horizon 50 years
    switch(msg->kind) {
    case toPFS:
      pfs();
      break;
    case toProg:
      state = Prog;
      report -> setUtility(0.6);
      scheduleAt(m3.randU(R::runif(0.0,1.0)) + now(), toAcD); // Caution: clock-reset
      break;
    case toAcD:
      state = AcD;
      report -> setUtility(0.0);
      break;
    case toEOS:
      ssim::Sim::stop_simulation();
      break;
    default:
      REprintf(\"Invalid kind of event: %i.\\n\", msg->kind);
      break;
    }
    if (id % 100000 == 0) Rcpp::checkUserInterrupt(); 
  }
  /**
      Exported function: Set up the report and process, run the simulations and return a report
  */
  //[[Rcpp::export]]
  List callSim_reduced(int n, List param, bool indivp) {
    Report report(n,indivp);
    report.setPartition(0.0,100.0,param(\"partitionBy\"));
    report.setDiscountRate(param(\"discountRate\"));    
    IDM person(param,&report);
    runSimulations(&person, n);
    Rcpp::List lst = report.asList();
    lst.push_back(param,\"param\");
    return lst;
  }")

simulations_reduced = function(n, param, simulator=callSim_reduced, indivp=TRUE) {
  object = simulator(n, param, indivp)
  if(!all(c("ut","costs","pt","events","prev") %in% names(object)))
    stop("simulator failed")
  stateT = c("PFS", "Prog", "AcD")
  eventT = c("toPFS", "toProg", "toAcD", "toEOS")
  for (name in c("ut","costs","pt","events","prev"))
    object[[name]] = transform(object[[name]], state=stateT[Key+1], Key=NULL)
  object$events = transform(object$events, event=eventT[event+1])
  class(object) = c("IDM","SummaryReport")
  object
}

##############################################################
##============================================================
## Complete model
##============================================================
##############################################################
sourceCpp(code="
  // BH is required for older versions of microsimulation
  // RcppArmadillo is required for newer versions of microsimulation
  // [[Rcpp::depends(BH)]]
  // [[Rcpp::depends(RcppArmadillo)]]
  // [[Rcpp::depends(microsimulation)]]
  #include <RcppArmadillo.h>
  #include <microsimulation.h>
  // #include <gsm.h>
  using namespace arma;
  using namespace ssim;
  enum state_t {PFS, Prog, AcD};
  enum event_t {toPFS, toProg, toAcD, toEOS};
  typedef ssim::SummaryReport<short,short> Report;
  /**
      Utility: Run a set of simulations for a single process
  */
  void runSimulations(ssim::cProcess* process, int n) {
    for (int i = 0; i < n; i++) {
      ssim::Sim::create_process(process);
      ssim::Sim::run_simulation();
      ssim::Sim::clear();
    }
  }
  /**
      Define a class for the process
  */
  class IDM : public ssim::cProcess
  {
  public:
    int id;
    state_t state;
    Rcpp::List param;
    Report *report;
    ssim::gsm m1; 
    ssim::gsm m2;
    ssim::gsm m3;
    int treat; // treat variable indicating 1=RFC or 0=FC
    // Cost data
    // RFC
    double c_RFC_rituximab;
    double c_RFC_admin_rituximab;
    double c_RFC_fludarabine;
    double c_RFC_admin_fludarabine;
    double c_RFC_cyclophosphamide;
    double c_RFC_admin_cyclophosphamide;
    double c_RFC_support_PFS;
    double c_RFC_bone_marrow_transplantation;
    double c_RFC_blood_transfusions;
    double c_RFC_support_Prog;
    double c_RFC_therapy;
    double c_RFC_bmt;
    double c_RFC_blood;
    
    // FC 
    double c_FC_rituximab;
    double c_FC_admin_rituximab;
    double c_FC_fludarabine;
    double c_FC_admin_fludarabine;
    double c_FC_cyclophosphamide;
    double c_FC_admin_cyclophosphamide;
    double c_FC_support_PFS;
    double c_FC_bone_marrow_transplantation;
    double c_FC_blood_transfusions;
    double c_FC_support_Prog;
    double c_FC_therapy;
    double c_FC_bmt;
    double c_FC_blood;
    
    IDM(Rcpp::List param, Report *report) : id(-1), param(param), report(report) {
      // read in from param...
      m1 = ssim::gsm(as<List>(param(\"m1\")));
      m2 = ssim::gsm(as<List>(param(\"m2\")));
      m3 = ssim::gsm(as<List>(param(\"m3\")));
      treat = param(\"treat\");
      // Cost data
      // RFC
      c_RFC_rituximab = param(\"c_RFC_rituximab\");
      c_RFC_admin_rituximab = param(\"c_RFC_admin_rituximab\");
      c_RFC_fludarabine = param(\"c_RFC_fludarabine\");
      c_RFC_admin_fludarabine = param(\"c_RFC_admin_fludarabine\");
      c_RFC_cyclophosphamide = param(\"c_RFC_cyclophosphamide\");
      c_RFC_admin_cyclophosphamide = param(\"c_RFC_admin_cyclophosphamide\");
      c_RFC_support_PFS = param(\"c_RFC_support_PFS\");
      c_RFC_bone_marrow_transplantation = param(\"c_RFC_bone_marrow_transplantation\");
      c_RFC_blood_transfusions = param(\"c_RFC_blood_transfusions\");
      c_RFC_support_Prog= param(\"c_RFC_support_Prog\");
      c_RFC_therapy = param(\"c_RFC_therapy\");
      c_RFC_bmt = param(\"c_RFC_bmt\");
      c_RFC_blood = param(\"c_RFC_blood\");

      // FC
      c_FC_rituximab = param(\"c_FC_rituximab\");
      c_FC_admin_rituximab = param(\"c_FC_admin_rituximab\");
      c_FC_fludarabine = param(\"c_FC_fludarabine\");
      c_FC_admin_fludarabine = param(\"c_FC_admin_fludarabine\");
      c_FC_cyclophosphamide = param(\"c_FC_cyclophosphamide\");
      c_FC_admin_cyclophosphamide = param(\"c_FC_admin_cyclophosphamide\");
      c_FC_support_PFS = param(\"c_FC_support_PFS\");
      c_FC_bone_marrow_transplantation = param(\"c_FC_bone_marrow_transplantation\");
      c_FC_blood_transfusions = param(\"c_FC_blood_transfusions\");
      c_FC_support_Prog = param(\"c_FC_support_Prog\");
      c_FC_therapy = param(\"c_FC_therapy\");
      c_FC_bmt = param(\"c_FC_bmt\");
      c_FC_blood = param(\"c_FC_blood\");

      }
    void pfs();
    void init(); // to be specified
    void handleMessage(const ssim::cMessage* msg); // to be specified
    void cancelEvents(); // utility function
  };
  void IDM::pfs() {
    state = PFS;
    report -> setUtility(0.8);
    if (treat == 1) { // RFC
    report -> addPointCost(state, c_RFC_rituximab);
    report -> addPointCost(state, c_RFC_admin_rituximab);
    report -> addPointCost(state, c_RFC_fludarabine);
    report -> addPointCost(state, c_RFC_admin_fludarabine);
    report -> addPointCost(state, c_RFC_cyclophosphamide);
    report -> addPointCost(state, c_RFC_admin_cyclophosphamide);
    report -> setCost(c_RFC_support_PFS); 
    report -> addPointCost(state, c_RFC_bone_marrow_transplantation);
    report -> addPointCost(state, c_RFC_blood_transfusions);
    report -> addPointCost(state, c_RFC_bmt);
    report -> addPointCost(state, c_RFC_blood);
    }
    else { // FC
    report -> addPointCost(state, c_FC_rituximab);
    report -> addPointCost(state, c_FC_admin_rituximab);
    report -> addPointCost(state, c_FC_fludarabine);
    report -> addPointCost(state, c_FC_admin_fludarabine);
    report -> addPointCost(state, c_FC_cyclophosphamide);
    report -> addPointCost(state, c_FC_admin_cyclophosphamide);
    report -> setCost(c_FC_support_PFS); 
    report -> addPointCost(state, c_FC_bone_marrow_transplantation);
    report -> addPointCost(state, c_FC_blood_transfusions);
    report -> addPointCost(state, c_FC_bmt);
    report -> addPointCost(state, c_FC_blood);
    }
    scheduleAt(m1.randU(R::runif(0.0,1.0),now()),toProg);
    scheduleAt(m2.randU(R::runif(0.0,1.0),now()),toAcD);
  }
  /**
      Initialise a simulation run for an individual
  */
  void IDM::init() {
    id++;
    pfs();
  }
  /**
      Handle receiving self-messages
  */
  void IDM::handleMessage(const ssim::cMessage* msg) {
    if (param(\"debug\")) Rprintf(\"id: %i, state: %i, kind: %i, previous: %f, now: %f\\n\",
                       id, state, msg->kind, this->previousEventTime, ssim::now());
    report->add(state, msg->kind, this->previousEventTime, ssim::now(), id);
    cancel_events();
    scheduleAt(15, toEOS); // End of study--Time horizon 15 years
    switch(msg->kind) {
    case toPFS:
      pfs();
      break;
    case toProg:
      state = Prog;
      report -> setUtility(0.6);
      if (treat == 1) { // RFC
        report -> setCost(c_RFC_support_Prog + c_RFC_therapy); 
      }
      else{ // FC
        report -> setCost(c_FC_support_Prog + c_FC_therapy); 
      }
      scheduleAt(m3.randU(R::runif(0.0,1.0)) + now(), toAcD);
      break;
    case toAcD:
    case toEOS:
      ssim::Sim::stop_simulation();
      break;
    default:
      REprintf(\"Invalid kind of event: %i.\\n\", msg->kind);
      break;
    }
    if (id % 100000 == 0) Rcpp::checkUserInterrupt(); 
  }
  /**
      Exported function: Set up the report and process, run the simulations and return a report
  */
  //[[Rcpp::export]]
  List callSim(int n, List param, bool indivp) {
    Report report(n,indivp);
    report.setPartition(0.0,100.0,param(\"partitionBy\"));
    report.setDiscountRate(param(\"discountRate\"));    
    IDM person(param,&report);
    runSimulations(&person, n);
    Rcpp::List lst = report.asList();
    lst.push_back(param,\"param\");
    return lst;
  }")

simulations = function(n, param, simulator=callSim, indivp=TRUE) {
  object = simulator(n, param, indivp)
  if(!all(c("ut","costs","pt","events","prev") %in% names(object)))
    stop("simulator failed")
  stateT = c("PFS", "Prog", "AcD")
  eventT = c("toPFS", "toProg", "toAcD", "toEOS")
  for (name in c("ut","costs","pt","events","prev"))
    object[[name]] = transform(object[[name]], state=stateT[Key+1], Key=NULL)
  object$events = transform(object$events, event=eventT[event+1])
  class(object) = c("IDM","SummaryReport")
  object
}

##############################################################
##============================================================
## Begin bootstrapping for probabilistic sensitivity analysis
##============================================================
##############################################################
## 1=RFC, 0=FC
treat_values <- c(0, 1) # treat = 0, 1

## Define bootstrap function
bootstrap <- function(iterations){
              PSA_results <- NULL
              
              for (iteration in iterations) {
                
              ## Fix seed 
              seed <- as.integer(as.numeric(Sys.time())) + sample.int(1000, 1)
              set.seed(seed)
              
              results <- list()  # Create an empty list to store the results
              results <- lapply(treat_values, function(treat_value){
                ## Bootstrap survival models
                bootstrap <- for (i in 1:3) {
                  if(i==1){
                    model_name <- paste0("m", i, "_fpm")
                  }
                  else if(i==2){
                    model_name <- paste0("m", i, "_fpm_ac")
                  }  
                  else if(i==3){
                    model_name <- paste0("m", i, "_semiMarkov_fpm_ac")
                  }
                  model <- get(model_name)
                  modified_model <- model
                  modified_model@fullcoef <- drop(rmvnorm(1, coef(model), vcov(model)))
                  
                  # Create m1star to m3star
                  assign(paste0(model_name, "_star"), modified_model, envir = .GlobalEnv)
                }
                
                ## Run the simulations_reduced to determine disLY_Prog_RFC and disLY_Prog_FC
                for (i in c(0, 1)) {
                  ndata <- data.frame(treat = i)
                  param <- list(treat = i,
                                partitionBy = 0.1,
                                discountRate = 0.035,
                                debug = FALSE,
                                ## Survival models
                                m1 = gsm_design(m1_fpm_star, newdata = ndata),
                                m2 = gsm_design(m2_fpm_ac_star, newdata = ndata),
                                m3 = gsm_design(m3_semiMarkov_fpm_ac_star, newdata = ndata))
                  
                  sim <- simulations_reduced(1e5, param = param, indivp = FALSE)
                  
                  if (i == 0) {
                    disLY_Prog_FC <- (sum(subset(sim$ut, state %in% c("Prog"))$utility)) / sim$n
                  } else{
                    disLY_Prog_RFC <- (sum(subset(sim$ut, state %in% c("Prog"))$utility)) / sim$n
                  }
                }
                ## Discounted life years (rate = 0.035) spent at the progression state
                disLY_Prog_RFC
                disLY_Prog_FC
                
                ## PSA for cost parameters
                ## Williams 2017 listed the following two variables cost_admin_oral and cost_admin_complex,
                ## but they did not include these two costs into PSA.
                # cost_admin_oral <- rpert(nruns, 174, 280, 482)      # Administration—deliver exclusively oral chemotherapy
                # cost_admin_complex <- rpert(nruns, 210, 430, 795)   # Administration—deliver complex chemotherapy, including
                # prolonged infusional treatment at first attendance
                cost_PFS_support <- rpert(1,14,28,42)                 # Monthly supportive care cost while in PFS
                temp1<- 5179/((disLY_Prog_RFC+disLY_Prog_FC)/2*12)
                cost_prog_support_dis<- rpert(1,0.5*temp1,temp1,1.5*temp1)  # Monthly supportive care and second-line and subsequent
                # therapy cost while in progressiona
                cost_bmt <- rpert(1, 34318.25, 47565.05, 54646.47)    # Bone marrow transplant
                cost_transfusion <- rpert(1, 173.84, 289.73, 405.62)  # Blood transfusion
                cost_blood_unit <- rpert(1, 96.67, 161.11, 225.26)    # 1 Unit of blood
                
                ## Run the complete model
                ndata <- data.frame(treat = treat_value)
                param <- list(treat = treat_value,
                              partitionBy = 0.1,
                              discountRate = 0.035,
                              debug = FALSE,
                              ## Survival models
                              m1 = gsm_design(m1_fpm_star, newdata = ndata),
                              m2 = gsm_design(m2_fpm_ac_star, newdata = ndata),
                              m3 = gsm_design(m3_semiMarkov_fpm_ac_star, newdata = ndata),
                              ## Cost data
                              ## RFC
                              ## PFS
                              c_RFC_rituximab = 10113,                  # Costs of rituximab                              £10113
                              c_RFC_admin_rituximab = 1224,             # Administration costs of rituximab                £1224
                              c_RFC_fludarabine = 2776,                 # Cost of fludarabine                              £2776
                              c_RFC_admin_fludarabine = 1109,           # Administration costs of fludarabine              £1109
                              c_RFC_cyclophosphamide = 21,              # Costs of cyclophosphamide                          £21
                              c_RFC_admin_cyclophosphamide = 1109,      # Administration costs of cyclophosphamide         £1109
                              c_RFC_support_PFS = cost_PFS_support*12,  # PSA: Yearly costs of supportive care in PFS      £?*12
                              c_RFC_bone_marrow_transplantation = 592,  # Cost of bone marrow transplantation               £592
                              c_RFC_blood_transfusions = 640,           # Cost of blood transfusions                        £640
                              c_RFC_bmt = cost_bmt*5/402,               # PSA:  Cost of bone marrow transplant
                              c_RFC_blood = (cost_transfusion*318+cost_blood_unit*1025)/402, # PSA: Cost of blood
                              ## Prog
                              c_RFC_support_Prog = 84*12,               # Yearly costs of supportive care in Progression  £84*12
                              c_RFC_therapy      = cost_prog_support_dis*12,
                              # PSA: Yearly (monthly*12)  cost of RFC's second-line and subsequent therapy £5179 weighted by the length of stay
                              
                              ## FC
                              ## PFS
                              c_FC_rituximab = 0,                      # Costs of rituximab                                   £0
                              c_FC_admin_rituximab = 0,                # Administration costs of rituximab                    £0
                              c_FC_fludarabine = 2790,                 # Cost of fludarabine                               £2790
                              c_FC_admin_fludarabine = 1115,           # Administration costs of fludarabine               £1115
                              c_FC_cyclophosphamide = 22,              # Costs of cyclophosphamide                           £22
                              c_FC_admin_cyclophosphamide = 1115,      # Administration costs of cyclophosphamide          £1115
                              c_FC_support_PFS = 28*12,                # Yearly costs of supportive care in PFS             £28
                              c_FC_bone_marrow_transplantation = 360,  # Cost of bone marrow transplantation                £360
                              c_FC_blood_transfusions = 507,           # Cost of blood transfusions                         £507
                              c_FC_bmt = cost_bmt*3/396,               # PSA:  Cost of bone marrow transplant
                              c_FC_blood = (cost_transfusion*269+cost_blood_unit*762)/396, # PSA: Cost of blood
                              ## Prog
                              c_FC_support_Prog = 84*12,               # Yearly costs of supportive care in Progression.  £84*12
                              c_FC_therapy   = cost_prog_support_dis*12
                              # PSA: Yearly (monthly*12) cost of FC's second-line and subsequent therapy £5179 weighted by the length of stay
                )
                sim <- simulations(1e5, param = param, indivp = FALSE)
                results[[treat_value+1]] <- sim
              })
                results_RFC <- results[[2]]
                disQALY_RFC <- (sum(results_RFC$ut$utility)) / results_RFC$n
                disCost_RFC <- (sum(results_RFC$costs$cost)) / results_RFC$n
                        
                results_FC <- results[[1]]  
                disQALY_FC <- (sum(results_FC$ut$utility)) / results_FC$n
                disCost_FC <- (sum(results_FC$costs$cost)) / results_FC$n
                      
                disQALY_inc <- disQALY_RFC - disQALY_FC
                disCost_inc <- disCost_RFC - disCost_FC
                ICER <- disCost_inc / disQALY_inc
                        
                ## Create a new row of data with the resample number
                new_row <- data.frame(
                            BootstrapNumber = iteration,
                            disCost_RFC = disCost_RFC,
                            disCost_FC = disCost_FC,
                            disQALY_RFC = disQALY_RFC,
                            disQALY_FC = disQALY_FC,
                            disCost_inc = disCost_inc,
                            disQALY_inc = disQALY_inc,
                            ICER)
                PSA_results <- rbind(PSA_results, new_row)
              }
              return(PSA_results)
}

bootstrap(1:2)

## Estimate 10 bootstrap times
system.time(mclapply(1:10, function(i) {    
  bootstrap(i)
}, mc.cores = 6))

output <- mclapply(1:1000, function(i) {    
          bootstrap(i)
          }, mc.cores = 6)

output_df <- do.call(rbind, output)

## Save results
## Prevent from changing the results. We put # here.                        
save(output_df, file = "../Output/07_fullmod_fpm_ac_PSA.RData")
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
# limitations under the License.orating relative survival extrapolation and multiple timescales for health technology assessment © 2023 by Chen EYT, Dickman PW, Clements MS is licensed under CC BY 4.0

