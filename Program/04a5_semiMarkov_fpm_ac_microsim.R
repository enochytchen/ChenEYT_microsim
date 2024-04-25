## Filename: 04a5_extrap_microsim_fpm_ac
## Purpose: Run semi-Markov (clock-reset) model using microsimulation package with
##          flexible parametric models within an all-cause survival framework
## Notes: Must have installed the lastest Rtools and R4.2.2+

## Set seed for coherency
set.seed(12345)

##############################################################
##============================================================
##Read packages
##============================================================
##############################################################
library(tidyverse)
library(rstpm2)
library(microsimulation)
library(Rcpp)

##############################################################
##============================================================
## Run the survival models
##============================================================
##############################################################
source("03a_survmod_fpm.R")

##############################################################
##============================================================
## The complete microsimulation model without cost data and no discount
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
    ssim::Rpexp background;
    ssim::gsm m1; 
    ssim::gsm m2;
    ssim::gsm m3;
    int treat; // treat variable indicating 1=RFC or 0=FC
    IDM(Rcpp::List param, Report *report) : id(-1), param(param), report(report), background() {
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
    scheduleAt(50, toEOS); // End of study--Time horizon 50 years
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
    object[[name]] = transform(object[[name]], state=stateT[Key+1], time=age, Key=NULL, age=NULL)
  object$events = transform(object$events, event=eventT[event+1])
  class(object) = c("IDM","SummaryReport")
  object
}

##############################################################
##============================================================
## Survival extrapolation (undiscounted life years)
##============================================================
##############################################################
results <- list()  # Create an empty list to store the results
treat_values <- c(0, 1) # treat = 0, 1

results <- lapply(treat_values, function(treat_value) {
    ndata <- data.frame(treat = treat_value)
    param <- list(treat = treat_value,
                  partitionBy = 0.1,
                  discountRate = 0,
                  debug = FALSE,
                  ## Survival models
                  m1 = gsm_design(m1_fpm, newdata = ndata),
                  m2 = gsm_design(m2_fpm_ac, newdata = ndata),
                  m3 = gsm_design(m3_semiMarkov_fpm_ac, newdata = ndata))
    
    sim <- simulations_reduced(1e5, param = param, indivp = FALSE)
    return(sim)
  })

## Assign results
results_FC <- results[[1]]
results_RFC <- results[[2]]

## Save results
saveRDS(results_FC, file = "../Data/04a5_semiMarkov_fpm_ac_microsim_FC.rds")
saveRDS(results_RFC, file = "../Data/04a5_semiMarkov_fpm_ac_microsim_RFC.rds")
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
