## Filename: 04a7_semiMarkov_williams_ac_microsim
## Purpose: Run semi-Markov (clock-reset) model using microsimulation package with
##          flexible parametric models within an all-cause survival framework
##          Models (m1: gompertz; m2: g-gamma; m3: gompertz) same as the base-case
##          analysis as Williams et al. 2017. 
##          Using all-cause survival framework
## Notes: Must have installed the lastest Rtools and R4.2.2+
## Reference: Williams C, Lewsey JD, Briggs AH, Mackay DF. 
##            Cost-effectiveness Analysis in R Using a Multi-state Modeling Survival Analysis Framework: A Tutorial. 
##            Med Decis Making. 2017 May;37(4):340–52.
##            Cost-effectiveness Analysis in R Using a Multi-state Modeling Survival Analysis Framework: A Tutorial © 2017 by Williams C. et al is licensed under CC BY 3.0. 

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
library(hesim)
library(data.table)

##############################################################
##============================================================
## Run the survival models
##============================================================
##############################################################
source("03a_survmod_williams.R")

##############################################################
##============================================================
## Hesim implementation 
## (m1: gompertz; m2: g-gamma; m3: gompertz) 
## Microsimulation within an all-cause survival framework
##============================================================
##############################################################
sourceCpp(code="
  // BH is required for older versions of microsimulation
  // RcppArmadillo is required for newer versions of microsimulation
  // [[Rcpp::depends(BH)]]
  // [[Rcpp::depends(hesim)]]
  // [[Rcpp::depends(RcppArmadillo)]]
  // [[Rcpp::depends(microsimulation)]]
  #include <RcppArmadillo.h>
  #include <microsimulation.h>
  #include <hesim.h>
  #include <hesim/ctstm/ctstm.h>
  using namespace arma;
  using namespace ssim;
  enum state_t {PFS, Prog, AcD};
  enum event_t {toPFS, toProg, toAcD, toEOS};
  typedef ssim::SummaryReport<short,short> Report;
  /**
      Utility: Random exponential using rate parameterisation
  */
  double rexpRate(double rate) { return R::rexp(1.0/rate); }
  // template<class T> double rexpRate(T rate) { return R::rexp(1.0/as<double>(rate)); }
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
  class IDMHesim : public ssim::cProcess
  {
  public:
    int id;
    state_t state;
    Rcpp::List param;
    Report *report;
    std::unique_ptr<hesim::ctstm::transmod> m1;  
    std::unique_ptr<hesim::ctstm::transmod> m2;  
    std::unique_ptr<hesim::ctstm::transmod> m3;  
    int treat; // treat variable indicating 1=RFC or 0=FC
    IDMHesim(Rcpp::List param, Report *report) : id(-1), param(param), report(report) {
       // get the rates from param
      treat = param(\"treat\");
      m1 = hesim::ctstm::transmod::create(as<Rcpp::Environment>(param(\"m1\")));
      m1 -> obs_index_.set_patient_index(0);
      m1 -> obs_index_.set_strategy_index(treat);
      m2 = hesim::ctstm::transmod::create(as<Rcpp::Environment>(param(\"m2\")));
      m2 -> obs_index_.set_patient_index(0);
      m2 -> obs_index_.set_strategy_index(treat);
      m3 = hesim::ctstm::transmod::create(as<Rcpp::Environment>(param(\"m3\")));
      m3 -> obs_index_.set_patient_index(0);
      m3 -> obs_index_.set_strategy_index(treat);
      }
    void pfs();
    void init(); // to be specified
    void handleMessage(const ssim::cMessage* msg); // to be specified
    void cancelEvents(); // utility function
  };
  void IDMHesim::pfs() {
    state = PFS;
    report -> setUtility(0.8);
    scheduleAt(m1->random(0,0) + now(), toProg); // do we need to change the sample?
    scheduleAt(m2->random(0,0) + now(), toAcD);}
  /**
      Initialise a simulation run for an individual
  */
  void IDMHesim::init() {
    id++;
    pfs();
  }
  /**
      Handle receiving self-messages
  */
  void IDMHesim::handleMessage(const ssim::cMessage* msg) {
    if (param(\"debug\")) Rprintf(\"id: %i, state: %i, kind: %i, previous: %f, now: %f\\n\",
                       id, state, msg->kind, this->previousEventTime, ssim::now());
    report->add(state, msg->kind, this->previousEventTime, ssim::now(), id);
    cancel_events();
    scheduleAt(50, toEOS); // End of study--Time horizon 15 years
    switch(msg->kind) {
    case toPFS:
      pfs();
      break;
    case toProg:
      state = Prog;
      report -> setUtility(0.6);
      scheduleAt(m3->random(0,0) + now(), toAcD);
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
  List callSimHesim(int n, List param, bool indivp) {
    Report report(n,indivp);
    report.setPartition(0.0,100.0,param(\"partitionBy\"));
    report.setDiscountRate(param(\"discountRate\"));    
    IDMHesim person(param,&report);
    runSimulations(&person, n);
    Rcpp::List lst = report.asList();
    lst.push_back(param,\"param\");
    return lst;
  }")

simulations = function(n, param, simulator=callSimHesim, indivp=TRUE) {
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
## Set parameters
##============================================================
##############################################################
results <- list()  # Create an empty list to store the results
treat_values <- c(0, 1) # treat = 0, 1

## Hesim transition models
prepare_input_data = function(strategies=NULL,
                              n_patients=NULL, patients=NULL,
                              tmat=NULL, states=NULL, state_names=NULL, ...) {
    stopifnot(!(is.null(n_patients) && is.null(patients)),
              !(is.null(tmat) && is.null(states) && is.null(state_names)),
              is.null(patients) || "patient_id" %in% names(patients))
    ## strategies
    if (is.null(strategies))
        strategies = data.frame(strategy_id = 1)
    ## patients
    if (is.null(patients))
        patients = data.frame(patient_id = 1:n_patients, ...)
    ## states
    if (is.null(states) && !is.null(tmat))
        states = data.frame( # Non-death health states
            state_id = 1:nrow(tmat),
            state_name = colnames(tmat)
        )
    if (is.null(states) && !is.null(state_names))
        states = data.frame( # Non-death health states
            state_id = 1:length(state_names),
            state_name = state_names
        )
    ## Organize data to read in Hesim
    hesim_dat = hesim_data(strategies = strategies,
                           patients = patients, 
                           states = states)
    transmod_data = hesim::expand(hesim_dat, by = c("strategies", "patients"))
    ## Add some default columns to X
    if (!("cons" %in% names(transmod_data)))
        transmod_data$cons = 1
    if (!("(Intercept)" %in% names(transmod_data)))
        transmod_data[["(Intercept)"]] = 1
    return(transmod_data)
}
tmat = matrix(c(NA,1, NA,NA), 2,2,TRUE)
colnames(tmat) = rownames(tmat) = c("Base", "Next")
transmod_data = prepare_input_data(n_patients = 1, tmat=tmat, 
                                   strategies = data.frame(strategy_id = 1:2), treat = 1)
transmod_data = transmod_data[strategy_id==1, treat := 0]
transmod_params1 <- params_surv_list(create_params(m1_gom,uncertainty="none"))
transmod1 <- create_IndivCtstmTrans(transmod_params1, 
                                     input_data = transmod_data,
                                     trans_mat = tmat,
                                     clock = "reset")
transmod_params2 <- params_surv_list(create_params(m2_gam,uncertainty="none"))
transmod2 <- create_IndivCtstmTrans(transmod_params2, 
                                     input_data = transmod_data,
                                     trans_mat = tmat,
                                     clock = "reset")
transmod_params3 <- params_surv_list(create_params(m3_semiMarkov_gom, uncertainty="none"))
transmod3 <- create_IndivCtstmTrans(transmod_params3, 
                                     input_data = transmod_data,
                                     trans_mat = tmat,
                                     clock = "reset")

##############################################################
##============================================================
## Run the microsimulation model
##============================================================
##############################################################
results<- lapply(treat_values, function(treat_value) {
                              ndata <- data.frame(treat = treat_value)
                              param <- list(treat = treat_value,
                                            partitionBy = 0.1,
                                            discountRate = 0.035,
                                            debug = FALSE,
                                            dxage = 61, 
                                            ## Survival models
                                            m1 = transmod1,
                                            m2 = transmod2,
                                            m3 = transmod3
                                            )
                              sim <- simulations(1e5, param = param, simulator=callSimHesim, indivp = FALSE)
                              return(sim)
})

## Assign results
results_FC <- results[[1]]
results_RFC <- results[[2]]

## Save results
saveRDS(results_FC, file = "../Data/04a7_semiMarkov_williams_ac_microsim_FC.rds")
saveRDS(results_RFC, file = "../Data/04a7_semiMarkov_williams_ac_microsim_RFC.rds")

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

