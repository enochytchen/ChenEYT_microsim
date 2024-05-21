## Filename: 05_fullmod_semiMarkov_fpm_rel_microsim.R
## Purpose: Run the complete microsimulation model same as "04a6_semiMarkov_fpm_rel_microsim.R"
##          semi-Markov (clock-reset) model using microsimulation package with
##          flexible parametric models within a relative survival framework
##          with costs
## Notes: Must have installed the lastest Rtools and R4.2.2+

## Set seed for coherency
set.seed(12345)
##############################################################
##============================================================
## Install/Read packages
##============================================================
##############################################################
library(tidyverse)
library(rstpm2)
library(microsimulation)
library(Rcpp)
library(flextable)

##############################################################
##============================================================
## Run the survival models
##============================================================
##############################################################
source("03a_survmod_fpm.R")

##############################################################
##============================================================
## Read population mortality files 
##============================================================
##############################################################
## Read the popmort file
popmort <- readRDS("../Data/01b_popmort.rds")
names(popmort)

## Transform sex variable to stay coherent with the data
popmort2 = transform(popmort, cohort=year-age)
popmort2 = arrange(popmort2, cohort, sex, age)

## Check the cohort data
table(popmort2$cohort)

## Define a function for drawing a piecewise constant rate
## Birth cohort approach 
## getRates_cohort(1, 1942) means male, dx 2003 and age 61 -> cohort = 1942
getRates_cohort = function(.sex, .cohort) { # .cohort = year-age
  subset(popmort2, sex==.sex & cohort==.cohort)
}

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
  enum state_t {PFS, Prog, ExcD, ExpD};
  enum event_t {toPFS, toProg, toExcD, toExpD, toEOS};
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
    double dxage;
     // Utility 
    double u_PFS;
    double u_Prog;
    int treat; // var_value variable indicating 1=RFC or 0=FC
    IDM(Rcpp::List param, Report *report) : id(-1), param(param), report(report), background() {
       // get the rates from param
      vec mu0 = as<vec>(param(\"mu0\"));
      vec ages = as<vec>(param(\"ages\"));
      background = ssim::Rpexp(mu0.memptr(), ages.memptr(), ages.size());
      // read in from param...
      m1 = ssim::gsm(as<List>(param(\"m1\")));
      m2 = ssim::gsm(as<List>(param(\"m2\")));
      m3 = ssim::gsm(as<List>(param(\"m3\")));
      dxage = param(\"dxage\");
      treat = param(\"treat\");
      // Utility 
      u_PFS = param(\"u_PFS\");
      u_Prog = param(\"u_Prog\");
      }
    void pfs();
    void init(); // to be specified
    void handleMessage(const ssim::cMessage* msg); 
    void cancelEvents(); 
  };
  void IDM::pfs() {
    state = PFS;
    report -> setUtility(u_PFS);
    scheduleAt(m1.randU(R::runif(0.0,1.0), now()), toProg);
    scheduleAt(m2.randU(R::runif(0.0,1.0), now()), toExcD);
    scheduleAt(background.rand(std::exp(-R::rexp(1.0)), dxage+now())-dxage, toExpD);
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
      report -> setUtility(u_Prog);
      scheduleAt(m3.randU(R::runif(0.0,1.0)) + now(), toExcD);
      scheduleAt(background.rand(std::exp(-R::rexp(1.0)), dxage+now())-dxage, toExpD);
      break;
    case toExcD:
    case toExpD:
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
  stateT = c("PFS", "Prog", "ExcD", "ExpD")
  eventT = c("toPFS", "toProg", "toExcD", "toExpD", "toEOS")
  for (name in c("ut","costs","pt","events","prev"))
    object[[name]] = transform(object[[name]], state=stateT[Key+1], time=age, Key=NULL, age=NULL)
  object$events = transform(object$events, event=eventT[event+1])
  class(object) = c("IDM","SummaryReport")
  object
  }

##############################################################
##============================================================
## The complete microsimulation model begins here
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
  enum state_t {PFS, Prog, ExcD, ExpD};
  enum event_t {toPFS, toProg, toExcD, toExpD, toEOS};
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
    double dxage;
    int treat; // treat variable indicating 1=RFC or 0=FC
    // Utility 
    double u_PFS;
    double u_Prog;
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

    IDM(Rcpp::List param, Report *report) : id(-1), param(param), report(report), background() {
       // get the rates from param
      vec mu0 = as<vec>(param(\"mu0\"));
      vec ages = as<vec>(param(\"ages\"));
      background = ssim::Rpexp(mu0.memptr(), ages.memptr(), ages.size());
      m1 = ssim::gsm(as<List>(param(\"m1\")));
      m2 = ssim::gsm(as<List>(param(\"m2\")));
      m3 = ssim::gsm(as<List>(param(\"m3\")));
      dxage = param(\"dxage\");
      treat = param(\"treat\");
      // Utility 
      u_PFS = param(\"u_PFS\");
      u_Prog = param(\"u_Prog\");
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
      }
    void pfs();
    void init(); // to be specified
    void handleMessage(const ssim::cMessage* msg); // to be specified
    void cancelEvents(); // utility function
  };
  void IDM::pfs() {
    state = PFS;
    report -> setUtility(u_PFS);
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
    }
    scheduleAt(m1.randU(R::runif(0.0,1.0), now()),toProg);
    scheduleAt(m2.randU(R::runif(0.0,1.0), now()),toExcD);
    scheduleAt(background.rand(std::exp(-R::rexp(1.0)), dxage+now())-dxage, toExpD);
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
      report -> setUtility(u_Prog);
      if (treat == 1) { // RFC
        report -> setCost(c_RFC_support_Prog + c_RFC_therapy); 
      }
      else{ // FC
        report -> setCost(c_FC_support_Prog + c_FC_therapy); 
      }
      scheduleAt(m3.randU(R::runif(0.0,1.0)) + now(), toExcD);
      scheduleAt(background.rand(std::exp(-R::rexp(1.0)), dxage+now())-dxage, toExpD);
      break;
    case toExcD:
    case toExpD:
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
  stateT = c("PFS", "Prog", "ExcD", "ExpD")
  eventT = c("toPFS", "toProg", "toExcD", "toExpD", "toEOS")
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

## Run the simulations_reduced to determine disLY_Prog_RFC and disLY_Prog_FC
for (i in c(0, 1)) {
        ndata <- data.frame(treat = i)
        param <- list(treat = i,
                      partitionBy = 0.1,
                      discountRate = 0.035,
                      debug = FALSE,
                      dxage = 61, 
                      ## Background rate
                      mu0 = getRates_cohort(3, 2005 - 61)$rate, # Draw a rate from sex 0/1, year 2003/2006 - age
                      ages = getRates_cohort(3, 2005 - 61)$age,
                      ## Survival models
                      m1 = gsm_design(m1_fpm, newdata = ndata),
                      m2 = gsm_design(m2_fpm_rel, newdata = ndata),
                      m3 = gsm_design(m3_semiMarkov_fpm_rel, newdata = ndata),
                      u_PFS = 1,    ## For estimating discounted LY
                      u_Prog = 1)   ## For estimating discounted LY

        sim <- simulations_reduced(1e6, param = param, indivp = FALSE)
        
        if (i == 0) { 
           disLY_FC <- sum(sim$ut$utility) / sim$n
           disLY_PFS_FC <- (sum(subset(sim$ut, state %in% c("PFS"))$utility)) / sim$n
           disLY_Prog_FC <- (sum(subset(sim$ut, state %in% c("Prog"))$utility)) / sim$n
        } else{
           disLY_RFC <- sum(sim$ut$utility) / sim$n
           disLY_PFS_RFC <- (sum(subset(sim$ut, state %in% c("PFS"))$utility)) / sim$n
           disLY_Prog_RFC <- (sum(subset(sim$ut, state %in% c("Prog"))$utility)) / sim$n
        }
      }
      
      ## Discounted life years (rate = 0.035) spent at the progression state
      disLY_Prog_RFC
      disLY_Prog_FC

results <- lapply(treat_values, function(treat_value) {      
      ## Run the complete model
      ndata <- data.frame(treat = treat_value)
      param <- list(treat = treat_value,
                    partitionBy = 0.1,
                    discountRate = 0.035,
                    debug = FALSE,
                    dxage = 61, 
                    ## Background rate
                    mu0 = getRates_cohort(3, 2005 - 61)$rate, # Draw a rate from sex 0/1, year 2003/2006 - age
                    ages = getRates_cohort(3, 2005 - 61)$age,
                    ## Survival models
                    m1 = gsm_design(m1_fpm, newdata = ndata),
                    m2 = gsm_design(m2_fpm_rel, newdata = ndata),
                    m3 = gsm_design(m3_semiMarkov_fpm_rel, newdata = ndata),
                    u_PFS = 0.8,
                    u_Prog = 0.6,
                    ## Cost data
                    ## RFC
                    ## PFS
                    c_RFC_rituximab = 10113,                  # Costs of rituximab                              £10113
                    c_RFC_admin_rituximab = 1224,             # Administration costs of rituximab                £1224
                    c_RFC_fludarabine = 2776,                 # Cost of fludarabine                              £2776
                    c_RFC_admin_fludarabine = 1109,           # Administration costs of fludarabine              £1109
                    c_RFC_cyclophosphamide = 21,              # Costs of cyclophosphamide                          £21
                    c_RFC_admin_cyclophosphamide = 1109,      # Administration costs of cyclophosphamide         £1109
                    c_RFC_support_PFS = 28*12,                # Yearly costs of supportive care in PFS          £28*12
                    c_RFC_bone_marrow_transplantation = 592,  # Cost of bone marrow transplantation               £592
                    c_RFC_blood_transfusions = 640,           # Cost of blood transfusions                        £640
                    ## Prog
                    c_RFC_support_Prog = 84*12,               # Yearly costs of supportive care in Progression  £84*12
                    c_RFC_therapy      = (5179/((disLY_Prog_RFC + disLY_Prog_FC)/2*12))*12,
                    # Yearly (monthly*12)  cost of RFC's second-line and subsequent therapy £5179 weighted by the length of stay
                    
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
                    ## Prog
                    c_FC_support_Prog = 84*12,               # Yearly costs of supportive care in Progression.  £84*12
                    c_FC_therapy   = (5179/((disLY_Prog_RFC + disLY_Prog_FC)/2*12))*12
                    # Yearly (monthly*12) cost of FC's second-line and subsequent therapy £5179 weighted by the length of stay
                    )
      sim <- simulations(1e6, param = param, indivp = FALSE)
      return(sim)
})

##############################################################
##============================================================
## Base-case analysis results
##============================================================
##############################################################

###########
### RFC ###
###########
results_RFC <- results[[2]]
results_RFC$param$treat     # Confirm it is treat = 1

## Undiscounted life years
undisLY_RFC <- (sum(results_RFC$pt$pt)) / results_RFC$n
## Undiscounted life years spent at PFS
undisLY_PFS_RFC <- (sum(subset(results_RFC$pt, state %in% c("PFS"))$pt)) / results_RFC$n
## Undiscounted life years spent at prog
undisLY_Prog_RFC <- (sum(subset(results_RFC$pt, state %in% c("Prog"))$pt)) / results_RFC$n

## Discounted life years
disLY_RFC
## Discounted life years spent at PFS
disLY_PFS_RFC
## Discounted life years spent at prog
disLY_Prog_RFC

## Discounted QALYs
disQALY_RFC <- (sum(results_RFC$ut$utility)) / results_RFC$n
## Discounted QALYs spent at PFS
disQALY_PFS_RFC <- (sum(subset(results_RFC$ut, state %in% c("PFS"))$utility)) / results_RFC$n
## Discounted QALYs spent at Prog
disQALY_Prog_RFC <- (sum(subset(results_RFC$ut, state %in% c("Prog"))$utility)) / results_RFC$n

## Discounted costs
disCost_RFC <- (sum(results_RFC$costs$cost)) / results_RFC$n
## Discounted costs spent at PFS
disCost_PFS_RFC <- (sum(subset(results_RFC$costs, state %in% c("PFS"))$cost)) / results_RFC$n
## Discounted costs spent at Prog
disCost_Prog_RFC <-  (sum(subset(results_RFC$costs, state %in% c("Prog"))$cost)) / results_RFC$n


##########
### FC ###
##########
results_FC <- results[[1]]  
results_FC$param$treat      # Confirm it is treat = 0

## Undiscounted life years
undisLY_FC <- (sum(results_FC$pt$pt) - sum(subset(results_FC$pt, state %in% c("ExpD", "ExcD"))$pt)) / results_FC$n
## Undiscounted life years spent at PFS
undisLY_PFS_FC <- (sum(subset(results_FC$pt, state %in% c("PFS"))$pt)) / results_FC$n
## Undiscounted life years spent at prog
undisLY_Prog_FC <- (sum(subset(results_FC$pt, state %in% c("Prog"))$pt)) / results_FC$n

## Discounted life years
disLY_FC 
## Discounted life years spent at PFS
disLY_PFS_FC
## Discounted life years spent at prog
disLY_Prog_FC

## Discounted QALYs
disQALY_FC <- (sum(results_FC$ut$utility)) / results_FC$n
## Discounted QALYs spent at PFS
disQALY_PFS_FC <- (sum(subset(results_FC$ut, state %in% c("PFS"))$utility)) / results_FC$n
## Discounted QALYs spent at Prog
disQALY_Prog_FC <- (sum(subset(results_FC$ut, state %in% c("Prog"))$utility)) / results_FC$n

## Discounted costs
disCost_FC <- (sum(results_FC$costs$cost)) / results_FC$n
## Discounted costs spent at PFS
disCost_PFS_FC <- (sum(subset(results_FC$costs, state %in% c("PFS"))$cost)) / results_FC$n
## Discounted costs spent at Prog
disCost_Prog_FC <-  (sum(subset(results_FC$costs, state %in% c("Prog"))$cost)) / results_FC$n

###################
### Incremental ###
###################
## Incremental life years
disLY_inc <- disLY_RFC - disLY_FC
## Incremental life years spent at PFS
disLY_PFS_inc <- disLY_PFS_RFC - disLY_PFS_FC
## Incremental life years spent at Prog
disLY_Prog_inc <- disLY_Prog_RFC - disLY_Prog_FC

## Incremental QALYs
disQALY_inc <- disQALY_RFC - disQALY_FC
## Incremental QALYs spent at PFS
disQALY_PFS_inc <- disQALY_PFS_RFC - disQALY_PFS_FC
## Incremental QALYs spent at Prog
disQALY_Prog_inc <- disQALY_Prog_RFC - disQALY_Prog_FC

## Incremental costs 
disCost_inc <- disCost_RFC - disCost_FC
## Incremental costs spent at PFS
disCost_PFS_inc <- disCost_PFS_RFC - disCost_PFS_FC
## Incremental costs spent at Prog
disCost_Prog_inc <- disCost_Prog_RFC - disCost_Prog_FC

## Cost per LY gained
disCost_per_LY <- disCost_inc / disLY_inc

## Cost per QALY gained
disCost_per_QALY <- disCost_inc / disQALY_inc

##############################################################
##============================================================
## Export table
##============================================================
##############################################################
## Round up
# Get a list of all variables in the global environment
global_vars <- ls()

# Filter and round variables as needed
for (var_name in global_vars) {
  var <- get(var_name)
  
  if (startsWith(var_name, "disLY") || startsWith(var_name, "disQALY")) {
    # Round to two decimal places
    var <-  format(round(var, digits = 2), nsmall = 2)
    
  } else if (startsWith(var_name, "disCost")) {
    # Round to zero decimal places
    var <- paste("£", format(round(var, digits = 0), nsmall = 0))
  }
  # Update the variable in the global environment
  assign(var_name, var)
}

## Arrange table
var_matrix <-  matrix(c("Mean life years", disLY_RFC, disLY_FC, disLY_inc,
                        "Mean life years PFS", disLY_PFS_RFC, disLY_PFS_FC, disLY_PFS_inc,
                        "Mean life years progression", disLY_Prog_RFC, disLY_Prog_FC, disLY_Prog_inc,
                        "Mean QALYs", disQALY_RFC, disQALY_FC, disQALY_inc,
                        "Mean QALYs PFS", disQALY_PFS_RFC, disQALY_PFS_FC, disQALY_PFS_inc,
                        "Mean QALYs progression", disQALY_Prog_RFC, disQALY_Prog_FC, disQALY_Prog_inc,
                        "Mean total costs", disCost_RFC, disCost_FC, disCost_inc,
                        "Mean costs of PFS", disCost_PFS_RFC, disCost_PFS_FC, disCost_PFS_inc,
                        "Mean costs of progression", disCost_Prog_RFC, disCost_Prog_FC, disCost_Prog_inc,
                        "Cost per life year gained", NA, NA, disCost_per_LY,
                        "Cost per QALY gained", NA, NA, disCost_per_QALY), nrow = 11, byrow = TRUE)
## Make data frame
df <- as.data.frame(var_matrix)

## Set colnames and rownames
colnames(df)[1:4] <- c("Parameter", "RFC", "FC", "Incremental")

## Export as docx
## Make a flextable and export it as docx 
ft <- flextable(df)  
ft <- autofit(ft)
ft <- align(ft, align = c("right"), i = NULL, j = 2:4)
ft <- add_footer_lines(ft,"FC, fludarabine and cyclophosphamide; RFC, rituximab, fludarabine, and cyclophosphamide; QALY, quality-adjusted life year; PFS, progression-free survival.")
ft <- set_caption(ft,  as_paragraph(as_b("Appendix C. Base case analysis results of the semi-Markov model using flexible parametric models within a relative survival framework (semi-Markov: FPMs, RSF).")))
ft

## Save results
## Prevent from changing the results. We put # here.
# save_as_docx(ft, path="../Output/05_fullmod_semiMarkov_fpm_rel_microsim.docx")

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

