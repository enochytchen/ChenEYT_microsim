## Filename: 07_rel_PSA
## Purpose: Run microsimulation model for probabilistic sensitivity analysis
##          Relative survival framework
## Caution: Must have installed the lastest Rtools and R4.2.2+

## Set the wd as where this R file is
if (require("rstudioapi") && isAvailable()) {
  original_wd <- getwd()  # Store the original working directory
  wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(wd)
}

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
library(mc2d) # for beta pert distribution
library(mvtnorm) # for bootstrapping parametric models 

##############################################################
##============================================================
## Run the other programs
##============================================================
##############################################################
source("02_msmcancer.R")
source("03_rel_survmod.R")

##############################################################
##============================================================
## Read popmort 
##============================================================
##############################################################
popmort <- readRDS("../Data/01b_popmort.rds")
names(popmort)

## Transform sex variable to stay coherent with the data
popmort2 = transform(popmort, cohort=year-age)
popmort2 = arrange(popmort2, cohort, sex, age)

## Check the cohort data
table(popmort2$cohort)

# Function for drawing a piecewise constant rate
# Birth cohort approach
# For example, male, dx 2003 and age 61 -> cohort = 1942
# getRates_cohort(1, 1942)
getRates_cohort = function(.sex, .cohort) { # .cohort = year-age
  subset(popmort2, sex==.sex & cohort==.cohort)$rate
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
  // #include <gsm.h>
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
    int treat; // var_value variable indicating 1=RFC or 0=FC
    IDM(Rcpp::List param, Report *report) : id(-1), param(param), report(report), background() {
       // get the rates from param
      vec mu0 = as<vec>(param(\"mu0\"));
      vec ages = arma::regspace(0.0,mu0.size()-1.0);
      background = ssim::Rpexp(mu0.memptr(), ages.memptr(), ages.size());
      // read in from param...
      m1 = ssim::gsm(as<List>(param(\"m1\")));
      m2 = ssim::gsm(as<List>(param(\"m2\")));
      m3 = ssim::gsm(as<List>(param(\"m3\")));
      dxage = param(\"dxage\");
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
    scheduleAt(m1.randU(R::runif(0.0,1.0),now()),toProg);
    scheduleAt(m2.randU(R::runif(0.0,1.0),now()),toExcD);
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
      report -> setUtility(0.6);
      scheduleAt(m3.randU(R::runif(0.0,1.0), now()), toExcD);
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
    
    IDM(Rcpp::List param, Report *report) : id(-1), param(param), report(report), background() {
       // get the rates from param
      vec mu0 = as<vec>(param(\"mu0\"));
      vec ages = arma::regspace(0.0,mu0.size()-1.0);
      background = ssim::Rpexp(mu0.memptr(), ages.memptr(), ages.size());
      // read in from param...
      m1 = ssim::gsm(as<List>(param(\"m1\")));
      m2 = ssim::gsm(as<List>(param(\"m2\")));
      m3 = ssim::gsm(as<List>(param(\"m3\")));
      dxage = param(\"dxage\");
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
    scheduleAt(m2.randU(R::runif(0.0,1.0),now()),toExcD);
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
      report -> setUtility(0.6);
      if (treat == 1) { // RFC
        report -> setCost(c_RFC_support_Prog + c_RFC_therapy); 
      }
      else{ // FC
        report -> setCost(c_FC_support_Prog + c_FC_therapy); 
      }
      scheduleAt(m3.randU(R::runif(0.0,1.0), now()), toExcD);
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
## Begin bootstrapping for probabilistic sensitivity analysis
##============================================================
##############################################################
## Number of bootstapping
m <- 1000

## 1=RFC, 0=FC
treat_values <- c(0, 1) # treat = 0, 1

## Create an empty dataset to save the results
PSA_results <- data.frame()

## Begin bootstrapping
for (iteration in 1:m) {
  results <- list()  # Create an empty list to store the results
  
  results <- lapply(treat_values, function(treat_value) {
  ## Fix seed 
  set.seed(as.numeric(Sys.time()))
  
  ## Bootstrap survival models
  bootstrap <- for (i in 1:3) {
    model_name <- paste0("m", i)
    model <- get(model_name)
    modified_model <- model
    modified_model@fullcoef <- drop(rmvnorm(1, coef(model), vcov(model)))
    
    # Create m1star to m3star
    assign(paste0(model_name, "star"), modified_model, envir = .GlobalEnv)
  }
  
  ## Run the simulations_reduced to determine disLY_Prog_RFC and disLY_Prog_FC
  for (i in c(0, 1)) {
        ndata <- data.frame(treat = i)
        param <- list(treat = i,
                      partitionBy = 0.1,
                      discountRate = 0.035,
                      debug = FALSE,
                      dxage = 61, 
                      ## Background rate
                      mu0 = getRates_cohort(3, 2005 - 61), # Draw a rate from sex 0/1, year 2003/2006 - age
                      ## Survival models
                      m1 = gsm_design(m1star, newdata = ndata),
                      m2 = gsm_design(m2star, newdata = ndata),
                      m3 = gsm_design(m3star, newdata = ndata))
        
        sim <- simulations_reduced(1e5, param = param, indivp = FALSE)

        if (i == 0) {
          disLY_Prog_FC <- (sum(subset(sim$pt, state %in% c("Prog"))$pt)) / sim$n
        } else{
          disLY_Prog_RFC <- (sum(subset(sim$pt, state %in% c("Prog"))$pt)) / sim$n
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
                dxage = 61, 
                ## Background rate
                mu0 = getRates_cohort(3, 2005 - 61), # Draw a rate from sex 0/1, year 2003/2006 - age
                ## Survival models
                m1 = gsm_design(m1star, newdata = ndata),
                m2 = gsm_design(m2star, newdata = ndata),
                m3 = gsm_design(m3star, newdata = ndata),
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
  return(sim)
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
      ICER
  )
  PSA_results <- rbind(PSA_results, new_row)
}

save(PSA_results, file = "../Data/07_rel_PSA.RData")
################################################################
setwd(original_wd)  # Reset to the original working directory
