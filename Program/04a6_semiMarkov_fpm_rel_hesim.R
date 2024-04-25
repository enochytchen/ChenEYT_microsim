## Filename: 04a6_semiMarkov_fpm_rel_hesim
## Purpose: Run semi-Markov (clock-reset) model using hesim package with
##          flexible parametric models within a relative survival framework
## Notes: Must have installed the lastest Rtools and R4.2.2+

library(data.table)
source("03a_survmod_fpm.R")
library(hesim)
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
  ## "hesim data" - invisibly EXPORTED
  hesim_dat <<- hesim_data(strategies = strategies,
                           patients = patients, 
                           states = states) 
  transmod_data = hesim::expand(hesim_dat, by = c("strategies", "patients"))
  ## Add some default columns to X
  if (!("cons" %in% names(transmod_data)))
    transmod_data[,cons := 1]
  if (!("(Intercept)" %in% names(transmod_data)))
    transmod_data[, "(Intercept)" := 1]
  return(transmod_data)
}
## for a full hesim model:
tmat = matrix(c(NA, 1, 2, 3,
                NA,NA, 4, 5,
                NA,NA,NA,NA,
                NA,NA,NA,NA), 4,4,TRUE)
colnames(tmat) = rownames(tmat) = c("PFS", "Prog", "ExpD", "ExcD")
n_patients = 1e5
as_dt_list <- function(x, n_samples=1)
  lapply(as.list(x), function(xi) if (is.numeric(xi))
    data.frame(cons=rep(xi,n_samples)) else xi)

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

mu0 = getRates_cohort(3, 2005 - 61)$rate # Draw a rate from sex 0/1, year 2003/2006 - age
ages = getRates_cohort(3, 2005 - 61)$age

background <- params_surv(coefs = as_dt_list(log(mu0)),
                          aux = list(time = ages),
                          dist = "pwexp")
transmod_data = prepare_input_data(n_patients = n_patients, tmat=tmat, 
                                   strategies = data.frame(strategy_id = 1:2)) # 1=RFC, 2=FC
transmod_data = transmod_data[,treat := 0][strategy_id==1, treat := 1]
transmod_params <- params_surv_list(create_params(m1_flexfpm_ac,uncertainty="none"),
                                    background,
                                    create_params(m2_flexfpm_rel,uncertainty="none"),
                                    background,
                                    create_params(m3_semiMarkov_flexfpm_rel,uncertainty="none"))
transmod <- create_IndivCtstmTrans(transmod_params, 
                                   input_data = transmod_data,
                                   trans_mat = tmat,
                                   clock = "mixt",
                                   transition_types=c("reset","age","reset","age","reset"),
                                   start_age = 61)
utility_tbl <- stateval_tbl(data.table(state_id = 1:2, est = c(0.8, 0.6)),
                            dist="fixed")
utilitymod <- create_StateVals(utility_tbl, hesim_data=hesim_dat, n=1)
## ## PFS
## c_PFS_pRFC = 17584  # Point cost of RFC arm at PFS state
## c_PFS_pFC  = 5909   # Point cost of FC arm at PFS state
## c_PFS_dRFC = 336    # Yearly cost of RFC arm at PFS state
## c_PFS_dFC  = 336    # Yearly cost of FC arm at PFS state
## ## Prog
## c_Prog_pRFC = 5179  # Point cost of RFC arm at Progression state
## c_Prog_pFC  = 5179  # Point cost of FC arm at Progression state
## c_Prog_dRFC = 1008  # Yearly cost of RFC arm at Progression state
## c_Prog_dFC  = 1008  # Yearly cost of FC arm at Progression state
point_costs_tbl <- stateval_tbl(data.table(state_id    = c(1,2,1,2),
                                           strategy_id = c(1,1,2,2),
                                           est = c(17584,5179,5909,5179)),
                                dist="fixed")
interval_costs_tbl <- stateval_tbl(data.table(state_id = 1:2, est = c(336,1008)),
                                   dist="fixed")
point_costs <- create_StateVals(point_costs_tbl, hesim_data=hesim_dat, n=1,
                                method="starting")
interval_costs <- create_StateVals(interval_costs_tbl, hesim_data=hesim_dat, n=1)
econmod <- IndivCtstm$new(trans_model = transmod,
                          utility_model = utilitymod,
                          cost_models = list(point_costs,interval_costs))
system.time(econmod$sim_disease(max_t = 15, max_age = 100))
head(econmod$disprog_)
econmod$sim_stateprobs(t = c(4,8,15))
econmod$stateprobs_
xtabs(prob~strategy_id+t,data=econmod$stateprobs_,subset=state_id %in% 1:2)
xtabs(prob~state_id+strategy_id,data=econmod$stateprobs_)
econmod$sim_qalys(dr = .035)
xtabs(qalys~strategy_id, data=econmod$qalys_)
econmod$sim_qalys(dr = 0)
econmod$qalys_
econmod$sim_costs(dr = 0.035)
econmod$costs_
xtabs(costs~strategy_id, data=econmod$costs_)

