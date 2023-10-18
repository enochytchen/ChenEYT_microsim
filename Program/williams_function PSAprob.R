##########################################################################################
##========================================================================================
## START OF PSAPROB FUNCTION
##========================================================================================
##########################################################################################

##========================================================================================
## THIS FUNCTION RETURNS STATE OCCUPANCY PROBABILITIES FOR EACH DRAW

## MODIFIBLE ARGUMENTS:
## ntrans       NUMBER OF TRANSITIONS IN THE MODELLING
## ncovs        NUMBER OF COVARIATE PARAMETERS IN THE MODEL (EXCLUDING INTERCEPT) FOR EACH
##              TRANSITION I.E. ONE FOR EACH BINARY AND CONTINUOUS VARIABLE AND K-1 FOR  
##              EACH CATEGORICAL VARIABLE, WHERE K= NO OF CATEGORIES
## covs         VARIABLE NAMES FOR THE COVARIATES
## coveval      VALUE AT WHICH TO EVALUATE EACH COVARIATE
## dist         DISTRIBUTION TO USE FOR EACH TRANSITION. OPTIONS ARE:
##              wei FOR WEIBULL, exp FOR EXPONENTIAL, gom FOR GOMPERTZ,
##              logl FOR LOGLOGISTIC,logn FOR LOGNORMAL AND gam FOR GENERALISED GAMMA.
##              IF FITTING THE SAME MODEL OVER THE OBSERVED AND EXTRAPOLATED PERIOD 
##              DISTRIBUTION WILL SPAN THE WHOLE TIME HORIZON. IF FITTING A DIFFERENT 
##              MODEL OVER THE OBSERVED PERIOD THAN FOR THE EXTRAPOLATION, DISTRIBUTION
##              WILL SPAN THE OBSERVED PERIOD ONLY.
## dist2        ONLY APPLICABLE IF FITTING A DIFFERENT MODEL OVER THE EXTRAPOLATION PERIOD
##              THAN THE OBSERVED PERIOD.DISTRIBUTION TO USE FOR EACH TRANSITION OVER THE
##              EXTRAPOLATION PERIOD. OPTIONS ARE:
##              wei FOR WEIBULL, exp FOR EXPONENTIAL, gom FOR GOMPERTZ,
##              logl FOR LOGLOGISTIC,logn FOR LOGNORMAL AND gam FOR GENERALISED GAMMA.
## timeseq      THE TIME POINTS TO USE FOR PREDICTIONS OVER THE OBSERVED PERIOD OF THE STUDY.
##              THE FIRST ARGUMENT OF seq SHOULD BE THE START TIME, THE SECOND ARGUMENT THE 
##              END TIME AND THE THIRD ARGUMENT THE TIME INCREMENT. 
## timeseq_ext  THE TIME POINTS TO USE FOR PREDICTIONS OVER THE EXTRAPOLATION PERIOD.
##              THE FIRST ARGUMENT OF seq SHOULD BE THE START TIME, THE SECOND ARGUMENT THE 
##              END TIME AND THE THIRD ARGUMENT THE TIME INCREMENT. 
## data         DATASET TO USE FOR MODELLING
## trans        TRANSITION MATRIX
## Markov       EITHER TRUE OR FALSE, FOR A MARKOV OR SEMI-MARKOV BASE CASE MODEL RESPECTIVELY
## nruns        NUMBER OF DRAWS FOR THE PROBABILISTIC SENSITIVITY ANALYSIS
## seedno       NUMBER TO USE TO SET THE RANDOM NUMBER GENERATOR SO THAT SIMULATIONS CAN
##              BE REPLICATED EXACTLY. IF NOT REQUIRED SET TO NULL
## M            NUMBER OF SIMULATIONS USED TO CALCULATE PROBABILITIES 
## addin        WHETHER TO SUBSEQUENTLY ADD BACK IN ANY MORE TIME POINTS. EITHER TRUE OR FALSE.
## addvalue     ONLY APPLICABLE IF addin = TRUE. ROW TO ADD BACK IN.
##========================================================================================

PSAprob<-function(ntrans=3, ncovs=c(1,1,2), 
               covs=rbind("treat", "treat", c("treat", "lemedian")) ,
               coveval=rbind(0,0,c(0,0)), 
               dist=cbind("wei", "wei", "wei"),
               dist2=cbind(NA, NA, NA),
               timeseq=seq(0,4,1/12),
               timeseq_ext=seq(49/12,15,1/12),
               data=msmcancer, trans=tmat, Markov=FALSE,
               nruns=1000,seedno=12345, M=100){

  #### set up required lists 
  models<-vector("list", ntrans)
  fmla<-vector("list", ntrans)
  covars<-vector("list", ntrans)
  datasub<-vector("list", ntrans)
  lp<-vector("list", ntrans)
  coeffs<-vector("list", ntrans)
  var<-vector("list", ntrans)
  tvar<-vector("list", ntrans)
  x<-vector("list", ntrans)
  xnam1<-vector("list", ntrans) 
  xnam2<-vector("list", ntrans) 
  xnam3<-vector("list", ntrans) 
  xnam<-vector("list", ntrans) 
  ses<-vector("list", ntrans) 
  CV<-vector("list", ntrans) 
  norm<-vector("list", ntrans) 
  gammapar<-vector("list", ntrans)
  lambdapar<-vector("list", ntrans)
  gammapar_ext<-vector("list", ntrans)
  lambdapar_ext<-vector("list", ntrans)
  kappapar<-vector("list", ntrans) 
  mupar<-vector("list", ntrans) 
  sigmapar<-vector("list", ntrans) 
  zpar<-vector("list", ntrans) 
  upar<-vector("list", ntrans) 
  value<-vector("list", ntrans)
  param<-vector("list", ntrans)
  temp<-vector("list", ntrans)
  H<-vector("list", ntrans)
  H_ext<-vector("list", ntrans)
  #### create the timepoints
  tt2<-timeseq_ext
  for (i in 1:ntrans) {
    
    if (is.na(dist2[i])==TRUE) tt<-c(timeseq,timeseq_ext)
    if (is.na(dist2[i])==FALSE) tt<-timeseq
  
  ##### set up models

    covars[[i]]<-covs[i,1:ncovs[i]]
    if (Markov==TRUE){
    fmla[[i]]<-as.formula(paste("Surv(Tstart,Tstop,status)~ ",paste(covars[[i]], collapse= "+"))) 
    }
    if (Markov==FALSE){
      fmla[[i]]<-as.formula(paste("Surv(time,status)~ ",paste(covars[[i]], collapse= "+"))) 
    }
    datasub[[i]]<-subset(data,trans==i) 
    
    if (dist[i]=="wei") {
      ##### base case models  
      coeffs[[i]]<-phreg(fmla[[i]],dist="weibull", data=datasub[[i]])$coeff 
      var[[i]]<-phreg(fmla[[i]],dist="weibull", data=datasub[[i]])$var
      
      ##### transform covariance matrices
      ##### transform covariates, log(scale), log(shape) to
      #####            covariates, intercept, log(shape) where 
      #####            intercept= -shape*log(scale)
      ses[[i]]=FALSE
      ##covariates
      xnam1[[i]] <- paste0("~x", 1:ncovs[i], collapse=", " )
      ##intercept
      xnam2[[i]]<-paste0("~ -exp(","x", ncovs[i]+2,")","*","x", ncovs[i]+1 )
      ##log(shape)
      xnam3[[i]]<-paste0("~x", ncovs[i]+2)
      xnam[[i]]<-paste0("list(", paste0(c(xnam1[[i]],xnam2[[i]],xnam3[[i]]),collapse=", " ), ")")
      tvar[[i]]<- deltamethod(eval(parse(text=xnam[[i]])), coeffs[[i]], var[[i]], ses=ses[[i]] )
      
      ##### Cholesky decomposition of transformed covariance matrix
      CV[[i]]<-chol(tvar[[i]])
      ##### start of simulations
      ##### draws from normal distribution
      set.seed(seedno+i)
      norm[[i]]<-matrix(data=rnorm((ncovs[i]+2)*nruns),ncol=ncovs[i]+2,nrow=nruns)
      ##### create random variables based on each row of the covariance matrix
      ##### for each draw
      for (j in 1:(ncovs[i]+2)){
        if (j!=(ncovs[i]+1)){
      param[[i]][[j]]<-apply(norm[[i]],1,function(x)
        coeffs[[i]][j]+ sum(x*CV[[i]][,j])) 
        }
        if (j==(ncovs[i]+1)){
          param[[i]][[j]]<-apply(norm[[i]],1,function(x)
            -exp(coeffs[[i]][(ncovs[i]+2)])*coeffs[[i]][(ncovs[i]+1)] + sum(x*CV[[i]][,j])) 
        }
      }
      ##### create gamma (shape) parameter from each draw
      gammapar[[i]]<-exp(param[[i]][[(ncovs[i]+2)]]) 
      ##### create lambda parameter from each draw 
      value[[i]]<-coveval[i,1:ncovs[i]]
      for (j in 1:ncovs[i]){
        temp[[i]][[j]]<-param[[i]][[j]]*value[[i]][j]
      }
      lambdapar[[i]]<-exp(param[[i]][[(ncovs[i]+1)]] + Reduce("+", temp[[i]]))
      ##### create the timepoints
      ##### create the cumulative hazard matrix
      H[[i]]<-matrix(ncol=nruns, nrow=length(tt))
      for (k in 1:nruns) {
        for (m in 1:length(tt)){
          H[[i]][m,k]<-lambdapar[[i]][k]*tt[m]^gammapar[[i]][k]
       }}  
    }
    
    if (is.na(dist2[i])==FALSE & dist2[i] =="wei") {
      ##### base case models  
      coeffs[[i]]<-phreg(fmla[[i]],dist="weibull", data=datasub[[i]])$coeff 
      var[[i]]<-phreg(fmla[[i]],dist="weibull", data=datasub[[i]])$var
      
      ##### transform covariance matrices
      ##### transform covariates, log(scale), log(shape) to
      #####            covariates, intercept, log(shape) where 
      #####            intercept= -shape*log(scale)
      ses[[i]]=FALSE
      ##covariates
      xnam1[[i]] <- paste0("~x", 1:ncovs[i], collapse=", " )
      ##intercept
      xnam2[[i]]<-paste0("~ -exp(","x", ncovs[i]+2,")","*","x", ncovs[i]+1 )
      ##log(shape)
      xnam3[[i]]<-paste0("~x", ncovs[i]+2)
      xnam[[i]]<-paste0("list(", paste0(c(xnam1[[i]],xnam2[[i]],xnam3[[i]]),collapse=", " ), ")")
      tvar[[i]]<- deltamethod(eval(parse(text=xnam[[i]])), coeffs[[i]], var[[i]], ses=ses[[i]] )
      
      ##### Cholesky decomposition of transformed covariance matrix
      CV[[i]]<-chol(tvar[[i]])
      ##### start of simulations
      ##### draws from normal distribution
      set.seed(seedno+i)
      norm[[i]]<-matrix(data=rnorm((ncovs[i]+2)*nruns),ncol=ncovs[i]+2,nrow=nruns)
      ##### create random variables based on each row of the covariance matrix
      ##### for each draw
      for (j in 1:(ncovs[i]+2)){
        if (j!=(ncovs[i]+1)){
          param[[i]][[j]]<-apply(norm[[i]],1,function(x)
            coeffs[[i]][j]+ sum(x*CV[[i]][,j])) 
        }
        if (j==(ncovs[i]+1)){
          param[[i]][[j]]<-apply(norm[[i]],1,function(x)
            -exp(coeffs[[i]][(ncovs[i]+2)])*coeffs[[i]][(ncovs[i]+1)] + sum(x*CV[[i]][,j])) 
        }
      }
      ##### create gamma (shape) parameter from each draw
      gammapar[[i]]<-exp(param[[i]][[(ncovs[i]+2)]]) 
      ##### create lambda parameter from each draw 
      value[[i]]<-coveval[i,1:ncovs[i]]
      for (j in 1:ncovs[i]){
        temp[[i]][[j]]<-param[[i]][[j]]*value[[i]][j]
      }
      lambdapar[[i]]<-exp(param[[i]][[(ncovs[i]+1)]] + Reduce("+", temp[[i]]))
      ##### create the timepoints
      ##### create the cumulative hazard matrix
      H_ext[[i]]<-matrix(ncol=nruns, nrow=length(tt2))
      for (k in 1:nruns) {
        for (m in 1:length(tt2)){
          H_ext[[i]][m,k]<-lambdapar[[i]][k]*tt2[m]^gammapar[[i]][k]
        }}  
    }
    
    if (dist[i]=="exp") {
      ##### base case models  
      coeffs[[i]]<-phreg(fmla[[i]],dist="weibull", shape=1, data=datasub[[i]])$coeff 
      var[[i]]<-phreg(fmla[[i]],dist="weibull", shape=1, data=datasub[[i]])$var
      
      ##### transform covariance matrices
      ##### transform covariates, log(scale) to
      #####            covariates, intercept,where 
      #####            intercept= -log(scale)
      ses[[i]]=FALSE
      ##covariates
      xnam1[[i]] <- paste0("~x", 1:ncovs[i], collapse=", " )
      ##intercept
      xnam2[[i]]<-paste0("~ -x", ncovs[i]+1 )
      xnam[[i]]<-paste0("list(", paste0(c(xnam1[[i]],xnam2[[i]]),collapse=", " ), ")")
      tvar[[i]]<- deltamethod(eval(parse(text=xnam[[i]])), coeffs[[i]], var[[i]], ses=ses[[i]] )
      
      ##### Cholesky decomposition of transformed covariance matrix
      CV[[i]]<-chol(tvar[[i]])
      ##### start of simulations
      ##### draws from normal distribution
      set.seed(seedno+i)
      norm[[i]]<-matrix(data=rnorm((ncovs[i]+1)*nruns),ncol=ncovs[i]+1,nrow=nruns)
      ##### create random variables based on each row of the covariance matrix
      ##### for each draw
      for (j in 1:(ncovs[i])){
        param[[i]][[j]]<-apply(norm[[i]],1,function(x)
          coeffs[[i]][j]+ sum(x*CV[[i]][,j])) 
      }
      if (j==(ncovs[i]+1)){
        param[[i]][[j]]<-apply(norm[[i]],1,function(x)
          -coeffs[[i]][j] + sum(x*CV[[i]][,j])) 
      }
    
    ##### create lambda parameter from each draw 
    value[[i]]<-coveval[i,1:ncovs[i]]
    for (j in 1:ncovs[i]){
      temp[[i]][[j]]<-param[[i]][[j]]*value[[i]][j]
    }
    lambdapar[[i]]<-exp(param[[i]][[(ncovs[i]+1)]] + Reduce("+", temp[[i]]))
    ##### create the cumulative hazard matrix
    H[[i]]<-matrix(ncol=nruns, nrow=length(tt))
    for (k in 1:nruns) {
      for (m in 1:length(tt)){
        H[[i]][m,k]<-lambdapar[[i]][k]*tt[m]
      }}  
    }
    
    if (is.na(dist2[i])==FALSE & dist2[i] =="exp") {
      ##### base case models  
      coeffs[[i]]<-phreg(fmla[[i]],dist="weibull", shape=1, data=datasub[[i]])$coeff 
      var[[i]]<-phreg(fmla[[i]],dist="weibull", shape=1, data=datasub[[i]])$var
      
      ##### transform covariance matrices
      ##### transform covariates, log(scale) to
      #####            covariates, intercept,where 
      #####            intercept= -log(scale)
      ses[[i]]=FALSE
      ##covariates
      xnam1[[i]] <- paste0("~x", 1:ncovs[i], collapse=", " )
      ##intercept
      xnam2[[i]]<-paste0("~ -x", ncovs[i]+1 )
      xnam[[i]]<-paste0("list(", paste0(c(xnam1[[i]],xnam2[[i]]),collapse=", " ), ")")
      tvar[[i]]<- deltamethod(eval(parse(text=xnam[[i]])), coeffs[[i]], var[[i]], ses=ses[[i]] )
      
      ##### Cholesky decomposition of transformed covariance matrix
      CV[[i]]<-chol(tvar[[i]])
      ##### start of simulations
      ##### draws from normal distribution
      set.seed(seedno+i)
      norm[[i]]<-matrix(data=rnorm((ncovs[i]+1)*nruns),ncol=ncovs[i]+1,nrow=nruns)
      ##### create random variables based on each row of the covariance matrix
      ##### for each draw
      for (j in 1:(ncovs[i])){
        param[[i]][[j]]<-apply(norm[[i]],1,function(x)
          coeffs[[i]][j]+ sum(x*CV[[i]][,j])) 
      }
      if (j==(ncovs[i]+1)){
        param[[i]][[j]]<-apply(norm[[i]],1,function(x)
          -coeffs[[i]][j] + sum(x*CV[[i]][,j])) 
      }
      
      ##### create lambda parameter from each draw 
      value[[i]]<-coveval[i,1:ncovs[i]]
      for (j in 1:ncovs[i]){
        temp[[i]][[j]]<-param[[i]][[j]]*value[[i]][j]
      }
      lambdapar[[i]]<-exp(param[[i]][[(ncovs[i]+1)]] + Reduce("+", temp[[i]]))
      ##### create the cumulative hazard matrix
      H_ext[[i]]<-matrix(ncol=nruns, nrow=length(tt2))
      for (k in 1:nruns) {
        for (m in 1:length(tt2)){
          H_ext[[i]][m,k]<-lambdapar[[i]][k]*tt2[m]
        }}  
    }

    
  if (dist[i]=="gom") {
    ##### base case models  
    coeffs[[i]]<-phreg(fmla[[i]],dist="gompertz",  param="rate", data=datasub[[i]])$coeff 
    var[[i]]<-phreg(fmla[[i]],dist="gompertz",  param="rate", data=datasub[[i]])$var
    #####transformation of covariance matrix not required 
    ##### Cholesky decomposition of covariance matrix
    CV[[i]]<-chol(var[[i]])
    ##### start of simulations
    ##### draws from normal distribution
    set.seed(seedno+i)
    norm[[i]]<-qnorm(matrix(data=runif((ncovs[i]+2)*nruns),ncol=ncovs[i]+2,nrow=nruns))
    ##### create random variables based on each row of the covariance matrix
    ##### for each draw
    for (j in 1:(ncovs[i]+2)){
      param[[i]][[j]]<-apply(norm[[i]],1,function(x)
        coeffs[[i]][j]+ sum(x*CV[[i]][,j])) 
    }
    ##### create gamma (shape) parameter from each draw
    gammapar[[i]]<-(param[[i]][[(ncovs[i]+1)]]) 
    ##### create lambda parameter from each draw 
    value[[i]]<-coveval[i,1:ncovs[i]]
    for (j in 1:ncovs[i]){
      temp[[i]][[j]]<-param[[i]][[j]]*value[[i]][j]
    }
    lambdapar[[i]]<-exp(param[[i]][[(ncovs[i]+2)]] + Reduce("+", temp[[i]]))
    ##### create the cumulative hazard matrix
    H[[i]]<-matrix(ncol=nruns, nrow=length(tt))
    for (k in 1:nruns) {
      for (m in 1:length(tt)){
        H[[i]][m,k]<-lambdapar[[i]][k]*(1/gammapar[[i]][k])*(exp(gammapar[[i]][k]*tt[m])-1)       
      }}  
  }
  
    if (is.na(dist2[i])==FALSE & dist2[i] =="gom") {
      ##### base case models  
      coeffs[[i]]<-phreg(fmla[[i]],dist="gompertz",  param="rate", data=datasub[[i]])$coeff 
      var[[i]]<-phreg(fmla[[i]],dist="gompertz",  param="rate", data=datasub[[i]])$var
      #####transformation of covariance matrix not required 
      ##### Cholesky decomposition of covariance matrix
      CV[[i]]<-chol(var[[i]])
      ##### start of simulations
      ##### draws from normal distribution
      set.seed(seedno+i)
      norm[[i]]<-qnorm(matrix(data=runif((ncovs[i]+2)*nruns),ncol=ncovs[i]+2,nrow=nruns))
      ##### create random variables based on each row of the covariance matrix
      ##### for each draw
      for (j in 1:(ncovs[i]+2)){
        param[[i]][[j]]<-apply(norm[[i]],1,function(x)
          coeffs[[i]][j]+ sum(x*CV[[i]][,j])) 
      }
      ##### create gamma (shape) parameter from each draw
      gammapar[[i]]<-(param[[i]][[(ncovs[i]+1)]]) 
      ##### create lambda parameter from each draw 
      value[[i]]<-coveval[i,1:ncovs[i]]
      for (j in 1:ncovs[i]){
        temp[[i]][[j]]<-param[[i]][[j]]*value[[i]][j]
      }
      lambdapar[[i]]<-exp(param[[i]][[(ncovs[i]+2)]] + Reduce("+", temp[[i]]))
      ##### create the cumulative hazard matrix
      H_ext[[i]]<-matrix(ncol=nruns, nrow=length(tt2))
      for (k in 1:nruns) {
        for (m in 1:length(tt2)){
          H_ext[[i]][m,k]<-lambdapar[[i]][k]*(1/gammapar[[i]][k])*(exp(gammapar[[i]][k]*tt2[m])-1)       
        }}  
    }
    
  if (dist[i]=="logl") {
    ##### base case models  
    coeffs[[i]]<-aftreg(fmla[[i]],dist="loglogistic", data=datasub[[i]])$coeff 
    var[[i]]<-aftreg(fmla[[i]],dist="loglogistic", data=datasub[[i]])$var
    
    ##### transform covariance matrices
    ##### transform covariates, log(scale), log(shape) to
    #####            -covariates, intercept, -log(shape) where 
    #####            intercept= log(scale)
    ses[[i]]=FALSE
    ##covariates
    xnam1[[i]] <- paste0("~ -x", 1:ncovs[i], collapse=", " )
    ##intercept
    xnam2[[i]]<-paste0("~x", ncovs[i]+1)
    ##log(shape)
    xnam3[[i]]<-paste0("~x", ncovs[i]+2)
    xnam[[i]]<-paste0("list(", paste0(c(xnam1[[i]],xnam2[[i]],xnam3[[i]]),collapse=", " ), ")")
    tvar[[i]]<- deltamethod(eval(parse(text=xnam[[i]])), coeffs[[i]], var[[i]], ses=ses[[i]] ) 
    ##### Cholesky decomposition of covariance matrix
    CV[[i]]<-chol(var[[i]])
    ##### start of simulations
    ##### draws from normal distribution
    set.seed(seedno+i)
    norm[[i]]<-matrix(data=rnorm((ncovs[i]+2)*nruns),ncol=ncovs[i]+2,nrow=nruns)
    ##### create random variables based on each row of the covariance matrix
    ##### for each draw
    for (j in 1:(ncovs[i]+2)){
      if (j!=(ncovs[i]+1)){
        param[[i]][[j]]<-apply(norm[[i]],1,function(x)
          -coeffs[[i]][j]+ sum(x*CV[[i]][,j])) 
      }
      if (j==(ncovs[i]+1)){
        param[[i]][[j]]<-apply(norm[[i]],1,function(x)
          coeffs[[i]][j]+ sum(x*CV[[i]][,j])) 
      }
    }
    ##### create gamma (shape) parameter from each draw
    gammapar[[i]]<-exp(param[[i]][[(ncovs[i]+2)]]) 
    ##### create lambda parameter from each draw 
    value[[i]]<-coveval[i,1:ncovs[i]]
    for (j in 1:ncovs[i]){
      temp[[i]][[j]]<-param[[i]][[j]]*value[[i]][j]
    }
    lambdapar[[i]]<-exp(param[[i]][[(ncovs[i]+1)]] + Reduce("+", temp[[i]]))
    ##### create the cumulative hazard matrix
    H[[i]]<-matrix(ncol=nruns, nrow=length(tt))
    for (k in 1:nruns) {
      for (m in 1:length(tt)){
        H[[i]][m,k]<--log(1/(1+lambdapar[[i]][k]*tt[m]^(1/gammapar[[i]][k])))
      }}  
  }
    
    
    if (is.na(dist2[i])==FALSE & dist2[i] =="logl") {
      ##### base case models  
      coeffs[[i]]<-aftreg(fmla[[i]],dist="loglogistic", data=datasub[[i]])$coeff 
      var[[i]]<-aftreg(fmla[[i]],dist="loglogistic", data=datasub[[i]])$var
      
      ##### transform covariance matrices
      ##### transform covariates, log(scale), log(shape) to
      #####            -covariates, intercept, -log(shape) where 
      #####            intercept= log(scale)
      ses[[i]]=FALSE
      ##covariates
      xnam1[[i]] <- paste0("~ -x", 1:ncovs[i], collapse=", " )
      ##intercept
      xnam2[[i]]<-paste0("~x", ncovs[i]+1)
      ##log(shape)
      xnam3[[i]]<-paste0("~x", ncovs[i]+2)
      xnam[[i]]<-paste0("list(", paste0(c(xnam1[[i]],xnam2[[i]],xnam3[[i]]),collapse=", " ), ")")
      tvar[[i]]<- deltamethod(eval(parse(text=xnam[[i]])), coeffs[[i]], var[[i]], ses=ses[[i]] ) 
      ##### Cholesky decomposition of covariance matrix
      CV[[i]]<-chol(var[[i]])
      ##### start of simulations
      ##### draws from normal distribution
      set.seed(seedno+i)
      norm[[i]]<-matrix(data=rnorm((ncovs[i]+2)*nruns),ncol=ncovs[i]+2,nrow=nruns)
      ##### create random variables based on each row of the covariance matrix
      ##### for each draw
      for (j in 1:(ncovs[i]+2)){
        if (j!=(ncovs[i]+1)){
          param[[i]][[j]]<-apply(norm[[i]],1,function(x)
            -coeffs[[i]][j]+ sum(x*CV[[i]][,j])) 
        }
        if (j==(ncovs[i]+1)){
          param[[i]][[j]]<-apply(norm[[i]],1,function(x)
            coeffs[[i]][j]+ sum(x*CV[[i]][,j])) 
        }
      }
      ##### create gamma (shape) parameter from each draw
      gammapar[[i]]<-exp(param[[i]][[(ncovs[i]+2)]]) 
      ##### create lambda parameter from each draw 
      value[[i]]<-coveval[i,1:ncovs[i]]
      for (j in 1:ncovs[i]){
        temp[[i]][[j]]<-param[[i]][[j]]*value[[i]][j]
      }
      lambdapar[[i]]<-exp(param[[i]][[(ncovs[i]+1)]] + Reduce("+", temp[[i]]))
      ##### create the cumulative hazard matrix
      H_ext[[i]]<-matrix(ncol=nruns, nrow=length(tt2))
      for (k in 1:nruns) {
        for (m in 1:length(tt2)){
          H_ext[[i]][m,k]<--log(1/(1+lambdapar[[i]][k]*tt2[m]^(1/gammapar[[i]][k])))
        }}  
    }    
    
  if (dist[i]=="logn") {
    ##### base case models  
    coeffs[[i]]<-aftreg(fmla[[i]],dist="lognormal", data=datasub[[i]])$coeff 
    var[[i]]<-aftreg(fmla[[i]],dist="lognormal", data=datasub[[i]])$var
    
    ##### transform covariance matrices
    ##### transform covariates, log(scale), log(shape) to
    #####            -covariates, intercept, -log(shape) where 
    #####            intercept= log(scale)
    ses[[i]]=FALSE
    ##covariates
    xnam1[[i]] <- paste0("~ -x", 1:ncovs[i], collapse=", " )
    ##intercept
    xnam2[[i]]<-paste0("~x", ncovs[i]+1)
    ##log(shape)
    xnam3[[i]]<-paste0("~x", ncovs[i]+2)
    xnam[[i]]<-paste0("list(", paste0(c(xnam1[[i]],xnam2[[i]],xnam3[[i]]),collapse=", " ), ")")
    tvar[[i]]<- deltamethod(eval(parse(text=xnam[[i]])), coeffs[[i]], var[[i]], ses=ses[[i]] )
    
    ##### Cholesky decomposition of transformed covariance matrix
    CV[[i]]<-chol(tvar[[i]])
    ##### start of simulations
    ##### draws from normal distribution
    set.seed(seedno+i)
    norm[[i]]<-matrix(data=rnorm((ncovs[i]+2)*nruns),ncol=ncovs[i]+2,nrow=nruns)
    ##### create random variables based on each row of the covariance matrix
    ##### for each draw
    for (j in 1:(ncovs[i]+2)){
      if (j!=(ncovs[i]+1)){
        param[[i]][[j]]<-apply(norm[[i]],1,function(x)
          -coeffs[[i]][j]+ sum(x*CV[[i]][,j])) 
      }
      if (j==(ncovs[i]+1)){
        param[[i]][[j]]<-apply(norm[[i]],1,function(x)
          coeffs[[i]][j]+ sum(x*CV[[i]][,j])) 
      }
    }
    ##### create sigma (shape) parameter from each draw
    sigmapar[[i]]<-exp(param[[i]][[(ncovs[i]+2)]]) 
    ##### create lambda parameter from each draw 
    value[[i]]<-coveval[i,1:ncovs[i]]
    for (j in 1:ncovs[i]){
      temp[[i]][[j]]<-param[[i]][[j]]*value[[i]][j]
    }
    mupar[[i]]<-exp(param[[i]][[(ncovs[i]+1)]] + Reduce("+", temp[[i]]))
    ##### create the cumulative hazard matrix
    H[[i]]<-matrix(ncol=nruns, nrow=length(tt))
    for (k in 1:nruns) {
      for (m in 1:length(tt)){
        H[[i]][m,k]<- -log(1-pnorm((log(tt[m]) -mupar[[i]][k])/sigmapar[[i]][k]))
      }}  
  }
    
    if (is.na(dist2[i])==FALSE & dist2[i] =="logn") {
      ##### base case models  
      coeffs[[i]]<-aftreg(fmla[[i]],dist="lognormal", data=datasub[[i]])$coeff 
      var[[i]]<-aftreg(fmla[[i]],dist="lognormal", data=datasub[[i]])$var
      
      ##### transform covariance matrices
      ##### transform covariates, log(scale), log(shape) to
      #####            -covariates, intercept, -log(shape) where 
      #####            intercept= log(scale)
      ses[[i]]=FALSE
      ##covariates
      xnam1[[i]] <- paste0("~ -x", 1:ncovs[i], collapse=", " )
      ##intercept
      xnam2[[i]]<-paste0("~x", ncovs[i]+1)
      ##log(shape)
      xnam3[[i]]<-paste0("~x", ncovs[i]+2)
      xnam[[i]]<-paste0("list(", paste0(c(xnam1[[i]],xnam2[[i]],xnam3[[i]]),collapse=", " ), ")")
      tvar[[i]]<- deltamethod(eval(parse(text=xnam[[i]])), coeffs[[i]], var[[i]], ses=ses[[i]] )
      
      ##### Cholesky decomposition of transformed covariance matrix
      CV[[i]]<-chol(tvar[[i]])
      ##### start of simulations
      ##### draws from normal distribution
      set.seed(seedno+i)
      norm[[i]]<-matrix(data=rnorm((ncovs[i]+2)*nruns),ncol=ncovs[i]+2,nrow=nruns)
      ##### create random variables based on each row of the covariance matrix
      ##### for each draw
      for (j in 1:(ncovs[i]+2)){
        if (j!=(ncovs[i]+1)){
          param[[i]][[j]]<-apply(norm[[i]],1,function(x)
            -coeffs[[i]][j]+ sum(x*CV[[i]][,j])) 
        }
        if (j==(ncovs[i]+1)){
          param[[i]][[j]]<-apply(norm[[i]],1,function(x)
            coeffs[[i]][j]+ sum(x*CV[[i]][,j])) 
        }
      }
      ##### create sigma (shape) parameter from each draw
      sigmapar[[i]]<-exp(param[[i]][[(ncovs[i]+2)]]) 
      ##### create lambda parameter from each draw 
      value[[i]]<-coveval[i,1:ncovs[i]]
      for (j in 1:ncovs[i]){
        temp[[i]][[j]]<-param[[i]][[j]]*value[[i]][j]
      }
      mupar[[i]]<-exp(param[[i]][[(ncovs[i]+1)]] + Reduce("+", temp[[i]]))
      ##### create the cumulative hazard matrix
      H_ext[[i]]<-matrix(ncol=nruns, nrow=length(tt2))
      for (k in 1:nruns) {
        for (m in 1:length(tt2)){
          H_ext[[i]][m,k]<- -log(1-pnorm((log(tt2[m]) -mupar[[i]][k])/sigmapar[[i]][k]))
        }}  
    }  
    
    
    
  if (dist[i]=="gam") {
    ##### base case models  
    coeffs[[i]]<-flexsurvreg(fmla[[i]],dist="gengamma",data=datasub[[i]])$res[,1]
    var[[i]]<-flexsurvreg(fmla[[i]],dist="gengamma", data=datasub[[i]])$cov
    
    ##### no transformation of covariance matrix required
    ##### Cholesky decomposition of covariance matrix
    CV[[i]]<-chol(var[[i]])
    ##### start of simulations
    ##### draws from normal distribution
    set.seed(seedno+i)
    norm[[i]]<-matrix(data=rnorm((ncovs[i]+3)*nruns),ncol=ncovs[i]+3,nrow=nruns)
    ##### create random variables based on each row of the covariance matrix
    ##### for each draw
    for (j in 1:(ncovs[i]+3)){
      if (j!=(ncovs[i]+1)){
        param[[i]][[j]]<-apply(norm[[i]],1,function(x)
          coeffs[[i]][j]+ sum(x*CV[[i]][,j])) 
      }
      if (j==(ncovs[i]+1)){
        param[[i]][[j]]<-apply(norm[[i]],1,function(x)
          log(coeffs[[i]][j])+ sum(x*CV[[i]][,j])) 
      } 
    }
    
    ##### create sigma parameter from each draw
    sigmapar[[i]]<-exp(param[[i]][[2]]) 
    ##### create kappa parameter for each draw
    kappapar[[i]]<-param[[i]][[3]] 
    ##### create gamma parameter from each draw
    gammapar[[i]]<-(abs(kappapar[[i]]))^(-2)
    ##### create mu parameter from each draw
    value[[i]]<-coveval[i,1:ncovs[i]]
    
    for (j in 1:ncovs[i]){
      temp[[i]][[j]]<-param[[i]][[j+3]]*value[[i]][j]
    }
    mupar[[i]]<-exp(param[[i]][[1]] + Reduce("+", temp[[i]]))
    ##### create the z parameter and u parameter from each draw 
    zpar[[i]]<-matrix(ncol=nruns, nrow=length(tt))
    upar[[i]]<-matrix(ncol=nruns, nrow=length(tt))
    for (k in 1:nruns) {
      for (m in 1:length(tt)){
        zpar[[i]][m,k]<- sign(kappapar[[i]][k])*((log(tt[m])-mupar[[i]][k])/sigmapar[[i]][k])
        upar[[i]][m,k]<-gammapar[[i]][k]*exp(abs(kappapar[[i]][k])*zpar[[i]][m,k])
      }}  
    
    #### create the cumulative hazard matrix
    H[[i]]<-matrix(ncol=nruns, nrow=length(tt))
    for (k in 1:nruns) {
      for (m in 1:length(tt)){
        if(kappapar[[i]][k]>0){
          H[[i]][m,k]<- -log(1-pgamma(upar[[i]][m,k],gammapar[[i]][k]))
        }
        if(kappapar[[i]][k]==0){
          H[[i]][m,k]<- -log(1-pnorm(zpar[[i]][m,k]))
        }
        if(kappapar[[i]][k]<0){
          H[[i]][m,k]<- -log(pgamma(upar[[i]][m,k],gammapar[[i]][k]))
        }
      }}  
  }

    if (is.na(dist2[i])==FALSE & dist2[i] =="gam") {
      ##### base case models  
      coeffs[[i]]<-flexsurvreg(fmla[[i]],dist="gengamma",data=datasub[[i]])$res[,1]
      var[[i]]<-flexsurvreg(fmla[[i]],dist="gengamma", data=datasub[[i]])$cov
      
      ##### no transformation of covariance matrix required
      ##### Cholesky decomposition of covariance matrix
      CV[[i]]<-chol(var[[i]])
      ##### start of simulations
      ##### draws from normal distribution
      set.seed(seedno+i)
      norm[[i]]<-matrix(data=rnorm((ncovs[i]+3)*nruns),ncol=ncovs[i]+3,nrow=nruns)
      ##### create random variables based on each row of the covariance matrix
      ##### for each draw
      for (j in 1:(ncovs[i]+3)){
        if (j!=(ncovs[i]+1)){
          param[[i]][[j]]<-apply(norm[[i]],1,function(x)
            coeffs[[i]][j]+ sum(x*CV[[i]][,j])) 
        }
        if (j==(ncovs[i]+1)){
          param[[i]][[j]]<-apply(norm[[i]],1,function(x)
            log(coeffs[[i]][j])+ sum(x*CV[[i]][,j])) 
        } 
      }
      
      ##### create sigma parameter from each draw
      sigmapar[[i]]<-exp(param[[i]][[2]]) 
      ##### create kappa parameter for each draw
      kappapar[[i]]<-param[[i]][[3]] 
      ##### create gamma parameter from each draw
      gammapar[[i]]<-(abs(kappapar[[i]]))^(-2)
      ##### create mu parameter from each draw
      value[[i]]<-coveval[i,1:ncovs[i]]
      
      for (j in 1:ncovs[i]){
        temp[[i]][[j]]<-param[[i]][[j+3]]*value[[i]][j]
      }
      mupar[[i]]<-exp(param[[i]][[1]] + Reduce("+", temp[[i]]))
      ##### create the z parameter and u parameter from each draw 
      zpar[[i]]<-matrix(ncol=nruns, nrow=length(tt2))
      upar[[i]]<-matrix(ncol=nruns, nrow=length(tt2))
      for (k in 1:nruns) {
        for (m in 1:length(tt2)){
          zpar[[i]][m,k]<- sign(kappapar[[i]][k])*((log(tt2[m])-mupar[[i]][k])/sigmapar[[i]][k])
          upar[[i]][m,k]<-gammapar[[i]][k]*exp(abs(kappapar[[i]][k])*zpar[[i]][m,k])
        }}  
      
      #### create the cumulative hazard matrix
      H_ext[[i]]<-matrix(ncol=nruns, nrow=length(tt2))
      for (k in 1:nruns) {
        for (m in 1:length(tt2)){
          if(kappapar[[i]][k]>0){
            H_ext[[i]][m,k]<- -log(1-pgamma(upar[[i]][m,k],gammapar[[i]][k]))
          }
          if(kappapar[[i]][k]==0){
            H_ext[[i]][m,k]<- -log(1-pnorm(zpar[[i]][m,k]))
          }
          if(kappapar[[i]][k]<0){
            H_ext[[i]][m,k]<- -log(pgamma(upar[[i]][m,k],gammapar[[i]][k]))
          }
        }}  
    }    
    
if (is.na(dist2[i])==FALSE) H[[i]]<-rbind(H[[i]], H_ext[[i]])
    
    
  }  
  tt<-c(timeseq,timeseq_ext)
    ##### create the cumulative hazard matrix continued
      Haz<-matrix(ncol=nruns, nrow=length(tt)*ntrans)
      Haz<-do.call("rbind", H)
      newtrans<-rep(1:ntrans,each=length(tt))
      time<-rep(tt,ntrans)
      
      nHaz<-vector("list", nruns)
      for (k in 1:nruns){
        nHaz[[k]]<-cbind(time=as.vector(time),Haz[,k]<-as.vector(Haz[,k])
                         ,trans=as.vector(newtrans))
        nHaz[[k]]<-as.data.frame(nHaz[[k]])
        names(nHaz[[k]])[2]<-"Haz"
      }
      ### function to count any differences of > 1 in consecutive hazards
      funtemp1<-function(x){
        length(which(diff(x$Haz)>1 & diff(x$trans)==0))  
      }
      ### function to identify position of any differences of > 1 in consecutive hazards
      funtemp2<-function(x){
        which(diff(x$Haz)>1 & diff(x$trans)==0)  
      }
      temp1<-sapply(X= nHaz, FUN=funtemp1)
      temp2<-sapply(X= nHaz, FUN=funtemp2)
      if(length(which(temp1!=0))!=0){
      Hazerror<-temp2[which(temp1==1)[1]]
      Hazerror<-as.numeric(Hazerror)
     
      ### adjust hazards that had differences of > 1 in consecutive hazards
      for (k in 1:nruns){
        if (temp1[k]==1 & temp2[[k]][1]==Hazerror & (nHaz[[k]]$Haz[(Hazerror+1)]-nHaz[[k]]$Haz[(Hazerror-1)])<2){
          nHaz[[k]]$Haz[Hazerror]<-(nHaz[[k]]$Haz[(Hazerror-1)]+ nHaz[[k]]$Haz[(Hazerror+1)])/2
        }
        else if (temp1[k]==1 & temp2[[k]][1]==Hazerror & (nHaz[[k]]$Haz[(Hazerror+1)]-nHaz[[k]]$Haz[(Hazerror-2)])<3 & (nHaz[[k]]$Haz[(Hazerror+1)]-nHaz[[k]]$Haz[(Hazerror-2)])>2){
          nHaz[[k]]$Haz[Hazerror]<- nHaz[[k]]$Haz[(Hazerror+1)]-(nHaz[[k]]$Haz[(Hazerror+1)]- nHaz[[k]]$Haz[(Hazerror-2)])/3
          nHaz[[k]]$Haz[(Hazerror-1)]<- nHaz[[k]]$Haz[(Hazerror+1)]-2*(nHaz[[k]]$Haz[(Hazerror+1)]- nHaz[[k]]$Haz[(Hazerror-2)])/3
        }
        
        else if (temp1[k]==1 & temp2[[k]][1]==Hazerror & (nHaz[[k]]$Haz[(Hazerror+1)]-nHaz[[k]]$Haz[(Hazerror-3)])<4 & (nHaz[[k]]$Haz[(Hazerror+1)]-nHaz[[k]]$Haz[(Hazerror-3)])>3){
          nHaz[[k]]$Haz[Hazerror]<- nHaz[[k]]$Haz[(Hazerror+1)]-(nHaz[[k]]$Haz[(Hazerror+1)]- nHaz[[k]]$Haz[(Hazerror-3)])/4
          nHaz[[k]]$Haz[(Hazerror-1)]<- nHaz[[k]]$Haz[(Hazerror+1)]-2*(nHaz[[k]]$Haz[(Hazerror+1)]- nHaz[[k]]$Haz[(Hazerror-3)])/4
          nHaz[[k]]$Haz[(Hazerror-2)]<- nHaz[[k]]$Haz[(Hazerror+1)]-3*(nHaz[[k]]$Haz[(Hazerror+1)]- nHaz[[k]]$Haz[(Hazerror-3)])/4
        }
        else if (temp1[k]==1 & temp2[[k]][1]==Hazerror & (nHaz[[k]]$Haz[(Hazerror+1)]-nHaz[[k]]$Haz[(Hazerror-4)])<5 & (nHaz[[k]]$Haz[(Hazerror+1)]-nHaz[[k]]$Haz[(Hazerror-4)])>4){
          nHaz[[k]]$Haz[Hazerror]<- nHaz[[k]]$Haz[(Hazerror+1)]-(nHaz[[k]]$Haz[(Hazerror+1)]- nHaz[[k]]$Haz[(Hazerror-4)])/5
          nHaz[[k]]$Haz[(Hazerror-1)]<- nHaz[[k]]$Haz[(Hazerror+1)]-2*(nHaz[[k]]$Haz[(Hazerror+1)]- nHaz[[k]]$Haz[(Hazerror-4)])/5
          nHaz[[k]]$Haz[(Hazerror-2)]<- nHaz[[k]]$Haz[(Hazerror+1)]-3*(nHaz[[k]]$Haz[(Hazerror+1)]- nHaz[[k]]$Haz[(Hazerror-4)])/5
          nHaz[[k]]$Haz[(Hazerror-3)]<- nHaz[[k]]$Haz[(Hazerror+1)]-4*(nHaz[[k]]$Haz[(Hazerror+1)]- nHaz[[k]]$Haz[(Hazerror-4)])/5
        }
        
      }
      }
      rm(Haz)
      if (Markov==TRUE){
      #### calculate state occupancy probabilities (Markov) 
        msf<-vector("list", nruns)
        msf<-vector("list", nruns)
        for (k in 1:nruns) {
          msf[[k]]<-list(nHaz[[k]],trans)
          msf[[k]]$trans<-trans
          msf[[k]]$Haz<-nHaz[[k]]
          class(msf[[k]])<-"msfit"
        }
        prob<- lapply(msf, function(x)
         try(probtrans(x, predt = 0,variance=FALSE), TRUE))
        err<-lapply(X=prob, FUN=is.error)
        errid<-which(err==1)
        prob<-prob[which(err==0)] 
      }
      if (Markov==FALSE){
        #### calculate state occupancy probabilities (semi-Markov)  
        set.seed(seedno)
          prob<- lapply(nHaz, function(x) 
          try(mssample(Haz=x, trans=trans,tvec=tt,M=M,clock="reset"), TRUE))
        err<-lapply(X=prob, FUN=is.error)
        errid<-which(err==1)
        prob<-prob[which(err==0)]
     }
     rm(nHaz)
      is.wholenumber <-
        function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
      
      for (k in 1:(length(which(err==0)))) {
        prob[[k]][,nrow(tmat)+2]<-is.wholenumber(prob[[k]][,1]*12)
        prob[[k]]<-prob[[k]][prob[[k]][,nrow(tmat)+2]==1, -(nrow(tmat)+2)]
      }
      if (length(errid)>0) return(list(errid,prob))
      if (length(errid)==0) return(prob)
}

##########################################################################################
##========================================================================================
## END OF PSAPROB FUNCTION
##========================================================================================
##########################################################################################
















