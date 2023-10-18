##########################################################################################
##========================================================================================
## START OF PSAMEANLY FUNCTION
##========================================================================================
##########################################################################################

##========================================================================================
## THIS FUNCTION RETURNS THE MEAN LIFE YEARS IN A GIVEN STATE FOR A GIVEN TREATMENT ARM. 
## IT CALCULATES THE MEAN ACROSS EACH DRAW IN THE PROBABILISTIC SENSITIVITY ANALYSIS

## MODIFIBLE ARGUMENTS:
## Markov         EITHER TRUE OR FALSE, FOR object BEING BASED ON A MARKOV OR SEMI-MARKOV 
##                MODEL RESPECTIVELY
## object         OBJECT RETURNED FROM THE PSAprob FUNCTION 
## state          STATE OF INTEREST
## instate        EITHER TRUE OR FALSE. TRUE GIVES THE MEAN LIFE YEARS SPENT IN THE STATE. 
##                FALSE GIVES THE MEAN LIFE YEARS OUTWITH THE STATE
## discounted     DISCOUNTED? EITHER TRUE OR FALSE.
## dis1yronwards  START DISCOUNTING FROM YEAR 1 INSTEAD OF TIME 0?  EITHER TRUE OR FALSE.
## rate           IF DISCOUNTED, RATE EXPRESSED AS A PROPORTION I.E. BETWEEN 0 AND 1

##========================================================================================

PSAmeanLY<-function(Markov=FALSE,object=SMPSAprobRFC, state=1, instate=TRUE,
                           discounted=TRUE, dis1yronwards=TRUE, rate=0.035){
  
  if (Markov==FALSE & dis1yronwards==FALSE) {
  
    if (instate==TRUE & discounted==TRUE) {
    meanLY<-mean(sapply(object, function(x) trapz(x$time,x[,state+1]/(1+rate)^x$time)))
    }
    if (instate==FALSE & discounted==TRUE) {
    meanLY<-mean(sapply(object, function(x) trapz(x$time,(1-x[,state+1])/(1+rate)^x$time)))
    }
    if (instate==TRUE & discounted==FALSE) {
    meanLY<-mean(sapply(object, function(x) trapz(x$time,x[,state+1])))
    }
    if (instate==FALSE & discounted==FALSE) {
    meanLY<-mean(sapply(object, function(x) trapz(x$time,1-x[,state+1])))
    }
  return(meanLY) 
  }
  
  if (Markov==TRUE  & dis1yronwards==FALSE) {
  
    if (instate==TRUE & discounted==TRUE) {
      temp<-sapply(object, 
        function(x) trapz(x[[1]]$time,x[[1]][,state+1]/((1+rate)^x[[1]]$time)))
    }
    if (instate==FALSE & discounted==TRUE) {
      temp<-sapply(object, 
        function(x) trapz(x[[1]]$time,(1-x[[1]][,state+1])/(1+rate)^x[[1]]$time))
    }
    if (instate==TRUE & discounted==FALSE) {
      temp<-sapply(object, 
        function(x) trapz(x[[1]]$time,x[[1]][,state+1]))
    }
    if (instate==FALSE & discounted==FALSE) {
      temp<-sapply(object, 
        function(x) trapz(x[[1]]$time,1-x[[1]][,state+1]))
    }
    #### remove any extreme outliers
    if (length(which(temp>1000|temp< -1000))!=0){
      temp<-temp[-which(temp>1000|temp< -1000)]
    }
    #### replace any negative values with zero
    temp<-replace(temp, temp<0, 0)
    meanLY<-mean(temp)
  return(meanLY)
  }  
  
  temp1a<- function(x){ 
    x$pstaterev<-x[,state+1]/(1+rate)^x$time
    for (i in 1: length(x$time)){
      if(x$time[i]<1) x$pstaterev[i]<-x[i,state+1]
    }
    trapz(x$time,x$pstaterev)
  }
  
  temp1b<- function(x){ 
    x$pstaterev<-(1-x[,state+1])/(1+rate)^x$time
    for (i in 1: length(x$time)){
      if(x$time[i]<1) x$pstaterev[i]<-1-x[i,state+1]
    }
    trapz(x$time,x$pstaterev)
  }
  
  
  if (Markov==FALSE & dis1yronwards==TRUE) {
    
    if (instate==TRUE & discounted==TRUE) {
      meanLY<-mean(sapply(object, temp1a))
    }
    if (instate==FALSE & discounted==TRUE) {
      meanLY<-mean(sapply(object, temp1b))
    }
    return(meanLY) 
  }
  
  temp2a<- function(x){ 
    x[[1]]$pstaterev<-x[[1]][,state+1]/(1+rate)^x[[1]]$time
    for (i in 1: length(x[[1]]$time)){
      if(x[[1]]$time[i]<1) x[[1]]$pstaterev[i]<-x[[1]][i,state+1]
    }
    trapz(x[[1]]$time,x[[1]]$pstaterev)
  }
  
  temp2b<- function(x){ 
    x[[1]]$pstaterev<-(1-x[[1]][,state+1])/(1+rate)^x[[1]]$time
    for (i in 1: length(x[[1]]$time)){
      if(x[[1]]$time[i]<1) x[[1]]$pstaterev[i]<-1-x[[1]][i,state+1]
    }
    trapz(x[[1]]$time,x[[1]]$pstaterev)
  }
  
  
  if (Markov==TRUE & dis1yronwards==TRUE) {
    
    if (instate==TRUE & discounted==TRUE) {
      meanLY<-mean(sapply(object, temp2a))
    }
    if (instate==FALSE & discounted==TRUE) {
      meanLY<-mean(sapply(object, temp2b))
    }
    return(meanLY) 
  }
}
  
 
##########################################################################################
##========================================================================================
## END OF PSAMEANLY FUNCTION
##========================================================================================
##########################################################################################


