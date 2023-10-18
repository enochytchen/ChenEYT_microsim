##########################################################################################
##========================================================================================
## START OF PSAQALY FUNCTION
##========================================================================================
##########################################################################################

##========================================================================================
## THIS FUNCTION RETURNS THE QALYs IN A GIVEN STATE FOR A GIVEN TREATMENT ARM FOR EACH
## DRAW OF THE PROBABILISTIC SENSITIVITY ANALYSIS. 

## MODIFIBLE ARGUMENTS:
## Markov       EITHER TRUE OR FALSE, FOR object BEING BASED ON A MARKOV OR SEMI-MARKOV 
##              MODEL RESPECTIVELY
## object       OBJECT RETURNED FROM THE PSAprob FUNCTION 
## state        STATE OF INTEREST
## discounted   DISCOUNTED? EITHER TRUE OR FALSE.
## dis1yronwards  START DISCOUNTING FROM YEAR 1 INSTEAD OF TIME 0?  EITHER TRUE OR FALSE.
## rate         IF DISCOUNTED, RATE EXPRESSED AS A PROPORTION I.E. BETWEEN 0 AND 1
## utility      FIXED UTILITY WEIGHT FOR THE STATE OF INTEREST

##========================================================================================
PSAQALY<-function(Markov=FALSE, object=SMPSAprobRFC, state=1, 
                           discounted=TRUE, dis1yronwards=TRUE, rate=0.035,
                                utility=utility1){
  if (Markov==FALSE & dis1yronwards==FALSE ) {
  
     if (discounted==TRUE) {
       QALY<-sapply(object, function(x) trapz(x$time,x[,state+1]/(1+rate)^x$time))*utility
     }
     if (discounted==FALSE) {
       QALY<-sapply(object, function(x) trapz(x$time,x[,state+1]))*utility
     }
  } 

  if (Markov==TRUE & dis1yronwards==FALSE) {
    
    if (discounted==TRUE) {
      QALY<-sapply(object, 
                   function(x) trapz(x[[1]]$time,x[[1]][,state+1]/(1+rate)^x[[1]]$time))*utility
    }
    if (discounted==FALSE) {
      QALY<-sapply(object, function(x) trapz(x[[1]]$time,x[[1]][,state+1]))*utility
    }
  }
  
  if (Markov==FALSE & dis1yronwards==TRUE) {
    
    temp1<- function(x, rate=0.035){ 
      x$pstaterev<-x[,state+1]/(1+rate)^x$time
      for (i in 1: length(x$time)){
        if(x$time[i]<1) x$pstaterev[i]<-x[i,state+1] 
        return(x)    
      }
    }
    temp2<-lapply(object,FUN=temp1)
    
      QALY<-sapply(temp2, function(x) trapz(x$time,x$pstaterev))*utility
      
    }

  if (Markov==TRUE & dis1yronwards==TRUE) {
    
    temp1<- function(x, rate=0.035){ 
      x[[1]]$pstaterev<-x[[1]][,state+1]/(1+rate)^x[[1]]$time
      for (i in 1: length(x[[1]]$time)){
        if(x[[1]]$time[i]<1) x[[1]]$pstaterev[i]<-x[[1]][i,state+1] 
        return(x)    
      }
    }
    temp2<-lapply(object,FUN=temp1)
    
    QALY<-sapply(temp2, function(x) trapz(x[[1]]$time,x[[1]]$pstaterev))*utility
    
  }
  
  return(QALY) 
  
  }
  
##########################################################################################
##========================================================================================
## END OF PSAQALY FUNCTION
##========================================================================================
##########################################################################################






