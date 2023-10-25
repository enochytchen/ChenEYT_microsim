## Filename: 02_msmcancer
## Purpose: Make the data provided in Williams et al. 2017 into multistate structure
## Notes:     The codes were modified from Williams et al. 2017's supplementary materials
## Reference: Williams C, Lewsey JD, Briggs AH, Mackay DF. 
##            Cost-effectiveness Analysis in R Using a Multi-state Modeling Survival Analysis Framework: A Tutorial. 
##            Med Decis Making. 2017 May;37(4):340–52.
##            Cost-effectiveness Analysis in R Using a Multi-state Modeling Survival Analysis Framework: A Tutorial © 2017 by Williams C. et al is licensed under CC BY 3.0. 

## Set the wd as where this R file is
if (require("rstudioapi") && isAvailable()) {
  original_wd <- getwd()  # Store the original working directory
  wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(wd)
}
##############################################################
##============================================================
## Read packages
##============================================================
##############################################################
library(mstate)
library(tidyverse)

##############################################################
##============================================================
## Read the digitised data in Williams et al.
##============================================================
##############################################################
data<-read.table("../Data/data.txt",sep="\t",skip=1) 
names(data)=c("patid","treat","prog","death", "prog_t","death_t",
              "progdeath", "progdeath_t")

## Data cleaning procedures in Williams et al. 2017
### variables in the dataset are:
### patid          patient identification number
### treat          treatment    0=FC, 1=RFC
### prog           progression  0=no, 1=yes
### death          death        0=no, 1=yes
### prog_t         time (in months) to progression or last known not to have a progression
### death_t        time (in months) to death or last known to be alive
### progdeath      progression or death 0=no, 1=yes
### progdeath_t    time (in months) to progression or death or last known to be alive  
###                without progression

### Convert times from months to years
data$prog_ty<-data$prog_t/12
data$death_ty<-data$death_t/12
data$progdeath_ty<-data$progdeath_t/12

### Plot K-M survival curves
mfit1 <- survfit(Surv(death_ty, death==1) ~ treat, data = data)
plot(mfit1, conf.int = FALSE, mark.time = FALSE, lty = 1, col = c("red","blue"),
     xlab = "Time", ylab = "Survival Probability",
     main = "Kaplan-Meier Survival Curves", xlim = c(0, 8))

### Make the follow-up only until 47 months
### There is a drastic drop for treat = 0 at the end of follow-up
data <- data %>%
        mutate(death = ifelse(death_t>=47, 0 , death))
mfit1 <- survfit(Surv(death_ty, death==1) ~ treat, data = data)
plot(mfit1, conf.int = FALSE, mark.time = FALSE, lty = 1, col = c("red","blue"),
     xlab = "Time", ylab = "Survival Probability",
     main = "Kaplan-Meier Survival Curves", xlim = c(0, 8))

#################################################################
##==============================================================
### Prepare relative survival modelling
### Merge the multistate-structured data with population mortality rates
##==============================================================
#################################################################
### Assign the patient characteristics from Hallek et al. 2010  
### age, year, and sex are required to merge with the appointed population
### mortality files
data <- data %>%
        mutate(age = 61, # Median age 61 years
               year = sample(2003:2006, size = nrow(.), replace = TRUE), # Enrolled during 2003-2006 
               sex  = 3, # To be coherent with the average mortality rate across sexes in the popmort file
               att_age = floor(age + death_ty),
               att_year = floor(year + death_ty))

### Merge with popmort file to obtain mortality rate
popmort <- readRDS("../Data/01b_popmort.rds")
names(popmort)
names(popmort) <- c("att_age", "att_year", "rate", "sex")
data <- merge(data, popmort, by = c("att_year", "att_age", "sex"), all = FALSE)

#################################################################
##==============================================================
##Create data subsets based on treatment arm
##==============================================================
#################################################################
## Treat 1: RFC; treat 0: FC
RFCdata<-subset(data, treat==1)
FCdata<-subset(data, treat==0)

## Generate competing risk indicator
## 0=did not die or have a progression, 1=had a progression, 2=died without progression
RFCdata$crstatus<-0
RFCdata$crstatus[RFCdata$prog==1]<-1
RFCdata$crstatus[RFCdata$prog==0 & RFCdata$death==1]<-2

FCdata$crstatus<-0
FCdata$crstatus[FCdata$prog==1]<-1
FCdata$crstatus[FCdata$prog==0 & FCdata$death==1]<-2

## Create data subset for those who had progression
RFCdataprog<-subset(RFCdata, prog==1)
RFCdataprog$time<-RFCdataprog$death_ty-RFCdataprog$prog_ty

## Illness-death model
tmat<- transMat(list(c(2,3), 3, c()), 
                names = c("progression-free", "progression", "death"))
tmat

## Create dataset that can be used for multi-state modelling
covs<-c("treat","prog_ty", "death_ty", "age", "sex", "year", "rate", "att_age", "att_year")
msmcancer<-msprep(time=c(NA, "prog_ty", "death_ty"),
                  status=c(NA, "prog", "death"), data=data, trans=tmat,keep=covs)

## Continue splitting the data by transition
## Create dataset specific to each transition
msmcancer1<-subset(msmcancer,trans==1) # Transition 1 PFS to prog
msmcancer2<-subset(msmcancer,trans==2) # Transition 2 PFS to death
msmcancer3<-subset(msmcancer,trans==3) # Transition 3 prog to death

## Create data subsets based on treatment arm
msmcancer1FC<-subset(msmcancer1,treat==0)
msmcancer2FC<-subset(msmcancer2,treat==0)
msmcancer3FC<-subset(msmcancer3,treat==0)
msmcancer1RFC<-subset(msmcancer1,treat==1)
msmcancer2RFC<-subset(msmcancer2,treat==1)
msmcancer3RFC<-subset(msmcancer3,treat==1)

################################################################
setwd(original_wd)  # Reset to the original working directory
################################################################
## A microsimulation model incorporating relative survival extrapolation and multiple timescales for health technology assessment © 2023 by Chen EYT, Dickman PW, Clements MS is licensed under CC BY 4.0
