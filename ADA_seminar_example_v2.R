
#read in ADA functions
#calculating probability weighted estimates, including balance checking diagnostics and KM Plots
# set  path to home directory, below only example
setwd('G:/My Drive/BBScausality_Code')
source('ADA_functions_v2.R')
library(survival)
library("survminer")
library("dplyr")
library("ggplot2")
library("Hmisc")
library("rms")
library("Matching")
library("psych")
library("haven")
library(pROC)

###################################################################################################
#Step 1: Read in your csv dataset
###################################################################################################
# data is restricted to relevant population i.e. patients in the safety population
# who did not die or were censored prior to the chosen ADA landmark (week 4, here)
# example data has non-missing ADA status for all patients in experimental arm at the landmark
# and no missing data for covariates

data.ldmrk = read.csv("example_data.csv")

#create additional covariates (grouped from existing ones) as needed for analyses
#in our case: re-group race
data.ldmrk$RACEGRP= ifelse(data.ldmrk$RACE=="WHITE", "WHITE", "OTHER")

#set study-specific covariates, list of names and corresponding formulas for weight approach
#BNLR is continuous, the rest is categorical (TOBHX has 3 levels, remaining ones only 2 levels)
study.cov = c("AGE65","BECOG","SEX","BNLR", "TOBHX","RACEGRP")

#create formula for modeling ADA-status as a function of baseline covariates, known to influence OS outcome
formula=as.formula(indicator~AGE65+BECOG+SEX+BNLR+TOBHX+RACEGRP)
#similar formula for modeling survival as a function of ADA status and same baseline covariates
#for adjusted model comparing ADA+ vs ADA- on treatment arm only 
sformula=as.formula(Surv(AVAL,CNSR==0)~indicator+AGE65+BECOG+SEX+BNLR+TOBHX+RACEGRP)

#create numerical variable reflecting ADA status at landmark (for easier programming in later functions performing 
# the principal stratum analysis)
## value of 3 would only occur in case of missing ADA status at landmark
data.ldmrk$ADALMnum = ifelse(data.ldmrk$ADALM4=="ADA- 4 weeks", 1,
                     ifelse(data.ldmrk$ADALM4=="ADA+ 4 weeks", 2,
                            ifelse(data.ldmrk$ADALM4=="CONTROL",0,3))) 

# set numerical variable reflecting experimental/control treatment arm (1= experimental, 0=control)
data.ldmrk$trt=ifelse(data.ldmrk$ACTARM=="EXPERIMENTAL", 1, 0)

#record sample sizes at landmark
samplesize=NULL
samplesize=rbind(samplesize,table(data.ldmrk$ADALMnum))

row.names(samplesize)=c("sample size")
colnames(samplesize)=c("control","ADA LM-","ADA LM+")

###################################################################################################
# Step 2: Analyze data using weighted placebo approach
###################################################################################################

##create dummy variables for each categorical variable, needed for model diagnostics###
## for feature count, it is important that you use the same naming convention as below i.e.
## for categorical variables with more than 2 categories, use same name upfront (e.g. TOBHX)
## then a "_" followed by the corresponding value (e.g. white, asian or other)
## for variables with only 2 categories, naming could also be different
data.ldmrk$BECOG = as.numeric(data.ldmrk$BECOG)
data.ldmrk$AGE65_old=ifelse(data.ldmrk$AGE65==">= 65",1,0)
data.ldmrk$SEX_male=ifelse(data.ldmrk$SEX=="M",1,0) 
data.ldmrk$TOBHX_previous=ifelse(data.ldmrk$TOBHX=="PREVIOUS",1,0)
data.ldmrk$TOBHX_current=ifelse(data.ldmrk$TOBHX=="CURRENT",1,0)
data.ldmrk$TOBHX_never=ifelse(data.ldmrk$TOBHX=="NEVER",1,0)
data.ldmrk$RACEGRP_white=ifelse(data.ldmrk$RACEGRP=="WHITE",1,0)


##vars is a vector with the names of the dichotomized variables, including the continuous variables
## (only BNLR in our case)
vars=c("BECOG", "AGE65_old", "SEX_male", "TOBHX_previous", "TOBHX_current", "TOBHX_never",
        "RACEGRP_white","BNLR")

########################################################################################################
## no study/endpoint/variable specifics from here onwards anymore
## using normal bootstrap for confidence interval calculation
########################################################################################################

# warnings converted into errors, there should be no warnings for the main call of the 
# weighted placebo function call below

# this option is also needed for correct count of warnings during bootstrap loops
# in case of not too many warnings, you may think of restricting your bootstrap samples only to those without warnings
# i.e. run bootstrap loop for more than 1000 samples and use the first 1000 samples without warning for the
# calculation of the standard error
# in case of too many warnings, you may need to investigate cause and consider using wild bootstrap instead

options(warn=2)

# set up format of report for write out of results (HR estimates + CIs)
# for ADA negative and ADA positive group vs their corresponding adjusted controls
result.report=matrix(NA,nrow=6,ncol=3)
row.names(result.report)=c("ADA- vs unadj control","ADA+ vs unadj control", "ADA- vs adj control","ADA+ vs adj control",
                           "ADA+ vs ADA- unadj", "ADA+ vs ADA- adj")
colnames(result.report)=c("HR estimate","95% Lower","95% Upper")

#set randomization ratio - to determine expected sample size for
# (unobserved counterfactual) ADA groups in control
#if 2:1 (2 active vs 1 control) rratio = 2;  if 1:1 rratio = 1;  if 1:2 (1 active vs 2 control) rratio=0.5
rratio=1

##unadjusted HRs against full control
d1=data.ldmrk[data.ldmrk$ADALMnum %in% c(0,1),]
fit1=coxph(Surv(AVAL,CNSR==0)~trt,ties="breslow",data=d1)
result.report[1,1]=exp(as.numeric(summary(fit1)$coefficients[1]))
result.report[1,2]=summary(fit1)$conf.int[3]
result.report[1,3]=summary(fit1)$conf.int[4]

d2=data.ldmrk[data.ldmrk$ADALMnum %in% c(0,2),]
fit2=coxph(Surv(AVAL,CNSR==0)~trt,ties="breslow",data=d2)
result.report[2,1]=exp(as.numeric(summary(fit2)$coefficients[1]))
result.report[2,2]=summary(fit2)$conf.int[3]
result.report[2,3]=summary(fit2)$conf.int[4]

##unadjusted and adjusted HRs of ADA+vs ADA- i.e. normal Cox
##with either only ADA as covariate or ADA+baseline
d3=data.ldmrk[data.ldmrk$ADALMnum %in% c(1,2),]
d3$indicator=d3$ADALMnum-1 # ADA+=1, ADA-=0
fit3=coxph(Surv(AVAL,CNSR==0)~indicator,ties="breslow",data=d3)
fit4=coxph(sformula,ties="breslow",data=d3)

result.report[5,1]=exp(as.numeric(summary(fit3)$coefficients[1]))
result.report[5,2]=summary(fit3)$conf.int[3]
result.report[5,3]=summary(fit3)$conf.int[4]

result.report[6,1]=exp(as.numeric(summary(fit4)$coefficients[1]))
result.report[6,2]=summary(fit4)$conf.int[1,3]
result.report[6,3]=summary(fit4)$conf.int[1,4]

# call weighted placebo approach function, for generation of HRs, KM-Plots and model diagnostics
loghrs=wpa(data_wpa=data.ldmrk,formula=formula,plot=TRUE,filename="KM",diagnosis = TRUE,
                    varlist=vars,randratio=rratio)

#write HR estimates into report, CIs are later calculated via bootstrapping
result.report[3:4,1]=exp(c(loghrs$loghr.neg,loghrs$loghr.pos))

######################################################################################
##### Step 3: Bootstrap for calculation of CIs for weighted approach  
######################################################################################

# split control arm, ADA-, and ADA+ into separate datasets (data0-data2)
# and get corresponding sample size (n0-n2)
n0=sum(as.numeric(data.ldmrk[,"ADALMnum"])==0)
n1=sum(as.numeric(data.ldmrk[,"ADALMnum"])==1)
n2=sum(as.numeric(data.ldmrk[,"ADALMnum"])==2)
data0=data.ldmrk[data.ldmrk$ADALMnum==0,]
data1=data.ldmrk[data.ldmrk$ADALMnum==1,]
data2=data.ldmrk[data.ldmrk$ADALMnum==2,]

# finally use Bootstrap=1000, you may set this first to only 10 bootstrap loops for quick run
# in case that you get a lot of warnings (mainly due to small sample size)
# check first if the problem is due to covariates with too few patients for a given level
# in small groups like ADA+ and whether you cannot meaningfully re-group
# another problem could also be extreme outlier values for continuous variables (data issues?)

Bootstrap=1000
#store log hazard ratios of all bootstrap samples here
loghrs.boot.all=matrix(NA,nrow=Bootstrap,ncol=2)

# set seed for bootstrap and initialize bootstrap count
set.seed(12345)
boot=1

# variables for tracking of warnings, resp. whether procedure was passed successfully
# during below bootstrap loop
warningnum.wpa=0

while (boot <=Bootstrap){
  data0.boot=data0[sample(x=1:n0,replace=TRUE,size=n0),]
  data1.boot=data1[sample(x=1:n1,replace=TRUE,size=n1),]
  data2.boot=data2[sample(x=1:n2,replace=TRUE,size=n2),]
  data.boot=rbind(data0.boot,data1.boot,data2.boot)
  #wpa bootstrap
  #check the fitting model
  fitdata=subset(data.boot, ADALMnum %in% c(1,2))
  fitdata$indicator=fitdata$ADALMnum-1
  fit.glm=try(glm(formula,family=binomial,data=fitdata),silent=TRUE)
  if (is.list(fit.glm)){
    q0=subset(data.boot, ADALMnum ==0)
    check1=try(predict(fit.glm,newdata=q0,type="response"),silent=TRUE)
    if (is.vector(check1)==TRUE){
      loghrs.boot=wpa(data_wpa=data.boot,formula=formula, plot=FALSE)
      loghrs.boot.all[boot,]=c(loghrs.boot$loghr.neg,loghrs.boot$loghr.pos)
      #print(boot)
      boot=boot+1
    }
    else (warningnum.wpa=warningnum.wpa+1)
  }
  else {warningnum.wpa=warningnum.wpa+1}
  
}

####################################################
##### Step 4: Write out results         ############
####################################################

# calculate standard error from bootstrap samples
# and write CIs into report (assuming normal distribution)
stderr.loghr.neg=sqrt(var(loghrs.boot.all[,1]))
stderr.loghr.pos=sqrt(var(loghrs.boot.all[,2]))

result.report[3,2]=exp(loghrs$loghr.neg-qnorm(0.975)*stderr.loghr.neg)
result.report[3,3]=exp(loghrs$loghr.neg+qnorm(0.975)*stderr.loghr.neg)
result.report[4,2]=exp(loghrs$loghr.pos-qnorm(0.975)*stderr.loghr.pos)
result.report[4,3]=exp(loghrs$loghr.pos+qnorm(0.975)*stderr.loghr.pos)  


#round numbers in report to 3 digits
result.report<-round(result.report,3)

## warnings from normal bootstrap
warningnum=c(Bootstrap,warningnum.wpa)
names(warningnum)=c("Number of Bootstrap samples", "Number of warnings")

#write out results
study = "Example study"
endpoint = "OS"
write.csv(samplesize,file=paste0("samplesize_", study, "_", endpoint, ".csv"))
write.csv(result.report,file=paste0("report_", study, "_", endpoint, ".csv"))
write.csv(warningnum,file=paste0("warningnum_", study, "_", endpoint, ".csv"))



