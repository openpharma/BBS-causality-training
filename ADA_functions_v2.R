#######################################################################
######### Supportive functions for R markdown of ADA example based
######### on training given at BBS training series 2021 
######### "Gentle introduction into causal thinking" Feb 2, 2021
#######################################################################


##################################################################
## A) Help functions for model diagnostics
##################################################################

#function to calculate standardized mean difference
# for balance checking
std.mean=function(data_var,weight,treatment){
  data0_var=subset(data_var,treatment==0)
  weight0=weight[treatment==0]
  data1_var=subset(data_var,treatment==1)
  weight1=weight[treatment==1]
  
  mean_std=NULL
  mean_std$var=colnames(data_var)
  for (i in 1:dim(data_var)[2]){
    Z.wmean0=sum(data0_var[,i]*weight0)/sum(weight0)
    Z.wmean1=sum(data1_var[,i]*weight1)/sum(weight1)
    if (Z.wmean0 ==Z.wmean1){
      Z.wmd = 0
    } 
    else {
      Z.ws0=sum(weight0)/(sum(weight0)^2 -sum(weight0^2))*sum(weight0*(data0_var[,i]-Z.wmean0)^2)
      Z.ws1=sum(weight1)/(sum(weight1)^2 -sum(weight1^2))*sum(weight1*(data1_var[,i]-Z.wmean1)^2)
      Z.wmd=(Z.wmean1-Z.wmean0)/sqrt((Z.ws0+Z.ws1)/2)
    }
    mean_std$wmd[i]=Z.wmd
  }
  return(mean_std)
}

#function to count the number of different critical variables in vector
#argument varlist contains the names of the variables
#variables with a "_" in the name and otherwise same should only be counted once
get_no_of_critvars=function(varlist){
  res=c()
  i=1
  ## iterate over variable names in vector
  while(i<=length(varlist)){
    if (is.na(varlist[i])) {i=i+1}
    else{
      #index of "_" in variable name, remove later characters
      ind = regexpr("_", varlist[i])[1]
      if (ind>0){
        res=cbind(res,substr(varlist[i],1,ind-1))
      }
      # no "_" in the name, so just keep variable
      if (ind<=0) {
        res=cbind(res,varlist[i])
      }
      i=i+1
    }
  }
  len=length(unique(as.vector(res)))
  return (len)
}

# function returning number of expected features with absolute standardized difference above threshold
# based on total number of features, sample size in active and control (as expected per randomization ratio)
# and threshold value (usually 0.1 and 0.25)
expected_fdiff <- function(n.feature, n.arm1, n.arm2, threshold){
  threshold.t <- threshold / sqrt( (1/n.arm1) + 1/(n.arm2))
  df <- n.arm1 + n.arm2 - 2
  stopifnot(threshold.t>=0)
  p <- (1-pt(threshold.t, df))*2
  res <- n.feature*p
  names(res) <- paste("Expected # features >", threshold)
  return (res)
}



##################################################################
## B) Help functions for control arm weighting approach 
##################################################################

#calculate adjusted HR in ADA+ and ADA- (assuming no missing ADA status at landmark, else see Kong, Heinzmann and Lauer 2021 https://arxiv.org/abs/2101.04263)
#when balance diagnostics are requested (diagnosis=TRUE), varlist and ratio both must be set

wpa=function(data_wpa,formula,plot=TRUE,filename="plot",diagnosis=FALSE,varlist=NULL, randratio=NULL){
  
  #set shorter label for KM plots
  data_wpa$ADALMC=ifelse(data_wpa$ADALMnum==1,"ADA-",
                      ifelse(data_wpa$ADALMnum==2, "ADA+",
                             ifelse(data_wpa$ADALMnum==0, "Control","")))
                      
  # separate treatment and ADA groups
  q1=subset(data_wpa, ADALMnum ==1)  # ADA-
  q2=subset(data_wpa, ADALMnum ==2)  # ADA+
  q0=subset(data_wpa, ADALMnum ==0)  # control
  
  # weight.pos = probability to be ADA-positive (i.e. 0 or 1 in treatment arm)
  q1$weight.pos=0
  q2$weight.pos=1 
  
  #model probability to be ADA+ based on baseline covariates in the treatment arm
  fitdata=rbind(q1,q2)
  fitdata$indicator=fitdata$ADALMnum-1  # Indicator for ADA neg is one, ADA pos is 2, change to Ada neg=0,ADA pos=1
  fit.glm=glm(formula,family=binomial,data=fitdata)
  
  #predict probability to be ADA+ in control based on above model
  q0$weight.pos=predict(fit.glm,newdata=q0,type="response")
  
  summary(q0$weight.pos)
  
  #combine control and ADA neg patients in one dataset
  data1=rbind(q0,q1)
  data1$weight.neg=1-data1$weight.pos  #probability to be ADA-negative is 1-probability to be ADA-positive
  data1$weight = data1$weight.neg
  #combine control and ADA pos. patients in one dataset
  data2=rbind(q0,q2)
  data2$weight.neg=1-data2$weight.pos
  data2$weight=data2$weight.pos
 
  fit1=coxph(Surv(AVAL,CNSR==0)~trt,weights=weight,ties="breslow",data=data1)
  loghr.neg=as.numeric(summary(fit1)$coefficients[1])
  
  fit2=coxph(Surv(AVAL,CNSR==0)~trt,weights=weight,ties="breslow",data=data2)
  loghr.pos=as.numeric(summary(fit2)$coefficients[1])
  
  
  if (plot){
    dataplot=list(data1,data2)
    # separate plots for ADA- and ADA+
    for(j in 1:2){
      grp_labelc=dataplot[[j]][dataplot[[j]]["trt"]==1,]$ADALMC[1]
      grp_label=dataplot[[j]][dataplot[[j]]["trt"]==1,]$ADALMnum[1]
      dataplot[[j]]$status=ifelse(dataplot[[j]]["trt"]==1, grp_labelc, paste(grp_labelc,'control'))
      
      fit_km=survfit(Surv(AVAL, CNSR==0)~status, weights=weight, data=dataplot[[j]])
      
      nam.plot<- ifelse(grp_label==1,"ADAneg","ADApos")
      png(paste(filename, '_', nam.plot,'.png',sep=''),width=720,height=480)
      print(ggsurvplot(fit_km,data = dataplot[[j]],  risk.table = F,   pval = F,          
                       conf.int = F,  axes.offset=F, xlim=c(0,max(dataplot[[j]]$AVAL)),break.time.by=3, xlab=paste('Time (', tolower(dataplot[[j]]$AVALU[1]), ')',sep=''), palette=c('red','blue')))
      dev.off()
    }
    
    #plot also only ADA groups 
    data3=rbind(q1,q2)
    data3$status=data3$ADALMC
    
    fit_km=survfit(Surv(AVAL, CNSR==0)~status, data=data3)
    
    png(paste(filename, '_ADAgrps.png',sep=''),width=720,height=480)
    print(ggsurvplot(fit_km,data = data3,  risk.table = F,   pval = F,          
                     conf.int = F,  axes.offset=F, xlim=c(0,max(data3$AVAL)),break.time.by=3, xlab=paste('Time (', tolower(data3$AVALU[1]), ')',sep=''), palette=c('red','dark red')))
    dev.off()
    
    #plot control, overall and both adjusted
    #3 times control, to see the split of the overall control arm in the 2 adjusted controls
    q0$weight.neg=1-q0$weight.pos
    q0_rep1=q0
    q0_rep2=q0
    q0_rep3=q0
    q0_rep1$status = "Overall control"
    q0_rep2$status = "Adj. control for ADA-"
    q0_rep3$status = "Adj. control for ADA+"
    q0_rep1$weight.all=1
    q0_rep2$weight.all=q0_rep2$weight.neg
    q0_rep3$weight.all=q0_rep3$weight.pos
    
    data4=rbind(q0_rep1,q0_rep2,q0_rep3)
    
    fit_km=survfit(Surv(AVAL, CNSR==0)~status, weights=weight.all, data=data4)
    
    png(paste(filename, '_Controlgrps.png',sep=''),width=720,height=480)
    print(ggsurvplot(fit_km,data = data4,  risk.table = F,   pval = F,          
                     conf.int = F,  axes.offset=F, xlim=c(0,max(data4$AVAL)),break.time.by=3, xlab=paste('Time (', tolower(data4$AVALU[1]), ')',sep=''), palette=c('blue','light blue', 'dark blue')))
    dev.off()
    
  }
  
  #all outputs for balance diagnostics i.e.
  #count of expected vs observed features
  #balance plots
  #absolute mean differences (against unadjusted and adjusted control)
  if (diagnosis){
    # get number of observed/reversed features with absolute adjusted mean difference > 0.1 resp. > 0.25
    wmd.crit=matrix(NA,nrow=2,ncol=9)
    row.names(wmd.crit)=c("ADA-","ADA+")
    colnames(wmd.crit)=c("Rand. ratio", "No. of covariates", "N",
                        "#Expected features>0.1", "#Obs. features>0.1","#Rev. features>0.1",
                         "#Expected features>0.25", "#Obs. features>0.25", "#Rev. features>0.25")
    
   
    #initialize (standard) threshold variables
    thresh1=0.1
    thresh2=0.25
    #number of different baseline covariates
    nvars=get_no_of_critvars(varlist)
    #number of ADA- and ADA+ patients in treatment arm
    n1=nrow(q1)
    n2=nrow(q2)
    
    #write out randomization ratio and number of covariates (same for all groups)
    wmd.crit[,1]=paste0("Ratio ", randratio, ":1")
    wmd.crit[,2]=nvars
    
    ######### analyze balance for ADA-
    data1_cov=data1[,(names(data1) %in% varlist)]
    # calculated unadjusted mean differences (weight = 1)
    wmd1_unadjusted=std.mean(data_var=data1_cov,weight=rep(1,dim(data1_cov)[1]),treatment=data1$trt)
    # calculated adjusted mean difference (weight = 1 in treatment but different in control)
    wmd1_adjusted=std.mean(data_var=data1_cov,weight=data1$weight,treatment=data1$trt)
    wmd1=NULL
    wmd1$var=rep(wmd1_adjusted$var,2)
    wmd1$absolute_std_diff=abs(c(wmd1_adjusted$wmd,wmd1_unadjusted$wmd))
    wmd1$type=c(rep("adjusted",length(wmd1_unadjusted$var)),rep("unadjusted",length(wmd1_unadjusted$var)))
    ## get number of adjusted variables exceeding 0.1 or 0.25 absolute difference
    wmd1_adjusted$critvar1=ifelse((abs(wmd1_adjusted$wmd)>thresh1), wmd1_adjusted$var, NA)
    wmd1_adjusted$critvar2=ifelse((abs(wmd1_adjusted$wmd)>thresh2), wmd1_adjusted$var, NA)
    # get number of adjusted variables exceeding 0.1 or 0.25 absolute difference and worse than non-adjusted
    wmd1_adjusted$critrev1=ifelse((abs(wmd1_adjusted$wmd)>thresh1 & abs(wmd1_adjusted$wmd)-abs(wmd1_unadjusted$wmd)>0), wmd1_adjusted$var, NA)
    wmd1_adjusted$critrev2=ifelse((abs(wmd1_adjusted$wmd)>thresh2 & abs(wmd1_adjusted$wmd)-abs(wmd1_unadjusted$wmd)>0), wmd1_adjusted$var, NA)
    # write out results from above, plus number of patients and expected feature differences (rounded to 1 digit)
    wmd1_crit = cbind(n1, round(expected_fdiff(nvars,n1,n1/randratio,thresh1), 1), get_no_of_critvars(wmd1_adjusted$critvar1), get_no_of_critvars(wmd1_adjusted$critrev1),
                      round(expected_fdiff(nvars,n1,n1/randratio,thresh2),1), get_no_of_critvars(wmd1_adjusted$critvar2), get_no_of_critvars(wmd1_adjusted$critrev2))
    wmd.crit[1,3:9] = as.vector(wmd1_crit)
    
    ######## now repeat same as above for ADA+
    data2_cov=data2[,(names(data2) %in% varlist)]
    wmd2_unadjusted=std.mean(data_var=data2_cov,weight=rep(1,dim(data2_cov)[1]),treatment=data2$trt)
    wmd2_adjusted=std.mean(data_var=data2_cov,weight=data2$weight,treatment=data2$trt)
    wmd2=NULL
    wmd2$var=rep(wmd2_adjusted$var,2)
    wmd2$absolute_std_diff=abs(c(wmd2_adjusted$wmd,wmd2_unadjusted$wmd))
    wmd2$type=c(rep("adjusted",length(wmd2_unadjusted$var)),rep("unadjusted",length(wmd2_unadjusted$var)))
    wmd2_adjusted$critvar1=ifelse((abs(wmd2_adjusted$wmd)>thresh1), wmd2_adjusted$var, NA)
    wmd2_adjusted$critvar2=ifelse((abs(wmd2_adjusted$wmd)>thresh2), wmd2_adjusted$var, NA)
    wmd2_adjusted$critrev1=ifelse((abs(wmd2_adjusted$wmd)>thresh1 & abs(wmd2_adjusted$wmd)-abs(wmd2_unadjusted$wmd)>0), wmd2_adjusted$var, NA)
    wmd2_adjusted$critrev2=ifelse((abs(wmd2_adjusted$wmd)>thresh2 & abs(wmd2_adjusted$wmd)-abs(wmd2_unadjusted$wmd)>0), wmd2_adjusted$var, NA)
    wmd2_crit = cbind(n2, round(expected_fdiff(nvars,n2,n2/randratio,thresh1),1), get_no_of_critvars(wmd2_adjusted$critvar1), get_no_of_critvars(wmd2_adjusted$critrev1),
                      round(expected_fdiff(nvars,n2,n2/randratio,thresh2),1), get_no_of_critvars(wmd2_adjusted$critvar2), get_no_of_critvars(wmd2_adjusted$critrev2))
    wmd.crit[2,3:9] = as.vector(wmd2_crit)
    
    #write out results for number of features exceeding 0.1 resp. 0.25 for adjusted mean difference
    write.csv(wmd.crit,"Feature_count.csv")
    
    ###################################################################################################
    ## save absolute differences##
    ## Remark: If numerical results are needed in addition to grphical display of ASMD, uncomment this part
      # 
      #wmd1.save=cbind(wmd1_unadjusted$var,wmd1_unadjusted$wmd,wmd1_adjusted$wmd)
      #wmd2.save=cbind(wmd2_unadjusted$var,wmd2_unadjusted$wmd,wmd2_adjusted$wmd)
      #  
      #colnames(wmd1.save)=c("var","unadjusted","adjusted")
      #colnames(wmd2.save)=c("var","unadjusted","adjusted")
      #
      #write.csv(wmd1.save,"ASMD_negative.csv")
      #write.csv(wmd2.save,"ASMD_positive.csv")
    ###################################################################################################

    
    
    ##plot for absolute difference##
    png('Balance_negative.png')
    plot1<-ggplot(as.data.frame(wmd1), aes(absolute_std_diff, var)) + geom_point(aes(color = type))+ scale_color_manual(values=c("blue","red"))+scale_x_continuous(breaks=seq(0,0.3,by=0.1))+expand_limits(x=0.3)+xlab("ASMD") + ggtitle("ADA-")+ylab("Covariate") +theme(
      axis.title.x = element_text(color="black", size=24, face="bold"),
      axis.title.y = element_text(color="black", size=24, face="bold"),
      text = element_text(size=20),
      plot.title = element_text(hjust=1),
    )
    print(plot1)
    dev.off()
    
    
    png('Balance_positive.png')
    plot2=ggplot(as.data.frame(wmd2), aes(absolute_std_diff, var)) +geom_point(aes(color = type))+scale_color_manual(values=c("blue","red"))+scale_x_continuous(breaks=seq(0,0.3,by=0.1))+expand_limits(x=0.3)+xlab("ASMD") + ggtitle("ADA+")+ylab("Covariate") +theme(
      axis.title.x = element_text(color="black", size=24, face="bold"),
      axis.title.y = element_text(color="black", size=24, face="bold"),
      text = element_text(size=20),
      plot.title = element_text(hjust=1),
    )
    print(plot2)
    dev.off()
    
    
  }
  
  return(list(loghr.pos=loghr.pos,loghr.neg=loghr.neg))
}

