library(rpart)
library(statmod)
library(survival)
library(randomForest)
library(party)
library(randomForestSRC)
library(MASS)
library(pec)
library(timereg)
library(prodlim)
library(rms)

library(survPresmooth)
library(MASS)

#Code written by M. Matabuena

# Compute weights (by default: Kaplan-Meier estimator)
# Choose between: random forest, Beran estimator, Cox model
pesos= function(X,time,indicadorcensura, metodo="kaplan"){
  time= as.vector(time)
  
  status2= as.numeric(indicadorcensura)
  status= as.factor(indicadorcensura)
  X= as.data.frame(X)
  n= dim(X)[1]
  p= dim(X)[2]
  nombres= paste("V",1:p,sep="")
  colnames(X)= nombres
  formula3= "Surv(time,status2)~V1"
  formula= "Hist(time,status)~V1"
  formula2= c("")
  if(p>1){
    for(i in 2:p){
      formula2= paste(formula2,paste("+",nombres[i],sep=""),sep="")
    }
  }  
  
  
  formula= as.formula(paste(formula,formula2,sep=""))
  
  formula3= as.formula(paste(formula3,formula2,sep=""))
  
  
  dat= data.frame(X,time,status,status2)
  
  ordentiempos= order(dat$time)
  dat <- dat[ordentiempos,]
  
  if(metodo=="kaplan"){
    W_kaplan=ipcw(formula,data=dat,method="marginal",times=sort(unique(dat$time)),subjectTimes=dat$time)
    W_kaplan=W_kaplan$IPCW.subjectTimes
    
    W_kaplan= W_kaplan
    return(W_kaplan)
  }
  
  
  
  if(metodo=="cox"){
    W_cox=ipcw(formula,data=dat,method="cox",times=sort(unique(dat$time)),subjectTimes=dat$time)
    W_cox= W_cox$IPCW.subjectTimes
    W_cox= W_cox
    return(W_cox)
    
  }
  
  if(metodo=="random_forest"){
    
    W_rfsrc=ipcw(formula, data=dat,method="rfsrc",times=sort(unique(dat$time)),subjectTimes=dat$time)
    W_rfsrc= W_rfsrc$IPCW.subjectTimes
    return(W_rfsrc)
  }
  
  if(metodo=="beran"){
    W_nonpar=ipcw(formula, data=dat,method="nonpar",times=sort(unique(dat$time)),subjectTimes=dat$time)
    W_nonpar= W_nonpar$IPCW.subjectTimes
    W_nonpar= W_nonpar
    return(W_nonpar)
  }
  
}
