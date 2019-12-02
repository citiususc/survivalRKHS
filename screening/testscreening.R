Rcpp::sourceCpp("screening.cpp")



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

#Alternative screening method 1 (see paper)
screeningcomparativa= function(x,delta,time,gamma){
pairwise.difference <- function(m){
  npairs <- choose( ncol(m), 2 )
  results <- matrix( NA, nc=npairs, nr=nrow(m) )
  cnames <- rep(NA, npairs)
  if(is.null(colnames(m))) colnames(m) <- paste("col", 1:ncol(m), sep="")
  k <- 1
  for(i in 1:ncol(m)){
    for(j in 1:ncol(m)){
      if(j <= i) next;
      results[ ,k] <- m[ ,i] - m[ ,j]
      cnames[k] <- paste(colnames(m)[ c(i, j) ], collapse=".vs.")
      k <- k + 1
    }
  }
  colnames(results) <- cnames
  rownames(results) <- rownames(m)
  return(results)
}


IPOD.cont= function(j,x,delta,time,gamma){
  #tau=quantile(time, prob=.9)
  
  N=nrow(x); Lambda=(3:round(log(N),0))
  out=array(0,dim=c(length(gamma), length(Lambda)))
  
  for( i in 1: length(Lambda)){##i
    R=Lambda[i] #of slicings, R=3,4,5,6,...
    q=c(quantile(x[,j],probs=c((1:(R-1))/R),na.rm = TRUE))
    
    ##create index for subgroup
    index=rep(1,nrow(x))
    for (r in 1:length(q)){
      ind=which(x[,j]>=q[r])
      index[ind]=r+1
    }
    
    time.pool=numeric(R)
    for(r in 1:R){##r
      time.pool[r]=max(time[index==r])
    }
    t=seq(min(time),min(time.pool),.1)
    h=c(0,diff(range(t))/10)
    
    temp.a=array(0,dim=c(length(t)-1,length(gamma),R))
    R.result=numeric(length(gamma))
    for(r in 1:R){##r
      sub.time=time[index==r]; sub.delta=delta[index==r]
      f=presmooth(sub.time, sub.delta, x.est=t, estimand = "f", bw.selec="fixed",fixed.bw =h)
      
      for(a in 1:length(gamma)){ ##a
        g=f$estimate^gamma[a]
        rec=(t[-1]-t[-length(t)])*g[-length(t)]
        temp.a[,a,r]<- cumsum(rec[1:length(rec)])
      } ##a
    }##r
    
    for( a in 1:length(gamma)){#a
      R.result[a]=max(apply(abs(pairwise.difference (temp.a[,a,])),1,max))
    }##a
    out[,i]=R.result
  }##i
  out=apply(out,1,sum)
  return(out)
}#function

### Discrete covariate
IPOD.disc= function(j,x,delta,time,gamma){
  R=sort(unique(x[,j]))
  time.pool=length(R)
  for(r in 1:length(R)){##r
    time.pool[r]=max(time[x[,j]==R[r]])
  }
  t=seq(min(time),min(time.pool),.1)
  h=c(0,diff(range(t))/10)
  out=numeric(length(gamma))
  temp.a=array(0,dim=c(length(t)-1,length(gamma),length(R)))
  for(r in 1:length(R)){##r
    sub.time=time[x[,j]==R[r]]; sub.delta=delta[x[,j]==R[r]]
    f=presmooth(sub.time, sub.delta, x.est=t, estimand = "f", bw.selec="fixed",fixed.bw =h)
    
    for(a in 1:length(gamma)){ ##a
      g=f$estimate^gamma[a]
      rec=(t[-1]-t[-length(t)])*g[-length(t)]
      temp.a[,a,r]<- cumsum(rec[1:length(rec)])
    }##a
  }##r
  for( a in 1:length(gamma)){#a
    out[a]=max(apply(abs(pairwise.difference (temp.a[,a,])),1,max))
  }##a
  return(out)
} 

## Wrapper function for continuous and discrete covariates
IPOD=function(x,delta,time,gamma){
  p = ncol(x)
  cep=numeric(p)
  one_model=function(j){
    if( length(unique(x[,j]))<5) {res=IPOD.disc(j,x,delta,time,gamma)
    }else  {res=IPOD.cont(j,x,delta,time,gamma)
    }
  }
  cep = sapply(1:p,one_model)
  return(cep)
}
return(IPOD(x,delta,time,gamma=gamma))
}

#Run our screening method
screeningRKHSdependencia= function(x,time,delta){
   p= dim(x)[2]
   n= dim(x)[1]
   x1= data.frame(x)
   
   pesoscalculados= pesos(x1,time,1-delta) 
   pesoscalculados= delta/pesoscalculados
   
   mediana= function(x){
     n= length(x)
     distancias= matrix(0,ncol=n,nrow=n)
     for(i in 1:n){
       for(j in 1:n){
         distancias[i,j]= abs(x[i]-x[j]) 
       }
     }
     
     distancias= as.vector(distancias)
     distancias= distancias[distancias>0]
     distancias= distancias^2
     return(sqrt(median(distancias)/2))
   }
   
   
   asociacion= numeric(p)
   
   for(i in 1:p){
     x1= x[,i]
     y1= time
     
     pesoscalculados1= pesos(x1,y1,1-delta)
     pesoscalculados1= delta/pesoscalculados1
     pesoscalculados2= pesoscalculados
     
     x1= x1[delta==1]
     y1= time[delta==1]

     pesoscalculados1= pesoscalculados1[delta==1]
     pesoscalculados2= pesoscalculados2[delta==1]
     
     
     n= length(x1)
     
     
     aaux2=HS1C(x1,y1,mediana(x1),mediana(y1),1,pesos1=rep(1,n),pesoscalculados2,pesoscalculados1)
    
     asociacion[i]= aaux2
   }
   return(asociacion)
   
 }
 
 
 #Simulated data example in order to testing 
 simul_dat_example2= function(N, p,rho, seed=100){
   if(!is.null(seed))
     set.seed(seed)
   mu0 = rep(2,p);
   
   sigMa0 = matrix(rho,p,p)
   
   diag(sigMa0)=1
   
   X = mvrnorm(N,mu0, sigMa0)
   eps = rnorm(N)
   # Calculate the response Y
   # Define the calculation funtion of Y from x for Case 1
   g1 = function(x)
   { return(0.5*x) }
   
   g2 = function(x,y)
   { return(0.3*x*y) }
   
   g3 = function(x,y)
   {return(0.2*(x+y)^2)}
   
   gY = function(x1,x2,x3,eps)
   {return(g1(x1)*g2(x2,x3)+g3(x3,x2)+eps)}
   
   Y = gY(X[,1],X[,2],X[,3],eps)
   
  CeY = abs(1.2+X[,1]+X[,2])+rnorm(N,0,0.5)
   
   active=1:3

   deltay = as.numeric(Y<=CeY); deltax = rep(1,N)
   obsY = apply(cbind(Y,CeY),1,min);
   
   delta = deltay[order(obsY)]
   print(table(delta))
   x = X[order(obsY),]
   time = obsY[order(obsY)]
   return(list(x=x,time=time,delta=delta,active=active))
 }
 
#Main code
p=200; N=300
dat=simul_dat_example2(N, p,rho=0.0,seed= floor(runif(1,0,100000)))
gamma=c(1,1.2) 
X=dat$x;time=dat$time;delta=dat$delta;active=dat$active
salida= screeningRKHSdependencia(X,as.vector(time),as.vector(delta))
salida= unlist(salida)
print("------------------------")
print("Ours")
print(rank(salida)[1:3])
print("------------------")
print("Alternative")
salida2=screeningcomparativa(X,delta,time,gamma)
print(apply(salida2,1,rank)[1:3,])





