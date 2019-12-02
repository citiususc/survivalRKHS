n=100 #Sample size
p=5 #Number of covariates

Rcpp::sourceCpp("censuracpp.cpp")

#Vector whose entries are the weights for each observation
#In the example below the weights correspond to non-censored data
#With censored data, use function pesos.R to obtain appropriate weights
pesos=rep(1,n) 
pesos= pesos/sum(pesos) #Normalize the weights.
wcens= pesos%*%t(pesos); 
wcens= wcens*n*(n-1)/2; #Matrix constructed with the weights. 

#Random sample of p incorrelated uniformly distributed variables
X<-matrix(0,nrow=n,ncol=p)

for (i in 1:p) {
  X[,i]<-runif(n,0,10)
}


#Create a dependent response
Y<-rep(0,n)
Y<-log(X[,3])
numberclose<-10 #Number of neighbour points 
totaltime<-3 #Number of cross validation iterations 

#dim(wcens) #Check it is nxn

res<-ModelFreeVS(X,Y,numberclose,totaltime,wcens)
print(res)

