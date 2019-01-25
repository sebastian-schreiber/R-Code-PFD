# This R code is for the paper 
# "When rarity has costs: coexistence under positive frequency-dependence and environmental stochasticity"
# by  Sebastian J. Schreiber, Masato Yamamichi and Sharon Y. Strauss
# that has been accepted for publication in Ecology. 

# This file provides the  base commands 
# for simulating and exploring the stochastic, two species competition model 
# Ni[t+1]=Ni[t]*(Ni[t]/(ai+Ni[t]+bi*Nj[t]))*(lambdai[t]/(1+alphai*Ni[t]+betai*Nj[t]))+si*Ni[t]
# where lambdai is the intrinsic growth rate of species i
# alphai/betai are the strength of intra/interspecific competition on species i
# ai and bi determine the strength of positive frequency dependence
# and si is the survival of species i
# Most of this code (as with the Ecology paper) only examine the model with ai=0


# required packages
require(rootSolve)
require(MASS)

# The following code solves for 
# equilibria given an initial guess
# input: list of parms lambda, a, b, alpha, beta and initial guess
# output: an equilibrium
equi=function(parms,guess){
  with(as.list(parms),{
    f1=function(x,y)lambda[1]*x^2/(x+b[1]*y)/(1+alpha[1]*x+beta[1]*y)+s[1]*x;
    f2=function(x,y)lambda[2]*y^2/(y+b[2]*x)/(1+alpha[2]*y+beta[2]*x)+s[2]*y;
    h=function(x){(f1(x[1],x[2])-x[1])^2+(f2(x[1],x[2])-x[2])^2}
    out=optim(guess,h)
    return(out)
  }
  )
}

# The following numerically approximates stable equilibria when they exist
# input: list of parms lambda, a, b, alpha, beta
# output: an equilibrium when it exists else NA
stable.equi=function(parms){
  with(as.list(parms),{
    f1=function(x,y)lambda[1]*x^2/(x+b[1]*y)/(1+alpha[1]*x+beta[1]*y)+s[1]*x;
    f2=function(x,y)lambda[2]*y^2/(y+b[2]*x)/(1+alpha[2]*y+beta[2]*x)+s[2]*y;
    L=10
    eq<-(lambda-1)/alpha
    c<-seq(0,1,length=L)
    xs<-c*eq[1]
    ys<-(1-c)*eq[2]
    Tf<-500
    for(t in 1:Tf){
      xs.new<-f1(xs,ys)
      ys.new<-f2(xs,ys)
      xs<-xs.new
      ys<-ys.new
    }
    if(length(which(xs*ys>0))>0){
      temp=min(which(xs*ys>0))
      return(cbind(xs[temp],ys[temp]))}
    else{return(NA)}
  }
  )
}


# this function simulates the stochastic model
# inputs: parameters include length Tf of simulation, mean value of lambdas, standard deviation of 
# log lambdas, autocorrelation rho.time, etc. and initial condition x0
# output: simulation as list of species densities, minimum z of species densities, and lambda values (t1 and t2)

um=function(parms,x0){
	reps=parms$reps
	Tf=parms$Tf
	s=parms$s
	b=parms$b
	alpha=parms$alpha
	beta=parms$beta
	mean.lambda=parms$lambda
	sd.log.lambda=parms$sd.log.lambda
	rho=parms$rho
	rho.time=parms$rho.time
	log.mean.lambda=log(mean.lambda)-sd.log.lambda^2/2
	Sigma=diag(sd.log.lambda)
	Sigma[1,2]=rho*sd.log.lambda[1]
	Sigma[2,1]=rho*sd.log.lambda[1]
	x=matrix(x0[1],Tf,reps)
	y=matrix(x0[2],Tf,reps)
	z=x
	env=mvrnorm(n=reps,mu=c(0,0),Sigma=Sigma)
	for(t in 2:Tf){
	  env=rho.time*env+sqrt(1-rho.time^2)*mvrnorm(n=reps,mu=c(0,0),Sigma=Sigma)
	  temp=exp(cbind(env[,1]+log.mean.lambda[1],env[,2]+log.mean.lambda[2]))
	  
	  x[t,]=x[t-1,]^2*temp[,1]/(x[t-1,]+b[1]*y[t-1,])/(1+alpha[1]*x[t-1,]+beta[1]*y[t-1,])+s[1]*x[t-1,];
	y[t,]=y[t-1,]^2*temp[,2]/(y[t-1,]+b[2]*x[t-1,])/(1+alpha[2]*y[t-1,]+beta[2]*x[t-1,])+s[2]*y[t-1,];
	z[t,]=pmin(x[t,],y[t,])
	}
	return(list(x=x,y=y,z=z))
}

# the following code
# computes the invasion growth rates for species 2 when b=0
# input: parameters
# output: estimate of invasion growth rate
igr=function(parms){
  s=parms$s # survival
  b=c(0,0) # remove PFD
  alpha=parms$alpha #intraspecific competition coeffs. 
  beta=parms$beta #interspecific competition coeffs.
  mean.lambda=parms$lambda # mean lambda values
  sd.log.lambda=parms$sd.log.lambda # sds of log lambda (assumed to be equal)
  rho=parms$rho # cross correlations
  rho.time=parms$rho.time # temporal autocorrelation
  log.mean.lambda=log(mean.lambda)-sd.log.lambda^2/2 # adjust log mean to get desired mean
  Sigma=diag(sd.log.lambda) # initiate covariance matrix with diagonal elements
  Sigma[1,2]=rho*sd.log.lambda[1] # off-diagonal elements of cov. matrix
  Sigma[2,1]=rho*sd.log.lambda[1]
  #
  # Burn in of one species "the resident"
  #
  x=100 # initial condition
  env=mvrnorm(n=1,mu=c(0,0),Sigma=Sigma)
  for(t in 1:1000){ 
    env=rho.time*env+sqrt(1-rho.time^2)*mvrnorm(n=1,mu=c(0,0),Sigma=Sigma)
    temp=exp(env+log.mean.lambda)
    x=x*temp[1]/(1+alpha[1]*x)+s[1]*x
  }
  #
  # calculate fitness of non-resident for a long time series
  #
  Tx=1000000 #1000000
  log.fitness=0
  for(t in 1:Tx){
    env=rho.time*env+sqrt(1-rho.time^2)*mvrnorm(n=1,mu=c(0,0),Sigma=Sigma)
    temp=exp(env+log.mean.lambda)
    x.old=x
    x=x*temp[1]/(1+alpha[1]*x)+s[1]*x;
    fitness=temp[2]/(1+beta[2]*x.old)+s[2];
    log.fitness=log.fitness+log(fitness)
  }
  return(log.fitness/Tx)
}
