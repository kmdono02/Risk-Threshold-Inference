Sim_data_logit <- function(Risk, beta1=-5, lambda=1, theta=Inf, n, p, p0, p1, seed_val){
  # Risk: disease risk
  # beta1: slope parameter
  # lambda: scale parameter for Dunning's scaled logit model
  # theta: truncation point (v_0), Inf reflects no truncation point
  # n: sample size
  # p0: immune response sampling rate for control group case-control study design
  # p1: immune response sampling rate for case group case-control study design
  # seed_val: seed for simualtion
  library(nleqslv)
  library(dplyr)
  # Setup model
  ## solve for beta0 using grid search (intercept for vaccine disease model) ##
  
  integrand <-function(s,beta0){	
    lin <- beta0 + beta1*s
    expit <- 1/(1+exp(-lin))
    ans <- lambda*expit*dgamma(s,4,1)
    ans
  }
  
  beta0s <- seq(-100,100,by=.1)
  nbeta <- length(beta0s)
  y <- beta0s*0
  
  for (ii in 1:nbeta)
    y[ii] <- integrate(integrand,0,theta,beta0=beta0s[ii])$value
  
  #plot(beta0s,y)
  
  mydiff <- abs(y-Risk)
  mymin <- min(mydiff)
  beta0 <- head(beta0s[mydiff==mymin],1)
  
  #####
  
  ## Find true threshold
  expit <- function(x){1/(1+exp(-x))}
  
  # expit1 is just expit when the distribution is truncated by theta.  If theta=Inf, they are the same.
  
  expit1 <- function(x){
    ifelse(0<=x & x<=theta, expit(beta0+beta1*x), 0)
  }
  
  # Need to find threshold that results in desired risk level, risk level term defined by ratio of 2 integrals.
  
  riskvacnum <- function(x){dgamma(x,4,1)*lambda*expit1(x)}
  riskvacden <- function(x){dgamma(x,4,1)}
  riskvac <- function(s){(integrate(riskvacnum,s,theta)$value)/(integrate(riskvacden,s,Inf)$value)}
  
  # Using grid search to find threshold
  v <- seq(0,8,by=.001)
  nv <- length(v)
  risks <- vector(, length=nv)
  for(i in 1:nv){
    risks[i]<-riskvac(v[i])
  }
  
  # Calculate true thresholds for risk level vector p
  truethres <- rep(NA, length=length(p))
    expit <- function(x){1/(1+exp(-x))}
    expit1 <- function(x){
      ifelse(0<=x & x<=theta, expit(beta0+beta1*x), 0)
    }
    
    v <- seq(0,8,by=.001)
    nv <- length(v)
    risks <- vector(, length=nv)
    for(i in 1:nv){
      risks[i]<-riskvac(v[i])
    }
    
  for(j in 1:length(p)){  
    mydiff2 <- abs(risks-p[j])
    mymin2 <- min(mydiff2)
    index <- which(mydiff2==mymin2)
    truethres[j] <- min(v[index])
  }
  
  ## Simulate data
  set.seed(seed_val)
  S <- rgamma(n,4,1)
  Y <- rbinom(n,1,lambda*expit1(S))
  vaccases <- sum(Y)
  
  ## Creating matricies for observed data under case control, incorporating resulting weights ##
  vacdata <- cbind(S, Y)
  vacy1 <- vacdata[Y==1, ]
  vacy0 <- vacdata[Y==0, ]
  
  v0obs <- sample(dim(vacy0)[1], trunc(p0*(dim(vacy0)[1])), replace=FALSE)
  v1obs <- sample(dim(vacy1)[1], trunc(p1*(dim(vacy1)[1])), replace=FALSE)
  
  # Creating subset for subjects observed and not observed respectively under case-cohort subsampling
  vacy0obs <- data.frame(vacy0[v0obs,]) %>% mutate(S_data=1)
  vacy1obs <- data.frame(vacy1[v1obs,]) %>% mutate(S_data=1)
  
  vacy0not <- data.frame(vacy0[-v0obs,]) %>% mutate(S_data=0)
  vacy1not <- data.frame(vacy1[-v1obs,]) %>% mutate(S_data=0)
  
  vacobs <- rbind(vacy0obs, vacy0not, vacy1obs, vacy1not)
  
  data_obs <- vacobs %>% mutate(S_obs=ifelse(S_data==1, S, NA)) %>% select(S_obs, Y, S)
  
  output <- list(data_obs, truethres, c(p0, p1))
  names(output)[1] <- "simmed_data"
  names(output)[2] <- "true_thresholds"
  names(output[[2]]) <- as.character(p)
  names(output)[3] <- "case_control_sampling_probs"
  names(output[[3]]) <- c("control_prob", "case_prob")
  return(output)
}