Sim_data_logistic <- function(Risk, beta1=2, t0=40, n, p, p0, p1, seed_val){
  # Risk: disease risk
  # beta1: slope parameter
  # t0: follow-up time for study
  # theta: truncation point (v_0), Inf reflects no truncation point
  # n: sample size
  # p0: immune response sampling rate for control group case-control study design
  # p1: immune response sampling rate for case group case-control study design
  # seed_val: seed for simualtion
  
  library(dplyr)
  library(nleqslv)
  # Setup model
  condft <- function(s,t,par0,par1,sigma,a){
    lin <- par0+par1*s
    prt <- t*exp(-lin)
    (prt/(1+prt))*dgamma(s,a,1)
  }
  
  objfn1 <- function(par0){
    integrate(condft, 0, Inf, t=t0, par0, par1=beta1, sigma, a=4)$value-Risk
  }
  
  sol <- nleqslv(beta1, objfn1)
  beta0 <- sol$x
  
  # Calculate true thresholds for risk level vector p
  truethres <- rep(NA, length=length(p))
  for(i in 1:length(p)){
    numthres <- function(s){
    integrate(condft, s, Inf, t=t0, par0=beta0, par1=beta1, sigma, a=4)$value
    }
    denthres <- function(s){
      1-pgamma(s, 4, 1)
    }
    
    thresprob <- function(s){
      (numthres(s)/denthres(s))-p[i]
    }
    
    truethres[i] <- uniroot(thresprob, c(0, 40))$root
  }
  
  ## Simulate data
  set.seed(seed_val)
  E <- rlogis(n, 0, 1)
  S0 <- rgamma(n, 2, 1)
  U <- runif(n, 0, 1)
  X <- rgamma(n, 3, 1)
  S <- X+(U*S0)
  T <- exp(beta0+beta1*S+E)
  CTexp <- rexp(n, rate=1/t0)
  CensFn <- function(x){
    min(t0, x)
  }
  CT <- unlist(lapply(CTexp, CensFn))
  
  ## Generating observed outcomes, censoring indicator (1->Not Censored), vector of weights ##
  TObs <- pmin(T, CT)
  Y <- as.numeric(T<CT)
  Cstatus <- ifelse(Y==0, 1, 0)
  vacdata <- data.frame(cbind(S,TObs,Cstatus))
  
  # Apply case-control design
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
  
  data_obs <- vacobs %>% mutate(S_obs=ifelse(S_data==1, S, NA)) %>% select(S_obs, TObs, Cstatus)
  
  output <- list(data_obs, truethres, c(p0, p1))
  names(output)[1] <- "simmed_data"
  names(output)[2] <- "true_thresholds"
  names(output[[2]]) <- as.character(p)
  names(output)[3] <- "case_control_sampling_probs"
  names(output[[3]]) <- c("control_prob", "case_prob")
  return(output)
}