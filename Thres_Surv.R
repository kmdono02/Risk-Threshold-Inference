Thres_Surv <- function(dataset, p, bts=2000, alpha=0.025, weighted=FALSE, include_naive=TRUE, 
                       rnd_digits=2){
  # dataset: data
  # p: vector of risk levels of interest
  # bts: number of bootstrap resamples
  # alpha: alpha level for confidence intervals (two sided)
  # weighted: logical value to indicate use of weighted methods, if false then all missing immune response values are removed from data
  # include_naive: logical value to calculate nonparametric estimator which ignors right censoring
  # rnd_digits: number of digits to round results to
  
  ## Check for valid arguments
  if(any(p<0) | any(p>1)){stop("Risk levels p all need to be between 0 and 1")}
  if(alpha>=1){stop("Alpha needs to be < 1; alpha/2 used for coverage probability")}
  if(!(weighted %in% c(TRUE,FALSE))){stop("Please choose a logical value for weighted (default is FALSE)")}
  if(!(include_naive %in% c(TRUE,FALSE))){stop("Please choose a logical value for include_naive (default is TRUE)")}
  if(is.numeric(bts)==FALSE|bts<1){stop("Please choose a numeric value for the number of bootstrap resamples (bts) which is >=1")}
  bts <- floor(bts) # forces boots to be an integer
  ## Adjust dataset based on weighted or non-weighted
  if(weighted==TRUE){
    p0_hat <- mean(!is.na(dataset[dataset[,3]==1,1]))
    p1_hat <- mean(!is.na(dataset[dataset[,3]==0,1]))
    dataset[,4] <- ifelse(is.na(dataset[,1])==0, 1, 0)
  }
  if(weighted==FALSE){
    dataset <- dataset[complete.cases(dataset),]
    p0_hat <- 1
    p1_hat <- 1
    dataset[,4] <- 1
  }
  ## Setting up weights ##
  IPW <- function(x){
    term <- ifelse(x==0, p0_hat, p1_hat)
    1/term
  }
  
  ## Sorting data matrix.  Data needs to be of columns in following order: 
  #  1) Obs. immune response (S) 2) Observed survival time 3) Censor Status (=1 if right censored)

  # Creating Y, Applying weights based on Y, datacomp is set of subjects with complete info
  Y <- 1-dataset[,3]
  Wghts <- unlist(lapply(Y, IPW))
  datawghted <- cbind(dataset[,-c(4)], Y, Wghts*dataset[,4])
  vaccases <- sum(Y)
  
  datacomp <- datawghted[datawghted[,5]!=0,]
  
  # Sorting data based on observed time
  datatsrt <- datacomp[order(datacomp[,2]),]
  StSt <- datatsrt[,1]
  TtSt <- datatsrt[,2]
  CtSt <- datatsrt[,3]
  YtSt <- datatsrt[,4]
  WtSt <- datatsrt[,5]
  
  # Sorting data based on immune response
  datassrt <- datacomp[order(datacomp[,1]),]
  SsSt <- datassrt[,1]
  TsSt <- datassrt[,2]
  CsSt <- datassrt[,3]
  YsSt <- datassrt[,4]
  WsSt <- datassrt[,5]
  
  ## KM Estimator
  
  # Using S indicies
  KMfn <- function(x){
    intx <- floor(x)
    Npt <- (-1*YtSt[StSt>=SsSt[x]]*WtSt[StSt>=SsSt[x]])+(-1*CtSt[StSt>=SsSt[x]]*WtSt[StSt>=SsSt[x]])
    Npt[1] <- sum(WtSt[StSt>=SsSt[x]])+Npt[1]
    Ndiff <- cumsum(Npt)+(CtSt[StSt>=SsSt[x]]*WtSt[StSt>=SsSt[x]])
    Nleft <- (WtSt[StSt>=SsSt[x]]*YtSt[StSt>=SsSt[x]])+Ndiff
    KMQ <- Ndiff/Nleft
    KMEst <- cumprod(KMQ)
    tail(KMEst[is.na(KMEst)==0])[1]
  }
  
  # To limit computation time, find s threshold that is root for KM-p and limit search to
  # highest immune response level in the infection (since rate of infection=0 for larger immune 
  # responses)
  
  maxsinf <- max(SsSt[YsSt==1])
  maxint <- which(SsSt==maxsinf)[length(which(SsSt==maxsinf))]

# Estimation of threshold for each level in p vector, nump=# of levels
nump <- length(p)
vKMest <- vector(,length=nump)
vthresEst <- vector(,length=nump)
    
  # Use uniroot in case where end points have opposite sign
for(i in 1:nump){  
  rtfn <- function(x){
    (1-KMfn(x))-p[i]
  }
  if(rtfn(1)>0 & rtfn(maxint)<0){
    rtm <- uniroot(rtfn, interval=c(1, maxint))$root
    SindKM <- ceiling(rtm)
    vKMest[i] <- SsSt[SindKM]
  }else if(rtfn(1)<=0){
    
  # If the risk < p at lowest immune response, lowest immune response is the estimate
    
    vKMest[i] <- SsSt[1]
  }else if(rtfn(1)>0 & rtfn(maxint)>0){
    
    # If risk > p at end points, can't use uniroot, optimize method used. Two cases
    # 1) if risk > p for all S in interval [1, maxint], risk=0 < p at S[maxint+1]
    # 2) if minimum of risk-p < 0 in [1, maxint] where endpt is spot where risk-p is minimized, 
    # can find estimate by using root finding process along [1, endpt]
    
    if(optimize(rtfn, c(1, maxint))$objective>0){
      if(maxint<length(SsSt)){
      vKMest[i] <- SsSt[maxint+1]
      }else{
        vKMest[i] <- Inf
      }
    }else{
      endpt <- floor(optimize(rtfn, c(1, maxint))$minimum)
      rtm2 <- uniroot(rtfn, interval=c(1, endpt))$root
      SindKM2 <- ceiling(rtm2)
      vKMest[i] <- SsSt[SindKM2]
    }
  }else if(rtfn(maxint)==0){
    vKMest[i] <- SsSt[maxint]
  }
}  
  #Naive estimator
  SsStRev <- rev(SsSt)
  
  Wrev <- rev(WsSt)
  YWprod <- YsSt*WsSt
  RevY <- rev(YWprod)
  num <- cumsum(RevY)
  num2 <- rev(num)
  SindW <- cumsum(rev(WsSt))
  den <- rev(SindW)
  prop <- num2/den

for(i in 1:nump){
  diff <- prop-p[i]
  if(any(diff<0)){
    index <- which(diff<0)[[1]]
    vthresEst[i] <- SsSt[index]
  }else{
    vthresEst[i] <- Inf
  }
}  
  # Bootstrapping
  vthresBt <- matrix(,nrow=bts, ncol=nump)
  vKMestBt <- matrix(,nrow=bts, ncol=nump)

  if(vaccases>0){
    for(k in 1:bts){
      SamMaster <- datawghted[sample(dim(datawghted)[1], replace=TRUE),]
      Resam <- SamMaster[SamMaster[,5]!=0, ]
      ResamSSrt <- Resam[order(Resam[,1]),]
      ResamTSrt <- Resam[order(Resam[,2]),]

      ## Sorting by S
      SBtSSt <- ResamSSrt[,1]
      TBtSSt <- ResamSSrt[,2]
      YBtSSt <- ResamSSrt[,4]
      WtsBtSSt <- ResamSSrt[,5]
      # In order to deal with weights and possible duplicates from re-sampling data,
      # replicating observations based on weight
      SBtRep <- rep(SBtSSt, WtsBtSSt)
      TBtRep <- rep(TBtSSt, WtsBtSSt)
      YBtRep <- rep(YBtSSt, WtsBtSSt)

      ## Constructing Naive Estimator
      # Finding S values for those whose Y=1
      RSXY1 <- SBtRep[YBtRep==1]
      # Finding S indicies/locations of subjects whose Y=1
      Ind2 <- match(RSXY1, RSXY1)
      # Calculating # of those infected with S>s for each s observed in infected,
      # in order from smallest to largest
      Ind3 <- length(RSXY1)-Ind2+1
      # Finding S indicies for those infected (Ind4) and not infected (Ind5)
      Ind4 <- which(YBtRep==1)
      Ind5 <- which(YBtRep==0)
      # Replace S indicies with those Ind3, # of those infected with S>s for each s observed
      # in infected
      Ind6 <- replace(YBtRep, Ind4, Ind3)
      numprt <- rev(cumsum(rev(YBtRep)))
      # num counts # of those infected with S>s for each s observed (for TOTAL WEIGHTED sample)
      num <- replace(Ind6, Ind5, numprt[Ind5])
      # Number of observed infections in weighted sample
      VaccasesBt <- sum(YBtSSt)

      # Counting total subjects with S>s for each observed s, smallest to largest in weighted sample
      Sind <- match(SBtRep, SBtRep)
      den <- length(SBtRep)-Sind+1
      # prop calculates estimator at each observed s in sample, finding min. s with estimate <p
      # if no such s exists, estimate is NA
      prop <- num/den
    for(j in 1:nump){
      diff <- prop-p[j]
      if(any(diff<0)){
        indexbt <- which(diff<0)[[1]]
        vthresBt[k,j] <- SBtRep[indexbt]
      }else{
        vthresBt[k,j] <- Inf
      }
    }
      ## Constructing KM estimator
      # Counts # of dupilcates for each S then sorting by T
      dups <- rle(SBtRep)$lengths
      BtDataUnique <- cbind(unique(ResamSSrt), dups)
      # Since weights were incorportaed by replicating obs. according to weights,
      # weights column can be removed
      ResamTSrt <- BtDataUnique[order(BtDataUnique[,2]),]
      SBtTSt <- ResamTSrt[,1]
      TBtTSt <- ResamTSrt[,2]
      CBtTSt <- ResamTSrt[,3]
      YBtTSt <- ResamTSrt[,4]
      WtsBtTSt <- ResamTSrt[,6]

      # Doing same adjusted estimation process as in original sample to bootstrap sample

      KMfnBt <- function(x){
        intx <- floor(x)
        Npt <- (-1*YBtTSt[SBtTSt>=SBtSSt[intx]]*WtsBtTSt[SBtTSt>=SBtSSt[intx]])+(-1*CBtTSt[SBtTSt>=SBtSSt[intx]]*WtsBtTSt[SBtTSt>=SBtSSt[intx]])
        Npt[1] <- sum(WtsBtTSt[SBtTSt>=SBtSSt[intx]])+Npt[1]
        Ndiff <- cumsum(Npt)+(CBtTSt[SBtTSt>=SBtSSt[intx]]*WtsBtTSt[SBtTSt>=SBtSSt[intx]])
        Nleft <- (WtsBtTSt[SBtTSt>=SBtSSt[intx]]*YBtTSt[SBtTSt>=SBtSSt[intx]])+Ndiff
        KMQ <- Ndiff/Nleft
        KMEst <- cumprod(KMQ)
        tail(KMEst[is.na(KMEst)==0])[1]
      }

    maxinfbt <- max(SBtSSt[YBtSSt==1])
    maxintbt <- which(SBtSSt==maxinfbt)[length(which(SBtSSt==maxinfbt))]
    
    for(j in 1:nump){
      rtfnBt <- function(x){
        (1-KMfnBt(x))-p[j]
      }

      if(rtfnBt(1)>0 & rtfnBt(maxintbt)<0){
        rtmBt <- uniroot(rtfnBt, interval=c(1, maxintbt))$root
        SindKMBt <- ceiling(rtmBt)
        vKMestBt[k,j] <- SBtSSt[SindKMBt]
      }else if(rtfnBt(1)<=0){

        # If the risk < p initially, that is the estimate

        vKMestBt[k,j] <- SBtSSt[1]
      }else if(rtfnBt(1)>0 & rtfnBt(maxintbt)>0){

        # If risk > p at end points, can't use uniroot, optimize method used as before

        if(optimize(rtfnBt, c(1, maxintbt))$objective>0){
          if(maxintbt<length(SBtSSt)){
            vKMestBt[k,j] <- SBtSSt[maxintbt+1]
          }else{
            vKMestBt[k,j] <- Inf
          }
        }else{
          endptbt <- floor(optimize(rtfnBt, c(1, maxintbt))$minimum)
          rtm2bt <- uniroot(rtfnBt, interval=c(1, endptbt))$root
          SindKM2Bt <- ceiling(rtm2bt)
          vKMestBt[k,j] <- SBtSSt[SindKM2Bt]
        }
        }else if(rtfnBt(maxintbt)==0){
          vKMestBt[k,j] <- SBtSSt[maxintbt]
        }
      }
    }
  }


    ## Constructing CI using quantiles
    # From naive estimator
  LLThres <- vector(,length=nump)
  ULThres <- vector(,length=nump)
  LLThresKM <- vector(,length=nump)
  ULThresKM <- vector(,length=nump)

  for(j in 1:nump){
    CLEndPts <- quantile(vthresBt[,j], c((alpha/2), 1-(alpha/2)))
    expand <- as.numeric(unlist(CLEndPts))
    LLThres[j] <- expand[1]
    ULThres[j] <- expand[2]

    # From KM Estimator
    CLEndPtsKM <- quantile(vKMestBt[,j], c((alpha/2), 1-(alpha/2)))
    expandKM <- as.numeric(unlist(CLEndPtsKM))
    LLThresKM[j] <- expandKM[1]
    ULThresKM[j] <- expandKM[2]
}
  # If no cases observed, no info. provided to construct CI so CI endpoints are labeled NA
  if(vaccases==0){
    LLThres[1:nump] <- NA
    ULThres[1:nump] <- NA
    LLThresKM[1:nump] <- NA
    ULThresKM[1:nump] <- NA
  }

  ## Adding CI width variables ##
  width <- vector(,length=nump)
  widthKM <- vector(,length=nump)

for(j in 1:nump){
  width[j] <- abs(ULThres[j]-LLThres[j])
  widthKM[j] <- abs(ULThresKM[j]-LLThresKM[j])
}

naive <- data.frame(matrix(nrow=length(p), ncol=6))
adjusted <- data.frame(matrix(nrow=length(p), ncol=6))
names(naive) <- c('Risk_Level','Method','Est', 'CI_LL', 'CI_UL', 'CI_Width')
names(adjusted) <- c('Risk_Level','Method','Est', 'CI_LL', 'CI_UL', 'CI_Width')

# Combining results into list for printing
for(j in 1:nump){
  naive[j,] <- c(p[j],NA,round(vthresEst[j], rnd_digits), round(LLThres[j], rnd_digits), 
                  round(ULThres[j], rnd_digits), round(width[j], rnd_digits))
  adjusted[j,] <- c(p[j],NA,round(vKMest[j], rnd_digits), round(LLThresKM[j], rnd_digits), 
                     round(ULThresKM[j], rnd_digits), round(widthKM[j], rnd_digits))
}
  naive$Method <- "Naive"
  adjusted$Method <- "Adjusted"

  if(include_naive==TRUE){
    results <- data.frame(rbind(adjusted,naive))
  }
  if(include_naive==FALSE){
    results <- data.frame(adjusted)
  }

return(results)
}
###### END FUNCTION ######