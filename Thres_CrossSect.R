Thres_CrossSect <- function(dataset, p, bts=2000, alpha=0.025, weighted=FALSE, rnd_digits=2){
  # dataset: data
  # p: vector of risk levels of interest
  # bts: number of bootstrap resamples
  # alpha: alpha level for confidence intervals (two sided)
  # weighted: logical value to indicate use of weighted methods, if false then all missing immune response values are removed from data
  # rnd_digits: number of digits to round results to

  # Note: Data needs to be of columns in following order: 
  #  1) Obs. immune response (S) 2) Disease status
  
  ## Check for valid arguments
  if(any(p==0)&weighted==TRUE){stop("For 0 risk level, only unweighted methods are implemented")}
  if(any(p<0) | any(p>1)){stop("Risk levels p all need to be between 0 and 1")}
  if(alpha>=1){stop("Alpha needs to be < 1; alpha/2 used for coverage probability")}
  if(!(weighted %in% c(TRUE,FALSE))){stop("Please choose a logical value for weighted (default is FALSE)")}
  if(is.numeric(bts)==FALSE|bts<1){stop("Please choose a numeric value for the number of bootstrap resamples (bts) which is >=1")}
  bts <- floor(bts) # forces boots to be an integer
  ## Adjust dataset based on weighted or non-weighted
  if(weighted==TRUE){
    p0_hat <- mean(!is.na(dataset[dataset[,2]==0,1]))
    p1_hat <- mean(!is.na(dataset[dataset[,2]==1,1]))
    dataset[,3] <- ifelse(is.na(dataset[,1])==0, 1, 0)
  }
  if(weighted==FALSE){
    dataset <- dataset[complete.cases(dataset),]
    p0_hat <- 1
    p1_hat <- 1
    dataset[,3] <- 1
  }
  ## Setting up weights ##
  IPW <- function(x){
    term <- ifelse(x==0, p0_hat, p1_hat)
    1/term
  }
  
  ## Sorting data matrix.  
  Wghts <- unlist(lapply(dataset$Y, IPW))
  datawghted <- cbind(dataset[,-c(3)], Wghts*dataset[,3])
  vaccases <- sum(datawghted$Y)
  
  datacomp <- datawghted[datawghted[,3]!=0,]
  # Creating nonparametric estimator
    vthresEst <- rep(NA,length=length(p))
    LLThres <- rep(NA,length=length(p))
    ULThres <- rep(NA,length=length(p))
    
    SSt <- datacomp[order(datacomp[,1]),1]
    SStRev <- rev(SSt)
    Ymatch <- datacomp[order(datacomp[,1]),2]
    Wmatch <- datacomp[order(datacomp[,1]),3]
    Wrev <- rev(Wmatch)
    YWprod <- Ymatch*Wmatch
    RevY <- rev(YWprod)
    num <- cumsum(RevY)
    num2 <- rev(num)
    SindW <- cumsum(rev(Wmatch))
    den <- rev(SindW)
    prop <- num2/den
    for(j in 1:length(p)){
      if(p[j]!=0){
      diff <- prop-p[j]
      if(vaccases>0){
        if(any(diff<0)){
          index <- which(diff<0)[[1]]
          vthresEst[j] <- SSt[index]
        }else{
          vthresEst[j] <- Inf
        }
      }else{vthresEst[j] <- NA}
      }
      if(p[j]==0){
        if(vaccases>1){
          maxs <- max(SSt[Ymatch==1])
          nexts <- max(SSt[Ymatch==1&SSt<maxs])
          SVst <- sort(SSt[Ymatch==1], decreasing=FALSE)
          sum <- 0
          for(k in 1:(vaccases-1)){
            sum<-sum+((k/vaccases)^(vaccases))*(SVst[k+1]-SVst[k])
          }
          vthresEst[j] <- maxs+sum
          ULThres[j] <- maxs+(((1-0.5*alpha)^(-1)-1)^(-1))*(maxs-nexts)
          LLThres[j] <- maxs+(((0.5*alpha)^(-1)-1)^(-1))*(maxs-nexts)
        }else{
          vthresEst[j] <- NA
          ULThres[j] <- NA
          LLThres[j] <- NA
          }
        }
    }
  ## Calc CI using bootstrap
    # Set-up empty objects for bootstrap results
    VaccasesBt <- rep(NA, length=bts)
    vthresBt <- matrix(,nrow=bts, ncol=length(p))
    width <- rep(NA,length=length(p))
    
    if(vaccases>0){
      if(any(p>0)){
      for (k in 1:bts){
        
        SamMaster <- datawghted[sample(dim(datawghted)[1], replace=TRUE),]
        Resam <- SamMaster[SamMaster[,3]!=0, ]
        ResamSrt <- Resam[order(Resam[,1]),]
        SBtSt <- ResamSrt[,1]
        YBt <- ResamSrt[,2]
        WtsBt <- ResamSrt[,3]
        SBtRep <- rep(SBtSt, WtsBt)
        YBtRep <- rep(YBt, WtsBt)

        RSXY1 <- SBtRep[YBtRep==1]
        Ind2 <- match(RSXY1, RSXY1)
        Ind3 <- length(RSXY1)-Ind2+1
        Ind4 <- which(YBtRep==1)
        Ind5 <- which(YBtRep==0)
        Ind6 <- replace(YBtRep, Ind4, Ind3)
        numprt <- rev(cumsum(rev(YBtRep)))
        num <- replace(Ind6, Ind5, numprt[Ind5])
        SBtStRev <- rev(SBtSt)
        VaccasesBt[k] <- sum(YBt)

        # Creating non-parametric estimator
        Sind <- match(SBtRep, SBtRep)
        den <- length(SBtRep)-Sind+1
        prop <- num/den
      for(j in 1:length(p)){
        if(p[j]>0){
          diff <- prop-p[j]
          if(any(diff<0)){
            indexbt <- which(diff<0)[[1]]
            vthresBt[k,j] <- SBtRep[indexbt]
          }else{
            vthresBt[k,j] <- Inf
          }
        }else{
          vthresBt[k,j] <- NA
        }  
      }  
    }
  }
      # CI using nonparametric estimator #
    for(j in which(p>0)){  
      CLEndPts <- quantile(vthresBt[,j], c((alpha/2), 1-(alpha/2)))
      expand <- as.numeric(unlist(CLEndPts))
      LLThres[j] <- expand[1]
      ULThres[j] <- expand[2]
      }
    }
    if(vaccases==0){
      LLThres[which(p>0)] <- NA
      ULThres[which(p>0)] <- NA
    }
    
    for(j in 1:length(p)){
      width[j] <- abs(ULThres[j]-LLThres[j])
    }
    
    results <- data.frame(matrix(nrow=length(p), ncol=5))
    names(results) <- c('Risk_Level','Est', 'CI_LL', 'CI_UL', 'CI_Width')

    # Combining results into list for printing
    for(j in 1:length(p)){
      results[j,] <- c(p[j],round(vthresEst[j], rnd_digits), round(LLThres[j], rnd_digits), 
                     round(ULThres[j], rnd_digits), round(width[j], rnd_digits))
    }
    
    return(results)
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  