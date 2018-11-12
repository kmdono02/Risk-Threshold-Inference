---
title: 'Threshold of Risk Inference: Examples with Simulated Data'
author: "Kevin Donovan"
date: "November 9, 2018"
output: 
  html_document:
    keep_md: TRUE
---





# Cross-sectional Data
## Positive risk levels
Let's consider the scenario where inference on risk thresholds is done with cross-sectional data collected at a single time point (as detailed in Sections 2 and 3 of the paper).  Let $Y$ denote disease status ($Y=1$ denoting disease, $Y=1$ denoting no disease) and $S$ denote immune response.  Suppose $S$ is observed for proportion $p_0$ of subjects with $Y=0$ and proportion $p_1$ of subjects with $Y=1$.  In the following simulated examples, sample size $n=12,500$ where $S$ is generated from a gamma distribution with mean and variance 4 and $Y|S$ is generated from a Bernoulli distribution with $\Pr(Y=1|S)=$exbit$(\beta_0+\beta_1S)$.  The slope $\beta_1$ was set to -5 and $\beta_0$ was determined based on the chosen disease risk $\Pr(Y=1)$; in the following simulations risks of 0.01 and 0.1 were used with $p_0=0.2$ and $p_1=1$.  The R function sim_data_logistic found in the GitHub repository simulates data under this model (based on choosen values of $n, \beta_1, t_0$, $p_0$, $p_1$, and disease risk).  

For disease risk 0.01, point estimates and 95$\%$ confidence intervals are computed of thresholds for risk levels 0.001, 0.003, 0.005, 0.007, and 0.009 and for disease risk 0.1, point estimates and 95$\%$ confidence intervals are computed of thresholds for risk levels 0.009, 0.01, 0.03, 0.05, 0.07, and 0.09.  The confidence intervals are computed using 2000 bootstrap resamples.  These calculations were done using the R function Thres_Surv, whose source code can be found in the GitHub repository, which implements the methods detailed for Section 4 for inference with time to event data.  These estimates and confidence intervals were computed using $\hat{v}_c^{w}$ as defined in Section 3 of the paper.  For each threshold, the point estimate and confidence interval using $\hat{v}_c^{w}$ are plotted with the corresponding risk level on x axis; the point estimate is colored red and the true threshold is colored blue.

The true thresholds under the model from which the data were simulated (models detailed previously) are also provided (function sim_data_logit calculates these for vector of choosen risk levels).


```r
source("sim_data_logit.R")
source("Thres_CrossSect.R")

## Low risk example
sim_data <- Sim_data_logit(Risk = 0.01, n=12500, p=c(0.001,0.003,0.005,0.007,0.009), p0=0.2, p1=1, seed_val = 012)

vac_data_analysis_LR <- sim_data$simmed_data %>% select(S_obs, Y) 
low_risk_cross_sec_results <- Thres_CrossSect(dataset = vac_data_analysis_LR, p=c(0.001,0.003,0.005,0.007,0.009), weighted=TRUE, bts=2000)
low_risk_cross_sec_results
```

```
##   Risk_Level  Est CI_LL CI_UL CI_Width
## 1      0.001 1.27  1.15  1.47     0.32
## 2      0.003 1.01  0.90  1.12     0.22
## 3      0.005 0.85  0.75  0.93     0.18
## 4      0.007 0.72  0.55  0.82     0.27
## 5      0.009 0.53  0.22  0.71     0.49
```

```r
# Now plot
plotCI(x=low_risk_cross_sec_results$Risk_Level, y=low_risk_cross_sec_results$Est, 
       ui=low_risk_cross_sec_results$CI_UL,li=low_risk_cross_sec_results$CI_LL, xlab="Risk Level", ylab="Threshold Estimate", col="red",scol="black", main="Threshold of risk inference on simulated data: disease risk=0.01\nCross sectional scenario", ylim=c(0.25, 1.75))
points(x=as.numeric(names(sim_data$true_thresholds)), y=sim_data$true_thresholds, col="blue")
```

![](Thres_SimData_Exs_files/figure-html/cross_sec_logit-1.png)<!-- -->

```r
## Medium risk example
sim_data <- Sim_data_logit(Risk = 0.1, n=12500, p=c(0.009, 0.01, 0.03, 0.05, 0.07, 0.09), p0=0.2, p1=1, seed_val = 012)

vac_data_analysis_MR <- sim_data$simmed_data %>% select(S_obs, Y) 
med_risk_cross_sec_results <- Thres_CrossSect(dataset = vac_data_analysis_MR, p=c(0.009, 0.01, 0.03, 0.05, 0.07, 0.09), weighted=TRUE, bts=2000)
med_risk_cross_sec_results
```

```
##   Risk_Level  Est CI_LL CI_UL CI_Width
## 1      0.009 1.99  1.94  2.04     0.11
## 2      0.010 1.96  1.92  2.02     0.11
## 3      0.030 1.63  1.60  1.68     0.08
## 4      0.050 1.43  1.37  1.48     0.10
## 5      0.070 1.22  1.13  1.29     0.16
## 6      0.090 0.90  0.68  1.03     0.35
```

```r
# Now plot
plotCI(x=med_risk_cross_sec_results$Risk_Level, y=med_risk_cross_sec_results$Est, 
       ui=med_risk_cross_sec_results$CI_UL, li=med_risk_cross_sec_results$CI_LL, xlab="Risk Level",
       ylab="Threshold Estimate",col="red",scol="black", main="Threshold of risk inference on simulated data: disease risk=0.1\nTime to event scenario")
points(x=as.numeric(names(sim_data$true_thresholds)), y=sim_data$true_thresholds, col="blue", sub="Cross sectional scenario")
```

![](Thres_SimData_Exs_files/figure-html/cross_sec_logit-2.png)<!-- -->

## Zero risk level
For the cross sectional data scenario, suppose $\Pr(Y=1|S)=$exbit$(\beta_0+\beta_1S)$ for $S<v_0$ and $\Pr(Y=1|S)=0$ otherwise (as described in Section 2 of the paper).  In the following simulated corss-sectional examples, everything else with the previous cross-sectional simulations is unchanged.  The R function sim_data_logit found in the GitHub repository simulates data under this model using the argument **theta** (denotes the value of $v_0$ to simulate data under).  The nonparametric estimator $\tilde{v}_0$ is computed, with the corresponding bootstrap confidence intervals as detailed in Section 2 of the paper.

The two simulated datasets used in the following examples, one simulated under a 0.01 disease risk and the other under 0.10 disease risk with $v_0$=2 for both scenarios, using the function sim_data_logit provided in the GitHub repository.


```r
source("sim_data_logit.R")
source("Thres_CrossSect.R")

## Low risk example
sim_data <- Sim_data_logit(Risk = 0.01, n=12500, theta=2, p0=1, p1=1, p=0, seed_val = 012)

vac_data_analysis_LR <- sim_data$simmed_data %>% select(S_obs, Y) 
low_risk_cross_sec_results <- Thres_CrossSect(dataset = vac_data_analysis_LR, p=0, weighted=FALSE, bts=2000)
low_risk_cross_sec_results
```

```
##   Risk_Level  Est CI_LL CI_UL CI_Width
## 1          0 1.96  1.91  7.24     5.33
```

```r
## Medium risk example
sim_data <- Sim_data_logit(Risk = 0.1, n=12500, theta=2, p0=0.2, p1=1, p=0,seed_val = 012)

vac_data_analysis_MR <- sim_data$simmed_data %>% select(S_obs, Y) 
med_risk_cross_sec_results <- Thres_CrossSect(dataset = vac_data_analysis_MR, p=0, weighted=FALSE, bts=2000)
med_risk_cross_sec_results
```

```
##   Risk_Level Est CI_LL CI_UL CI_Width
## 1          0   2  1.99  2.17     0.17
```

# Time to Event Data
Now, let's consider the scenario where inference on risk thresholds is done with time to event data (as described in Section 4 of the paper).  Let time $t_0$ denote end of follow-up ($t_0>0$), $T$ denote survival time, $S$ denote immune response, $C$ denote censoring status ($C=1$ denoting right censoring at or before time $t_0$, and $C=0$ denoting observed disease before time $t_0$).  Suppose $S$ is observed for proportion $p_0$ of subjects with $C=1$ and proportion $p_1$ of subjects with $C=0$.  In the following simulated examples, sample size $n=12,500$ where $S$ is generated from a gamma distribution with mean and variance 4 and $T|S$ is generated by log$(T)=\beta_0+\beta_1S+\epsilon$ where $\epsilon$ is generated from a logistic distribution.  The slope $\beta_1$ sas set to 2 and $\beta_0$ was determined based on the chosen disease risk $\Pr(T<t_0)$; in the following simulations risks of 0.01 and 0.1 were used with $p_0=0.2$ and $p_1=1$.  The R function sim_data_logistic found in the GitHub repository simulates data under this model (based on choosen values of $n, \beta_1, t_0$, $p_0$, $p_1$, and disease risk).  

For disease risk 0.01, point estimates and 95$\%$ confidence intervals are computed of thresholds for risk levels 0.001, 0.003, 0.005, 0.007, and 0.009 and for disease risk 0.1, point estimates and 95$\%$ confidence intervals are computed of thresholds for risk levels 0.009, 0.01, 0.03, 0.05, 0.07, and 0.09.  The confidence intervals are computed using 2000 bootstrap resamples.  These calculations were done using the R function Thres_Surv, whose source code can be found in the GitHub repository, which implements the methods detailed for Section 4 for inference with time to event data.  These estimates and confidence intervals were computed using $\hat{v}_c^{w}$ (referred to as the "naive" method in the R function Thres_Surv) and $\hat{v}_c^{wk}$ ("adjusted" method), as defined in Section 4 of the paper.  For each threshold, the point estimate and confidence interval using $\hat{v}_c^{wk}$ are plotted with the corresponding risk level on x axis; the point estimate is colored red and the true threshold is colored blue.

The true thresholds under the model from which the data were simulated (models detailed previously) are also provided (function sim_data_logistic calculates these for vector of choosen risk levels).


```r
source("sim_data_logistic.R")
source("Thres_Surv.R")

## Low risk example
sim_data <- Sim_data_logistic(Risk = 0.01, n=12500, p=c(0.001,0.003,0.005,0.007,0.009), p0=0.2, p1=1, seed_val = 012)

vac_data_analysis_LR <- sim_data$simmed_data %>% select(S_obs, TObs, Cstatus) 
low_risk_tte_results <- Thres_Surv(dataset = vac_data_analysis_LR, p=c(0.001,0.003,0.005,0.007,0.009), weighted=TRUE, bts=2000)
low_risk_tte_results
```

```
##    Risk_Level   Method  Est CI_LL CI_UL CI_Width
## 1       0.001 Adjusted 2.07  1.91  3.05     1.15
## 2       0.003 Adjusted 1.79  1.34  1.98     0.64
## 3       0.005 Adjusted 1.26  0.94  1.76     0.83
## 4       0.007 Adjusted 0.96  0.51  1.28     0.77
## 5       0.009 Adjusted 0.67  0.34  1.03     0.69
## 6       0.001    Naive 2.00  1.79  2.21     0.42
## 7       0.003    Naive 1.39  1.02  1.79     0.76
## 8       0.005    Naive 0.89  0.34  1.19     0.85
## 9       0.007    Naive 0.34  0.34  0.77     0.43
## 10      0.009    Naive 0.34  0.34  0.46     0.12
```

```r
adjusted_results_low <- low_risk_tte_results %>% filter(Method=="Adjusted")

# Now plot
plotCI(x=adjusted_results_low$Risk_Level, y=adjusted_results_low$Est, 
       ui=adjusted_results_low$CI_UL,li=adjusted_results_low$CI_LL, xlab="Risk Level", ylab="Threshold Estimate",
       col="red",scol="black", main="Threshold of risk inference on simulated data: disease risk=0.01\nTime to event scenario")
points(x=as.numeric(names(sim_data$true_thresholds)), y=sim_data$true_thresholds, col="blue")
```

![](Thres_SimData_Exs_files/figure-html/time_to_event_logistic-1.png)<!-- -->

```r
## Medium risk example
sim_data <- Sim_data_logistic(Risk = 0.1, n=12500, p=c(0.009, 0.01, 0.03, 0.05, 0.07, 0.09), p0=0.2, p1=1, seed_val = 012)

vac_data_analysis_MR <- sim_data$simmed_data %>% select(S_obs, TObs, Cstatus) 
med_risk_tte_results <- Thres_Surv(dataset = vac_data_analysis_MR, p=c(0.009, 0.01, 0.03, 0.05, 0.07, 0.09), 
                                   weighted=TRUE, bts=2000)
med_risk_tte_results
```

```
##    Risk_Level   Method  Est CI_LL CI_UL CI_Width
## 1       0.009 Adjusted 2.94  2.75  3.08     0.33
## 2       0.010 Adjusted 2.88  2.71  3.06     0.35
## 3       0.030 Adjusted 2.14  2.02  2.26     0.24
## 4       0.050 Adjusted 1.76  1.64  1.90     0.26
## 5       0.070 Adjusted 1.41  1.28  1.53     0.25
## 6       0.090 Adjusted 1.01  0.62  1.19     0.56
## 7       0.009    Naive 2.71  2.52  2.86     0.35
## 8       0.010    Naive 2.63  2.47  2.77     0.29
## 9       0.030    Naive 1.85  1.75  1.95     0.20
## 10      0.050    Naive 1.33  1.20  1.44     0.24
## 11      0.070    Naive 0.51  0.27  0.90     0.63
## 12      0.090    Naive 0.27  0.27  0.34     0.07
```

```r
adjusted_results_med <- med_risk_tte_results %>% filter(Method=="Adjusted")

# Now plot
plotCI(x=adjusted_results_med$Risk_Level, y=adjusted_results_med$Est, 
       ui=adjusted_results_med$CI_UL, li=adjusted_results_med$CI_LL, xlab="Risk Level", ylab="Threshold Estimate",col="red",scol="black", main="Threshold of risk inference on simulated data: disease risk=0.1\nTime to event scenario")
points(x=as.numeric(names(sim_data$true_thresholds)), y=sim_data$true_thresholds, col="blue", sub="Time to event scenario")
```

![](Thres_SimData_Exs_files/figure-html/time_to_event_logistic-2.png)<!-- -->
