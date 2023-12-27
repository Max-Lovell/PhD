# Adaptation in R by Max Lovell of the following MATLAB code by Brian Maniscalco: 
  #www.columbia.edu/~bsm2105/type2sdt/fit_meta_d_MLE.m

library(NlcOptim)
#example from original code for comparison:
  #nR_S1 <- c(100,50,20,10,5,1) ; nR_S2 <- c(3,7,8,12,27,89)
  #fit_meta_d_MLE(nR_S1,nR_S2)

#nR_S1=nr_sx_adj[i,grep('S1',colnames(nr_sx))];nR_S2=nr_sx_adj[i,grep('S2',colnames(nr_sx))]
fit_meta_d_MLE <- function(nR_S1, nR_S2){ # nR_S1=nr_sx[i,grep('S1',colnames(nr_sx))]+.5; nR_S2=nr_sx[i,grep('S2',colnames(nr_sx))]+.5
  
  ## Maximum Likelihood Estimator
  fit_meta_d_logL <- function(parameters){ # parameters<-guess
    #setup
    meta_d1  <- parameters[1]
    t2c1     <- parameters[-1]
    S1mu <- -meta_d1/2
    S2mu <-  meta_d1/2
    S1sd <- 1;  S2sd <- S1sd/s
    S1mu <- S1mu - eval(parse(text=constant_criterion))
    S2mu <- S2mu - eval(parse(text=constant_criterion))
    t1c1 <- 0
    
    #likelihood uses multinomial model: Each discrete type 2 outcome y,s,r is 
    #an event with a fixed probability that occurred a certain number of times. 
    #The probability of the entire set of type 2 outcomes across all trials is then 
    #proportional to the product of the probability of each individual type 2 outcome
    
    #Total new type-1 hit,miss,FA,CR considering new d' values
    C_area_rS1 <- pnorm(t1c1,S1mu,S1sd)
    I_area_rS1 <- pnorm(t1c1,S2mu,S2sd)
    C_area_rS2 <- 1-pnorm(t1c1,S2mu,S2sd)
    I_area_rS2 <- 1-pnorm(t1c1,S1mu,S1sd)
    
    #new t2 c values
    t2c1x <- c(-Inf, t2c1[1:(nRatings-1)],t1c1,t2c1[nRatings:length(t2c1)],Inf)
    #hit,miss,FA,CR for each confidence level
    prC_rS1<-c(); prI_rS1<-c(); prC_rS2<-c(); prI_rS2<-c()
    for(i in 1:nRatings){
      prC_rS1[i] <- (pnorm(t2c1x[i+1],S1mu,S1sd) - pnorm(t2c1x[i],S1mu,S1sd)) / C_area_rS1
      prI_rS1[i] <- (pnorm(t2c1x[i+1],S2mu,S2sd) - pnorm(t2c1x[i],S2mu,S2sd)) / I_area_rS1
      prC_rS2[i] <- ((1-pnorm(t2c1x[nRatings+i],S2mu,S2sd)) - (1-pnorm(t2c1x[nRatings+i+1],S2mu,S2sd))) / C_area_rS2
      prI_rS2[i] <- ((1-pnorm(t2c1x[nRatings+i],S1mu,S1sd)) - (1-pnorm(t2c1x[nRatings+i+1],S1mu,S1sd))) / I_area_rS2
    }
    
    #t2 resp count - count each time each response type actually occurred in the dataset
    nC_rS1 <- nR_S1[1:nRatings]
    nC_rS2 <- nR_S2[(nRatings+1):length(nR_S2)]
    nI_rS1 <- nR_S2[1:nRatings]
    nI_rS2 <- nR_S1[(nRatings+1):length(nR_S1)]
    
    #likelihood of ocurring
    logL <- sum(nC_rS1*log(prC_rS1)) + sum(nI_rS1*log(prI_rS1)) + sum(nC_rS2*log(prC_rS2)) + sum(nI_rS2*log(prI_rS2))
    if(is.nan(logL)){ logL=-Inf }
    logL <- -logL
  }
  
  
  # vars for use in MLE fit
  s <- 1
  nRatings <- length(nR_S1)/2
  nCriteria <- 2*nRatings-1
  constant_criterion <- 'meta_d1 * (t1c1 / d1)'
  
  #Linear constraints
  A <- matrix(0, nCriteria-2, nCriteria)
  A[cbind(1:nrow(A), 2:(nCriteria-1))] <- 1
  A[cbind(1:nrow(A), 3:nCriteria)] <- -1
  b <- rep(-0.00001,nCriteria-2)
  LB <- c(-10, rep(c(-20,0),each=(nCriteria-1)/2))
  UB <- c(10, rep(c(0,20),each=(nCriteria-1)/2))
  
  # guess
  ratingHR <- c()
  ratingFAR <- c()
  for(c in 2:(nRatings*2)){
    ratingHR[c-1] <- sum(nR_S2[c:length(nR_S2)]) / sum(nR_S2)
    ratingFAR[c-1] <- sum(nR_S1[c:length(nR_S1)]) / sum(nR_S1)
  }
  d1 <- qnorm(ratingHR[nRatings])-qnorm(ratingFAR[nRatings])
  meta_d1 <- d1
  c1 <- -.5*(qnorm(ratingHR)+qnorm(ratingFAR))
  t1c1 <- c1[nRatings]
  t2c1 <- c1[-nRatings]
  guess <- c(meta_d1, t2c1-eval(parse(text=constant_criterion)))
  
  #fit
  xf <- suppressWarnings(solnl(guess, objfun=fit_meta_d_logL, A=A, B=b, lb=LB, ub=UB,maxnFun=100000))
  meta_d1 <- xf$par[1]
  t2c1 <- xf$par[-1]+eval(parse(text=constant_criterion))
  logL <- -xf$fn
  
  #observed t2 HR/FAR
  I_nR_rS2 <- nR_S1[(nRatings+1):length(nR_S1)]
  I_nR_rS1 <- nR_S2[nRatings:1]
  C_nR_rS2 <- nR_S2[(nRatings+1):length(nR_S2)]
  C_nR_rS1 <- nR_S1[nRatings:1]
  obs_FAR2_rS2<-c(); obs_HR2_rS2<-c(); obs_FAR2_rS1<-c(); obs_HR2_rS1<-c();
  for(i in 2:nRatings){
    obs_FAR2_rS2[i-1] <- sum(I_nR_rS2[i:length(I_nR_rS2)])/sum(I_nR_rS2)
    obs_HR2_rS2[i-1]  <- sum(C_nR_rS2[i:length(C_nR_rS2)])/sum(C_nR_rS2)
    obs_FAR2_rS1[i-1] <- sum(I_nR_rS1[i:length(I_nR_rS1)])/sum(I_nR_rS1)
    obs_HR2_rS1[i-1]  <- sum(C_nR_rS1[i:length(C_nR_rS1)])/sum(C_nR_rS1)
  }
  S1mu <- -meta_d1/2
  S2mu <-  meta_d1/2
  S1sd <- 1;  S2sd <- S1sd/s
  mt1c1 <- eval(parse(text=constant_criterion))
  C_area_rS2 <- 1-pnorm(mt1c1,S2mu,S2sd)
  I_area_rS2 <- 1-pnorm(mt1c1,S1mu,S1sd)
  C_area_rS1 <- pnorm(mt1c1,S1mu,S1sd)
  I_area_rS1 <- pnorm(mt1c1,S2mu,S2sd)
  
  est_FAR2_rS2<-c(); est_HR2_rS2<-c(); est_FAR2_rS1<-c(); est_HR2_rS1<-c();
  for(i in 1:(nRatings-1)){
    t2c1_lower <- t2c1[nRatings-i]
    t2c1_upper <- t2c1[nRatings-1+i]
    
    I_FAR_area_rS2 <- 1-pnorm(t2c1_upper,S1mu,S1sd)
    C_HR_area_rS2  <- 1-pnorm(t2c1_upper,S2mu,S2sd)
    I_FAR_area_rS1 <- pnorm(t2c1_lower,S2mu,S2sd)
    C_HR_area_rS1  <- pnorm(t2c1_lower,S1mu,S1sd)

    est_FAR2_rS2[i] <- I_FAR_area_rS2 / I_area_rS2
    est_HR2_rS2[i]  <- C_HR_area_rS2 / C_area_rS2
    est_FAR2_rS1[i] <- I_FAR_area_rS1 / I_area_rS1
    est_HR2_rS1[i]  <- C_HR_area_rS1 / C_area_rS1
  }
  
  # Output
  fit <- c()
  fit$da        <- sqrt(2/(1+s^2))*s*d1
  fit$s         <- s
  fit$meta_da   <- sqrt(2/(1+s^2))*s*meta_d1
  fit$M_diff    <- fit$meta_da-fit$da
  fit$M_ratio   <- fit$meta_da/fit$da
  fit$meta_ca   <- (sqrt(2)*s/sqrt(1+s^2))*eval(parse(text=constant_criterion))
  t2ca          <- (sqrt(2)*s/sqrt(1+s^2))*t2c1
  fit$t2ca_rS1  <- t2ca[1:(nRatings-1)]
  fit$t2ca_rS2  <- t2ca[nRatings:length(t2ca)]
  fit$S1units$d1        <- d1
  fit$S1units$meta_d1   <- meta_d1
  fit$S1units$s         <- s
  fit$S1units$meta_c1   <- eval(parse(text=constant_criterion))
  fit$S1units$t2c1_rS1  <- t2c1[1:(nRatings-1)]
  fit$S1units$t2c1_rS2  <- t2c1[nRatings:length(t2c1)]
  fit$logL         <- logL
  fit$est_HR2_rS1  <- est_HR2_rS1
  fit$obs_HR2_rS1  <- obs_HR2_rS1
  fit$est_FAR2_rS1 <- est_FAR2_rS1
  fit$obs_FAR2_rS1 <- obs_FAR2_rS1
  fit$est_HR2_rS2  <- est_HR2_rS2
  fit$obs_HR2_rS2  <- obs_HR2_rS2
  fit$est_FAR2_rS2 <- est_FAR2_rS2
  fit$obs_FAR2_rS2 <- obs_FAR2_rS2
  return(fit)
}
