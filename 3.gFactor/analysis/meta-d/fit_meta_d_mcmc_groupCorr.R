# Adaptation in R by Max Lovell of the following MATLAB code by Steve Fleming: 
  # https://github.com/metacoglab/HMeta-d/blob/master/Matlab/fit_meta_d_mcmc_groupCorr.m
library(rjags)

# Example----------
#example data from extract_corr_sample.m, which follows exampleFit_corr.m
  #setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  #corr_data <- read.csv("example_corr_data.csv",header=F)
  #nResponse <- (ncol(corr_data)-1)/4
  #names(corr_data) <- c("task_n",
  #                      paste(rep(c("R1_S1","R2_S1","R1_S2","R2_S2"),each=nResponse), 
  #                            rep(c(nResponse:1,1:nResponse),2), 
  #                            sep="_")
  #                      )
  #
  #nR_S1_raw <- corr_data[,1:(nResponse*2+1)]
  #nR_S2_raw <- corr_data[,c(1,(nResponse*2+2):ncol(corr_data))]
  #nR_S1 <- split(nR_S1_raw[,-1], nR_S1_raw$task_n)
  #nR_S2 <- split(nR_S2_raw[,-1], nR_S2_raw$task_n)
#
  #mcmc_params <- list(
  #  response_conditional = 0, # response-conditional meta-d? 0=no,1=both simultaneously, 2=alternating
  #  rc_alternating = c(0,1,2,2), # if r_c=2, sets response types. 0=non-rc, 1=rS_1, 2=rS_2
  #  estimate_dprime = 0, # also estimate dprime in same model?
  #  nchains = 3, # How Many Chains?
  #  nadapt = 1000, # How many adaptation samples? rjags default is 1000
  #  nburnin = 1000, # How Many Burn-in Samples?
  #  nsamples = 10000,  # How Many Recorded Samples? 30000 recommended
  #  nthin = 1, # How Often is a Sample Recorded?
  #  dic = 0 # note DIC is computationally costly
  #)
  #
  #mcmc_params$response_conditional <- 0
  #corr <- fit_meta_d_mcmc_groupCorr(nR_S1, nR_S2,mcmc_params)
  #mcmc_params$response_conditional <- 1
  #both_rc <- fit_meta_d_mcmc_groupCorr(nR_S1, nR_S2,mcmc_params)
  #mcmc_params$response_conditional <- 2
  #mcmc_params$rc_alternating <- c(0,1,2,2)
  #mixed_rc <- fit_meta_d_mcmc_groupCorr(nR_S1, nR_S2,mcmc_params)
  #
  #for(fit in list(corr,mixed_rc)){
  #  print(fit$summary$group_level)
  #  print(fit$summary$correlations)
  #  print(fit$summary$subject_level_means)
  #  prsfs <- fit$mcmc$Rhat$psrf[,2]
  #  print(prsfs[prsfs>1.1]) #convergence is lower for alternating RC
  #}
  
# JAGS model creators###############
#Standard
createJAGSModel <- function(Ntask){ #dots,breath,
  data_mod<-c(); l_multi<-c(); group_level<-c(); rhos<-c(); vcov_mat<-c();
  #number of rhos is different
  Nrho <- factorial(Ntask)/(factorial(Ntask-2)*factorial(2))
  rhos <- paste0(
  '\n\tfor(r in 1:',Nrho,'){',
  '\n\t\trho[r] ~ dunif(-1,1)',
  '\n\t}')
  #create trial number variable strings
  for(t in 1:Ntask){ #t<-1 ; https://rstudio-pubs-static.s3.amazonaws.com/272658_ae4d482c86514674be17042c852ebbfc.html
    count_name <- paste0("counts",t)
    data_mod <- c(data_mod, paste0(
      "\n\t\tH[s,",t,"] <- sum(",count_name,"[s,(nratings[",t,"]*3+1):(nratings[",t,"]*4)])",
      "\n\t\tM[s,",t,"] <- sum(",count_name,"[s,(nratings[",t,"]*2+1):(nratings[",t,"]*3)])",
      "\n\t\tFA[s,",t,"] <- sum(",count_name,"[s,(nratings[",t,"]+1):(nratings[",t,"]*2)])",
      "\n\t\tCR[s,",t,"] <- sum(",count_name,"[s,1:nratings[",t,"]])",
      "\n"))
    
    l_multi <- c(l_multi,paste0(
      "\n\t\t",count_name,"[s,1:nratings[",t,"]] ~ dmulti(prT[s,1:nratings[",t,"],",t,"],CR[s,",t,"])",
      "\n\t\t",count_name,"[s,(nratings[",t,"]+1):(nratings[",t,"]*2)] ~ dmulti(prT[s,(nratings[",t,"]+1):(nratings[",t,"]*2),",t,"],FA[s,",t,"])",
      "\n\t\t",count_name,"[s,(nratings[",t,"]*2+1):(nratings[",t,"]*3)] ~ dmulti(prT[s,(nratings[",t,"]*2+1):(nratings[",t,"]*3),",t,"],M[s,",t,"])",
      "\n\t\t",count_name,"[s,(nratings[",t,"]*3+1):(nratings[",t,"]*4)] ~ dmulti(prT[s,(nratings[",t,"]*3+1):(nratings[",t,"]*4),",t,"],H[s,",t,"])",
      "\n"))
    
    #variance-covariance matrix
    unique_rhos <- t(combn(1:Ntask, 2))
    for(tt in 1:Ntask){ #tt<-3
      if(t==tt){
        vcov_mat <- c(vcov_mat, paste0("\n\tT[",t,",",tt,"] <- sigma_logMratio[",t,"]^2"))
      } else {
        rho <- which((unique_rhos[,1]==t & unique_rhos[,2]==tt) | (unique_rhos[,1]==tt & unique_rhos[,2]==t))
        vcov_mat <- c(vcov_mat, paste0("\n\tT[",t,",",tt,"] <- rho[",rho,"]*sigma_logMratio[",t,"]*sigma_logMratio[",tt,"]"))
      }
    }
  }
  
  #concatenate with rest of model
  jags_mod <- c(
    "data{",
    "\n\tfor (s in 1:nsubj) {",
    data_mod,
    "\n\t}",
    "\n}",
    "\n",
    "\nmodel {",
    "\n\tfor (s in 1:nsubj) {",
    l_multi,
    "\n\t\tfor (task in 1:ntask) {",
    "\n\t\t\tmu[s,task] <- Mratio[s,task]*d1[s,task]",
    "\n\t\t\tS2mu[s,task] <- mu[s,task]/2",
    "\n\t\t\tS1mu[s,task] <- -mu[s,task]/2",
    "\n",
    "\n\t\t\t# Calculate normalisation constants",
    "\n\t\t\tC_area_rS1[s,task] <- phi(c1[s,task] - S1mu[s,task])",
    "\n\t\t\tI_area_rS1[s,task] <- phi(c1[s,task] - S2mu[s,task])",
    "\n\t\t\tC_area_rS2[s,task] <- 1-phi(c1[s,task] - S2mu[s,task])",
    "\n\t\t\tI_area_rS2[s,task] <- 1-phi(c1[s,task] - S1mu[s,task])",
    "\n",
    "\n\t\t\t# Get nC_rS1 probs",
    "\n\t\t\tpr[s,1,task] <- phi(cS1[s,1,task] - S1mu[s,task])/C_area_rS1[s,task]",
    "\n\t\t\tfor (k in 1:(nratings[task]-2)) {",
    "\n\t\t\t\tpr[s,k+1,task] <- (phi(cS1[s,k+1,task] - S1mu[s,task])-phi(cS1[s,k,task] - S1mu[s,task]))/C_area_rS1[s,task]",
    "\n\t\t\t}",
    "\n\t\t\tpr[s,nratings[task],task] <- (phi(c1[s,task] - S1mu[s,task])-phi(cS1[s,nratings[task]-1,task] - S1mu[s,task]))/C_area_rS1[s,task]",
    "\n",
    "\n\t\t\t# Get nI_rS2 probs",
    "\n\t\t\tpr[s,nratings[task]+1,task] <- ((1-phi(c1[s,task] - S1mu[s,task]))-(1-phi(cS2[s,1,task] - S1mu[s,task])))/I_area_rS2[s,task]",
    "\n\t\t\tfor (k in 1:(nratings[task]-2)) {",
    "\n\t\t\t\tpr[s,nratings[task]+1+k,task] <- ((1-phi(cS2[s,k,task] - S1mu[s,task]))-(1-phi(cS2[s,k+1,task] - S1mu[s,task])))/I_area_rS2[s,task]",
    "\n\t\t\t}",
    "\n\t\t\tpr[s,nratings[task]*2,task] <- (1-phi(cS2[s,nratings[task]-1,task] - S1mu[s,task]))/I_area_rS2[s,task]",
    "\n\t\t\t",
    "\n\t\t\t# Get nI_rS1 probs",
    "\n\t\t\tpr[s,(nratings[task]*2)+1, task] <- phi(cS1[s,1,task] - S2mu[s,task])/I_area_rS1[s,task]",
    "\n\t\t\tfor (k in 1:(nratings[task]-2)) {",
    "\n\t\t\t\tpr[s,nratings[task]*2+1+k,task] <- (phi(cS1[s,k+1,task] - S2mu[s,task])-phi(cS1[s,k,task] - S2mu[s,task]))/I_area_rS1[s,task]",
    "\n\t\t\t}",
    "\n\t\t\tpr[s,nratings[task]*3,task] <- (phi(c1[s,task] - S2mu[s,task])-phi(cS1[s,nratings[task]-1,task] - S2mu[s,task]))/I_area_rS1[s,task]",
    "\n\t\t\t",
    "\n\t\t\t# Get nC_rS2 probs",
    "\n\t\t\tpr[s,(nratings[task]*3)+1,task] <- ((1-phi(c1[s,task] - S2mu[s,task]))-(1-phi(cS2[s,1,task] - S2mu[s,task])))/C_area_rS2[s,task]",
    "\n\t\t\tfor (k in 1:(nratings[task]-2)) {",
    "\n\t\t\t\tpr[s,nratings[task]*3+1+k,task] <- ((1-phi(cS2[s,k,task] - S2mu[s,task]))-(1-phi(cS2[s,k+1,task] - S2mu[s,task])))/C_area_rS2[s,task]",
    "\n\t\t\t}",
    "\n\t\t\tpr[s,nratings[task]*4,task] <- (1-phi(cS2[s,nratings[task]-1,task] - S2mu[s,task]))/C_area_rS2[s,task]",
    "\n\t\t\t",
    "\n\t\t\t# Avoid underflow of probabilities",
    "\n\t\t\tfor (i in 1:(nratings[task]*4)) {",
    "\n\t\t\t\tprT[s,i,task] <- ifelse(pr[s,i,task] < Tol, Tol, pr[s,i,task])",
    "\n\t\t\t}",
    "\n",
    "\n\t\t\t# Specify ordered prior on criteria (bounded above and below by Type 1 c)",
    "\n\t\t\tfor (j in 1:(nratings[task]-1)) {",
    "\n\t\t\t\tcS1_raw[s,j,task] ~ dnorm(-mu_c2[task], lambda_c2[task]) T(,c1[s,task])",
    "\n\t\t\t\tcS2_raw[s,j,task] ~ dnorm(mu_c2[task], lambda_c2[task]) T(c1[s,task],)",
    "\n\t\t\t}",
    "\n\t\t\tcS1[s,1:(nratings[task]-1),task] <- sort(cS1_raw[s,1:(nratings[task]-1),task])",
    "\n\t\t\tcS2[s,1:(nratings[task]-1),task] <- sort(cS2_raw[s,1:(nratings[task]-1),task])",
    "\n",
    "\n\t\t\tMratio[s,task] <- exp(logMratio[s,task])",
    "\n\t\t}",
    "\n\t\t# Draw log(M)'s from bi/multivariate Gaussian",
    "\n\t\tlogMratio[s,1:ntask] ~ dmnorm.vcov(mu_logMratio[], T[,]) #dmnorm(mu_logMratio[], TI[,])",
    "\n\t}",
    "\n\t",
    "\n\t#hyperpriors",
    "\n\tfor(t in 1:ntask){",
    "\n\t\tmu_c2[t] ~ dnorm(0, 0.01) #note paper uses N(M,SD), JAGS uses N(M,Tau), where T=1/(SD^2)",
    "\n\t\tsigma_c2[t] ~ dnorm(0, 0.01) I(0, )",
    "\n\t\tlambda_c2[t] <- pow(sigma_c2[t], -2)",
    "\n",
    "\n\t\tmu_logMratio[t] ~ dnorm(0, 1)",
    "\n\t\tlambda_logMratio[t] ~ dgamma(0.001,0.001) ",
    "\n\t\tsigma_logMratio[t] <- 1/sqrt(lambda_logMratio[t])",
    "\n\t}",
    "\n",rhos,
    "\n",vcov_mat,
    "\n}")

  cat(jags_mod, file=paste0(getSrcDirectory(function(x){x}),"/Bayes_metad_group_corr_R.txt"),append=F)
  collapsed_mod <- paste(jags_mod,collapse="")
  return(collapsed_mod)
}

#Simultaenous response-conditional model
createJAGSModel_rc <- function(Ntask){ #dots,breath,
  data_mod<-c(); l_multi<-c(); group_level<-c(); vcov_mat<-c();
    
  #create trial number variable strings
  for(t in 1:Ntask){ #t<-1 ; https://rstudio-pubs-static.s3.amazonaws.com/272658_ae4d482c86514674be17042c852ebbfc.html
    count_name <- paste0("counts",t)
    data_mod <- c(data_mod, paste0(
      "\n\t\tH[s,",t,"] <- sum(",count_name,"[s,(nratings[",t,"]*3+1):(nratings[",t,"]*4)])",
      "\n\t\tM[s,",t,"] <- sum(",count_name,"[s,(nratings[",t,"]*2+1):(nratings[",t,"]*3)])",
      "\n\t\tFA[s,",t,"] <- sum(",count_name,"[s,(nratings[",t,"]+1):(nratings[",t,"]*2)])",
      "\n\t\tCR[s,",t,"] <- sum(",count_name,"[s,1:nratings[",t,"]])",
      "\n"))
    
    l_multi <- c(l_multi,paste0(
      "\n\t\t",count_name,"[s,1:nratings[",t,"]] ~ dmulti(prT[s,1:nratings[",t,"],",t,"],CR[s,",t,"])",
      "\n\t\t",count_name,"[s,(nratings[",t,"]+1):(nratings[",t,"]*2)] ~ dmulti(prT[s,(nratings[",t,"]+1):(nratings[",t,"]*2),",t,"],FA[s,",t,"])",
      "\n\t\t",count_name,"[s,(nratings[",t,"]*2+1):(nratings[",t,"]*3)] ~ dmulti(prT[s,(nratings[",t,"]*2+1):(nratings[",t,"]*3),",t,"],M[s,",t,"])",
      "\n\t\t",count_name,"[s,(nratings[",t,"]*3+1):(nratings[",t,"]*4)] ~ dmulti(prT[s,(nratings[",t,"]*3+1):(nratings[",t,"]*4),",t,"],H[s,",t,"])",
      "\n"))
  }
  
  #variance-covariance matrix
  unique_rhos <- t(combn(1:Ntask, 2))
  for(s in 1:2){
    #number of rhos is different
    Nrho <- factorial(Ntask)/(factorial(Ntask-2)*factorial(2))
    #rhos
    vcov_mat <- c(vcov_mat,paste0(
      "\n\tfor(r in 1:",Nrho,"){",
      "\n\t\trho_rS",s,"[r] ~ dunif(-1,1)",
      "\n\t}\n"))
    
    for(t in 1:Ntask){
      for(tt in 1:Ntask){ #tt<-3
        if(t==tt){
          vcov_mat <- c(vcov_mat, paste0("\n\tT_rS",s,"[",t,",",tt,"] <- sigma_logMratio_rS",s,"[",t,"]^2"))
        } else {
          rho <- which((unique_rhos[,1]==t & unique_rhos[,2]==tt) | (unique_rhos[,1]==tt & unique_rhos[,2]==t))
          vcov_mat <- c(vcov_mat, paste0("\n\tT_rS",s,"[",t,",",tt,"] <- rho_rS",s,"[",rho,"]*sigma_logMratio_rS",s,"[",t,"]*sigma_logMratio_rS",s,"[",tt,"]"))
        }
      }
    }
    vcov_mat <- c(vcov_mat,"\n")
  }
  
  #concatenate with rest of model
  jags_mod <- c(
    "data{",
    "\n\tfor (s in 1:nsubj) {",
    data_mod,
    "\n\t}",
    "\n}",
    "\n",
    "\nmodel {",
    "\n\tfor (s in 1:nsubj) {",
    l_multi,
    "\n\t\tfor (task in 1:ntask) {",
    "\n\t\t\tmu_rS1[s,task] <- Mratio_rS1[s,task]*d1[s,task]",
    "\n\t\t\tS2mu_rS1[s,task] <- mu_rS1[s,task]/2",
    "\n\t\t\tS1mu_rS1[s,task] <- -mu_rS1[s,task]/2",
    "\n\t\t\tmu_rS2[s,task] <- Mratio_rS2[s,task]*d1[s,task]",
    "\n\t\t\tS2mu_rS2[s,task] <- mu_rS2[s,task]/2",
    "\n\t\t\tS1mu_rS2[s,task] <- -mu_rS2[s,task]/2",
    "\n",
    "\n\t\t\t# Calculate normalisation constants",
    "\n\t\t\tC_area_rS1[s,task] <- phi(c1[s,task] - S1mu_rS1[s,task])",
    "\n\t\t\tI_area_rS1[s,task] <- phi(c1[s,task] - S2mu_rS1[s,task])",
    "\n\t\t\tC_area_rS2[s,task] <- 1-phi(c1[s,task] - S2mu_rS2[s,task])",
    "\n\t\t\tI_area_rS2[s,task] <- 1-phi(c1[s,task] - S1mu_rS2[s,task])",
    "\n",
    "\n\t\t\t# Get nC_rS1 probs",
    "\n\t\t\tpr[s,1,task] <- phi(cS1[s,1,task] - S1mu_rS1[s,task])/C_area_rS1[s,task]",
    "\n\t\t\tfor (k in 1:(nratings[task]-2)) {",
    "\n\t\t\t\tpr[s,k+1,task] <- (phi(cS1[s,k+1,task] - S1mu_rS1[s,task])-phi(cS1[s,k,task] - S1mu_rS1[s,task]))/C_area_rS1[s,task]",
    "\n\t\t\t}",
    "\n\t\t\tpr[s,nratings[task],task] <- (phi(c1[s,task] - S1mu_rS1[s,task])-phi(cS1[s,nratings[task]-1,task] - S1mu_rS1[s,task]))/C_area_rS1[s,task]",
    "\n",
    "\n\t\t\t# Get nI_rS2 probs",
    "\n\t\t\tpr[s,nratings[task]+1,task] <- ((1-phi(c1[s,task] - S1mu_rS2[s,task]))-(1-phi(cS2[s,1,task] - S1mu_rS2[s,task])))/I_area_rS2[s,task]",
    "\n\t\t\tfor (k in 1:(nratings[task]-2)) {",
    "\n\t\t\t\tpr[s,nratings[task]+1+k,task] <- ((1-phi(cS2[s,k,task] - S1mu_rS2[s,task]))-(1-phi(cS2[s,k+1,task] - S1mu_rS2[s,task])))/I_area_rS2[s,task]",
    "\n\t\t\t}",
    "\n\t\t\tpr[s,nratings[task]*2,task] <- (1-phi(cS2[s,nratings[task]-1,task] - S1mu_rS2[s,task]))/I_area_rS2[s,task]",
    "\n\t\t\t",
    "\n\t\t\t# Get nI_rS1 probs",
    "\n\t\t\tpr[s,(nratings[task]*2)+1, task] <- phi(cS1[s,1,task] - S2mu_rS1[s,task])/I_area_rS1[s,task]",
    "\n\t\t\tfor (k in 1:(nratings[task]-2)) {",
    "\n\t\t\t\tpr[s,nratings[task]*2+1+k,task] <- (phi(cS1[s,k+1,task] - S2mu_rS1[s,task])-phi(cS1[s,k,task] - S2mu_rS1[s,task]))/I_area_rS1[s,task]",
    "\n\t\t\t}",
    "\n\t\t\tpr[s,nratings[task]*3,task] <- (phi(c1[s,task] - S2mu_rS1[s,task])-phi(cS1[s,nratings[task]-1,task] - S2mu_rS1[s,task]))/I_area_rS1[s,task]",
    "\n\t\t\t",
    "\n\t\t\t# Get nC_rS2 probs",
    "\n\t\t\tpr[s,(nratings[task]*3)+1,task] <- ((1-phi(c1[s,task] - S2mu_rS2[s,task]))-(1-phi(cS2[s,1,task] - S2mu_rS2[s,task])))/C_area_rS2[s,task]",
    "\n\t\t\tfor (k in 1:(nratings[task]-2)) {",
    "\n\t\t\t\tpr[s,nratings[task]*3+1+k,task] <- ((1-phi(cS2[s,k,task] - S2mu_rS2[s,task]))-(1-phi(cS2[s,k+1,task] - S2mu_rS2[s,task])))/C_area_rS2[s,task]",
    "\n\t\t\t}",
    "\n\t\t\tpr[s,nratings[task]*4,task] <- (1-phi(cS2[s,nratings[task]-1,task] - S2mu_rS2[s,task]))/C_area_rS2[s,task]",
    "\n\t\t\t",
    "\n\t\t\t# Avoid underflow of probabilities",
    "\n\t\t\tfor (i in 1:(nratings[task]*4)) {",
    "\n\t\t\t\tprT[s,i,task] <- ifelse(pr[s,i,task] < Tol, Tol, pr[s,i,task])",
    "\n\t\t\t}",
    "\n",
    "\n\t\t\t# Specify ordered prior on criteria (bounded above and below by Type 1 c)",
    "\n\t\t\tfor (j in 1:(nratings[task]-1)) {",
    "\n\t\t\t\tcS1_raw[s,j,task] ~ dnorm(-mu_c2[task], lambda_c2[task]) T(,c1[s,task])",
    "\n\t\t\t\tcS2_raw[s,j,task] ~ dnorm(mu_c2[task], lambda_c2[task]) T(c1[s,task],)",
    "\n\t\t\t}",
    "\n\t\t\tcS1[s,1:(nratings[task]-1),task] <- sort(cS1_raw[s,1:(nratings[task]-1),task])",
    "\n\t\t\tcS2[s,1:(nratings[task]-1),task] <- sort(cS2_raw[s,1:(nratings[task]-1),task])",
    "\n",
    "\n\t\t\tMratio_rS1[s,task] <- exp(logMratio_rS1[s,task])",
    "\n\t\t\tMratio_rS2[s,task] <- exp(logMratio_rS2[s,task])",
    "\n\t\t}",
    "\n\t\t# Draw log(M)'s from bi/multivariate Gaussian",
    "\n\t\tlogMratio_rS1[s,1:ntask] ~ dmnorm.vcov(mu_logMratio_rS1[], T_rS1[,]) #dmnorm(mu_logMratio[], TI[,])",
    "\n\t\tlogMratio_rS2[s,1:ntask] ~ dmnorm.vcov(mu_logMratio_rS2[], T_rS2[,]) #dmnorm(mu_logMratio[], TI[,])",
    "\n\t}",
    "\n\t",
    "\n\t#hyperpriors",
    "\n\tfor(t in 1:ntask){",
    "\n\t\tmu_c2[t] ~ dnorm(0, 0.01) #note paper uses N(M,SD), JAGS uses N(M,Tau), where T=1/(SD^2)",
    "\n\t\tsigma_c2[t] ~ dnorm(0, 0.01) I(0, )",
    "\n\t\tlambda_c2[t] <- pow(sigma_c2[t], -2)",
    "\n",
    "\n\t\tmu_logMratio_rS1[t] ~ dnorm(0, 1)",
    "\n\t\tlambda_logMratio_rS1[t] ~ dgamma(0.001,0.001) ",
    "\n\t\tsigma_logMratio_rS1[t] <- 1/sqrt(lambda_logMratio_rS1[t])",
    "\n\t\tmu_logMratio_rS2[t] ~ dnorm(0, 1)",
    "\n\t\tlambda_logMratio_rS2[t] ~ dgamma(0.001,0.001) ",
    "\n\t\tsigma_logMratio_rS2[t] <- 1/sqrt(lambda_logMratio_rS2[t])",
    "\n\t}",
    "\n",vcov_mat,
    "\n}")
  
  cat(jags_mod, file=paste0(getSrcDirectory(function(x){x}),"/Bayes_metad_rc_group_corr_R.txt"),append=F)
  collapsed_mod <- paste(jags_mod,collapse="")  #this is stored in the final model
  return(collapsed_mod)
}

# Alternating Response-conditional models
createJAGSModel_alt_rc <- function(Ntask,resp_conds=c(0,2,2,2)){ #Ntask<-4
  data_mod<-c(); l_multi<-c(); norm_const<-c(); n_rSP<-c(); underflow<-c(); rhos<-c(); vcov_mat<-c();
  
  #variance-covariance matrix
  Nrho <- factorial(Ntask)/(factorial(Ntask-2)*factorial(2)) #rho for each combination
  rhos <- paste0(
    "\n\tfor(r in 1:",Nrho,"){",
    "\n\t\trho[r] ~ dunif(-1,1)",
    "\n\t}")
  
  #create trial number variable strings
  for(t in 1:Ntask){ #t<-1 ; https://rstudio-pubs-static.s3.amazonaws.com/272658_ae4d482c86514674be17042c852ebbfc.html
    count_name <- paste0("counts",t)
    rc <- resp_conds[t]
      #intro <- c(intro,paste0("#Task ",t," Model rS_",ifelse(rc==0,"1&2",rs)," ----------------"))
    if(rc==0 || rc==1){
      data_mod <- paste0(data_mod,paste0(
        "\n\t\tM[s,",t,"] <- sum(",count_name,"[s,(nratings[",t,"]*2+1):(nratings[",t,"]*3)])",
        "\n\t\tCR[s,",t,"] <- sum(",count_name,"[s,1:nratings[",t,"]])\n"))
      l_multi <- paste0(l_multi,paste0(
        "\n\t\t",count_name,"[s,1:nratings[",t,"]] ~ dmulti(prT[s,1:nratings[",t,"],",t,"],CR[s,",t,"])",
        "\n\t\t",count_name,"[s,(nratings[",t,"]*2+1):(nratings[",t,"]*3)] ~ dmulti(prT[s,(nratings[",t,"]*2+1):(nratings[",t,"]*3),",t,"],M[s,",t,"])"))
      norm_const <- paste0(norm_const,paste0(
        "\n\t\tC_area_rS1[s,",t,"] <- phi(c1[s,",t,"] - S1mu[s,",t,"])",
        "\n\t\tI_area_rS1[s,",t,"] <- phi(c1[s,",t,"] - S2mu[s,",t,"])"))
      n_rSP <- paste0(n_rSP,paste0(
        "\n\t\t#Get nC_rS1 probs task ",t," rS_",ifelse(rc==0,"1&2",rc),
        "\n\t\tpr[s,1,",t,"] <- phi(cS1[s,1,",t,"] - S1mu[s,",t,"])/C_area_rS1[s,",t,"]",
        "\n\t\tfor (k in 1:(nratings[",t,"]-2)) {",
        "\n\t\t\tpr[s,k+1,",t,"] <- (phi(cS1[s,k+1,",t,"] - S1mu[s,",t,"])-phi(cS1[s,k,",t,"] - S1mu[s,",t,"]))/C_area_rS1[s,",t,"]",
        "\n\t\t}",
        "\n\t\tpr[s,nratings[",t,"],",t,"] <- (phi(c1[s,",t,"] - S1mu[s,",t,"])-phi(cS1[s,nratings[",t,"]-1,",t,"] - S1mu[s,",t,"]))/C_area_rS1[s,",t,"]",
        "\n",
        "\n\t\t#Get nI_rS1 probs task ",t," rS_",ifelse(rc==0,"1&2",rc),
        "\n\t\tpr[s,(nratings[",t,"]*2)+1, ",t,"] <- phi(cS1[s,1,",t,"] - S2mu[s,",t,"])/I_area_rS1[s,",t,"]",
        "\n\t\tfor (k in 1:(nratings[",t,"]-2)) {",
        "\n\t\t\tpr[s,nratings[",t,"]*2+1+k,",t,"] <- (phi(cS1[s,k+1,",t,"] - S2mu[s,",t,"])-phi(cS1[s,k,",t,"] - S2mu[s,",t,"]))/I_area_rS1[s,",t,"]",
        "\n\t\t}",
        "\n\t\tpr[s,nratings[",t,"]*3,",t,"] <- (phi(c1[s,",t,"] - S2mu[s,",t,"])-phi(cS1[s,nratings[",t,"]-1,",t,"] - S2mu[s,",t,"]))/I_area_rS1[s,",t,"]",
        "\n"))
      underflow <- paste0(underflow,paste0(
        "\n\t\tfor(i in c(1:nratings[",t,"], (nratings[",t,"]*2+1):(nratings[",t,"]*3))){",
        "\n\t\t\tprT[s,i,",t,"] <- ifelse(pr[s,i,",t,"] < Tol, Tol, pr[s,i,",t,"])",
        "\n\t\t}"))
    } 
    if(rc==0 || rc==2){
      data_mod <- paste0(data_mod, paste0(
        "\n\t\tH[s,",t,"] <- sum(",count_name,"[s,(nratings[",t,"]*3+1):(nratings[",t,"]*4)])",
        "\n\t\tFA[s,",t,"] <- sum(",count_name,"[s,(nratings[",t,"]+1):(nratings[",t,"]*2)])\n"))
      l_multi <- paste0(l_multi,paste0(
        "\n\t\t",count_name,"[s,(nratings[",t,"]+1):(nratings[",t,"]*2)] ~ dmulti(prT[s,(nratings[",t,"]+1):(nratings[",t,"]*2),",t,"],FA[s,",t,"])",
        "\n\t\t",count_name,"[s,(nratings[",t,"]*3+1):(nratings[",t,"]*4)] ~ dmulti(prT[s,(nratings[",t,"]*3+1):(nratings[",t,"]*4),",t,"],H[s,",t,"])"))
      norm_const <- paste0(norm_const,paste0(
        "\n\t\tC_area_rS2[s,",t,"] <- 1-phi(c1[s,",t,"] - S2mu[s,",t,"])",
        "\n\t\tI_area_rS2[s,",t,"] <- 1-phi(c1[s,",t,"] - S1mu[s,",t,"])"))
      n_rSP <- paste0(n_rSP,paste0(
        "\n\t\t#Get nI_rS2 probs task ",t," rS_",ifelse(rc==0,"1&2",rc),
        "\n\t\tpr[s,nratings[",t,"]+1,",t,"] <- ((1-phi(c1[s,",t,"] - S1mu[s,",t,"]))-(1-phi(cS2[s,1,",t,"] - S1mu[s,",t,"])))/I_area_rS2[s,",t,"]",
        "\n\t\tfor (k in 1:(nratings[",t,"]-2)) {",
        "\n\t\t\tpr[s,nratings[",t,"]+1+k,",t,"] <- ((1-phi(cS2[s,k,",t,"] - S1mu[s,",t,"]))-(1-phi(cS2[s,k+1,",t,"] - S1mu[s,",t,"])))/I_area_rS2[s,",t,"]",
        "\n\t\t}",
        "\n\t\tpr[s,nratings[",t,"]*2,",t,"] <- (1-phi(cS2[s,nratings[",t,"]-1,",t,"] - S1mu[s,",t,"]))/I_area_rS2[s,",t,"]",
        "\n",
        "\n\t\t#Get nC_rS2 probs task ",t," rS_",ifelse(rc==0,"1&2",rc),
        "\n\t\tpr[s,(nratings[",t,"]*3)+1,",t,"] <- ((1-phi(c1[s,",t,"] - S2mu[s,",t,"]))-(1-phi(cS2[s,1,",t,"] - S2mu[s,",t,"])))/C_area_rS2[s,",t,"]",
        "\n\t\tfor (k in 1:(nratings[",t,"]-2)) {",
        "\n\t\t\tpr[s,nratings[",t,"]*3+1+k,",t,"] <- ((1-phi(cS2[s,k,",t,"] - S2mu[s,",t,"]))-(1-phi(cS2[s,k+1,",t,"] - S2mu[s,",t,"])))/C_area_rS2[s,",t,"]",
        "\n\t\t}",
        "\n\t\tpr[s,nratings[",t,"]*4,",t,"] <- (1-phi(cS2[s,nratings[",t,"]-1,",t,"] - S2mu[s,",t,"]))/C_area_rS2[s,",t,"]",
        "\n"))
      underflow <- paste0(underflow,paste0(
        "\n\t\tfor(i in c((nratings[",t,"]+1):(nratings[",t,"]*2), (nratings[",t,"]*3+1):(nratings[",t,"]*4))){",
        "\n\t\t\tprT[s,i,",t,"] <- ifelse(pr[s,i,",t,"] < Tol, Tol, pr[s,i,",t,"])",
        "\n\t\t}"))
    }
    

    unique_rhos <- t(combn(1:Ntask, 2))
    for(tt in 1:Ntask){ #tt<-3
      if(t==tt){
        vcov_mat <- c(vcov_mat, paste0("\n\tT[",t,",",tt,"] <- sigma_logMratio[",t,"]^2")) #use sigma_logMratio directly if using dmnorm.vcov()
      } else {
        rho <- which((unique_rhos[,1]==t & unique_rhos[,2]==tt) | (unique_rhos[,1]==tt & unique_rhos[,2]==t))
        vcov_mat <- c(vcov_mat, paste0("\n\tT[",t,",",tt,"] <- rho[",rho,"]*sigma_logMratio[",t,"]*sigma_logMratio[",tt,"]"))
      }
    }
  }
  
  #concatenate with rest of model
  jags_mod <- c(
    "\ndata{",
    "\n\tfor(s in 1:nsubj){",
    data_mod,
    "\n\t}",
    "\n}",
    "\n",
    "\nmodel {",
    "\n\tfor(s in 1:nsubj){",
    "\n\t\t## TYPE 2 SDT MODEL (META-D)",
    "\n\t\t# Multinomial likelihood for response counts ordered as c(nR_S1,nR_S2)",
    l_multi,
    "\n",
    "\n\t\tfor (task in 1:ntask) {",
    "\n\t\t\t# Means of SDT distributions",
    "\n\t\t\tmu[s,task] <- Mratio[s,task]*d1[s,task]",
    "\n\t\t\tS2mu[s,task] <- mu[s,task]/2",
    "\n\t\t\tS1mu[s,task] <- -mu[s,task]/2",
    "\n",   
    "\n\t\t\t# Specify ordered prior on criteria (bounded above and below by Type 1 c)",
    "\n\t\t\tfor (j in 1:(nratings[task]-1)) {",
    "\n\t\t\t\tcS1_raw[s,j,task] ~ dnorm(-mu_c2[task], lambda_c2[task]) T(,c1[s,task])",
    "\n\t\t\t\tcS2_raw[s,j,task] ~ dnorm(mu_c2[task], lambda_c2[task]) T(c1[s,task],)",
    "\n\t\t\t}",
    "\n\t\t\tcS1[s,1:(nratings[task]-1),task] <- sort(cS1_raw[s,1:(nratings[task]-1),task])",
    "\n\t\t\tcS2[s,1:(nratings[task]-1),task] <- sort(cS2_raw[s,1:(nratings[task]-1),task])",
    "\n",
    "\n\t\t\tMratio[s,task] <- exp(logMratio[s,task])",
    "\n\t\t}",
    "\n",
    "\n\t\t# Calculate normalisation constants",
    norm_const,
    "\n",
    n_rSP,
    "\n\t\t# Avoid underflow of probabilities",
    underflow,
    "\n",
    "\n\t\t# Draw log(M)'s from bi/multivariate Gaussian",
    "\n\t\tlogMratio[s,1:ntask] ~ dmnorm.vcov(mu_logMratio[], T[,]) #dmnorm(mu_logMratio[], TI[,])",
    "\n\t}",
    "\n",
    "\n\t#hyperpriors",
    "\n\tfor(t in 1:ntask){",
    "\n\t\tmu_c2[t] ~ dnorm(0, 0.01) #note paper uses N(M,SD), JAGS uses N(M,Tau), where T=1/(SD^2)",
    "\n\t\tsigma_c2[t] ~ dnorm(0, 0.01) I(0, )",
    "\n\t\tlambda_c2[t] <- pow(sigma_c2[t], -2)",
    "\n\t\tmu_logMratio[t] ~ dnorm(0, 1)",
    "\n\t\tlambda_logMratio[t] ~ dgamma(0.001,0.001)",
    "\n\t\tsigma_logMratio[t] <- 1/sqrt(lambda_logMratio[t])",
    "\n\t}",
    "\n",rhos, 
    "\n",vcov_mat,
    "\n}")
  cat(jags_mod, file=paste0(getSrcDirectory(function(x){x}),"/Bayes_metad_alt_rc_group_corr_R.txt"),append=F)
  collapsed_mod <- paste(jags_mod,collapse="")
  return(collapsed_mod)
}

# H-metad --------------
fit_meta_d_mcmc_groupCorr <- function(nR_S1, nR_S2, mcmc_params, fncdf=pnorm, fninv=qnorm){ # fncdf=pnorm; fninv=qnorm; nR_S1<-c_split[c(1,3),6];nR_S2<-c_split[c(2,4),6] ;nR_S1 <- c_nrsxs[1,];nR_S2 <- c_nrsxs[2,]
  gc() #can get quite heavy
  require('rjags')
  require('parallel')
  Ntask <- length(nR_S1)
  ## Sampling
  if(missing('mcmc_params')){
    mcmc_params <- list(
      response_conditional = 0, # response-conditional meta-d? 0=no,1=both simultaneously, 2=alternating
      rc_alternating = c(0,1,2,1,2,2), # if r_c=2, sets response types. 0=non-rc, 1=rS_1, 2=rS_2
      estimate_dprime = 0, # also estimate dprime in same model?
      nchains = min(c(detectCores()-1,max(Ntask,3))), # How Many Chains?
      #suggestion on numbers: https://stackoverflow.com/a/62732851/7705626
      nadapt = Ntask*1000, # How many adaptation samples? rjags default is 1000
      nburnin = Ntask*1000, # How Many Burn-in Samples? 
      nsamples = Ntask*10000, # How Many Recorded Samples? 30000 recommended
      nthin = Ntask, #floor(Ntask/2), # How Often is a Sample Recorded? worth upping this...
      dic = 0, # Not needed as for comparing separte models. # note: computationally expensive. uses 2 chains. >1 provides n.iter for an extra DIC.samples() run e.g. use 10000
      rhat = 1,  # note:computationally expensive, and best done with divergent starting values. consider across 3 chains: observed values, 0, random (or min or max?)
      parallel = 1, # note: running in parallel increases memory useage, if the model doesn't fit you might need to run sequentially
      monitorparams = 0,  # select monitors, 0=all params. =='rho' can speed things up, and subj level estimates take up a lot of memory. do a practice run to see what is monitored by default.
      saveallsamples = 0 # These take up a large amount of computer memory
    )
  }

  Nsubj <- nrow(nR_S1[[1]])
  Nratings <- c()
  counts<-c() #counts <- lapply(1:Ntask, matrix,data=NA,nrow=Nsubj,ncol=Nratings*4)
  d1 <- matrix(NA,Nsubj,Ntask)
  c1 <- matrix(NA,Nsubj,Ntask)
  nTot <- matrix(NA,Nsubj,Ntask)
  datastruct <- list()
  #create datastruct and add in list elements for each counts matrix.
  for(task in 1:Ntask){ #task<-1
    # Get type 1 SDT parameter values
    nrat <- ncol(nR_S1[[task]])/2
    Nratings[task] <- nrat
    counts[[task]] <- as.matrix(cbind(nR_S1[[task]],nR_S2[[task]])) #https://sourceforge.net/p/mcmc-jags/discussion/610037/thread/16f90079/
    datastruct[[paste0('counts',task)]] <- counts[[task]]
    adj_f <- 1/(nrat*2) # Adjust to ensure non-zero counts for type 1 d' point estimate (not necessary if estimating d' inside JAGS)
    nR_S1_adj <- nR_S1[[task]] + adj_f
    nR_S2_adj <- nR_S2[[task]] + adj_f
    nTot[,task] <- apply(counts[[task]],1,sum)
    
    ratingHR<-matrix(NA,Nsubj,(nrat*2)-1)
    ratingFAR<-matrix(NA,Nsubj,(nrat*2)-1)
    for(c in 2:(nrat*2)){
      ratingHR[,c-1] <- apply(nR_S2_adj[c:ncol(nR_S2_adj)],1,sum) / apply(nR_S2_adj,1,sum)
      ratingFAR[,c-1] <- apply(nR_S1_adj[c:ncol(nR_S1_adj)],1,sum) / apply(nR_S1_adj,1,sum)
    }
    d1[,task] <- fninv(ratingHR[,nrat]) - fninv(ratingFAR[,nrat])
    c1[,task] <- -0.5*(fninv(ratingHR[,nrat]) + fninv(ratingFAR[,nrat]))
  } #note to run AM's example here needs data$nratings <- c(4,4); data$ntask <- 2

  # Assign variables to the observed nodes note: "Any numeric objects in data corresponding to node arrays used in file are taken to represent the values of observed nodes in the model" https://www.rdocumentation.org/packages/rjags/versions/4-13/topics/jags.model
  datastruct <- c(datastruct, list('d1'=d1, 'c1'=c1, 'nsubj'=Nsubj, 'ntask'=Ntask, 'nratings'=Nratings, 'Tol'=1e-05))
  
  #Trace monitors: d1 and c1 can't be used here as these values are static and not touched by JAGS. The following vars from group H-metad all appear fine to monitor. See 'Find Trace Parameters' section below
  switch(as.character(mcmc_params$response_conditional),
    '0' = {
      jags_txt <- createJAGSModel(Ntask)
      model_file <- 'Bayes_metad_group_corr_R.txt'
      monitorparams <- c('mu_logMratio', 'sigma_logMratio', 'rho', 'Mratio', 'mu_c2','sigma_c2','cS1','cS2') # Monitors all JAGS modelled parameters
    },
    '1' = {
      jags_txt <- createJAGSModel_rc(Ntask)
      model_file <- 'Bayes_metad_rc_group_corr_R.txt'
      monitorparams <- c('mu_logMratio_rS1','mu_logMratio_rS2','sigma_logMratio_rS1','sigma_logMratio_rS2','rho_rS1','rho_rS2','Mratio_rS1','Mratio_rS2','mu_c2','sigma_c2','cS1','cS2')
      },
    '2' = {
      jags_txt <- createJAGSModel_alt_rc(Ntask,mcmc_params$rc_alternating)
      model_file <- 'Bayes_metad_alt_rc_group_corr_R.txt'
      monitorparams <- c('mu_logMratio', 'sigma_logMratio', 'rho', 'Mratio', 'mu_c2','sigma_c2','cS1','cS2') # Monitors all JAGS modelled parameters
      }
  )
  if(mcmc_params$monitorparams!=0){
    monitorparams <- mcmc_params$monitorparams
  }
  #load DIC monitoring
  if(mcmc_params$dic>0){ #https://sourceforge.net/p/mcmc-jags/discussion/610037/thread/ea46dc43/
    monitorparams <- c(monitorparams,'deviance') #this is for the trace monitor 'deviance' https://cchecastaldo.github.io/BayesianShortCourse/content/lectures/JAGSPrimerMCMCVis.pdf
  }
  mod_path <- file.path(getSrcDirectory(function(x){x}),model_file)# mod_path<-file.path('meta-d/',model_file)
  message(paste('monitored parameters:',paste(monitorparams,collapse=', ')))
  message(paste0('model file: ',mod_path))
  
  
  
  start_time <- Sys.time() #timer
  # Use JAGS to Sample
  if(mcmc_params$parallel==0){
    load.module('dic');load.module('glm');load.module('lecuyer')
    jags_mod <- jags.model(mod_path, datastruct, n.chains=mcmc_params$nchains,n.adapt=mcmc_params$nadapt) #n.adapt=1000 # https://cran.r-project.org/web/packages/rjags/rjags.pdf
    #note if adaptation incomplete: https://stats.stackexchange.com/questions/429483/adaptation-incomplete-in-rjags-how-bad-is-it
    update(jags_mod, n.iter=mcmc_params$nburnin)
    mod_samp <- coda.samples(jags_mod, monitorparams, n.iter=mcmc_params$nsamples, thin=mcmc_params$nthin)
    #organised as: mod_samp[chain_n][sample,stat] with a stat for each task (and each subj OR at group level)
    unload.module('dic');unload.module('glm');unload.module('lecuyer')
    
  } else if(mcmc_params$parallel==1){
    time_est <- round(((mcmc_params$nadapt+mcmc_params$nburnin+mcmc_params$nsamples)*Nsubj*((Ntask/2)-.5))/60000,2)
    message(paste0("No progress bar available in parallel. VERY rough time estimate: ",time_est,' mins'))
    #cl <- makePSOCKcluster(mcmc_params$nchains,methods=F)
    jinits <- function(){list(.RNG.name = 'lecuyer::RngStream',
                              .RNG.seed = round(1e+06*runif(1)))}
    cl <- makeCluster(mcmc_params$nchains,methods=F, type="PSOCK")
    JAGSmod <- function(seed){
      set.seed(seed)
      jags_mod <- jags.model(mod_path,datastruct,inits=jinits,n.adapt=mcmc_params$nadapt) #n.adapt=1000 # https://cran.r-project.org/web/packages/rjags/rjags.pdf
      update(jags_mod, n.iter=mcmc_params$nburnin)
      mod_samp <- coda.samples(jags_mod, monitorparams, n.iter=mcmc_params$nsamples, thin=mcmc_params$nthin)
      return(mod_samp)
    }
    clusterExport(cl,c('mod_path','datastruct','mcmc_params','monitorparams','jinits'),envir=environment()) ##'JAGSmod','jags.model','coda.samples','load.module',
    clusterEvalQ(cl,{library(rjags);load.module('lecuyer');load.module('glm');load.module('dic')})
    mod_samp <- parLapply(cl,1:mcmc_params$nchains,JAGSmod)
    stopCluster(cl)
    mod_samp <- as.mcmc.list(lapply(mod_samp,as.mcmc))
  }
  message('Sampling complete. Now packaging samples...')
  print(Sys.time()-start_time)

  
  
  ## PACKAGE SAMPLES ## this is done first as the samples can require a lot of memory and is best deleted ASAP
  fit <- c()
  # Note: plots will do each participant and task, need to do group-level only.
  #message('Plotting samples...')
  #fit$mcmc$plot <- tryCatch(plot(mod_samp),exception=function(e){message(e)})
  fit$mcmc$samples <- c()
  #parameters to extract (flexible to monitorparams and r.s.)
  t_names <- paste("Task",1:Ntask)
  param_n <- table(gsub("(.*)(\\[.*)","\\1",colnames(mod_samp[[1]]))) 
  group_level <- names(param_n[param_n==Ntask]) #group level estimates e.g. mu & sigma for log(M-ratio) & c2, so long as there's more than Ntask participants this'll be fine
  subj_level <- names(param_n)[grep('^Mratio',names(param_n))] #names(param_n[param_n==Nsubj*Ntask])
  cSx <- names(param_n)[grep('^cS',names(param_n))] #names(param_n[param_n==sum(Nratings-1)*Nsubj])
  
  #function to extract smaller sample list
  statSamples <- function(stat,n){
    stat_samples <- sapply(1:n,function(x){
      c(sapply(mod_samp,function(y){
        y[,grep(paste0(stat,'\\[',x,'\\]'),colnames(y))
        ]}))
    }) # i.e.  mod_samp[[1]][,grep(paste0('rho','\\[',1,'\\]'),colnames(mod_samp[[1]]))]
    
    if(n==Ntask){ colnames(stat_samples)<-t_names}
    return(stat_samples)
    #subject level samples not dealt with by this, but are organised: [sample,"Mratio[subj,task]"]
    #e.g.: mod_samp[[1]][,grep(paste0('^Mratio','\\[[0-9]*,',1,'\\]'),colnames(mod_samp[[1]]))]
  }

  if(length(group_level)>0){
    #group level samples
    group_samples <- lapply(group_level,statSamples,Ntask)
    names(group_samples) <- group_level
    fit$mcmc$samples <- group_samples
  }
  
  #Correlation samples
  rho_idx <- t(combn(1:Ntask, 2)) #names(rhos) <- paste0("Tasks ",rho_idx[,1],"&",rho_idx[,2])
  Nrho <- factorial(Ntask)/(factorial(Ntask-2)*factorial(2))
  if(Nrho>1 && 'rho'%in%monitorparams){ #with 1 estimate the apply functions don't work
    rho_samples <- statSamples('rho',Nrho)
    colnames(rho_samples) <- paste0("Tasks ",rho_idx[,1],"&",rho_idx[,2])
    fit$mcmc$samples$rho <- rho_samples
  }
  HDIs <- lapply(fit$mcmc$samples,function(x){
    t(apply(x,2,function(y){HPDinterval(as.mcmc(y))}))
  })
  HDIs <- lapply(HDIs,function(x){colnames(x)<-c('low','high');x})
  
  if(Nrho==1 && 'rho'%in%monitorparams){ #otherwise apply functions above will throw errors
    rho_samples <- c(sapply(mod_samp,function(x){x[,grep('rho',colnames(x))]}))
    HDI_rho <- HPDinterval(as.mcmc(rho_samples))
    fit$mcmc$samples$rho <- rho_samples
    HDIs <- c(HDIs,list(rho=HDI_rho))
  }
  fit$mcmc$HDI <- HDIs

  if(mcmc_params$saveallsamples==1){ #NOT recommended
    fit$mcmc$samples$all <- mod_samp
  }
  
  ## EXTRACT SUMMARY#
  message('extracting model window....') #y'know I'm not sure this does anything, as burn-in and adaptation aren't saved anyway
  mod_window <- window(mod_samp)  #https://stackoverflow.com/questions/63507610/regarding-a-warning-message-in-jags-for-r 
  rm(mod_samp)
  if(mcmc_params$rhat==1){
    message('running Rhat...')
    require('stableGR')
    rgx_rhat <- paste0("^(",paste(c('d1','c1','Mratio','cS1','cS2'),collapse="|"),")\\[") # seems like subject levels stats have too much info. note: deviance shoud maybe be removed here?
    mod_win_cut <- lapply(mod_window,function(x){x[,-grep(rgx_rhat,colnames(x))]}) 
    fit$mcmc$rhat$full <- tryCatch({ gelman.diag(mod_win_cut) #PRSF and MPRSF This is not the improved Rhat: https://arxiv.org/pdf/1903.08008.pdf
    }, error = function(e) { message(e); print('try gelman.diag() on mcmc$samples directly') }) #PRSF and MPRSF This is not the improved Rhat: https://arxiv.org/pdf/1903.08008.pdf
    fit$mcmc$rhat$stable  <- tryCatch({ stable.GR(mod_win_cut) #PRSF and MPRSF This is not the improved Rhat: https://arxiv.org/pdf/1903.08008.pdf
    }, error = function(e) { message(e); print('try stable.GR() on mcmc$samples directly') }) #PRSF and MPRSF This is not the improved Rhat: https://arxiv.org/pdf/1903.08008.pdf
    rm(mod_win_cut)
    # gelman.diag n_parameters https://stackoverflow.com/questions/57501259/coda-gelman-diag-error-in-chol-defaultw-the-leading-minor-of-order-nn-is
    # PSRF should be <1.1 https://www.stata.com/features/overview/gelman-rubin-convergence-diagnostic ; MPSRF too https://jbds.isdsa.org/index.php/jbds/article/view/45
    # https://www.rdocumentation.org/packages/ggmcmc/versions/1.5.1.1/topics/ggs
  }
  message('Summarising model window...')
  mod_sum <- summary(mod_window)
  message('Packaging summary...')
  rm(mod_window)
  #package summary
  #stats individually, as with MATLAB toolbox
  fit$d1 <- d1
  fit$c1 <- c1
  #outputs flexibile to new trace monitors, r.s. analyses, and task numbers
    # mcmc_list <- mod_sum

  if(length(group_level)>0){
    fit$group_level <- sapply(group_level, function(x){ mod_sum$statistics[grep(x,rownames(mod_sum$statistics)),"Mean"] })
    rownames(fit$group_level) <- t_names
  }
  #subject level estimates e.g. M-ratio on each task
  if(length(subj_level)>0){
    for(s in subj_level){
      x_data <- mod_sum$statistics[grep(paste0('^',s,'\\['),rownames(mod_sum$statistics)),"Mean"] 
      fit[[s]] <- matrix(x_data, ncol=Ntask, dimnames=list(c(),t_names))
    }
  }
  #type-2 criteria, e.g. cS1[subj,rating,task]
  if(length(cSx)>0){
    cSx_names <- paste0(rep(cSx,each=sum(Nratings-1)), unlist(lapply(1:Ntask,function(x){ paste(paste0("c",seq(Nratings[x]-1)),paste0("t",x),sep="_")}))) #x<-1 ; rep(cSx,each=Nratings[x]-1),
    cSx_data <- mod_sum$statistics[grep(paste(cSx,collapse="|"),rownames(mod_sum$statistics)),"Mean"] #cs1[p,c,t]: increases participant, then conf, then task, then does cS2
    fit$cSx <- matrix(cSx_data, ncol=length(cSx_names),dimnames=list(c(),cSx_names))
  }
  
  #correlations (list and matrix format)
  fit$rho <- c()
  rhos <- mod_sum$statistics[grep("rho",rownames(mod_sum$statistics))]
  fit$rho$list <- rhos
  rho_mat <- matrix(1,Ntask,Ntask, dimnames=list(t_names,t_names))
  for(x in 1:length(rhos)){
    rho_mat[rho_idx[x,2],rho_idx[x,1]] <- rhos[x]
    rho_mat[rho_idx[x,1],rho_idx[x,2]] <- rhos[x]
  }
  fit$rho$mat <- rho_mat
  rm(mod_sum)
  #details on mcmc fit
  if(mcmc_params$dic>1){
    print('fitting DIC model...')
    jags_mod <- jags.model(mod_path, datastruct, n.chains=2,n.adapt=mcmc_params$nadapt) #n.adapt=1000 # https://cran.r-project.org/web/packages/rjags/rjags.pdf
    fit$mcmc$dic <- dic.samples(jags_mod, mcmc_params$dic) #popt harsher for complex models but also unstable
    print(Sys.time()-start_time)
  }
  
  # gelman.diag n_parameters https://stackoverflow.com/questions/57501259/coda-gelman-diag-error-in-chol-defaultw-the-leading-minor-of-order-nn-is
    # PSRF should be <1.1 https://www.stata.com/features/overview/gelman-rubin-convergence-diagnostic ; MPSRF too https://jbds.isdsa.org/index.php/jbds/article/view/45
  #https://www.rdocumentation.org/packages/ggmcmc/versions/1.5.1.1/topics/ggs
  fit$mcmc$params <- mcmc_params
  fit$mcmc$model <- jags_txt
  fit$mcmc$time_taken <- print(Sys.time()-start_time)
  gc()
  return(fit)
}

# Find trace monitors--------------
# note an error gets thrown by coda::gelman.diag(): Error in chol.default(W) : the leading minor of order 7 is not positive definite
#'You have derived quantities in you coda output that are functions of parameters you estimate'from https://cchecastaldo.github.io/BayesianShortCourse/content/lectures/JAGSPrimerMCMCVis.pdf
# Code below is from: https://stackoverflow.com/questions/57501259/coda-gelman-diag-error-in-chol-defaultw-the-leading-minor-of-order-nn-is
# Through testing, the only issue is monitoring d1 and c1 if not changed in JAGS after being passed in through R.

my.gelman.diag <- function(x, confidence = 0.95, transform = FALSE, autoburnin = FALSE, multivariate = TRUE){
  x <- as.mcmc.list(x)
  if (nchain(x) < 2)
    stop("You need at least two chains")
  if (autoburnin && start(x) < end(x)/2)
    x <- window(x, start = end(x)/2 + 1)
  Niter <- niter(x)
  Nchain <- nchain(x)
  Nvar <- nvar(x)
  xnames <- varnames(x)
  if (transform)
    x <- gelman.transform(x)
  x <- lapply(x, as.matrix)
  S2 <- array(sapply(x, var, simplify = TRUE),
              dim = c(Nvar, Nvar, Nchain)
  )
  W <- apply(S2, c(1, 2), mean)
  xbar <- matrix(sapply(x, apply, 2, mean, simplify = TRUE),
                 nrow = Nvar, ncol = Nchain)
  B <- Niter * var(t(xbar))
  if (Nvar > 1 && multivariate) {  #ph-edits 
    # CW <- chol(W)
    #    #This is W^-1*B.
    # emax <- eigen(
    #  backsolve(CW, t(backsolve(CW, B, transpose = TRUE)), transpose = TRUE),
    # symmetric = TRUE, only.values = TRUE)$values[1]
    emax <- 1
    mpsrf <- sqrt((1 - 1/Niter) + (1 + 1/Nvar) * emax/Niter)
  }  else {
    mpsrf <- NULL
  }
  
  w <- diag(W)
  b <- diag(B)
  s2 <- matrix(apply(S2, 3, diag), nrow = Nvar, ncol = Nchain)
  muhat <- apply(xbar, 1, mean)
  var.w <- apply(s2, 1, var)/Nchain
  var.b <- (2 * b^2)/(Nchain - 1)
  cov.wb <- (Niter/Nchain) * diag(var(t(s2), t(xbar^2)) - 2 *
                                    muhat * var(t(s2), t(xbar)))
  V <- (Niter - 1) * w/Niter + (1 + 1/Nchain) * b/Niter
  var.V <- ((Niter - 1)^2 * var.w + (1 + 1/Nchain)^2 * var.b +
              2 * (Niter - 1) * (1 + 1/Nchain) * cov.wb)/Niter^2
  df.V <- (2 * V^2)/var.V
  df.adj <- (df.V + 3)/(df.V + 1)
  B.df <- Nchain - 1
  W.df <- (2 * w^2)/var.w
  R2.fixed <- (Niter - 1)/Niter
  R2.random <- (1 + 1/Nchain) * (1/Niter) * (b/w)
  R2.estimate <- R2.fixed + R2.random
  R2.upper <- R2.fixed + qf((1 + confidence)/2, B.df, W.df) *
    R2.random
  psrf <- cbind(sqrt(df.adj * R2.estimate), sqrt(df.adj * R2.upper))
  dimnames(psrf) <- list(xnames, c("Point est.", "Upper C.I."))
  out <- list(psrf = psrf, mpsrf = mpsrf, B = B, W = W) #added ph
  class(out) <- "gelman.diag"
  return( out )
}

##mod_samp = output of coda.samples() function above
#l <- my.gelman.diag(window(mod_samp, start=mcmc_params$nburnin+1001))
#W <- l$W  #Within-sequence variance: d1 and c1 appear to be the issue here.
#rownames(W) <- rownames(l$psrf)
#colnames(W) <- rownames(l$psrf)
#B <- l$B  #Between-sequence variance
##Eigenvalues for d1 and c1 are fine?
#evals.W <- eigen(W, only.values = TRUE)$values
#names(evals.W) <- rownames(l$psrf)
#min.evals <- min(evals.W)
##library(matrixNormal) #install.packages('matrixNormal')
#W.tol <- matrixNormal::is.positive.definite(W, tol = 1e-18) #not positive-definite
#eval <- eigen(solve(W)%*%B, only.values = TRUE)$values[1] #singularity
#dc1 <- grep("^d1|c1",rownames(l$psrf))
#W.tol_short <- matrixNormal::is.positive.definite(W[-dc1,-dc1], tol = 1e-18)
