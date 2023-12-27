# Adaptation in R by Max Lovell of the following MATLAB code by Steve Fleming: 
  #https://github.com/metacoglab/HMeta-d/blob/master/Matlab/fit_meta_d_mcmc.m

#library(rjags)

#example from original code for comparison:
  # nR_S1 <- c(100,50,20,10,5,1)
  # nR_S2 <- c(3,7,8,12,27,89)
  # fit <- fit_meta_d_mcmc(nR_S1,nR_S2)
   
  # fits <- c()
  # for(i in 1:10){
  #   abc <- proc.time()
  #   fit <- fit_meta_d_mcmc(nR_S1,nR_S2)
  #   print(proc.time() - abc)
  #   fits <- c(fits, fit$meta_d)
  # }
  # fits
  
fit_meta_d_mcmc <- function(nR_S1, nR_S2, mcmc_params, fncdf=pnorm, fninv=qnorm){ 
  # fncdf=pnorm; fninv=qnorm; nR_S1<-c_nrsxs[[1,1]][1,];nR_S2<-c_nrsxs[[2,1]][1,]
  gc()
  ##SETUP
  require('rjags')
  require('parallel')
  
  # Optional parameters
  if(missing(mcmc_params)){
    mcmc_params <- list(
      response_conditional = 0, # response-conditional?
      estimate_dprime = 0, # also estimate dprime in same model?
      nchains = 3, # How Many Chains?
      nadapt = 1000,
      nburnin = 1000, # How Many Burn-in Samples?
      nsamples = 10000, # How Many Recorded Samples?
      nthin = 1, # How Often is a Sample Recorded?
      dic = 0, # note: computationally expensive. >1 provides n.iter for an extra DIC.samples() run e.g. use 10000
      rhat = 1, # note: computationally expensive. rhat best done with divergent starting values. consider across 3 chains: observed values, 0, random (or min or max?)
      parallel = 1,
      monitorparams = 0, # select monitors, 0=all params. ==c('d1','meta_d') can speed things up, and might be necessary if using rhat. do a practice run to see what is monitored by default.
      saveallsamples = 0,
      estimate_mratio = 0
    )
  }
  
  #extract counts
  counts <- c(nR_S1, nR_S2)
  nRatings <- length(nR_S1)/2
  adj_f <- 1/length(nR_S1) #change to .5?
  nR_S1_adj <- nR_S1 + adj_f
  nR_S2_adj <- nR_S2 + adj_f
  
  # Observed
  ratingHR<-c(); ratingFAR<-c();
  for(c in 2:(nRatings*2)){
    ratingHR[c-1] <- sum(nR_S2_adj[c:length(nR_S2_adj)]) / sum(nR_S2_adj)
    ratingFAR[c-1] <- sum(nR_S1_adj[c:length(nR_S1_adj)]) / sum(nR_S1_adj)
  }
  d1 <- fninv(ratingHR[nRatings])-fninv(ratingFAR[nRatings])
  c1 <- -.5*(fninv(ratingHR[nRatings])+fninv(ratingFAR[nRatings]))
  
  # Assign variables to the observed nodes
  datastruct <- list('counts'=counts, 'nratings'=nRatings, 'Tol'=1e-05)
  if(!mcmc_params$estimate_dprime){ 
    datastruct <- c(list('d1'=d1,'c1'= c1), datastruct)
  }
  # Select model file and parameters to monitor
  if(mcmc_params$response_conditional){
    model_file <- 'JAGS models/Bayes_metad_rc_R.txt'
    monitorparams <- c('meta_d_rS1','meta_d_rS2')
  } else {
    model_file <- 'JAGS models/Bayes_metad_R.txt'
    monitorparams <- c('meta_d')
  }
  
  if(mcmc_params$estimate_mratio){
    model_file <- 'JAGS models/Bayes_mratio_R'
    monitorparams <- c('m_ratio')
  }
  #model_file <- 'Bayes_mratio_R.txt'
  #monitorparams <- c('meta_d','Mratio','logMratio')
  monitorparams <- c(monitorparams,'d1','c1','cS1','cS2')
  if(all(mcmc_params$monitorparams!=0)){
    monitorparams <- mcmc_params$monitorparams
  }
  if(mcmc_params$dic>0){
    monitorparams <- c(monitorparams,'deviance') #note need jags.samples not coda.samples if monitoring 'pD' aswell
  }
  
  #Run model
  #do.call(file.remove, list(list.files(file.path(pwd(),'tmpjags'), full.names = TRUE)))
  #for reproduciability set seeds before separate chains: https://stackoverflow.com/questions/45872938/how-can-i-set-a-random-seed-using-the-jags-function
  require(rjags)
  #https://sourceforge.net/p/mcmc-jags/discussion/610037/thread/ea46dc43/
  mod_path <- file.path(getSrcDirectory(function(x){x}),model_file) # mod_path<-paste0('meta-d/',model_file)
  
  start_time <- Sys.time() #timer
  if(mcmc_params$parallel==0){
    set.seed(1)
    load.module('glm');load.module('lecuyer');load.module('dic')
    jags_mod <- jags.model(mod_path, datastruct, n.chains=mcmc_params$nchains, quiet=T) #n.adapt=1000 # https://cran.r-project.org/web/packages/rjags/rjags.pdf
    update(jags_mod, n.iter=mcmc_params$nburnin)
    mod_samp <- coda.samples(jags_mod, monitorparams, n.iter=mcmc_params$nsamples, thin=mcmc_params$nthin)
    unload.module('glm');unload.module('lecuyer');unload.module('dic')
    
  } else if(mcmc_params$parallel==1){
    message('running parallel chains, no progress bar available...')
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
    clusterExport(cl,c('mod_path','datastruct','mcmc_params','monitorparams','jinits'),envir=environment()) ##add vars to each cl env 'JAGSmod','jags.model','coda.samples','load.module',
    clusterEvalQ(cl, { #run code in each cl env
      library(rjags)
      load.module('lecuyer')
      load.module('glm')
      load.module('dic')
    })
    mod_samp <- parLapply(cl,1:mcmc_params$nchains,JAGSmod)
    stopCluster(cl)
    mod_samp <- as.mcmc.list(lapply(mod_samp,as.mcmc))
  }

  mod_window <- window(mod_samp)#, start=mcmc_params$nburnin+1001)
  mod_sum <- summary(mod_window) #accounts for default 1000 adaptation iterations #https://stackoverflow.com/questions/63507610/regarding-a-warning-message-in-jags-for-r 
  print(Sys.time()-start_time)
  
  # Package Monitored Stats
  stats <- mod_sum$statistics[!grepl('cS',rownames(mod_sum$statistics)),"Mean"]
  fit <- lapply(split(stats, names(stats)), unname)
  if(any(monitorparams%in%c('cS1','cS2'))){
    fit$t2ca_rS1  <- mod_sum$statistics[grep('cS1',rownames(mod_sum$statistics)),"Mean"] #c("cS1[1]","cS1[2]")
    fit$t2ca_rS2  <- mod_sum$statistics[grep('cS2',rownames(mod_sum$statistics)),"Mean"]
  }
  
  # Calculate M-ratio
  d1_samples <- unlist(lapply(mod_samp,function(x){x[,grep('d1',colnames(x))]})) #only does something if 'estimate d prime' is on, otherwise d1 same as in datastruct
  d1_samples[d1_samples==0] <- 0.0001 #avoid /0
  metad <- monitorparams[grep('meta_d', monitorparams)]
  smu<-c()
  for(p in metad){ #p<-metad[1]
    meta_d_samples <- unlist(lapply(mod_samp,function(x){x[,grep(p,colnames(x))]}))
    fit[[sub('meta_d','M_ratio',p)]] <- mean(meta_d_samples/d1_samples)
    fit[[sub('meta_d','M_diff',p)]]  <- mean(meta_d_samples-d1_samples)
    smu[[sub('meta_d','S1mu',p)]] <- -as.numeric(fit[[p]])/2
    smu[[sub('meta_d','S2mu',p)]] <- as.numeric(fit[[p]])/2
  }

  # R.S. Pseudo-T1 counts
  pt1 <- matrix(NA,8,nRatings-1)
  rownames(pt1) <- paste(rep(c('obs','est'),each=4),c('FAR2','HR2'),rep(c('rS2','rS1'),each=2),sep="_")
  # observed data
  I_nR_rS2 <- nR_S1[(nRatings+1):length(nR_S1)]
  I_nR_rS1 <- nR_S2[nRatings:1]
  C_nR_rS2 <- nR_S2[(nRatings+1):length(nR_S2)]
  C_nR_rS1 <- nR_S1[nRatings:1]
  for(i in 2:nRatings){
    pt1['obs_FAR2_rS2',i-1] <- sum(I_nR_rS2[i:nRatings])/sum(I_nR_rS2)
    pt1['obs_HR2_rS2',i-1]  <- sum(C_nR_rS2[i:nRatings])/sum(C_nR_rS2)
    pt1['obs_FAR2_rS1',i-1] <- sum(I_nR_rS1[i:nRatings])/sum(I_nR_rS1)
    pt1['obs_HR2_rS1',i-1]  <- sum(C_nR_rS1[i:nRatings])/sum(C_nR_rS1)
  }
  
  # model estimated counts
  if(all(c('c1','cS1','cS2')%in%monitorparams)&& any(grepl('meta_d',monitorparams)) ){
    s<-1;S1sd<-1; S2sd<-S1sd/s
    S2mu_rS2 <- smu[names(smu)%in%c("S2mu_rS2","S2mu")][[1]]
    S1mu_rS2 <- smu[names(smu)%in%c("S1mu_rS2","S1mu")][[1]]
    S1mu_rS1 <- smu[names(smu)%in%c("S1mu_rS1","S1mu")][[1]]
    S2mu_rS1 <- smu[names(smu)%in%c("S2mu_rS1","S2mu")][[1]]
    
    C_area_rS2 <- 1-fncdf(fit$c1,S2mu_rS2,S2sd)
    I_area_rS2 <- 1-fncdf(fit$c1,S1mu_rS2,S1sd)
    C_area_rS1 <- fncdf(fit$c1,S1mu_rS1,S1sd)
    I_area_rS1 <- fncdf(fit$c1,S2mu_rS1,S2sd)
    t2c1 <- c(fit$t2ca_rS1,fit$t2ca_rS2)
    for(i in 1:(nRatings-1)){
      t2c1_lower <- t2c1[nRatings-i]
      t2c1_upper <- t2c1[nRatings-1+i]
      
      pt1["est_FAR2_rS2",i] <- 1-fncdf(t2c1_upper,S1mu_rS2,S1sd) / I_area_rS2
      pt1["est_HR2_rS2",i]  <- 1-fncdf(t2c1_upper,S2mu_rS2,S2sd) / C_area_rS2
      pt1["est_FAR2_rS1",i] <- fncdf(t2c1_lower,S2mu_rS1,S2sd) / I_area_rS1
      pt1["est_HR2_rS1",i]  <- fncdf(t2c1_lower,S1mu_rS1,S1sd) / C_area_rS1
    }
  }
  fit$pseudo_t1 <- pt1

  fit$mcmc<-c()
  if(mcmc_params$dic>1){
    message('fitting DIC model...') #https://sourceforge.net/p/mcmc-jags/discussion/610037/thread/ea46dc43/
    start_time <- Sys.time() #timer
    jags_mod <- jags.model(mod_path, datastruct, n.chains=2,n.adapt=mcmc_params$nadapt) #n.adapt=1000 # https://cran.r-project.org/web/packages/rjags/rjags.pdf
    fit$mcmc$dic <- dic.samples(jags_mod, mcmc_params$dic) #popt harsher for complex models but also unstable
    message(Sys.time()-start_time)
  }
  
  if(mcmc_params$rhat==1){
    require('stableGR')
    message('calculating Rhat...')
    #cS1 and 2 can throw issues here, although unsure if they should be kept in or not?
    if(!mcmc_params$estimate_dprime){ rhatparams <- "^(d1|c1|cS)" 
    } else { rhatparams <- "^cS" }
    
    mod_win_cut <- lapply(mod_window,function(x){x[,!grepl(rhatparams,colnames(x)),drop=F]}) #(colnames(x)%in%c('d1','c1')) Note: change this for 'estimate_dprime=1'
    fit$mcmc$rhat$full <- gelman.diag(mod_win_cut) #PRSF and MPRSF This is not the improved Rhat: https://arxiv.org/pdf/1903.08008.pdf
    fit$mcmc$rhat$stable <- stable.GR(mod_win_cut)
    # gelman.diag n_parameters https://stackoverflow.com/questions/57501259/coda-gelman-diag-error-in-chol-defaultw-the-leading-minor-of-order-nn-is
    # PSRF should be <1.1 https://www.stata.com/features/overview/gelman-rubin-convergence-diagnostic ; MPSRF too https://jbds.isdsa.org/index.php/jbds/article/view/45
    #https://www.rdocumentation.org/packages/ggmcmc/versions/1.5.1.1/topics/ggs
  }
  if(mcmc_params$saveallsamples==1){ #NOT recommended
    fit$mcmc$samples <- mod_samp
  }
  fit$mcmc$params <- mcmc_params
  rm(mod_samp)
  rm(mod_window)
  return(fit)
}

