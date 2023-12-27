fit_meta_d_mcmc_group <- function(nR_S1, nR_S2, mcmc_params, fncdf=pnorm, fninv=qnorm){ 
  # fncdf=pnorm; fninv=qnorm; nR_S1<-c_nrsxs[[1,1]]; nR_S2<-c_nrsxs[[2,1]]
  
  ## SETUP ##
  require('rjags')
  require('parallel')
  
  Nsubj <- nrow(nR_S1)
  if(missing('mcmc_params')){
    mcmc_params <- list(
      response_conditional = 0, # response-conditional meta-d? 0=no,1=both simultaneously, 2=alternating
      estimate_dprime = 0, # also estimate dprime in same model?
      nchains = 3, # How Many Chains?
      #suggestion on numbers: https://stackoverflow.com/a/62732851/7705626
      nadapt = max(c(1000,Nsubj*10)), # How many adaptation samples? rjags default is 1000
      nburnin = max(c(1000,Nsubj*10)), # How Many Burn-in Samples? 
      nsamples = max(c(10000,Nsubj*100)), # How Many Recorded Samples? 30000 recommended
      nthin = 1, #floor(Ntask/2), # How Often is a Sample Recorded? worth upping this...
      dic = 0, # Not needed as for comparing separte models. # note: computationally expensive. uses 2 chains. >1 provides n.iter for an extra DIC.samples() run e.g. use 10000
      rhat = 1,  # note:computationally expensive, and best done with divergent starting values. consider across 3 chains: observed values, 0, random (or min or max?)
      parallel = 1, # note: running in parallel increases memory useage, if the model doesn't fit you might need to run sequentially
      monitorparams = 0,  # select monitors, 0=all params. =='rho' can speed things up, and subj level estimates take up a lot of memory. do a practice run to see what is monitored by default.
      saveallsamples = 0 # These take up a large amount of computer memory
    )
  }
  
  
  ## OBSERVED DATA ###
  Nratings <- ncol(nR_S1)/2
  counts <- as.matrix(cbind(nR_S1,nR_S2)) #https://sourceforge.net/p/mcmc-jags/discussion/610037/thread/16f90079/
  nTot <- apply(counts,1,sum)

  # Get type 1 SDT parameter values
  adj_f <- 1/(Nratings*2) # Adjust to ensure non-zero counts for type 1 d' point estimate (not necessary if estimating d' inside JAGS)
  nR_S1_adj <- nR_S1 + adj_f
  nR_S2_adj <- nR_S2 + adj_f
  ratingHR<-matrix(NA,Nsubj,(Nratings*2)-1)
  ratingFAR<-matrix(NA,Nsubj,(Nratings*2)-1)
  for(c in 2:(Nratings*2)){
    ratingHR[,c-1] <- apply(nR_S2_adj[c:ncol(nR_S2_adj)],1,sum) / apply(nR_S2_adj,1,sum)
    ratingFAR[,c-1] <- apply(nR_S1_adj[c:ncol(nR_S1_adj)],1,sum) / apply(nR_S1_adj,1,sum)
  }
  d1 <- fninv(ratingHR[,Nratings]) - fninv(ratingFAR[,Nratings])
  c1 <- -0.5*(fninv(ratingHR[,Nratings]) + fninv(ratingFAR[,Nratings]))
  datastruct <- list('counts'=counts, 'nratings'=Nratings, 'nsubj'=Nsubj, 'Tol'=1e-05)
  
  ## MODEL & MONITORS ##
  monitorparams <- c('mu_c2','sigma_c2','cS1','cS2')
  if(mcmc_params$estimate_dprime){ 
    monitorparams <- c(monitorparams,'mu_d1','sigma_d1','mu_c','sigma_c')
  } else {
    datastruct <- c(list('d1'=d1,'c1'= c1), datastruct)
  }
  # Select model file and parameters to monitor
  if(mcmc_params$response_conditional){
    monitorparams <- c(monitorparams,'mu_logMratio_rS1','mu_logMratio_rS2','sigma_logMratio_rS1','sigma_logMratio_rS2','Mratio_rS1','Mratio_rS2')
    if(mcmc_params$estimate_dprime){
      model_file <- 'JAGS models/Bayes_metad_rc_group_R.txt'
    } else {
      model_file <- 'JAGS models/Bayes_metad_rc_group_nodp_R.txt'
    }
  } else {
    monitorparams <- c(monitorparams,'mu_logMratio','sigma_logMratio','Mratio')
    if(mcmc_params$estimate_dprime){
      model_file <- 'JAGS models/Bayes_metad_group_R.txt'
    } else {
      model_file <- 'JAGS models/Bayes_metad_group_nodp_R.txt'
    }
  }
  #change monitored params for speed
  #if(mcmc_params$monitorparams!=0){ monitorparams <- mcmc_params$monitorparams }
  #load DIC monitoring #https://sourceforge.net/p/mcmc-jags/discussion/610037/thread/ea46dc43/
  #this is for the trace monitor 'deviance' https://cchecastaldo.github.io/BayesianShortCourse/content/lectures/JAGSPrimerMCMCVis.pdf
  #if(mcmc_params$dic>0){ monitorparams <- c(monitorparams,'deviance') }

  mod_path <- file.path(getSrcDirectory(function(x){x}),model_file)
  
  ## FIT MODEL ##
  start_time <- Sys.time() #timer
  if(mcmc_params$parallel==0){
    load.module('dic');load.module('glm');load.module('lecuyer')
    jags_mod <- jags.model(mod_path, datastruct, n.chains=mcmc_params$nchains,n.adapt=mcmc_params$nadapt) #n.adapt=1000 # https://cran.r-project.org/web/packages/rjags/rjags.pdf
    #note if adaptation incomplete: https://stats.stackexchange.com/questions/429483/adaptation-incomplete-in-rjags-how-bad-is-it
    update(jags_mod, n.iter=mcmc_params$nburnin)
    mod_samp <- coda.samples(jags_mod, monitorparams, n.iter=mcmc_params$nsamples, thin=mcmc_params$nthin)
    #organised as: mod_samp[chain_n][sample,stat] with a stat for each task (and each subj OR at group level)
    unload.module('dic');unload.module('glm');unload.module('lecuyer')
    
  } else if(mcmc_params$parallel==1){
    time_est <- round(((mcmc_params$nadapt+mcmc_params$nburnin+mcmc_params$nsamples)*Nsubj)/60000,2)
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
    on.exit(stopCluster(cl))
    mod_samp <- as.mcmc.list(lapply(mod_samp,as.mcmc))
  }
  print(Sys.time()-start_time)
  message('Sampling complete. Now packaging samples...')
  
  ## PACKAGE MODEL ##
  fit <- c()
  fit$samples  <- do.call(rbind,lapply(mod_samp,function(x){ x[,grep('mu_logMratio',colnames(x)),drop=F] }))
  fit$HDI <- HPDinterval(as.mcmc(exp(fit$samples)))#apply(fit$samples,2,function(x){ HPDinterval(as.mcmc(x)) })
  if(mcmc_params$rhat==1){
    require('stableGR')
    message('calculating Rhat...')
    if(!mcmc_params$estimate_dprime){ rhatparams <- "^(d1|c1|cS)" 
    } else { rhatparams <- "^cS" }    #cS1 and 2 can throw issues here, although unsure if they should be kept in or not?
    mod_cut <- lapply(mod_samp,function(x){x[,!grepl(rhatparams,colnames(x))]})
    fit$rhat$full <- gelman.diag(mod_cut)
    fit$rhat$stable <- stable.GR(mod_cut)
  }
  message('Calculating model summary...')
  print(Sys.time()-start_time)
  mod_sum <- summary(mod_samp)
  rm(mod_samp)
  fit$group_means <- mod_sum$statistics[!grepl('^(Mratio|cS1|cS2)',rownames(mod_sum$statistics)),'Mean']
  fit$Mratio <- mod_sum$statistics[grepl('^Mratio',rownames(mod_sum$statistics)),'Mean']
  fit$observed_d1 <- d1
  fit$observed_c1 <- c1
  fit$time_taken <- print(Sys.time()-start_time)
  return(fit)
}