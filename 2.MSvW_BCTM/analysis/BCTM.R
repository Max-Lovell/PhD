# Analysis code for MSvW BCT-M study by Max Lovell #
# Setup ----
#make working dir same as script file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#load environment if code already run
if(file.exists('BCTM.RData')){
  load('BCTM.RData')
}
# Packages ----
list.of.packages <- c("jsonlite", # read current study's raw data
                      "osfr",     # download schmidt's data
                      "archive",  # unpack compressed files from previous studies
                      'rmatio',   # read carpenter's data from matlab
                      'openxlsx', # read file on previous scale measure effects
                      "bfrr",     # Bayes Factors and RR search
                      "viridis"   # colour palette
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

source("fit_meta_d_mcmc.R")

#might also need these:
#install.packages("devtools")
#devtools::install_github("debruine/bfrr")
# Helpers -----
se <- function(x,na.rm=F){ 
  if(na.rm){x<-x[!is.na(x)]}
  sd(x)/sqrt(length(x)) 
}
ci <- function(x){ se(x)*qnorm(.975) }

#Short pivot function
widen <- function(dataset,id_cols,names_from,values_from){
  dataset <- data.frame(dataset)
  #select id cols
  id_comb <- do.call(paste, dataset[id_cols])
  id_idx <- which(!duplicated(id_comb))
  short <- dataset[id_idx,id_cols]
  #create new cols
  pivot_rows <- do.call(paste, dataset[names_from])
  pivot_dat <- sapply(dataset[values_from],tapply,pivot_rows,function(y) y)
  pivot_bind <- do.call(cbind,pivot_dat)
  colnames(pivot_bind) <- gsub(" ","",paste(rep(colnames(pivot_dat),each=length(rownames(pivot_dat))),rownames(pivot_dat),sep="_"))
  #combine
  cbind(short,pivot_bind)
}

pnames <- c("TMS"="TMS-D", "SOW_F"="Observe 1st", "SOW_S"="Observe 2nd",
            "PHQ"="PHQ-8", "GAD"="GAD-7", "RRS"="RRS",
            "RRS_D"="RRS Depression", "RRS_B"="RRS Brooding", "RRS_R"="RRS Rumination",
            "WBSI"="WBSI", "WBSI_sup"="WBSI Supression", "WBSI_int"="WBSI Intrusion",
            "acc"="Accuracy", "dp"="d'", "c"="c", "adj_dp"="Adj. d'",
            "dp2"="Type-2 d'", "c2"="Type-2 c", "adj_dp2"="Adj. type-2 d'",
            "adj_meta_d"="Adj. meta-d'", "H_meta_d"="Hmeta-d", "H_M_ratio"="HM-ratio",
            "ln_M_ratio"="log(M-ratio)","dp2rs"="Type-2 d' (rs)")

round2 <- function(x, digits = 2) {  # Function to always round 0.5 up
  posneg <- sign(x)
  z <- abs(x) * 10^digits
  z <- z + 0.5
  z <- trunc(z)
  z <- z / 10^digits
  z * posneg
}

adj <- function(dataset,outcome,covariates,condition_var='condition_pre'){#dataset<-data_c;outcome<-'dp_diff';covariates<-'c_diff'
  eq <- paste0(outcome,"~",condition_var,"+",paste(covariates,collapse='+'))
  anc <- summary(lm(eq,data=dataset))
  anc$coefficients[2,1:2]
}

# SDT functions -------------------
NRSX <- function(dataset,resps=c('up','down'),br_acc=T){
  #pesudoT1 for meta-d'
  S1_R1 <- table(dataset[dataset$stim==resps[1]&dataset$resp==resps[1],c("ID","conf")])
  S1_R2 <- table(dataset[dataset$stim==resps[1]&dataset$resp==resps[2],c("ID","conf")])
  S2_R1 <- table(dataset[dataset$stim==resps[2]&dataset$resp==resps[1],c("ID","conf")])
  S2_R2 <- table(dataset[dataset$stim==resps[2]&dataset$resp==resps[2],c("ID","conf")])
  #reverse
  S1_R1 <- S1_R1[,ncol(S1_R1):1] 
  S2_R1 <- S2_R1[,ncol(S2_R1):1]
  #create dataset
  nr_sx <- cbind(S1_R1,S1_R2,S2_R1,S2_R2)
  nconf <- ncol(nr_sx)/4
  colnames(nr_sx) <- c(paste(rep(c("S1_R1","S1_R2","S2_R1","S2_R2"),each=nconf),rep(c(nconf:1,1:nconf),2),sep="_"))
  return(nr_sx)
}

type1 <- function(nr_sx_adj,nconf){
  ## TYPE-1 ##
  HR <- apply(nr_sx_adj[,(nconf*3+1):(nconf*4)],1,sum)/apply(nr_sx_adj[,(nconf*2+1):(nconf*4)],1,sum)
  FAR <- apply(nr_sx_adj[,(nconf+1):(nconf*2)],1,sum)/apply(nr_sx_adj[,1:(nconf*2)],1,sum)
  dp <- qnorm(HR) - qnorm(FAR)
  c <- -.5*(qnorm(HR)+qnorm(FAR))
  data.frame(dp,c)
}

type2 <- function(nr_sx_adj,nconf){
  ## TYPE-2 ##
  C <- apply(nr_sx_adj[,c(1:nconf,(nconf*3+1):(nconf*4))],1,sum) #correct
  I <- apply(nr_sx_adj[,c((nconf+1):(nconf*3))],1,sum) #incorrect
  t2 <- data.frame(nr_sx_adj[,F])
  for(x in 2:nconf){
    t2[,paste0("Hc",x)] <- apply(nr_sx_adj[,c(1:(nconf+1-x),(nconf*3+x):(nconf*4))],1,sum)/C
    t2[,paste0("Fc",x)] <- apply(nr_sx_adj[,c((nconf+x):(nconf*2),(nconf*2+1):(nconf*3+1-x))],1,sum)/I
  }
  t2_z <- apply(t2,2,qnorm)
  
  if(nconf>2){
    dp2 <- apply(t2_z,1,function(x){
      t2_lm <- lm(x[grep("H",names(x))]~x[grep("F",names(x))]) #;plot(x[grep("FA",names(x))],x[grep("H",names(x))]); abline(t2_lm)
      return(sqrt(2/(1+(t2_lm$coefficients[2]^2)))*t2_lm$coefficients[1]) #see p. 62-63: Hautus, Macmillan & Creelman (2021) Detection Theory: A User's Guide 3rd ed.
    })
    c2 <- -.5*(t2_z[,paste0("Hc",ceiling(nconf/2))]-t2_z[,paste0("Fc",ceiling(nconf/2))])
  } else {
    dp2 <- t2_z[,"Hc2"] - t2_z[,"Fc2"]
    c2 <- -.5*(t2_z[,"Hc2"] + t2_z[,"Fc2"])
  }
  
  data.frame(dp2,c2)
}

type2RS <- function(nr_sx_adj,nconf){
  ## TYPE-2 RS ##
  C <- apply(nr_sx_adj[,(nconf*3+1):(nconf*4)],1,sum) #correct
  I <- apply(nr_sx_adj[,(nconf+1):(nconf*2)],1,sum) #incorrect
  t2 <- data.frame(nr_sx_adj[,F])
  for(x in 2:nconf){
    t2[,paste0("Hc",x)] <- apply(nr_sx_adj[,(nconf*3+x):(nconf*4),drop=F],1,sum)/C
    t2[,paste0("Fc",x)] <- apply(nr_sx_adj[,(nconf+x):(nconf*2),drop=F],1,sum)/I
  }
  t2_z <- apply(t2,2,qnorm)
  
  if(nconf>2){
    dp2rs <- apply(t2_z,1,function(x){
      t2_lm <- lm(x[grep("H",names(x))]~x[grep("F",names(x))]) #;plot(x[grep("FA",names(x))],x[grep("H",names(x))]); abline(t2_lm)
      return(sqrt(2/(1+(t2_lm$coefficients[2]^2)))*t2_lm$coefficients[1]) #see p. 62-63: Hautus, Macmillan & Creelman (2021) Detection Theory: A User's Guide 3rd ed.
    })
    c2rs <- -.5*(t2_z[,paste0("Hc",ceiling(nconf/2))]-t2_z[,paste0("Fc",ceiling(nconf/2))])
  } else {
    dp2rs <- t2_z[,"Hc2"]-t2_z[,"Fc2"]
    c2rs <- -.5*(t2_z[,"Hc2"]+t2_z[,"Fc2"])
  }
  
  data.frame(dp2rs,c2rs)
}

Hmetad <- function(nr_sx){
  ## META-D' ##
  output <- data.frame(nr_sx[,F])
  mcmc_params <- list(response_conditional=0,estimate_dprime=0,nchains=3,
                      nadapt=1000,nburnin=1000,nsamples=10000,nthin=1,
                      dic=0,rhat=0,parallel=1,monitorparams=0,saveallsamples=0,
                      estimate_mratio=0)
  for(i in 1:nrow(nr_sx)){ # i<-1
    fit <- fit_meta_d_mcmc(nr_sx[i,grep('S1',colnames(nr_sx))],nr_sx[i,grep('S2',colnames(nr_sx))],mcmc_params)
    output[i,c("d1","c1","H_meta_d","H_M_ratio")] <- fit[c("d1","c1","meta_d","M_ratio")]
  }
  output$ln_M_ratio <- log(output$H_M_ratio)
  return(output)
}

HmetadRS <- function(nr_sx){
  ## META-D' ##
  output <- data.frame(nr_sx[,F])
  mcmc_params <- list(response_conditional=1,estimate_dprime=0,nchains=3,
                      nadapt=1000,nburnin=1000,nsamples=10000,nthin=1,
                      dic=0,rhat=0,parallel=1,monitorparams=0,saveallsamples=0,
                      estimate_mratio=0)
  for(i in 1:nrow(nr_sx)){ # i<-1
    fit <- fit_meta_d_mcmc(nr_sx[i,grep('S1',colnames(nr_sx))],nr_sx[i,grep('S2',colnames(nr_sx))],mcmc_params)
    output[i,c("H_meta_d_rS1","H_M_ratio_rS1","H_meta_d_rS2","H_M_ratio_rS2")] <- fit[c('meta_d_rS1','meta_d_rS2','M_ratio_rS1','M_ratio_rS2')]
  }
  output$ln_M_ratio_rS1 <- log(output$H_M_ratio_rS1)
  output$ln_M_ratio_rS2 <- log(output$H_M_ratio_rS2)
  return(output)
}


SDT <- function(dataset,resps=c('up','down'),br_acc=T){ # dataset<-schmidt;resps=c(-1,1)
  nr_sx <- NRSX(dataset,resps)
  nconf <- ncol(nr_sx)/4
  #adjust for 0 counts
  nr_sx_adj <- nr_sx+.5 #from Hmeta-d is adjusted by (1/(n_rat*2)) instead - should standardise this?
  t1 <- type1(nr_sx_adj,nconf)
  t2 <- type2(nr_sx_adj,nconf)
  t2_rs <- type2RS(nr_sx_adj,nconf)
  metad <- Hmetad(nr_sx)
  metad_rs <- HmetadRS(nr_sx)
  #combine
  output <- data.frame(t1,t2,t2_rs,metad,metad_rs) #dataframe for easier column assignment below
  #breath accuracy
  if(br_acc){
    output$n_breath <- table(dataset$ID)
    output$acc <- (tapply(dataset$stim==resps[2]&dataset$resp==resps[2],dataset$ID,sum)/tapply(dataset$resp==resps[2],dataset$ID,sum))
  }
  #Add ID and time
  output$ID <- rownames(nr_sx)
  output$time <- tapply(dataset$time,dataset$ID,function(x){x[1]})
  return(output)
}

#################################################################### Priors ----
# Schmidt et al. (2019) -----------
downloadSchmidt <- function(){
  #Download
  dir.create('priors', showWarnings = FALSE)
  dir.create('priors/schmidt', showWarnings = FALSE)
  dwnld <- osf_download(osf_ls_files(osf_retrieve_node("https://osf.io/trwq3")), path='priors/schmidt', conflicts = "overwrite")
  archive_extract('priors/schmidt/Data MetaMFN.rar','priors/schmidt')
}

compileSchmidt <- function(){ #download schmidt (2019) & combine
  #read
  sch_files <- list.files('priors/schmidt/Data/R/',".txt$",full.names=T)
  sch_list <- lapply(sch_files, read.delim, header=TRUE,sep='')   #note: participant '42  ' is listed with two spaces
  names(sch_list) <- gsub(".*/(.*).txt","\\1", sch_files)
  return(sch_list)
}

#clean Schmidt's data
cleanSchmidt <- function(sch_list){ # sch_list<-sch_raw
  #clean
  sch_list[grep('F10', names(sch_list))] <- lapply(sch_list[grep('F10', names(sch_list))],function(x){x[x$TrialType=='Bloc',]})
  sch_short <- lapply(sch_list,function(x){x[colnames(x) %in% c('suj','sujet','reponse','resp','acc','confidence','conf','rt')]})
  sch_short <- lapply(sch_short,function(x){colnames(x)<-gsub('([a-z]{1,2}).*','\\1',colnames(x));x}) #get first two letters of var names to allow combining
  sch_short <- lapply(1:length(sch_short),function(x){sch_short[[x]]$cond_time <- names(sch_list)[x];sch_short[[x]]}) #get the condition and time variables
  #bind
  sch <- do.call(rbind,sch_short)
  #recode
  sch$re <- ifelse(sch$re%in%c('q',-1),-1,1)
  sch$stim <- ifelse(sch$ac,sch$re,sch$re*-1)
  sch$co <- factor(round2(sch$co,0)) #or floor?
  sch[,c('cond','task','time')] <- do.call(rbind,strsplit(tolower(sch$cond_time),'_'))   #F10 is gabor, F14 is memory
  sch[,'cond'] <- ifelse(sch[,'cond']=='in','exp','control')
  sch[,'task'] <- ifelse(sch[,'task']=='f10','perceptual','memory')
  #remove rt in line with original study: perceptual task: exclude <100ms >2000ms (5.2% of trials), memory task: >100ms & >10000ms (2% of trials) - %'s don't seem to be true?
  sch <- sch[!sch$rt<100 | (sch$task=='perceptual' & sch$rt>2000) | (sch$task=='memory' & sch$rt>10000),]
  
  #rename
  sch <- sch[,c('su','cond','time','task','stim','re','co')] #subset
  names(sch) <- c('ID','cond','time','task','stim','resp','conf')
  #incomplete cases
  id_tab <- table(sch$ID)
  sch <- sch[sch$ID %in% names(id_tab[id_tab==1000]),]
  sch$ID <- factor(sch$ID)
  return(sch)
}

if(!file.exists('priors/schmidt/Data MetaMFN.rar')){
  downloadSchmidt()
}

sch_raw <- compileSchmidt()
schmidt <- cleanSchmidt(sch_raw)

# Carpenter et al. (2019) -----------
downloadCarpenter <- function(){
  dir.create('priors', showWarnings = FALSE)
  dir.create('priors/carpenter', showWarnings = FALSE)
  download.file('https://github.com/metacoglab/CarpenterMetaTraining/archive/refs/heads/master.zip','priors/carpenter/carpenter.zip')
  archive_extract('priors/carpenter/carpenter.zip','priors/carpenter')
}

compileCarpentner <- function(){
  #welcome to the worst structured dataset to ever exist
  message('Loading massive MATLAB file...')
  carp_mat <- read.mat('priors/carpenter/CarpenterMetaTraining-master/Data/results.mat')
  message('done!')
  # hideous one-shot to compile dataset
  carp_list <- lapply(names(carp_mat$data),function(sub_name){ #sub_name<-names(carp_mat$data)[1]
    sub_dat <- carp_mat$data[[sub_name]]
    session <- lapply(names(sub_dat),function(ses_name){ #ses_name <- names(sub_dat)[1]
      ses_dat <- sub_dat[[ses_name]]
      list_vars <- c('trial_type','domain','stimType','trialType','stimCond')
      cbind(do.call(cbind,ses_dat[list_vars]),
            do.call(cbind,lapply(ses_dat[!names(ses_dat)%in%list_vars],function(x){ x[[1]] })),
            session=ses_name, id=sub_name)
    })
    do.call(rbind,session)
  })
  #bind participants
  carp_df <- do.call(rbind,carp_list)
  #stop being list-cols 
  carp_df <- apply(carp_df,2,unlist) #else you get: "Error in order(y) : unimplemented type 'list' in 'orderVector1'"
  
  #cleaning function - call here to not clutter global with above dataset
  cleanCarpenter(carp_df)
}

#clean carpenter's data
cleanCarpenter <- function(carp){ # carp<-carp_df
  carp <- carp[carp[,'key_press']!=-1,] #missing data?
  
  #get condition data
  session_2 <- carp[carp[,'session']=='session_02',]
  conditions <- session_2[!duplicated(session_2[,'id']),c('stimCond','feedbackCond','id')]
  
  #subset dataset
  carp <- carp[carp[,'session'] %in% c('session_01','session_10'),] # select pre and post-test
  carp <- carp[!carp[,'trialType']%in%c("confidence_rating_feedback","memorize_regular",NaN) & carp[,'trial_type']!="fullscreen",] #unwanted rows
  carp <- carp[,c("id",'session',"domain","stimType","trialType","trialNum","key_press","correct")] #subset of cols
  
  #make sure each response is followed by a confidence rating
  row_cat <- apply(carp[,c('id','session','domain','stimType','trialNum')],1,paste,collapse="")
  cat_tab <- table(row_cat)
  carp <- carp[row_cat %in% names(cat_tab[cat_tab==2]),]
  
  #put confidence ratings on same row as response
  carp[carp[,'trialType']=="trial","trialType"] <- "response" #trialType: perception 'response'== memory 'trials'
  confidence <- carp[carp[,'trialType']=="confidence_rating","key_press"]
  carp <- cbind(carp[carp[,'trialType']=="response",],conf=confidence)
  
  #recode key presses
  #carp <- carp[carp[,'key_press']!=-1 & carp[,'conf']!=-1,] #Missing data??
  carp[,'conf'] <- as.numeric(carp[,'conf'])-48 #JavaScript Key codes
  carp[,'key_press'] <- ifelse(carp[,'key_press']==79,-1,1) # technically 79=o,80=p
  
  #time var recode
  carp[,'session'] <- ifelse(carp[,'session']=='session_01',"pre","post")
  
  #add stim and condition vars
  reps <- table(unlist(carp[,'id']))
  carp <- cbind(carp, 
                stim = ifelse(as.numeric(carp[,'correct']),carp[,'key_press'],as.numeric(carp[,'key_press'])*-1), # stim displayed
                stim_cond = rep(unlist(conditions[,'stimCond']),reps),
                feedback_cond = rep(unlist(conditions[,'feedbackCond']),reps))
  
  #recode condition
  carp[,'feedback_cond'] <- ifelse(carp[,'feedback_cond']==2,"exp","control")
  
  #subset/rename vars
  carp <- carp[,c("id",'stim_cond','feedback_cond',"session","domain","stimType","stim","key_press","conf")]
  colnames(carp) <- c("ID",'stim_cond','feedback_cond',"time","domain","stim_type","stim","resp","conf")
  
  #add factors
  carp <- data.frame(carp)
  carp$ID <- factor(carp$ID)
  carp$conf <- factor(carp$conf)
  return(carp)
}

if(!file.exists('priors/carpenter/carpenter.zip')){
  downloadCarpenter()
}
carp <- compileCarpentner()

# Type-2 SDT priors --------------
priorSDT <- function(dataset,condition_var){ #dataset<-schmidt
  #schmidt <- schmidt[schmidt$task=='perceptual',] #select relevant task
  #carp <- carp[carp$domain=='perception',] #select relevant tasks
  pre <- SDT(dataset[dataset$time=='pre',],c(-1,1),F)
  post <- SDT(dataset[dataset$time=='post',],c(-1,1),F)
  id_time <- which(colnames(post)%in%c('ID','time'))
  diff <- mapply(function(x,y){x-y},post[,-id_time],pre[,-id_time])
  cbind(diff,cond=dataset[!duplicated(dataset$ID),condition_var])
}

priorEst <- function(prev_diff){ # prev_diff<-sch_carp
  prev_diff <- data.frame(prev_diff)
  prev_diff[,-ncol(prev_diff)] <- apply(prev_diff[,-ncol(prev_diff)],2,as.numeric)
  prev_diff <- prev_diff[!is.nan(prev_diff$ln_M_ratio),] #remove one row with a negative H_M_ratio
  #analyse
  prev_m <- apply(prev_diff[-ncol(prev_diff)],2,tapply,prev_diff$cond,mean)
  prev_se <- apply(prev_diff[-ncol(prev_diff)],2,tapply,prev_diff$cond,se)
  #create H1 table
  cond_int <- cbind(m=prev_m['exp',]-prev_m['control',],
                    se=apply(prev_se**2,2,sum)**.5)
  H1_sdt <- rbind(cond_int, 
        adj_dp = adj(prev_diff,'dp','c','cond'), #adjusted type-1 d'
        adj_dp2 = adj(prev_diff,'dp2',c('c2','dp','c'),'cond'), #adjusted type-2 d'
        adj_meta_d = adj(prev_diff,'H_meta_d','dp','cond')) #adjusted Hmeta-d
  round2(H1_sdt)
}

#sch_SDT <- priorSDT(schmidt,'cond')
#carp_sub <- carp[carp$stim_type==carp$stim_cond,] #select trained stimuli
#carp_SDT <- priorSDT(carp_sub,'feedback_cond')
#sch_carp <- rbind(sch_SDT,carp_SDT)
#saveRDS(sch_carp,'sch_carp.rds')
sch_carp <- readRDS('sch_carp.rds')
H1_sdt <- priorEst(sch_carp)
# Estimate SDT sample sizes -------------
sampleSize <- function(n_guess,div=1){
  prior_n <- round2(1/mean(1/table(sch_carp[,'cond'])),0) #harmonic mean
  H1s <- cbind(H1_sdt,est_n=NA,n_BF=NA)
  for(v in 1:length(n_guess)){ #v<-1
    v <- n_guess[v]
    n_group <- v/3
    m_est <- abs(H1s[names(v),"m"])/div #/2 if taking as a plausible maximum
    SE_est  <- H1s[names(v),"se"]*sqrt(prior_n/n_group)
    BF <- bfrr(sample_mean=0,sample_se=SE_est,
               sample_df=n_group-1,model="normal", mean=0, sd=m_est, tail=1,
               criterion=3,rr_interval=list(mean=c(0,0),sd=c(0,5)),precision=.01)
    H1s[names(v),c('est_n','n_BF')] <- c(v,BF$BF)
  }
  H1s <- H1s[names(n_guess),]
  rownames(H1s) <- pnames[rownames(H1s)]
  round2(H1s)
}

n_guess <- c("H_meta_d"=70,"H_M_ratio"=100,"ln_M_ratio"=80,"adj_dp2"=110,"adj_meta_d"=60)
n_est <- sampleSize(n_guess)

n_guess_half <- c("H_meta_d"=280,"H_M_ratio"=390,"ln_M_ratio"=320,"adj_dp2"=420,"adj_meta_d"=240)
n_est_half <- sampleSize(n_guess_half,2)

n_est[,'n_BF'] <- n_est_half[,'est_n']
colnames(n_est) <- c('mean','SE','','est_n_half')

#USE AS PLAUSIBLE MAXIMUM AS 240 IS DOABLE.
write.csv(n_est,'csv_files/estimated_sample_sizes.csv')

# Rowland et al. (2019) ----
downloadRowland <- function(){
  dir.create('priors', showWarnings = FALSE)
  dir.create('priors/rowland', showWarnings = FALSE)
  dwnld <- osf_download(osf_ls_files(osf_retrieve_node("https://osf.io/3gsh6")), path='priors/rowland', conflicts = "overwrite")
}

cleanRowland <- function(){
  rowland <- read.csv("priors/rowland/lab_osf.csv")
  rowland <- rowland[rowland$time%in%c(0,6),c('id','time','group','bct_acc')]
  rowland <- rowland[!rowland$id%in%rowland$id[is.na(rowland$bct_acc)],]
  rw <- widen(rowland,c('id','group'),'time','bct_acc')
  rw$acc_diff <- rw$bct_acc_6-rw$bct_acc_0
  rw$group <- ifelse(rw$group,'exp','control')
  return(rw)
}

analyseRowland <- function(rowland){
  ms <- tapply(rowland$acc_diff,rowland$group,mean)
  ses <- tapply(rowland$acc_diff,rowland$group,se)
  c(m=ms['exp']-ms['control'], se=sum(ses**2)**.5)
}

if(!file.exists('priors/rowland/lab_osf.csv')){
  downloadRowland()
}
rowland <- cleanRowland()

row_int <- analyseRowland(rowland)

# BCT Accuracy prior -----
prevBCT <- function(){
  prevs <- read.csv("priors/prior_BCT.csv")
  prevs <- prevs[prevs$Use=="experimental",!names(prevs)%in%c('Use','Citation','Link','Data','Notes')]
  # polsinelli decimals to %
  prevs[prevs$Study=="polsinelli_2017",grep('pre|post',names(prevs))] <- prevs[prevs$Study=="polsinelli_2017",grep('pre|post',names(prevs))]*100
  # convert cell SDs to SES
  prevs[prevs$error_type=='SD',grep('exp_err_',names(prevs))] <- prevs[prevs$error_type=='SD',grep('exp_err_',names(prevs))]/prevs[prevs$error_type=='SD',grep('exp_n',names(prevs))]
  prevs[prevs$error_type=='SD',grep('con_err_',names(prevs))] <- prevs[prevs$error_type=='SD',grep('con_err_',names(prevs))]/prevs[prevs$error_type=='SD',grep('con_n',names(prevs))]
  prevs[prevs$error_type=='SD','error_type'] <- 'SE'
  # diffs
  prevs[is.na(prevs$exp_m_diff),c('exp_m_diff','con_m_diff')] <- prevs[is.na(prevs$exp_m_diff),c('exp_m_post','con_m_post')]-prevs[is.na(prevs$exp_m_diff),c('exp_m_pre','con_m_pre')]
  prevs[is.na(prevs$int_m_diff),'int_m_diff'] <- prevs[is.na(prevs$int_m_diff),'exp_m_diff']-prevs[is.na(prevs$int_m_diff),'con_m_diff']
  #interaction SEs from F or T values
  prevs[is.na(prevs$int_se_diff),"int_se_diff"] <- (sqrt(prevs[,"int_ForT"])/prevs[,"int_m_diff"])[is.na(prevs$int_se_diff)]
  # combine errors
  prevs[is.na(prevs$exp_err_diff),"exp_err_diff"] <- (prevs$exp_err_pre^2 + prevs$exp_err_post^2)[is.na(prevs$exp_err_diff)]**.5
  prevs[is.na(prevs$con_err_diff),"con_err_diff"] <- (prevs$con_err_pre^2 + prevs$con_err_post^2)[is.na(prevs$con_err_diff)]**.5
  prevs[is.na(prevs$int_se_diff),"int_se_diff"] <- (prevs$exp_err_diff^2 + prevs$con_err_diff^2)[is.na(prevs$int_se_diff)]**.5
  
  return(prevs)
}

prev_bct <- prevBCT()
prev_bct[prev_bct$Study=='rowland_2019',c('int_m_diff','int_se_diff')] <- row_int
#add to H1s
#harm_acc <- round2(1/mean(1/c(prev_bct$exp_n,prev_bct$con_n),na.rm=T))
bct_int <- c(m=mean(prev_bct$int_m_diff), se=sqrt(sum(prev_bct$int_se_diff^2)))



# Clapper et al. (2020) -----
downloadClapper <- function(){
  dir.create('priors', showWarnings = FALSE)
  dir.create('priors/clapper', showWarnings = FALSE)
  clap_url <- list(
    e1_br='https://scholarworks.lib.csusb.edu/cgi/viewcontent.cgi?filename=3&article=1016&context=psychology-publications&type=additional',
    e1_qs='https://scholarworks.lib.csusb.edu/cgi/viewcontent.cgi?filename=2&article=1016&context=psychology-publications&type=additional',
    e2_br='https://scholarworks.lib.csusb.edu/cgi/viewcontent.cgi?filename=1&article=1016&context=psychology-publications&type=additional',
    e2_qs='https://scholarworks.lib.csusb.edu/cgi/viewcontent.cgi?filename=0&article=1016&context=psychology-publications&type=additional')
  cl <- lapply(1:length(clap_url),function(x){ download.file(clap_url[[x]],paste0('priors/clapper/clapper_',names(clap_url[x]),'.xlsx'),mode="wb") })
}

compileClapper <- function(){ #download schmidt (2019) & combine
  clap_files <- list.files('priors/clapper/',".xlsx$",full.names=T)
  clap_lst <- lapply(clap_files,read.xlsx)
  #add/rename vars
  clap_lst <- lapply(clap_lst,function(x){colnames(x)<-gsub('_','.',tolower(colnames(x)));x})
  e_names <- gsub('.*_(e[1-2])_.*','\\1',clap_files)
  clap_lst <- lapply(1:length(clap_lst),function(x){ clap_lst[[x]]$experiment <- e_names[x];clap_lst[[x]]})
  clap_br <- rbind(clap_lst[[1]],clap_lst[[3]])
  clap_qs <- merge(clap_lst[[2]],clap_lst[[4]],by=colnames(clap_lst[[4]])[colnames(clap_lst[[4]])%in%colnames(clap_lst[[2]])],all=T)
  list(br=clap_br,qs=clap_qs)
}

if(!file.exists('priors/clapper/clapper_e1_br.xlsx')){
  downloadClapper()
}
clap_raw <- compileClapper()

breathClapper <- function(cbr){ #cbr<-clap_raw$br
  #clean
  #n <- table(cbr[!duplicated(cbr$participt),'experiment']) #note experiment 1 missing a lot of people?
  counts <- cbind(hit = cbr$count==9&cbr$resp.n,
                  miss = cbr$count==9&cbr$resp.b,
                  FA = cbr$count!=9&cbr$resp.n,
                  CR = cbr$count!=9&cbr$resp.b)
  counts <- counts+.5
  sums <- apply(counts,2,tapply,cbr$participant,sum)
  rates <- cbind(HR = sums[,'hit']/(sums[,'hit']+sums[,'miss']),
                 FAR = sums[,'FA']/(sums[,'FA']+sums[,'CR']))
  qrates <- qnorm(rates)
  perf <- cbind(acc=sums[,'hit']/(sums[,'hit']+sums[,'FA']),
                dp=qrates[,'HR']-qrates[,'FAR'],
                c=apply(qrates,1,sum)*-.5)
}

clap_br <- breathClapper(clap_raw$br)

adjDp <- function(cbr){ # cbr<-clap_br
  cbr <- data.frame(cbr)
  cbr$acc <- cbr$acc*100
  slope <- round2(summary(lm(dp~c+acc, data=cbr))$coefficients)[3,1]
  #no_c <- round2(summary(lm(dp~acc, data=cbr))$coefficients)
  round2(bct_int['m']*slope/2)
}

clap_int <- adjDp(clap_br)

# Scale Measures -----
#get previous scale measures
prevScales <- function(){
  prev <- openxlsx::read.xlsx("priors/prior_scales.xlsx", sheet="scales_data")
  prev <- prev[,!(names(prev) %in% c("alt_control","sd"))]
  downfill <- apply(prev,2,function(x){ rep(x[!is.na(x)], diff(c(which(!is.na(x)),length(x)+1))) })
  wide <- widen(downfill,id_cols=c('study','measure','construct','scale_range','num_q'),names_from= c('condition','time'),values_from='m')
  #Remove some scales: SPANE already subtracted one subscale from another. #BITe and BRS already averaged
  wide <- wide[!wide$measure %in% c("Scale of Positive and Negative Experience - affect balance","Brief Irritability Test","Brief Resilience Scale"),]
  # numeric conversion
  wide[,grep('Pre|Post|num_q',colnames(wide))] <- lapply(wide[,grep('Pre|Post|num_q',colnames(wide))],as.numeric)
  #get scale length from range
  wide$scale_length <- sapply(strsplit(wide$scale_range, '-'),function(x){x<-as.numeric(x);x[2]-x[1]})
  #hypothesis direction
  wide$hypothesis_direction <- ifelse(wide$construct%in%c('Well-being','Mindfulness'),1,-1)
  return(wide)
}

prev_scales <- prevScales()

scaleH1 <- function(prev){ # prev<-prev_scales
  #average on 4-point scale
  prev[,grep('Pre|Post',colnames(prev))] <- sapply(prev[,grep('Pre|Post',colnames(prev))],function(x){x/prev$num_q/prev$scale_length*4})
  #mean interaction
  prev$interaction <- (prev$m_ExperimentalPost-prev$m_ExperimentalPre)-(prev$m_ControlPost-prev$m_ControlPre)
  prev$interaction_h <- prev$interaction*prev$hypothesis_direction
  round2(mean(prev$interaction_h))
}

scales_H1 <- scaleH1(prev_scales)

# H1s table -----------
H1s <- rbind(H1_sdt,#cbind(H1_sdt[,'m']/2,H1_sdt[,'se']), 
             scales=c(scales_H1,NA),
             acc=c(bct_int[1]/2,bct_int[2]))

H1s['adj_dp',] <- c(clap_int,NA)
H1s['acc',] <- H1s['acc',]/100
H1s <- round2(H1s)
#saveRDS(H1s,'rds_files/final_priors.rds')
#write.csv(H1s,'csv_files/final_priors.csv')
H1s <- readRDS('rds_files/final_priors.rds')

###################################################################### Data ----
# Task Data ------------
readBreath <- function(){
  task_files <- list.files("task_data/")
  task <- c()
  for(file_n in task_files){ #file_n<-vjs_files[10]
    file <- fromJSON(paste0("task_data/",file_n))
    if(length(file)==0){print(f);next} #print names of incorrect files
    file$ID <- sub("\\D+","",file_n)
    file$time <- sub(".*(pre|post).*","\\1",file_n)
    task <- rbind(task,file)
  }
  task$ID <- factor(task$ID, levels=unique(task$ID)) #ensures they appear in the table() output
  
  return(task)
}

br_raw <- readBreath()

cleanBreath <- function(raw_data){ # raw_data <- br_raw
  #change responses
  raw_data[,c("br","conf")] <- apply(raw_data[,c("br","conf")],2,function(x){ tolower(gsub("Arrow","",x)) })
  raw_data$conf <- gsub("left",1,raw_data$conf)
  raw_data$conf <- gsub("right",2,raw_data$conf)
  
  ## reset deletes everything up to last down press
  downs <- c(0,which(grepl("down",raw_data$br)|!duplicated(raw_data$ID,fromLast=T)))
  resets <- which(grepl("reset",raw_data$br))
  reset_seq <- unlist(lapply(resets,function(x){
    dlx <- downs[downs<x]
    (dlx[length(dlx)]+1):x
  })) #check: raw_data$reset <- F; raw_data[reset_seq,"reset"] <- T
  raw_data <- raw_data[-c(reset_seq,which(raw_data$conf=='no_conf')),] #delete reset data
  
  # create count variable
  seq_ends <- diff(c(0,which(grepl("down",raw_data$br)|!duplicated(raw_data$ID,fromLast=T))))
  raw_data$count <- unlist(lapply(seq_ends,seq))
  raw_data$sequence <- unlist(tapply(grepl("down",raw_data$br),raw_data$ID,function(x){ #index of down presses split by ID
    x[length(x)]<-T #end of sequence always 'True' even if 'up' press (i.e. if timer runs out, etc.)
    rep(1:sum(x),diff(c(0,which(x)))) # number of presses between TRUEs (starting at 0), with sequence number repeated the relevant number of times
  }))
  
  #correct and incorrect responses
  raw_data$target <- ifelse(raw_data$count==9,'down','up')
  raw_data$correct <- raw_data$br==raw_data$target
  
  #format names
  names(raw_data) <- c('resp','resp_rt','conf','conf_rt','ID','time','count','sequence','stim','correct')
  raw_data <- raw_data[,c('ID','count','sequence','stim','resp','resp_rt','correct','conf','conf_rt','time')]
  return(raw_data)
}

br_cl <- cleanBreath(br_raw)
write.csv(br_cl,"csv_files/raw_BCT_data.csv")
#tester <- br_cl[br_cl$ID%in%unique(br_cl$ID)[1:3],];tester$ID<-droplevels(tester$ID);
#br_SDT <- SDT(br_cl)
#write.csv(br_SDT,'csv_files/BCT_performance.csv')
#saveRDS(br_SDT,'rds_files/BCT_performance.rds')
br_SDT <- readRDS('rds_files/BCT_performance.rds')
SDT_vars <- colnames(br_SDT)[!colnames(br_SDT)%in%c("ID","time",'n_breath')]

# Survey Data -----------------
#Note: actual raw survey data available on request - no time to anonymise that atm!

getSurvey <- function(){
  survey_files <- list.files("./survey_data")
  all_rgxs <- paste("^breath",c("pre",paste0("[w|m]",1:10),"post"),"",sep="_")
  #merge files
  for(file_rgx in all_rgxs){ #file_rgx<-all_rgxs[2] #loop through the files for pre, post, and each day separately # file_rgx<-all_rgxs[1]
    day_files <- survey_files[grep(file_rgx,survey_files)] #read in first dataset to add the others to
    day_data <- read.csv(paste0("survey_data/",day_files[1]))
    day_data <- day_data[3:nrow(day_data),]
    for(file_name in 2:length(day_files)){ #loop through the rest and add the rows # file_name<-2
      file_name <- day_files[file_name]
      day_file <- read.csv(paste0("survey_data/",file_name))
      day_file <- day_file[3:nrow(day_file),!names(day_file)%in%c("ismobile","px_cm")]
      day_data <- rbind(day_data,day_file) #add file to rest of data for that day
    }
    suffix <- gsub('\\^|(breath)|_|(\\[w\\|m\\])','',file_rgx)
    names(day_data)[!grepl("Email",names(day_data))] <- paste(names(day_data)[!grepl("Email",names(day_data))],suffix,sep="_")
    #bit of data cleaning, fix variable names, etc
    if(file_rgx=="^breath_pre_"){
      #remove unfinished responses
      day_data <- day_data[day_data$Progress_pre==100 & day_data$Duration..in.seconds._pre>0 & day_data$Finished_pre==1 & day_data$DistributionChannel_pre!="preview",]
      br_pre <- br_SDT[br_SDT$time=='pre',]
      names(br_pre) <- paste(names(br_pre),suffix,sep="_")
      names(day_data)[names(day_data)=='email_pre'] <- 'Email'
      day_data <- day_data[!day_data$Email%in%day_data$Email[duplicated(day_data$Email)],]
      survey <- merge(day_data,br_pre, by.x="Random_ID_pre", by.y="ID_pre", all.x=F, all.y=F) #add task data #note technical error means some task data incomplete - remove
    } else {
      if(file_rgx=="^breath_post_"){#fixing some badly named vars in qualtrics:
        # remove unfinished responses
        day_data <- day_data[day_data$Progress_post==100 & day_data$Duration..in.seconds._post>0 & day_data$Finished_post==1 & day_data$DistributionChannel_post!="preview",]
        #rename
        colnames(day_data) <- gsub('^Q1(_[1-7]_post)$','TMS\\1', colnames(day_data))
        colnames(day_data) <- gsub('^Q1(_[1-7]).1(_post)$','GAD\\1\\2', colnames(day_data))
        colnames(day_data) <- gsub('PHQ.8','PHQ', colnames(day_data))
        colnames(day_data) <- gsub('Q16','RRS', colnames(day_data))
        colnames(day_data) <- gsub('Q17','WBSI', colnames(day_data))
        colnames(day_data) <- gsub('Q16','RRS', colnames(day_data))
        colnames(day_data) <- gsub('Q32','time_meditating', colnames(day_data))
        colnames(day_data) <- gsub('Q33','course_content', colnames(day_data))
        colnames(day_data) <- gsub('Q34','useful', colnames(day_data))
        colnames(day_data) <- gsub('Q22','authentic', colnames(day_data))
        colnames(day_data) <- gsub('Q35','feedback', colnames(day_data))

        #add breath data
        br_post <- br_SDT[br_SDT$time=='post',]
        names(br_post) <- paste(names(br_post),suffix,sep="_")
        day_data <- merge(day_data,br_post, by.x="Random_ID_post", by.y="ID_post", all.x=F, all.y=F) #add breath task data
      } else if(grepl('([1-9]|10)',file_rgx)){ #rename focus and yesterdays activities questions
        colnames(day_data) <- sub('Q2','activities',colnames(day_data))
        colnames(day_data) <- sub('Q1_1','focus',colnames(day_data))
        colnames(day_data) <- sub('focus_1_1','focus_1',colnames(day_data))
        #clean before merge (stops some duplications)
        day_data <- day_data[day_data[,grep('Progress',colnames(day_data))]==100,]
        day_data <- day_data[day_data[,grep('Finished',colnames(day_data))]==1,]
        day_data <- day_data[day_data[,grep('Duration..in.seconds.',colnames(day_data))]>0,]
        day_data <- day_data[day_data[,grep('DistributionChannel',colnames(day_data))]!="preview",]
      }
      names(day_data)[names(day_data)=='RecipientEmail'] <- 'Email'
      #add to the other days' data by email
      survey <- merge(survey,day_data,by="Email",all=T) 
    }
  }
  
  #rename exp vars
  names(survey)[grep("^exp_GAD_([8-9]|[0-5]{2})",names(survey))] <- paste("exp_PHQ",seq(1:8),"1",sep="_")
  
  #Remove duplicate emails: people who completed the same survey more than once.
  survey$Email <- tolower(survey$Email)
  dups <- survey$Email[duplicated(survey$Email)]
  print(paste0('n with duplicated emails: ',length(unique(dups))))
  survey <- survey[!survey$Email %in% dups,]
  rownames(survey) <- NULL #reset row numbers
  
  #delete identifying information
  nonymous <- c('Location','UserAgent','meta_info','IPAddress','Email','email','prolific_id','SONA_ID')
  survey <- survey[,!grepl(paste(nonymous,collapse='|'),colnames(survey))]
  
  return(survey)
}

survey <- getSurvey()
saveRDS(survey,'rds_files/anon_raw_data.rds')
#survey <- readRDS('rds_files/anon_raw_data.rds')
write.csv(survey,'csv_files/anon_raw_data.csv')

cleanSurvey <- function(survey){
  #date cleanup
  dates <- data.frame(date_pre=as.Date(survey$EndDate_pre),date_post=as.Date(survey$EndDate_post)) #put these aside
  dates$time_taken <- dates$date_post - dates$date_pre
  
  #drop unneeded variables
  drops <- c("Status","IPAddress","Progress","Finished","RecordedDate","ResponseId",
             "Name","RecipientEmail","ExternalReference", "Location","DistributionChannel","UserLanguage",
             "Email.validation","Prolific.ID","SC0","PROLIFIC_PID","always1",
             "Date","Duration","Consent","Resolution","Click","timer",".Submit","BCT.M")
  
  survey <- survey[!grepl(paste(drops, collapse = "|"),names(survey))]
  #add dates back in
  survey[,names(dates)] <- dates
  
  #recoding
  survey$Gender_pre[survey$Gender_pre == 3] <- "male"
  survey$Gender_pre[survey$Gender_pre == 4] <- "female"
  survey$Gender_pre[survey$Gender_pre == 5] <- "other"
  survey$Gender_pre[survey$Gender_pre == 6] <- "opt_out"
  
  #type coversion
  non_numeric <- grepl("^(TMS|GAD|PHQ|SOW|RRS|WBSI)_[0-9]{1,2}_(pre|post)$|^exp_",names(survey))
  survey[,non_numeric] <- sapply(survey[,non_numeric],as.numeric)
  survey <- survey[!is.na(survey$condition_pre),] #missing data?
  survey$condition_pre <- factor(survey$condition_pre)
  rownames(survey) <- survey$Random_ID_pre
  
  #tag dropouts & per-protocol
  survey$per_protocol <- (!is.na(survey$Random_ID_post)) & (survey$condition_pre=='control'| (survey$condition_pre!='control' & !is.na(survey$focus_10)))
  survey$dropout <- is.na(survey$Random_ID_post)
  survey$returned <- survey$condition_pre!='control' & is.na(survey$focus_10) & !is.na(survey$Random_ID_post)
  
  return(survey)
}

br_data <- cleanSurvey(survey)
# Scale Re-numbering ---------------------------
#place on 0-max scale, rather than 1-max
scale_recode <- function(dataset, scale_name, q_nums) {#dataset<-dots_data; scale_name<-"TMS";q_nums<-c(1:7)
  q_nums <- as.character(q_nums)
  q_names <- c(paste0(scale_name,"_",q_nums,"_pre"), paste0(scale_name,"_",q_nums,"_post"))
  dataset[,q_names] <- sapply(dataset[,q_names],as.numeric)
  dataset[,q_names] <- dataset[,q_names] - 1
  return(dataset)
}

br_data <- scale_recode(br_data,"TMS",c(1:7)) #  TMS: https://sci-hub.yncjkj.com/10.1891/0889-8391.23.3.185.
br_data <- scale_recode(br_data,"GAD",c(1:7))
br_data <- scale_recode(br_data,"PHQ",c(1:8)) # PHQ-8: https://asset-pdf.scinapse.io/prod/1985329916/1985329916.pdf

# Scale and Subscale items list ---------------------------  
scales_list <- list(TMS=c(1:7), GAD=c(1:7), PHQ=c(1:8), SOW_F=seq(1,14,2), SOW_S=seq(2,14,2),
                    RRS=c(1:22), RRS_D=c(1:4,6,8,9,14,17:19,22), RRS_B=c(5,10,13,15,16), RRS_R=c(7,11,12,20,21),
                    WBSI=c(1:15), WBSI_sup=c(13,1,8,10,11,14), WBSI_int=c(2,3,4,5,6,7,9,12,15))

expandScales <- function(scale_name, q_nums){ #scale_name <- "SOW_1";q_nums<-seq(1,14,2)
  scale_name <- sub("_.*", "", scale_name) #for subscales
  q_nums <- as.character(q_nums)
  pre_scale <- paste(scale_name,q_nums,"pre",sep="_")
  post_scale <- paste(scale_name,q_nums,"post",sep="_")
  scale <- list(pre = pre_scale, post = post_scale)
  return(scale)
}

scales <- Map(expandScales, names(scales_list), scales_list)

for(i in 1:length(scales)){
  sc_names <- c(scales[[i]]$pre,scales[[i]]$post)
  br_data[,sc_names] <- sapply(br_data[,sc_names], as.numeric)
}

# Mean Scores and difference over time ---------------------------
meanDiffs <- function(dataset){ # dataset<-br_data
  #scales
  pre <- as.data.frame(sapply(scales,function(x){apply(dataset[x$pre],1,mean)}))
  post <- as.data.frame(sapply(scales,function(x){apply(dataset[x$post],1,mean)}))
  diff <- post-pre
  names(pre) <- paste0(names(pre),"_pre")
  names(post) <- paste0(names(post),"_post")
  names(diff) <- paste0(names(diff),"_diff")
  dataset <- cbind(dataset,pre,post,diff)
  
  #SDT
  for(SDT_var in SDT_vars){ #SDT_var<-SDT_vars[1]
    dataset[,paste0(SDT_var,"_diff")] <- dataset[,paste0(SDT_var,"_post")]-dataset[,paste0(SDT_var,"_pre")]
  }
  
  #expectancies
  for(expec in  c("TMS","GAD","PHQ")){
    col_names <- paste0("^exp_",expec)
    dataset[,paste0("exp_",expec)] <- apply(dataset[,grep(col_names,names(dataset))],1,mean)
  }
  
  return(dataset)
}

#br_data <- meanDiffs(br_data)
#saveRDS(br_data,"rds_files/final_dataset.rds")
#write.csv(br_data,"csv_files/final_dataset.csv")
br_data <- readRDS("rds_files/final_dataset.rds")

############################################################## Pre-analysis ----
# Demographics -------
demos <- c("Age mean" = round2(mean(as.numeric(br_data$Age_pre))),
           "Age range" = paste(range(as.numeric(br_data$Age_pre)),collapse=", "),
           "gender" = table(br_data$Gender_pre))

# Time taken -------
timeTaken <- function(dataset){ # dataset<-br_data
  tt_d <- dataset[!is.na(dataset$TMS_diff) & dataset$condition_pre!='control',c('TMS_diff','time_taken','condition_pre')]
  tt_d$condition_pre <-  droplevels(tt_d$condition_pre)
  #H1 calculated = regressing the TMS against time taken to completion
  tt_lm <- summary(lm('TMS_diff~time_taken',tt_d))$coefficients[2,1:2]
  #dividing by the raw difference in TMS scores between interventions by the recorded slope.
  ms <- round2(sapply(tt_d[,c('time_taken','TMS_diff')],tapply,tt_d$condition_pre,mean))
  ses <- round2(sapply(tt_d[,c('time_taken','TMS_diff')],tapply,tt_d$condition_pre,se))
  time_diff <- round2(rbind(m=apply(ms,2,function(x){x['mental']-x['world']}),
                     se=apply(ses**2,2,sum)**.5))
  tt_H1 <- round2(time_diff['m','TMS_diff']/tt_lm['Estimate'])
  #BFs
  BF <- bfrr(sample_mean=abs(time_diff['m',"time_taken"]),sample_se=time_diff['se',"time_taken"],sample_df=nrow(tt_d)-1,
             model="normal", mean=0, sd=abs(tt_H1), tail=1, criterion=3,
             rr_interval=list(mean=c(0,0),sd=c(0,0)),precision=.01)
  #outputs
  #plot(tt_d$time_taken,tt_d$TMS_diff,main='TMS change by time taken')
  outs <- round2(c(m=time_diff['m',],
                   se=time_diff['se',],
                   lm=tt_lm, H1=tt_H1,
                   BF=BF$BF,RR=BF$RR$sd),2)
  list(main_output=outs,mean=ms,se=ses)
}

time_taken <- timeTaken(br_data)


# Feedback -------
#How well did you stick to the course content?
content <- function(dataset){
  #No idea why the numbering ended up this way on qualtrics
  #content_q <- c('1'='I meditated significantly more than the course required',
  #               '5'='I did the course almost exactly',
  #               '3'='I meditated significantly less than the course asked for',
  #               '4'='I did not meditate or do mindfulness during this course')
  content <- dataset[,c("condition_pre","course_content_post")]
  content$course_content_post <- as.numeric(content$course_content_post)
  content <- content[!is.na(content$course_content_post),]
  content[content$course_content_post == 1,'course_content_post'] <- 'more'
  content[content$course_content_post == 5,'course_content_post'] <- 'exact'
  content[content$course_content_post == 3,'course_content_post'] <- 'less'
  content[content$course_content_post == 4,'course_content_post'] <- 'none'
  cont_t <- table(content)[c('mental','world'),] #controls fucking around?
  exact_more <- (cont_t[,'exact'] + cont_t[,'more'])/apply(cont_t,1,sum)*100
  cbind(cont_t,'Exact or More'=round2(exact_more))
}

cont <- content(br_data)

#Do you feel this course featured authentic mindfulness meditation? 4: any other comments?
useful <- function(dataset){
  #useful_q <- c('1'='I found this course very useful',
  #              '2'='I did not find it useful',
  #              '3'='I disliked this course',
  #              '4'='Other (please elaborate if you can!)')
  #qual
  use_qual <- dataset[!is.na(dataset$useful_4_TEXT_post),c("condition_pre","useful_post","useful_4_TEXT_post")]
  #quant
  useful <- dataset[!is.na(dataset$useful_post),c("condition_pre","useful_post")]
  useful$useful_post <- gsub(",[1-4]","",useful$useful_post)
  useful[useful$useful_post=='1', 'useful_post'] <- 'very'
  useful[useful$useful_post=='2', 'useful_post'] <- 'not_useful'
  useful[useful$useful_post=='3', 'useful_post'] <- 'dislike'
  useful[useful$useful_post=='4', 'useful_post'] <- 'other'
  use_t <- table(useful[,c("condition_pre","useful_post"),])[c('mental','world'),-1]
  cbind(use_t,'Very'=round2(use_t[,'very']/apply(use_t,1,sum)*100))
}

usef <- useful(br_data)

#Do you feel this course featured authentic mindfulness meditation?
authentic <- function(dataset){
  #qual
  auth_qual <- dataset[!is.na(dataset$authentic_3_TEXT_post),c("condition_pre","authentic_post","authentic_3_TEXT_post")]
  #quant
  auth <- dataset[!is.na(dataset$authentic_post) & dataset$authentic_post!='',c("condition_pre","authentic_post")]#,"authentic_3_TEXT_post")]
  auth$authentic_post <- gsub(",3","",auth$authentic_post) #remove extra comments
  auth[auth$authentic_post=='1', 'authentic_post'] <- 'yes'
  auth[auth$authentic_post=='2', 'authentic_post'] <- 'no'
  auth[auth$authentic_post=='3', 'authentic_post'] <- 'comments'
  auth_t <- table(auth)[c('world','mental'),c('no','yes')]
  cbind(auth_t,Yes=round2(auth_t[,'yes']/apply(auth_t,1,sum)*100))
}

auth <- authentic(br_data)

#Please let us know if you have any other comments or feedback on the course? + gen qual
qual <- function(dataset){
  feedback <- dataset[dataset$condition_pre!='control',c("condition_pre","feedback_post","authentic_post","authentic_3_TEXT_post","useful_post","useful_4_TEXT_post")]
  write.csv(feedback,'csv_files/feedback.csv')
}

quali <- qual(br_data)

# Focus/Activities -------
engagement <- function(dataset){ # dataset<-br_data
  #prep
  foc_ac <- dataset[grep('condition_pre|focus|activities|Random_ID_post',colnames(dataset))]
  foc_ac <- foc_ac[foc_ac$condition_pre!='control',]
  foc_ac$condition_pre <- droplevels(foc_ac$condition_pre)
  foc_ac <- foc_ac[apply(is.na(foc_ac),1,sum)==0,]
  #activities
  foc_ac$activ <- apply(foc_ac[,grep('activities',colnames(foc_ac))]==1,1,sum)>=7
  #focus
  foc_ac$foc <- apply(foc_ac[,grep('focus',colnames(foc_ac))],1,function(x){mean(as.numeric(x))})
  foc_ac$foc <- floor(foc_ac$foc/10) * 10 #round down to nearest 10
  foc_ac$foc_70 <- foc_ac$foc>=70
  #percentages >7
  output <- cbind(focus=tapply(foc_ac$foc_70,foc_ac$condition_pre,sum),
                  activities=tapply(foc_ac$activ,foc_ac$condition_pre,sum))
  output <- round2(apply(output,2,function(x){x/table(foc_ac$condition_pre)*100}))
  #odds ratio
  foc_or1 <- fisher.test(table(foc_ac$foc_70,foc_ac$condition_pre))
  foc_or2 <- paste0('OR=',round2(foc_or1$estimate),', CI95%=[',paste(round2(foc_or1$conf.int),collapse=", "),']')
  act_or1 <- fisher.test(table(foc_ac$activ,foc_ac$condition_pre))
  act_or2 <- paste0('OR=',round2(act_or1$estimate),', CI95%=[',paste(round2(act_or1$conf.int),collapse=", "),']')
  
  list(percents=output,foc=foc_or2,act=act_or2)
}

engage <- engagement(br_data)


# Feedback + Engagement table -------------
fdbk <- cbind('Content: Exact or More'=cont[,ncol(cont)],
                        'Useful: Very'=usef[,ncol(usef)],
                        'Authentic: Yes'=auth[,ncol(auth)],
                        'Focus >70%'=engage[,'focus'],
                        '7+ Activities'=engage[,'activities'])
fdbk_per <- apply(round2(fdbk,0), 2, function(x){paste0(x,"%")})
rownames(fdbk_per) <- rownames(fdbk)
write.csv(fdbk_per,'csv_files/engagement.csv')
# Exclusions ---------------------------
outliers <- function(var,cutoff=3){  # var<-br_data$H_M_ratio_post
  pos <- mean(var,na.rm=T)+(3*sd(var,na.rm=T))
  neg <- mean(var,na.rm=T)-(3*sd(var,na.rm=T))
  (var<pos & var>neg) | is.na(var)
} #not used, but +/-3SD cutoff on HMratio gets the same ps

exclude <- function(dataset,score0=T,rate=F){ # dataset<-br_data
  if(score0){
    #remove due to score of 0
    low_score <- which(dataset$dp_pre==0|dataset$dp_post==0|is.na(dataset$dp_pre)|is.nan(dataset$dp_diff)|
                      dataset$acc_pre==0|dataset$acc_post==0|is.na(dataset$acc_pre)|is.nan(dataset$acc_diff))
    print(length(low_score))
    #note same as = all(low_scores==which(dataset$acc_pre==0|is.nan(dataset$acc_pre))) #VIEW: dataset[excluded,grep('acc|dp|H_M',colnames(dataset))]
    if(length(low_score)>0){
      dataset <- dataset[-low_score,]
    }
  }
  
  if(rate){
    #remove due to breath rate #dataset[,c('Random_ID_pre','n_breath_pre','n_breath_post')]
    br_range <- c(12-5,20+5)*15 #+/-5bpm buffer around 12-20bpm (Sapra et al., 2022) for 15 mins
    br_rate <- which(dataset$n_breath_pre<br_range[1]|dataset$n_breath_pre>br_range[2]|
                       dataset$n_breath_post<br_range[1]|dataset$n_breath_post>br_range[2])
    print(length(br_rate))
    if(length(br_rate)>0){
      dataset <- dataset[-br_rate,]
    }
  }
  
  dataset$condition_pre <- droplevels(dataset$condition_pre)
  print(table(dataset$condition_pre))
  return(dataset)
}

# br_data <- readRDS("rds_files/final_dataset.rds")
og_n <- nrow(br_data)
br_data <- exclude(br_data,score0=T,rate=F) #rate gets all the people from score anyway

# Dropouts -------
dropoutsTab <- function(dataset){ # dataset <- br_data
  drops <- t(sapply(dataset[,c('per_protocol','returned','dropout')],tapply,dataset$condition_pre,sum))
  drops <- rbind(drops,total=apply(drops,2,sum))
  drops <- cbind(drops,total=apply(drops,1,sum))
  dimnames(drops) <- lapply(dimnames(drops),tools::toTitleCase)
  rownames(drops) <- gsub('_','-',rownames(drops))
  return(drops)
}
  
dropouts <- dropoutsTab(br_data)
write.csv(dropouts,'csv_files/dropouts.csv')

################################################################## Analysis ----
# Descriptives ----------
shortlist <- c("TMS","SOW_F","SOW_S","PHQ","GAD","RRS","WBSI","acc","dp","H_M_ratio")

groupMeans <- function(dataset){ # dataset<-br_data
  dataset <- dataset[!is.na(dataset$TMS_diff),]
  dataset$condition_pre <- droplevels(dataset$condition_pre)
  all_v <- gsub('_diff','',names(dataset)[grep('_diff',names(dataset))])
  out <- c()
  for(i in c("pre","post","diff")){ # i<-"pre"
    time_data <- dataset[paste0(all_v,"_",i)]
    conds <- t(sapply(time_data,tapply,dataset$condition_pre,function(x){ 
      paste0(round2(mean(x))," (",round2(se(x)),")")
    }))
    colnames(conds) <- paste0(colnames(conds),"_",i)
    rownames(conds) <- all_v
    out <- cbind(out,conds)
  }
  tab <- out[,c(grep("control",colnames(out)),grep("world",colnames(out)),grep("mental",colnames(out)))]
  
  #pretty print
  tab <- tab[shortlist,-grep("diff",colnames(tab))]
  colnames(tab) <- paste(rep(c("Waitlist","World",'Mental'),each=2),c("Pre","Post"),sep=" ")
  rownames(tab) <- pnames[shortlist]
  return(tab)
}

group_means <- groupMeans(br_data)
write.csv(group_means,'csv_files/group_means.csv')

# Expectancies ---------
expectancies <- function(dataset){ # dataset<-br_data
  #setup
  exp_data <- dataset[!is.na(dataset$exp_TMS),c('condition_pre','exp_TMS','exp_GAD','exp_PHQ','exp_2SO_1_1','exp_WBSI_1_1','exp_RRS_1_1')]
  exp_data$condition_pre <- droplevels(exp_data$condition_pre)
  #condition means
  exp_m <- sapply(exp_data[,-1],function(x){
    c(mean=tapply(x,exp_data$condition_pre,mean),
      se=tapply(x,exp_data$condition_pre,se))
  })
  exp_m <- round2(exp_m)
  #names(exp_cond) <- paste(levels(exp_data$condition_pre),rep(c("m","se"),each=2),sep="_")
  exps <- round2(cbind(MS_m=exp_m['mean.mental',],MS_se=exp_m['se.mental',],
                       W_m=exp_m['mean.world',],W_se=exp_m['se.world',],
                       m=exp_m['mean.mental',]-exp_m['mean.world',],
                       se=sqrt(exp_m['se.mental',]^2+exp_m['se.world',]^2),
                       BF=NA,RRl=NA,RRu=NA))
  rownames(exps) <- gsub("exp_([1-9A-Z]{3,4})(_1_1)?","\\1",rownames(exps))
  #BF
  for(i in 1:nrow(exps)){ #i<-1
    #TMS:0-4;GAD,PHQ:0-3;SOW:1-5; RRS:1-4; WBSI:1-5
    BF <- bfrr(sample_mean=abs(exps[i,'m']),sample_se=exps[i,'se'],sample_df=nrow(exp_data)-1,
               model="normal", mean=.2, sd=.1, tail=2, criterion=3,
               rr_interval=list(mean=c(.2,.2),sd=c(0,5)),precision=.01)
    exps[i,c('BF','RRl','RRu')] <- c(round2(BF$BF),BF$RR$sd)
  }
  
  #table cleanup
  exps_cl <- exps
  #RRs
  exps_cl[exps_cl[,'BF']>=.33|exps_cl[,'BF']<3,'RRl'] <- 0
  exps_cl[exps_cl[,'RRl']==.01,'RRl'] <- 0
  # Format 2 d.p.
  exps_cl <- format(exps_cl, nsmall=2,trim=T) # ensure 2 d.p. trailing 0s
  #Mean SE concat
  exps_cl[,'MS_m'] <- paste0(exps_cl[,'MS_m']," (",exps_cl[,'MS_se'],")")
  exps_cl[,'W_m'] <- paste0(exps_cl[,'W_m']," (",exps_cl[,'W_se'],")")
  exps_cl[,'m'] <- paste0(exps_cl[,'m']," (",exps_cl[,'se'],")")
  #RR concat
  exps_cl[,'RRl'] <- paste0("[",exps_cl[,'RRl'],", ",exps_cl[,'RRu'],"]")
  #select cols
  exps_cl <- exps_cl[,c('MS_m','W_m','m','BF','RRl')]
  #rename
  colnames(exps_cl) <- c("MS Mean (SE)","W Mean (SE)","Mean diff. (SE)","BF","RR")
  rownames(exps_cl) <- pnames[rownames(exps_cl)]
  return(exps_cl)
}

expects <- expectancies(br_data)
write.csv(expects,"csv_files/expectancies.csv")
# Analysis -----------
shortlist_group <-c("TMS","SOW_F","SOW_S","PHQ","GAD","RRS","WBSI","adj_dp2","adj_meta_d","H_M_ratio")

interactions <- function(dataset){ # dataset <- br_data
  dataset <- dataset[!is.na(dataset$TMS_diff),]
  diffs <- dataset[grepl("condition|diff",names(dataset))]
  colnames(diffs) <- gsub("_diff",'',colnames(diffs))
  #get means and SEs
  ms <- t(sapply(diffs[-1],function(x){tapply(x,diffs$condition_pre,mean)}))
  ses <- t(sapply(diffs[-1],function(x){tapply(x,diffs$condition_pre,se)}))
  ses_sq <- ses**2
  #comparisons listing function
  compareGroup <- function(comp){ #comp <- comparisons[[1]]
    #means and ses
    tc <- cbind(m=ms[,comp[1]]-ms[,comp[2]], 
                se=(ses_sq[,comp[1]]+ses_sq[,comp[2]])**.5)
    #Covariate adjustment
    data_c <- diffs[diffs$condition_pre%in%comp[1:2],]
    data_c$condition_pre <- droplevels(data_c$condition_pre)
    tc <- rbind(tc, 
                'adj_dp' = adj(data_c,'dp','c'), #adjusted type-1 d'
                'adj_dp2' = adj(data_c,'dp2',c('c2','dp','c')), #adjusted type-2 d'
                'adj_meta_d'= adj(data_c,'H_meta_d',c('dp','c'))) #adjusted Hmeta-d
    #add cols and names
    tc <- cbind(tc,bf=NA,rrl=NA,rru=NA)
    colnames(tc) <- paste0(comp[3],'_',colnames(tc))
    return(tc)
  }
  #comparisons
  comparisons <- list(c("mental","world",name="mvw"),c("world","control",name="wvwl"),c("mental","control",name="mvwl"))
  comps <- do.call(cbind,lapply(comparisons,compareGroup))
  #shortlist
  #comps <- comps[shortlist_group,]
  list(n=nrow(diffs),comps=comps) #allows n to be pulled through
}

BFs <- function(dataset){ # dataset <- br_data
  ints <- interactions(dataset)
  comps <- cbind(H1=NA,round2(ints$comps))
  comps <- comps[!rownames(comps)%in%c('c','c1','c2','ln_M_ratio'),] # c is too small to run - log(M) has NAs
  # Bayes Factors
  for(comp in c('mvw','wvwl','mvwl')){ #comp<-'mvw'
    for(i in rownames(comps)){ #i<-rownames(comps)[1] 
      #print(paste0(comp,': ',rownames(comps)[i]))
      sample_mean <- comps[i,paste0(comp,'_m')]
      if(grepl('GAD|PHQ|RRS|WBSI',i)){ sample_mean <- sample_mean*-1 } # Hypothesised decrease
      #priors
      if(i%in%rownames(H1s)){
        comps[i,'H1'] <- H1s[i,'m']
      } else { comps[i,'H1'] <- H1s['scales','m'] }
      #calc BFs
      BF <- bfrr(sample_mean=sample_mean,sample_se=comps[i,paste0(comp,'_se')],
                 sample_df=ints$n-1,model="normal", mean=0, sd=comps[i,'H1'], tail=1,
                 criterion=3,rr_interval=list(mean=c(0,0),sd=c(0,5)),precision=.01)
      #add to tabe
      comps[i,paste(comp,c('bf','rrl','rru'),sep='_')] <- c(round2(BF$BF),BF$RR$sd)
    }
  }
  round2(comps)
}

res_pp <- BFs(br_data)
############################################################## Presentation ----
# Format results ------
shortlist_group <-c("TMS","SOW_F","SOW_S","PHQ","GAD","RRS","WBSI","acc","adj_dp","adj_dp2","adj_meta_d","H_M_ratio")

formatResults <- function(tab){ # tab<-res_pp
  max_rr <- tab[,grep("rru",colnames(tab))]==5
  tab[,colnames(max_rr)][max_rr] <- '>5'
  min_rr <- tab[,grep("rrl",colnames(tab))]==.01
  tab[,colnames(min_rr)][min_rr] <- '0'
  diffs <- cbind(
    'H1'= tab[,"H1"],
    'MS-World Mean'= paste0(tab[,"mvw_m"]," (",tab[,"mvw_se"],")"),
    'MS-World BF'= tab[,"mvw_bf"],
    'MS-World RR'= paste0("[",tab[,"mvw_rrl"],", ",tab[,"mvw_rru"],"]"),
    
    'World-Waitlist Mean'= paste0(tab[,"wvwl_m"]," (",tab[,"wvwl_se"],")"),
    'World-Waitlist BF'= tab[,"wvwl_bf"],
    'World-Waitlist RR'= paste0("[",tab[,"wvwl_rrl"],", ",tab[,"wvwl_rru"],"]"),
    
    'MS-Waitlist Mean'= paste0(tab[,"mvwl_m"]," (",tab[,"mvwl_se"],")"),
    'MS-Waitlist BF'= tab[,"mvwl_bf"],
    'MS-Waitlist RR'= paste0("[",tab[,"mvwl_rrl"],", ",tab[,"mvwl_rru"],"]")
  )
  short <- diffs[shortlist_group,]
  rownames(short) <- pnames[shortlist_group]
  return(short)
}

res_pp_f <- formatResults(res_pp)
write.csv(res_pp_f,'csv_files/PP_results.csv')

# Bar graphs -----------
dir.create('graphs', showWarnings = FALSE)
shorts <- unique(c(shortlist,shortlist_group))
colours <- turbo(length(shorts))
names(colours) <- shorts
pconds <- c(control='Waitlist',world='World',mental='Mental States')

# plot diffs
plotDiff <- function(dataset){ # dataset<-dots_data
  #Setup
  dataset <- dataset[!is.na(dataset$TMS_diff),]
  diffs <- dataset[grepl("diff",names(dataset))]
  colnames(diffs) <- gsub("_diff",'',colnames(diffs))
  #get means and CIs
  m <- t(sapply(diffs,function(x){tapply(x,dataset$condition_pre,mean)}))[shortlist,names(pconds)]
  cis <- t(sapply(diffs,function(x){tapply(x,dataset$condition_pre,ci)}))[shortlist,names(pconds)]
  #formatting
  colnames(m) <- pconds[colnames(m)]
  par(mar=rep(2,4))
  #plot
  plt <- barplot(m,col=colours[shortlist],beside=T,ylim=c(min(m-cis),max(m+cis)),main='Change in group means over time (+/- 95% CI)')
  arrows(plt, m-cis, plt, m+cis,angle=90,code=3,length=.05)
  legend("topright",legend=pnames[shortlist],col=colours[shortlist],pch=15,cex=0.7)
}

jpeg('graphs/group_means.jpeg', height=150, width=200,units='mm',quality=200,res=300)
plotDiff(br_data)
dev.off()

plotRes <- function(res,main='Decomposed Time*Condition Interactions (+/- 95% CI)'){ # res<-res_pp
  #res<-res_pp
  short <- res[shortlist_group,]
  m <- short[,grep('_m',colnames(short))]
  cis <- short[,grep('_se',colnames(short))]*qnorm(.975)
  #format
  pcomps <- c(mvw_m='Mental States-World', wvwl_m='World-Waitlist', mvwl_m='Mental States-Waitlist')
  colnames(m) <- pcomps[colnames(m)]
  #colours, desaturate exploratory analyses
  cols <- rep(colours[shortlist_group],3)
  desat <- function(cols, sat=0.5) { #https://stackoverflow.com/a/26322593/7705626
    x <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
    hsv(x[1,], x[2,], x[3,])
  }
  cols[(length(shortlist_group)*2+1):length(cols)] <- desat(cols[(length(shortlist_group)*2+1):length(cols)],.5)
  #arr_cols <- c(rep('black',length(shortlist_group)*2),rep('darkgrey',length(shortlist_group)))
  #plot
  par(mar=rep(3,4),xpd=T)
  plt <- barplot(m,col=cols,beside=T,ylim=c(min(m-cis),max(m+cis)),main=main)
  arrows(plt, m-cis, plt, m+cis,angle=90,code=3,length=.05)
  legend('topright',legend=pnames[rownames(m)],col=colours[shortlist_group],pch=15,cex=0.7, inset=c(-0.05, 0))
}

jpeg('graphs/interactions.jpeg', height=150, width=200,units='mm',quality=200,res=300)
plotRes(res_pp,'Time*Condition Interactions (Per-Protocol)')
dev.off()

################################################################## Save env ----
save.image(file='BCTM.RData')