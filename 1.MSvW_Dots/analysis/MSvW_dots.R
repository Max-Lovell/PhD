# Analaysis code for MSvW Dots study by Max Lovell #
# Setup ----
#make working dir same as script file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#load environment if code already run
if(file.exists('MSvW.RData')){
  load('MSvW.RData')
}
# Packages ----
list.of.packages <- c("jsonlite", # read current study's raw data
                      "osfr",     # download schmidt's data
                      "archive",  # unpack compressed files from previous studies
                      "foreign",  # Read Schmidt's data
                      'rmatio',   # read carpenter's data from matlab
                      'openxlsx', # read file on previous scale measure effects
                      "mice",     # multiple imputation
                      "bfrr",     # Bayes Factors and RR search
                      "viridis"   # colour palette
                      )

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

source("fit_meta_d_MLE.R")
source("fit_meta_d_mcmc.R")

#might also need these:
#install.packages("devtools")
#devtools::install_github("debruine/bfrr")
# Helpers -----
#names
pnames <- c("TMS"="TMS-D", "SOW"="SOW", "SOW_F"="SOW 1st", "SOW_S"="SOW 2nd",
            "PHQ"="PHQ-8", "GAD"="GAD-7", "RRS"="RRS",
            "RRS_D"="RRS Depression", "RRS_B"="RRS Brooding", "RRS_R"="RRS Rumination",
            "WBSI"="WBSI", "WBSI_sup"="WBSI Supression", "WBSI_int"="WBSI Intrusion",
            "acc"="Accuracy", "dp"="d'", "c"="c", "adj_dp"="Adj. d'",
            "adj_meta_d"="Adj. meta-d'", "H_meta_d"="Hmeta-d", "H_M_ratio"="HM-ratio",
            "ln_M_ratio"="log(M-ratio)")
shortlist <- c("TMS","SOW_F","SOW_S","PHQ","GAD","RRS","WBSI","dp","H_M_ratio")
shortlist_group <-c("TMS","SOW_F","SOW_S","PHQ","GAD","RRS","WBSI","adj_meta_d","H_M_ratio")

#errors
se <- function(x,na.rm=F){ 
  if(na.rm){x<-x[!is.na(x)]}
  sd(x)/sqrt(length(x)) 
  }
ci <- function(x){ se(x)*qnorm(.975) }

#Short pivot function
widen <- function(dataset,id_cols,names_from,values_from){ # dataset=mfn;id_cols=c("Id","Grupo_");names_from="Session__";values_from=c('dprime','mdDiff',"mdRatio",'ln_mdRatio')
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
# SDT function ----------------
SDT <- function(dataset,resps=c('l','r'),ML=F){ # dataset<-dots_task 
  #Experimentor error: confidence scale of 0-100 for some subjects
  print(paste0('N lost to conf scale size: ',length(unique(dataset[dataset$conf %in% (0:10*10),"ID"]))))
  dataset <- dataset[!dataset$conf %in% (0:10*10),]
  dataset$ID <- droplevels(dataset$ID)
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
  #adjust for 0 counts
  nr_sx_adj <- nr_sx+.5 #from Hmeta-d is adjusted by (1/(n_rat*2)) instead - should standardise this?
  #t1 dp and c
  colnames(nr_sx_adj)
  FAR <- apply(nr_sx_adj[,(nconf+1):(nconf*2)],1,sum) / apply(nr_sx_adj[,1:(nconf*2)],1,sum)
  HR <- apply(nr_sx_adj[,(nconf*3+1):(nconf*4)],1,sum) / apply(nr_sx_adj[,(nconf*2+1):(nconf*4)],1,sum)
  dp <- qnorm(HR) - qnorm(FAR)
  c <- -.5*(qnorm(HR)+qnorm(FAR))
  #combine
  output <- data.frame(dp=dp,c=c) #dataframe for easier column assignment below
  #Meta-d'
  mcmc_params <- list(response_conditional=0,estimate_dprime=0,nchains=3,
                      nadapt=1000,nburnin=1000,nsamples=10000,nthin=1,
                      dic=0,rhat=0,parallel=1,monitorparams=0,saveallsamples=0,
                      estimate_mratio=0)
  if(ML) {
    for(i in 1:nrow(nr_sx)){ # i<-1
      output[i,c("ML_meta_d","ML_M_ratio")] <- fit_meta_d_MLE(nr_sx_adj[i,grep('S1',colnames(nr_sx))],nr_sx_adj[i,grep('S2',colnames(nr_sx))])[c('meta_da','M_ratio')]
    }
    output$ln_M_ratio <- log(output$ML_M_ratio)
  } else {
    for(i in 1:nrow(nr_sx)){ # i<-1
      fit <- fit_meta_d_mcmc(nr_sx[i,grep('S1',colnames(nr_sx))],nr_sx[i,grep('S2',colnames(nr_sx))],mcmc_params)
      output[i,c("d1","c1","H_meta_d","H_M_ratio")] <- fit[c("d1","c1","meta_d","M_ratio")]
    }
    output$ln_M_ratio <- log(output$H_M_ratio)
  }

  #Add ID and time
  output$ID <- rownames(nr_sx)
  return(output)
}

#################################################################### Priors ----
# Schmidt et al. (2019) -----------

# Note - see Breath study for better approach to this. For now, H1s must be in line with Registered Report.
downloadSchmidt <- function(){
  #Download
  dir.create('priors', showWarnings = FALSE)
  dir.create('priors/schmidt', showWarnings = FALSE)
  dwnld <- osf_download(osf_ls_files(osf_retrieve_node("https://osf.io/trwq3")), path='priors/schmidt', conflicts = "overwrite")
  archive_extract('priors/schmidt/Data MetaMFN.rar','priors/schmidt')
}

if(!file.exists('priors/schmidt/Data MetaMFN.rar')){
  downloadSchmidt()
}

oldAdj <- function(dataset,metad,dp,condition){ #dataset=carpenter_raw;metad='diff_metad';dp='diff_d';condition='condition'
  lm_out <- lm(paste0(metad,'~',condition,'+',dp), data=dataset)
  lm_sum <- summary(lm_out)
  exp_adj  <- lm_sum$coefficients[1,1] + lm_sum$coefficients[2,1] + (lm_sum$coefficients[3,1]*mean(dataset[,dp]))
  cont_adj  <- lm_sum$coefficients[1,1] + (lm_sum$coefficients[3,1]*mean(dataset[,dp]))
  exp_adj-cont_adj
}

schmidtDelta <- function(){
  DeltaBases <- read.spss("priors/schmidt/Data/SPSS/DeltaBases.sav", to.data.frame=TRUE)
  #change vars
  delta <- DeltaBases[DeltaBases$Task_ == "gab",c('Id','Grupo_','deltaDprime','deltaMeta_Dprime','mdRatio_Delta21')]
  #pre-post vars for log(M-ratio)
  mfn <- read.spss("priors/schmidt/Data/SPSS/MFN_.sav", to.data.frame=TRUE)
  mfn <- mfn[mfn$Task_=="gab",c("Id","Grupo_","Session__",'dprime','metaDprime',"mdRatio")]
  #create log(m-ratio)
  mfn$ln_mdRatio <- log(mfn$mdRatio)
  #widen dataset
  mfn_w <- widen(mfn,c("Id","Grupo_"),"Session__",c('dprime',"metaDprime",'ln_mdRatio'))
  #difference vars
  mfn_w[,gsub('post','diff',names(mfn_w)[grep('post',names(mfn_w))])] <- mfn_w[,grep('post',names(mfn_w))]-mfn_w[,grep('pre',names(mfn_w))]
  #merge log(M-ratio). Note some people with missing data on other vars was included for registration so has to stay....
  mfn_s <- mfn_w[,c('Id','Grupo_','ln_mdRatio_diff')]
  sch_del <- merge(delta,mfn_s,by=c('Id','Grupo_'),all=T) 
  #mean diffs
  sch_del[,c('deltaMeta_Dprime','mdRatio_Delta21','ln_mdRatio_diff')] <- lapply(sch_del[,c('deltaMeta_Dprime','mdRatio_Delta21','ln_mdRatio_diff')],as.numeric)
  ms <- sapply(sch_del[,c('deltaMeta_Dprime','mdRatio_Delta21','ln_mdRatio_diff')],tapply,sch_del$Grupo_,mean,na.rm=T)
  ses <- sapply(sch_del[,c('deltaMeta_Dprime','mdRatio_Delta21','ln_mdRatio_diff')],tapply,sch_del$Grupo_,se,na.rm=T)
  #interactions
  ints <- rbind(mean=ms['in',]-ms['ex',], se=apply(ses**2,2,sum)**.5)
  adj_metad <- oldAdj(mfn_w,'metaDprime_diff','dprime_diff','Grupo_')
  cbind(ints,adj_metad=c(adj_metad,ints['se','deltaMeta_Dprime']))
}

sch_ints <- schmidtDelta()
sch_se_adj <- sch_ints['se',]/(sqrt(27)/sqrt(61))

# Carpenter et al. (2019) ------

downloadCarpenter <- function(){
  dir.create('priors', showWarnings = FALSE)
  dir.create('priors/carpenter', showWarnings = FALSE)
  download.file('https://github.com/metacoglab/CarpenterMetaTraining/archive/refs/heads/master.zip','priors/carpenter/carpenter.zip')
  archive_extract('priors/carpenter/carpenter.zip','priors/carpenter')
}

#note MATLAB files to extract this data are in the carpenter folder
#if(!file.exists('priors/carpenter/carpenter.zip')){
#  downloadCarpenter()
#}

carpenterDiffs <- function(){
  #adjusted meta-d' means
  carpenter_raw <- read.csv("priors/carpenter/carpenter_raw.csv") #can calculate this with relevant matlab file
  #dprime
  carpenter_raw$diff_d <- carpenter_raw$s10_d-carpenter_raw$s1_d  
  carpenter_raw$diff_metad <- carpenter_raw$s10_metad-carpenter_raw$s1_metad
  #mratio
  carpenter_raw$s1_mratio <- carpenter_raw$s1_metad/carpenter_raw$s1_d
  carpenter_raw$s10_mratio <- carpenter_raw$s10_metad/carpenter_raw$s10_d
  carpenter_raw$diff_mratio <- carpenter_raw$s10_mratio-carpenter_raw$s1_mratio
  #logmratio
  carpenter_raw$s1_ln_mratio <- log(carpenter_raw$s1_mratio)
  carpenter_raw$s10_ln_mratio <- log(carpenter_raw$s10_mratio)
  carpenter_raw$diff_ln_mratio <- carpenter_raw$s10_ln_mratio-carpenter_raw$s1_ln_mratio
  #means and ses
  ms <- sapply(carpenter_raw[,c('diff_metad','diff_mratio','diff_ln_mratio')],tapply,carpenter_raw$condition,mean)
  ses <- sapply(carpenter_raw[,c('diff_metad','diff_mratio','diff_ln_mratio')],tapply,carpenter_raw$condition,se)
  #interactions
  ints <- rbind(mean=ms['exp',]-ms['control',], se=apply(ses**2,2,sum)**.5)
  adj_metad <- oldAdj(carpenter_raw,'diff_metad','diff_d','condition')
  cbind(ints,adj_metad=c(adj_metad,ints['se','diff_metad']))
}

carp_ints <- carpenterDiffs()

# SDT Priors ------

both_ints <- round2(cbind(
  mean = apply(rbind(sch_ints['mean',],carp_ints['mean',]),2,mean),
  se_paper = sch_se_adj,
  se_actual = apply(rbind(sch_se_adj,carp_ints['se',])**2,2,sum)**.5)
  )
rownames(both_ints) <- c('H_meta_d','H_M_ratio','ln_M_ratio','adj_meta_d')

# Estimate SDT sample sizes -------------
sampleSize <- function(n_guess){
  prev_n <- c(schmidt_control = 13,
              schmidt_exp = 14,
              carpenter_control = 32,
              carpenter_exp = 29)
  prior_n <- round2(1/mean(1/prev_n),0) #harmonic mean
  H1s <- cbind(both_ints[,c('mean','se_paper')],est_n=NA,n_BF=NA)
  for(v in 1:length(n_guess)){ #v<-1
    v <- n_guess[v]
    n_group <- v/3
    m_est <- abs(H1s[names(v),"mean"])#/2 #if taking as a plausible maximum
    SE_est  <- H1s[names(v),"se_paper"]*sqrt(prior_n/n_group)
    BF <- bfrr(sample_mean=0,sample_se=SE_est,
               sample_df=n_group-1,model="normal", mean=0, sd=m_est, tail=1,
               criterion=3,rr_interval=list(mean=c(0,0),sd=c(0,5)),precision=.01)
    H1s[names(v),c('est_n','n_BF')] <- c(v,BF$BF)
  }
  H1s <- H1s[names(n_guess),]
  return(H1s)
}

n_guess <- c("H_meta_d"=110,"H_M_ratio"=300,"ln_M_ratio"=330,"adj_meta_d"=160)
n_est <- round2(sampleSize(n_guess))
write.csv(n_est,'csv_files/estimated_sample_sizes.csv')

# Scale Measures -----
#get previous scale measures
prevScales <- function(){
  prev <- openxlsx::read.xlsx("priors/MSvW_prior_effects.xlsx", sheet = "scales_data")
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
# H1s Table -----
H1s <- cbind(m=round2(both_ints[,'mean'],1),se=both_ints[,'se_paper'])
H1s <- rbind(H1s,scales=c(scales_H1,NA))     

################################################################### Compile ----
# Dots Task -------------------
dotsData <- function(){ #note 96831 and 77712 are corrupted, hence warnings
  task_files <- list.files("./task_data")
  ##Vanilla JS files
  vjs_files <- task_files[grep("_p(re|ost).json$",task_files)]
  task <- c()
  for(file_n in vjs_files){ #file_n<-vjs_files[10]
    file <- fromJSON(paste0("task_data/",file_n))
    if(nrow(file)<194){next}
    file$ID <- sub("\\D+","",file_n)
    task <- rbind(task,file)
  }
  #simplify
  task <- task[task$confidence!=-1,c("ID","target","response","confidence")] #delete practice trials
  task[,c("target","response")] <- apply(task[,c("target","response")],2,function(x){substr(x,1,1)})
  names(task)<-c("ID","stim","resp","conf")
  
  ##JSPsych files
  jsp_files <- task_files[grep("_short.csv$",task_files)]
  #options(warn=2)
  for(file_n in jsp_files){ # file_n<-jsp_files[1]#"dots_77712_short.csv"   
    file <- read.table(paste0("task_data/",file_n),header=TRUE,stringsAsFactors=F,sep=",")
    file <- file[file$trial_type!="practice",c("trial_type","response","correct_response")]
    row.names(file) <- NULL
    if(nrow(file)/2<168){next} #note 96831 and 77712 are corrupted, hence warnings
    file[!!as.numeric(row.names(file))%%2,"confidence"] <- file[!as.numeric(row.names(file))%%2,"response"]
    file <- file[file$trial_type=="experimental",c("response","correct_response","confidence")]
    file$ID <- gsub("[a-z_.]","\\1",file_n)
    file <- file[,c("ID","correct_response","response","confidence")]
    names(file) <- c("ID","stim","resp","conf")
    file$stim <- ifelse(file$stim=="e","l","r")
    file$resp <- ifelse(file$resp=="e","l","r")
    task <- rbind(task,file)
  }
  
  task$ID <- factor(task$ID, levels=unique(task$ID)) #ensures they appear in the table() output
  
  return(task)
}

dots_task <- dotsData()
#dots_SDT <- SDT(dots_task)
#saveRDS(dots_SDT,'rds_files/dots_SDT.rds')
#write.csv(dots_SDT,'csv_files/dots_SDT.csv')
dots_SDT <- readRDS('rds_files/dots_SDT.rds')
SDT_vars <- colnames(dots_SDT)[-ncol(dots_SDT)]

#note: quite a few negative meta-d' values: dots_metad[dots_metad$meta_d<0,]

# Load/Merge Survey Data -----------------
#before hand in would be good to do this - will be easier to anonymise emails and can shorten code
#pre <- read.csv('./survey_data_redownload/survey_pre.csv')
#survey2 <- merge(pre,dots_SDT, by.x="Random_ID", by.y="ID", all.x=F, all.y=F) #add task data #note technical error means some task data incomplete - remove
#survey$email_valid[duplicated(survey$email_valid)]

surveyData <- function(){ #file_rgx<-all_rgxs[2]
  #NOTE RAW DATA AVAILABLE ON REQUEST
  #system('pscp -pw [PASSWORD] [USER EMAIL]:[SERVER PATH]*.csv "[LOCAL PATH]"')
  survey_files <- list.files("./survey_data")
  all_rgxs <- paste("^metacog",c("pre",paste0("[w|m]",1:10),"post"),"",sep="_")
  #merge files
  for(file_rgx in all_rgxs){ #file_rgx<-all_rgxs[1] #loop through the files for pre, post, and each day separately # file_rgx<-all_rgxs[1]
    day_files <- survey_files[grep(file_rgx,survey_files)] #read in first dataset to add the others to
    day_data <- read.csv(paste0("survey_data/",day_files[1]))
    day_data <- day_data[3:nrow(day_data),] #remove internal qualtrics headers
    for(file_name in 2:length(day_files)){ #loop through the rest and add the rows
      file_name <- day_files[file_name]
      day_file <- read.csv(paste0("survey_data/",file_name))
      day_file <- day_file[3:nrow(day_file),!names(day_file)%in%c("ismobile","px_cm")] #added these to exp later
      day_data <- rbind(day_data,day_file) #add file to rest of data for that day
    }
    suffix <- gsub('\\^|(metacog)|_|(\\[w\\|m\\])','',file_rgx) #add to the end of variables to tag as relevant day
    names(day_data)[!grepl("Email",names(day_data))] <- paste(names(day_data)[!grepl("Email",names(day_data))],suffix,sep="_")
    #bit of data cleaning, fix variable names, etc
    if(file_rgx=="^metacog_pre_"){
      print(paste0('original pre-test sample size: ',nrow(day_data)))
      dots_pre <- dots_SDT
      names(dots_pre) <- paste(names(dots_pre),suffix,sep="_")
      names(day_data)[names(day_data)=='email_pre'] <- 'Email'
      #pre-clean
      print(paste0('n unnfinished responses: ',sum(day_data$Progress_pre!=100 | day_data$Duration..in.seconds._pre<=0 | day_data$Finished_pre!=1 | day_data$DistributionChannel_pre=="preview")))
      day_data <- day_data[day_data$Progress_pre==100 & day_data$Duration..in.seconds._pre>0 & day_data$Finished_pre==1 & day_data$DistributionChannel_pre!="preview",]
      print(paste0('n with duplicated emails: ',sum(day_data$Email%in%day_data$Email[duplicated(day_data$Email)])))
      day_data <- day_data[!day_data$Email%in%day_data$Email[duplicated(day_data$Email)],] #THIS MIGHT NOT BE GOOD IF DATA WAS EVER DOWNLOADED TWICE
      print(paste0('n pre-test before task data merge: ',nrow(day_data)))
      survey <- merge(day_data,dots_pre, by.x="Random_ID_pre", by.y="ID_pre", all.x=F, all.y=F) #add task data #note technical error means some task data incomplete - remove
      print(paste0('n pre-test after task data merge: ',nrow(survey)))
    } else {
      if(file_rgx=="^metacog_post_"){#fixing some badly named vars in qualtrics:
        day_data <- day_data[day_data$Progress_post==100 & day_data$Duration..in.seconds._post>0 & day_data$Finished_post==1 & day_data$DistributionChannel_post!="preview",]
        dots_post <- dots_SDT
        names(dots_post) <- paste(names(dots_post),suffix,sep="_")
        day_data <- merge(day_data,dots_post, by.x="Random_ID_post", by.y="ID_post", all.x=F, all.y=F) #add breath task data
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

      #merge
      survey <- merge(survey,day_data,by="Email",all=T) #add to the other days' data by email
    }
  }
  
  #Remove duplicate emails: people who completed the same survey more than once.
  survey$Email <- tolower(survey$Email)
  dups <- survey$Email[duplicated(survey$Email)]
  print(paste0('n with duplicated emails: ',length(unique(dups))))
  survey <- survey[!survey$Email %in% dups,]
  rownames(survey) <- NULL #reset row numbers

  #delete identifying information
  nonymous <- c('Location','UserAgent','meta_info','IPAddress','Email','email','prolific_id')
  survey <- survey[,!grepl(paste(nonymous,collapse='|'),colnames(survey))]
  return(survey)
}

#survey <- surveyData()
#saveRDS(survey,'rds_files/anonymised_data.rds')
#write.csv(survey,'csv_files/anonymised_data.csv')
survey <- readRDS('rds_files/anonymised_data.rds')

# Clean dataset  -----------------

cleanSurvey <- function(survey){ # write.csv(survey,'csv_files/survey.csv')
  survey <- survey[!is.na(survey$Random_ID_pre),] # might be where pre-test task data wasn't recieved?
  #convert to numeric
  non_numeric <- grepl("^(TMS|GAD|PHQ|X2SO|RRS|WBSI)_[0-9]{1,2}_(pre|post)$|^exp_",names(survey))
  survey[,non_numeric] <- sapply(survey[,non_numeric],as.numeric)
  
  #Some people missing one TMS question
  print(paste0('n missing one TMS question: ',sum((!is.na(survey$Random_ID_post) & (is.na(survey$TMS_4_post)|is.na(survey$TMS_5_post))))))
  survey <- survey[!(!is.na(survey$Random_ID_post) & (is.na(survey$TMS_4_post)|is.na(survey$TMS_5_post))),]

  #Gender
  survey$gender_pre[survey$gender_pre == 3] <- "male"
  survey$gender_pre[survey$gender_pre == 4] <- "female"
  survey$gender_pre[survey$gender_pre == 5] <- "other"
  survey$gender_pre[survey$gender_pre == 6] <- "opt_out"
  
  #RENAMING: PHQ expectancies & Orders of Sensory Observation Scale
  names(survey)[grep("^exp_GAD_([8-9]|[0-5]{2})",names(survey))] <- paste("exp_PHQ",seq(1:8),"1",sep="_")
  names(survey) <- gsub('X2SO','SOW',names(survey))

  #date cleanup
  dates <- data.frame(date_pre=as.Date(survey$EndDate_pre),date_post=as.Date(survey$EndDate_post)) #put these aside
  dates$time_taken <- dates$date_post - dates$date_pre
  
  #drop unneeded variables
  drops <- c("Status","IPAddress","Progress","Finished","RecordedDate","ResponseId",
             "Name","RecipientEmail","ExternalReference", "Location","DistributionChannel","UserLanguage",
             "Email.validation","Prolific.ID","SC0","PROLIFIC_PID","always1",
             "Date","Duration","Consent","Resolution","Click","timer",".Submit")
  survey <- survey[!grepl(paste(drops, collapse = "|"),names(survey))]
  
  #add date back in
  survey[,names(dates)] <- dates
  
  #condition variable
  print(paste0('n missing condition: ',sum(is.na(survey$condition_pre))))
  survey <- survey[!is.na(survey$condition_pre),]
  survey$condition_pre <- factor(survey$condition_pre)
  
  #tag dropouts & per-protocol
  survey$per_protocol <- (!is.na(survey$Random_ID_post)) & (survey$condition_pre=='control'| (survey$condition_pre!='control' & !is.na(survey$focus_10)))
  survey$dropout <- is.na(survey$Random_ID_post)
  survey$returned <- survey$condition_pre!='control' & is.na(survey$focus_10) & !is.na(survey$Random_ID_post)
  
  #rownames
  rownames(survey) <- survey$Random_ID_pre
  
  return(survey)
}

#dots_data <- cleanSurvey(survey)
#saveRDS(dots_data,'rds_files/dots_data.rds')
#write.csv(dots_data,'csv_files/dots_data.csv')
dots_data <- readRDS('rds_files/dots_data.rds')
#################################################################### Create ----
# Scale Re-numbering ---------------------------
#place on 0-max scale, rather than 1-max
scale_recode <- function(dataset, scale_name, q_nums) {#dataset<-dots_data; scale_name<-"TMS";q_nums<-c(1:7)
  q_nums <- as.character(q_nums)
  q_names <- c(paste0(scale_name,"_",q_nums,"_pre"), paste0(scale_name,"_",q_nums,"_post"))
  dataset[,q_names] <- sapply(dataset[,q_names],as.numeric)
  dataset[,q_names] <- dataset[,q_names] - 1
  return(dataset)
}

dots_data <- scale_recode(dots_data,"TMS",c(1:7)) #  TMS: https://sci-hub.yncjkj.com/10.1891/0889-8391.23.3.185.
dots_data <- scale_recode(dots_data,"GAD",c(1:7))
dots_data <- scale_recode(dots_data,"PHQ",c(1:8)) # PHQ-8: https://asset-pdf.scinapse.io/prod/1985329916/1985329916.pdf

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
  dots_data[,sc_names] <- sapply(dots_data[,sc_names], as.numeric)
}

# Mean Scores and difference over time ----
meanDiffs <- function(dataset){
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

dots_data <- meanDiffs(dots_data)
saveRDS(dots_data,'rds_files/dots_data.rds')
write.csv(dots_data,"csv_files/dots_data.csv")
dots_data <- readRDS('rds_files/dots_data.rds')
# Exclusions --------------------------
excluded <- which(dots_data$dp_pre==0|dots_data$dp_post==0|
                    is.na(dots_data$dp_pre)|is.nan(dots_data$dp_diff))

print(paste0("N excluded on d'=0: ",length(excluded)))
#view: dots_data[,grep('dp|H_M',colnames(dots_data))]
if(length(excluded)>0){
  dots_data <- dots_data[-excluded,]
}
dots_data$condition_pre <- factor(dots_data$condition_pre)
################################################################### Explore ----
# Dropouts -------
dropoutsTab <- function(dataset){ # dataset <- dots_data
  drops <- t(sapply(dataset[,c('per_protocol','returned','dropout')],tapply,dataset$condition_pre,sum))
  drops <- rbind(drops,total=apply(drops,2,sum))
  drops <- cbind(drops,total=apply(drops,1,sum))
  dimnames(drops) <- lapply(dimnames(drops),tools::toTitleCase)
  rownames(drops) <- gsub('_','-',rownames(drops))
  return(drops)
}

dropouts <- dropoutsTab(dots_data)
write.csv(dropouts,'csv_files/dropouts.csv')

# Demographics -------
#Demo
demos <- c("Age mean" = round2(mean(as.numeric(dots_data$age_pre))),
  "Age range" = paste(range(as.numeric(dots_data$age_pre)),collapse=", "),
  "gender" = table(dots_data$gender_pre))

#Time meditating
time_med <- dots_data[!is.na(dots_data$time_meditating_post) & dots_data$time_meditating_post!="" & dots_data$condition_pre!='control',
                      c("condition_pre","time_meditating_post")]
write.csv(time_med,'csv_files/time_med.csv') #could code these?

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

cont <- content(dots_data)

#Do you 
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


usef <- useful(dots_data)

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

auth <- authentic(dots_data)

quant_feedback <- cbind('Content: Exact or More'=cont[,ncol(cont)],
      'Useful: Very'=usef[,ncol(usef)],
      'Authentic: Yes'=auth[,ncol(auth)])
write.csv(quant_feedback,'csv_files/quant_feedback.csv')

#Please let us know if you have any other comments or feedback on the course? + gen qual
qual <- function(dataset){ #dataset<-dots_data
  feedback <- dataset[dataset$condition_pre!='control',c("condition_pre","feedback_post","authentic_post","authentic_3_TEXT_post","useful_post","useful_4_TEXT_post")]
  write.csv(feedback,'csv_files/feedback.csv')
  return(feedback)
}

quali <- qual(dots_data)

# Focus/Activities -------
engagement <- function(dataset){
  #prep
  foc_ac <- dataset[grep('condition_pre|focus|activities|Random_ID_post',colnames(dataset))]
  foc_ac <- foc_ac[foc_ac$condition_pre!='control',]
  foc_ac$condition_pre <- droplevels(factor(foc_ac$condition_pre))
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
  round2(apply(output,2,function(x){x/table(foc_ac$condition_pre)*100}))
}

engage <- engagement(dots_data)

################################################################### Analyse ----
# Multiple Imputation -----
predictorMatrix <- function(dataset){
  # change the matrix of predictor-outcome analyses
  scale_rgx <- paste0('^(',paste(c("SOW",names(scales_list)),collapse='|'),')_[0-9]{1,2}_p(re|ost)$')
  item_names <- colnames(dataset)[grep(scale_rgx,colnames(dataset))]
  single_names <- paste(rep(c(SDT_vars,names(scales_list)),each=2),c('pre','post'),sep='_')
  pred_names <- c('condition_pre',item_names,single_names)
  pred_mat <- matrix(0,length(pred_names),length(pred_names),dimnames=list(pred_names,pred_names)) #make a predictor matrix
  
  post_items <- item_names[grep('_post',item_names)]
  
  #Scale predictors
  for(i in post_items){ # i <-post_items[1] #rows are the things being predicted
    sc <- sub("(.*)_.*_.*","\\1",i)
    items <- colnames(pred_mat)[grep(paste0(sc,"_[0-9]{1,2}_pre"),colnames(pred_mat))]
    subscales <- paste0(names(scales_list)[!grepl(paste0(sc,"|^(WBSI|RRS)$"),names(scales_list))],"_pre")
    pred_mat[i,c(items,subscales,'condition_pre')] <- 1
  }
  
  #STD predictors
  subscales <- names(scales_list)[!grepl("^(WBSI|RRS)$",names(scales_list))]
  for(i in SDT_vars){
    pred_mat[paste0(i,"_post"),paste0(c("condition",i,subscales),"_pre")] <- 1
  }# write.csv(pred_mat,'csv_files/pred_mat.csv')
  
  return(pred_mat)
}

pred_mat <- predictorMatrix(dots_data)

methodsList <- function(dataset,pred_mat){
  scale_rgx <- paste0('^(',paste(c("SOW",names(scales_list)),collapse='|'),')_[0-9]{1,2}_post$')
  post_items <- colnames(dataset)[grep(scale_rgx,colnames(dataset))]
  
  #list of methods
  meths <- rep('',ncol(pred_mat))
  names(meths) <- colnames(pred_mat)
  meths[c(post_items,paste0(SDT_vars,"_post"))] <- 'norm'
  return(meths)
}

meths <- methodsList(dots_data,pred_mat)

#imputation
impute <- function(dataset,meths,pred_mat,dropouts){
  iters <- round((dropouts['Dropout',]/dropouts['Total',])['Total']*100)
  mice(dots_data[colnames(pred_mat)],m=100,maxit=iters,meth=meths,pred=pred_mat,print=F,seed=1)
}

#message('imputing data')
#start_time <- Sys.time()
#imp <- impute(dots_data,meths,pred_mat,dropouts)
#comps <- mice::complete(imp,"all")
#dots_mi <- lapply(comps, meanDiffs)
#print(Sys.time()-start_time)
#imp$loggedEvents
#saveRDS(dots_mi,'rds_files/dots_mi.rds')
dots_mi <- readRDS('rds_files/dots_mi.rds')

# Descriptives ----------
groupMeans <- function(dataset){ # dataset<-dots_data[dots_data$per_protocol,]
  dataset <- dataset[!is.na(dataset$TMS_diff),]
  dataset$condition_pre <- droplevels(dataset$condition_pre)
  all_v <- gsub('_diff','',names(dataset)[grep('_diff',names(dataset))])
  out <- c()
  for(i in c("pre","post")){ # i<-"pre"
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
  colnames(tab) <- paste(rep(c("Waitlist","World",'Mental States'),each=2),c("Pre","Post"),sep=" ")
  tab <- tab[shortlist,]
  rownames(tab) <- pnames[rownames(tab)]
  return(tab)
}

group_means <- groupMeans(dots_data[dots_data$per_protocol,])
write.csv(group_means,'csv_files/group_means.csv')

# d' difference ----
dpDiff <- function(dataset){ # dataset<-dots_data[dots_data$per_protocol,]
  dataset <- data.frame(dataset)
  #dp
  dp_bfs <- round2(rbind(
    dp = tapply(dataset$dp_diff,dataset$condition_pre,mean),
    dp_se = tapply(dataset$dp_diff,dataset$condition_pre,se),
    metad = tapply(dataset$H_meta_d_diff,dataset$condition_pre,mean),
    metad_se = tapply(dataset$H_meta_d_diff,dataset$condition_pre,se)
    ))
  dp_bfs <- dp_bfs[,c('mental','world')]
  #interactions
  dp_bfs <- cbind(dp_bfs,diff=dp_bfs[,'mental']-dp_bfs[,'world'])
  dp_bfs['dp_se','diff'] <- (dp_bfs['dp_se','mental']**2+dp_bfs['dp_se','world']**2)**.5
  dp_bfs <- round2(dp_bfs)

  #bfs
  BF <- bfrr(sample_mean=dp_bfs['dp','diff'],sample_se=dp_bfs['dp_se','diff'],
             sample_df=nrow(dataset)-1,model="normal", mean=0, sd=dp_bfs['metad','diff'], tail=1,
             criterion=3,rr_interval=list(mean=c(0,0),sd=c(0,5)),precision=.01)
  
  #format
  BFs <- round2(c(BF=BF$BF,RR=BF$RR$sd))
  dp_bfs <- format(dp_bfs,nsmall=2)
  dp_bfs <- gsub(' ','',dp_bfs)
  dp_out <- rbind("d'" = paste0(dp_bfs['dp',],' (',dp_bfs['dp_se',],')'),
                  "meta-d'" = paste0(dp_bfs['metad',],' (',dp_bfs['metad_se',],')'))
  colnames(dp_out) <- c('World','Mental States','Difference')
  return(list(dp_meatd=dp_out,BF=BFs))
}

dp_difference <- dpDiff(dots_data[dots_data$per_protocol,])
write.csv(dp_difference$dp_meatd,'csv_files/dp_difference.csv')

# Time taken -------

timeTaken <- function(dataset){ # dataset<-dots_data[dots_data$per_protocol,]
  
  tt_d <- dataset[!is.na(dataset$TMS_diff) & dataset$condition_pre!='control',c('TMS_diff','time_taken','condition_pre')]
  tt_d$condition_pre <-  droplevels(tt_d$condition_pre)
  #H1 calculated = regressing the TMS against time taken to completion
  tt_lm <- round2(summary(lm('TMS_diff~time_taken',tt_d))$coefficients[2,1:2])
  #dividing by the raw difference in TMS scores between interventions by the recorded slope.
  ms <- round2(sapply(tt_d[,c('time_taken','TMS_diff')],tapply,tt_d$condition_pre,mean))
  ses <- round2(sapply(tt_d[,c('time_taken','TMS_diff')],tapply,tt_d$condition_pre,se))
  time_diff <- round2(rbind(m=apply(ms,2,function(x){x['mental']-x['world']}),
                            se=apply(ses**2,2,sum)**.5))
  tt_H1 <- abs(round2(time_diff['m','TMS_diff']/tt_lm['Estimate']))
  #BFs
  BF <- bfrr(sample_mean=abs(time_diff['m',"time_taken"]),sample_se=time_diff['se',"time_taken"],sample_df=nrow(tt_d)-1,
             model="normal", mean=0, sd=abs(tt_H1), tail=1, criterion=3,
             rr_interval=list(mean=c(0,0),sd=c(-20,20)),precision=1)
  
  
  #formatting
  #plot(tt_d$time_taken,tt_d$TMS_diff,main='TMS change by time taken')
  ms <- format(ms,nsmall=2,trim=T)
  ms <- gsub(' ','',ms)
  time_rnd <- paste0(ms,' (',ses,')')
  time_out <- matrix(time_rnd,2,2,dimnames=list(c('Mental States','World'),c('Time Taken','TMS-D')))
  time_diff <- format(time_diff,nsmall=2,trim=T)
  time_diff <- gsub(' ','',time_diff)
  diff_rnd <- paste0(time_diff['m',],' (',time_diff['se',],')')
  time_out <- rbind(time_out,Difference=diff_rnd)
  time_out <- t(time_out)
  #BFs
  lm_out <- round2(c(lm=tt_lm[1], H1=tt_H1, BF=BF$BF))
  RRs <- ifelse(BF$RR$sd==20,'>20',BF$RR$sd)
  RRs <- ifelse(BF$RR$sd==-20,'<-20',RRs)
  RRs[1] <- ifelse(BF$BF<=.33|BF$BF>=.3,RRs[1],'0')
  lm_out <- c(lm_out,RR=paste0('[',paste(RRs,collapse=', '),']'))
  names(lm_out) <- c('slope','H1','BF','RR')
  return(list(desc=time_out,analysis=lm_out))
}

time_taken <- timeTaken(dots_data[dots_data$per_protocol,])
write.csv(time_taken$desc,'csv_files/time_taken.csv')

# Expectations ---------

expectancies <- function(dataset){ # dataset<-dots_data[dots_data$per_protocol,]
  #setup
  exp_data <- dataset[!is.na(dataset$exp_TMS),c('condition_pre','exp_TMS','exp_GAD','exp_PHQ','exp_SOW_1_1','exp_WBSI_1_1','exp_RRS_1_1')]
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
  rownames(exps) <- gsub("exp_([A-Z]{3,4})(_1_1)?","\\1",rownames(exps))
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

expects <- expectancies(dots_data)
write.csv(expects,"csv_files/expectancies.csv")
# Analysis -----------
interactions <- function(dataset){ # dataset <- dots_data
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
  list(n=nrow(diffs),comps=comps) #allows n to be pulled through
}

BFs <- function(dataset){ # dataset <- dots_data
  ints <- interactions(dataset)
  comps <- cbind(H1=NA,round2(ints$comps))
  comps <- comps[!rownames(comps)%in%c('c','c1','c2','ln_M_ratio'),] # c is too small to run - log(M) has NAs
  # Bayes Factors
  for(comp in c('mvw','wvwl','mvwl')){ #comp<-'mvw'
    for(i in rownames(comps)){ #i<-1 
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

#Per-protocol
res_pp <- BFs(dots_data[dots_data$per_protocol,])
#saveRDS(res_pp,'rds_files/res_pp.rds')
res_pp <- readRDS('rds_files/res_pp.rds')

#Multiple Imputation
#res_mi <- lapply(dots_mi,BFs)
#saveRDS(res_mi,'rds_files/res_mi.rds')
res_mi <- readRDS('rds_files/res_mi.rds')
res_red <- round2(Reduce("+", res_mi)/length(res_mi))

################################################################### Present ----
# Format results ------
formatResults <- function(tab){ # tab<-res_pp
  tab_f <- format(tab, nsmall=2, trim=T) # ensure 2 d.p. trailing 0s
  #format rrs
  max_rr <- tab[,grep("rru",colnames(tab))]==5
  tab_f[,colnames(max_rr)][max_rr] <- '>5.00' # too high
  min_rr <- tab[,grep("rrl",colnames(tab))]==.01
  tab_f[,colnames(min_rr)][min_rr] <- '0.00' # too low
  bfs <- tab[,grep("bf",colnames(tab))]
  insen <- bfs>=.33 & bfs<3
  tab_f[,colnames(min_rr)][insen] <- '0.00' # insensitive
  
  ## insensitive
  #bind to dataset
  diffs <- data.frame(
    'H1'= tab_f[,"H1"],
    'MS-World Mean'= paste0(tab_f[,"mvw_m"]," (",tab_f[,"mvw_se"],")"),
    'MS-World BF'= tab_f[,"mvw_bf"],
    'MS-World RR'= paste0("[",tab_f[,"mvw_rrl"],", ",tab_f[,"mvw_rru"],"]"),
    
    'World-Waitlist Mean'= paste0(tab_f[,"wvwl_m"]," (",tab_f[,"wvwl_se"],")"),
    'World-Waitlist BF'= tab_f[,"wvwl_bf"],
    'World-Waitlist RR'= paste0("[",tab_f[,"wvwl_rrl"],", ",tab_f[,"wvwl_rru"],"]"),
    
    'MS-Waitlist Mean'= paste0(tab_f[,"mvwl_m"]," (",tab_f[,"mvwl_se"],")"),
    'MS-Waitlist BF'= tab_f[,"mvwl_bf"],
    'MS-Waitlist RR'= paste0("[",tab_f[,"mvwl_rrl"],", ",tab_f[,"mvwl_rru"],"]")
  )
  diffs <- apply(diffs,2,function(x){paste0(x,'\t')}) #hacky approach to force excel to allow trailing zeroes
  rownames(diffs) <- rownames(tab_f) #get back rownames
  short <- diffs[shortlist_group,]
  rownames(short) <- pnames[rownames(short)]
  return(short)
}

res_pp_f <- formatResults(res_pp)
write.csv(res_pp_f,'csv_files/PP_results.csv')

res_mi_f <- formatResults(res_red)
write.csv(res_mi_f,'csv_files/MI_results.csv')

# Bar graphs -----------
shorts <- unique(c(shortlist,shortlist_group))
colours <- turbo(length(shorts))
names(colours) <- shorts
pconds <- c(control='Waitlist',world='World',mental='Mental States')

# plot diffs
plotDiff <- function(dataset){ # dataset<-dots_data[dots_data$per_protocol,]
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
plotDiff(dots_data[dots_data$per_protocol,])
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
  legend('bottomright',legend=pnames[rownames(m)],col=colours[shortlist_group],pch=15,cex=0.7, inset=c(-0.05, 0))
}

jpeg('graphs/interactions_pp.jpeg', height=150, width=200,units='mm',quality=200,res=300)
plotRes(res_pp,'Time*Condition Interactions (Per-Protocol)')
dev.off()

jpeg('graphs/interactions_mi.jpeg', height=150, width=200,units='mm',quality=200,res=300)
plotRes(res_red,'Time*Condition Interactions (Multiple Imputated)')
dev.off()

####################################################################### Fin ----
# End timer & save -----
save.image(file='MSvW.RData')