# Analysis code for G-factor study by Max Lovell #
# Setup ----
# note any time consuming analyses have their outputs stored in rds files.
#make working dir same as script file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#load environment if code already run
if(file.exists('gfactor.RData')){
  #load('gfactor.RData')
}
start_time <- Sys.time()
# Packages ----
list.of.packages <- c('jsonlite',  # ::fromJSON() for reding in raw data
                      'rjags',     # run Hmeta-d models (JAGS must be installed outside R)
                      'parallel',  # run JAGS in parallel
                      'stableGR',  # ::stable.GR() for improved Rhat
                      'beepr',     # ::beep() for alerts when code is done
                      'readr',     # ::write_excel_csv() for writing Rho character to csv
                      'devtools',  # ::install_github() as bfrr not working
                      'bfrr',      # ::bfrr() Bayes Factors and Robustness regions
                      'psych',     # ::corr.test() for non-H correlations
                      'NlcOptim',  # ::solnl() for Maxmimum Likelihood Estimation
                      'viridis',    # colour palette
                      'uaparserjs' # parse user agents
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
if('bfrr' %in% new.packages){ devtools::install_github("debruine/bfrr") }
lapply(list.of.packages, require, character.only = TRUE)

source("meta-d/fit_meta_d_MLE.R")
source("meta-d/fit_meta_d_mcmc.R")
source("meta-d/fit_meta_d_mcmc_group.R")
source("meta-d/fit_meta_d_mcmc_groupCorr.R")

# Helpers -----
idx <- t(combn(6,2))
# names
pnames <- c("TMS"="TMS-D", "SOW_F"="Observe 1st", "SOW_S"="Observe 2nd",
            "acc"="Accuracy", "dp"="d'", "c"="c", "adj_dp"="Adj. d'",
            "dp2"="Type-2 d'", "c2"="Type-2 c", "adj_dp2"="Adj. type-2 d'",
            "adj_meta_d"="Adj. meta-d'", "H_meta_d"="Hmeta-d", "H_M_ratio"="HM-ratio",
            "ln_M_ratio"="log(M-ratio)")

#Variables for correlations
var_names <- c('M_ratio','M_ratio_rS1','M_ratio_rS2')
tnames <- c("Dots","Gabor YN","Gabor 2AFC","Span YN","Span 2AFC","Breath")
AFC <- c('dots','gabor_2AFC','span_2AFC') #2AFC tasks
YN <- c('gabor_YN','span_YN','breath')
  
#Fisher's r<=>z transforms
r2z <- function(r){ .5*(log((1+r)/(1-r))) }
z2r <- function(z){ (exp(2*z)-1)/(exp(2*z)+1) }

# s-b adjust correlations
sb <- function(x){(2*x)/(1+x)}
# error calculation
se <- function(x,na.rm=F){ 
  if(na.rm){x<-x[!is.na(x)]}
  sd(x)/sqrt(length(x)) 
}
ci <- function(x){ se(x)*qnorm(.975) }

# Round function to stop R even rounding as we only need 2 d.p.
  # https://www.r-bloggers.com/2023/04/rounding-in-r-common-data-wrangling-frustrations-and-workarounds-in-r-julia-and-python
  # https://stackoverflow.com/questions/12688717/round-up-from-5/12688836#12688836
round2 <- function(x, digits = 2) {  # Function to always round 0.5 up
  posneg <- sign(x)
  z <- abs(x) * 10^digits
  z <- z + 0.5
  z <- trunc(z)
  z <- z / 10^digits
  z * posneg
}

###################################################################### Read ----
# Read data ----
# note putty must be installed for pscp to work
# replace square brackets: system('pscp -pw [password] [server space root]:[data folder]/*.json "C:\\Users\\[folder to download to]"')

readTasks <- function(tasks){
  all_files <- list.files("data/") #list downloaded files
  tasks_raw <- list()
  for(t in task_names){ #t<-task_names[4]
    task_files <- all_files[grep(t,all_files)] #list files for each task
    task <- c()
    for(f in task_files){ #f<-task_files[1] ;#loop through files
      file <- fromJSON(paste0("data/",f)) # Read in raw datafile
      if(length(file)==0){print(f) #print names of incorrect files
        file.remove(paste0("data/",f)) #remove incorrect files
      } else {
        file$ID <- sub("\\D+","",f) # get ID
        if(t=='survey'){ file <- data.frame(file) }
        task <- rbind(task,file) # bind subj data to matrix
      } 
    }
    write.csv(data.frame(apply(task,2,as.character)),paste0('csv_files/',t,".csv"))#convert to character so commas are correct in CSV
    tasks_raw <- c(tasks_raw,list(task)) # combine tasks into list
  }
  names(tasks_raw) <- task_names
  return(tasks_raw)
}

task_names <- c("survey","dots","gabor","span","breath")
#tasks_raw <- readTasks(task_names)
#saveRDS(tasks_raw,'rds_files/tasks_raw.rds')
tasks_raw <- readRDS('rds_files/tasks_raw.rds')
# Dropouts ----
all_files <- list.files("data/") #list downloaded files
id_count <- table(gsub('\\D+','',all_files)) #count number of files for each person
#note if(length(file)==0){print(f);next} above removes empty datasets from the folder
valid_ids <- names(id_count[id_count==5]) # participants with full data sets

dropouts <- cbind( #count different number of dropouts
  student_n = sum(nchar(names(id_count))==5),
  known_n = sum(nchar(names(id_count))==6),
  dropout_n = sum(id_count!=5),
  protocol_n = sum(id_count==5),
  dropout_rate = (sum(id_count!=5)/length(id_count))*100
)

# Check User Agents, Debugging ------
#user agents and ids
#srv <- tasks_raw$survey[,c('ID','user_agent')]
#uas <- data.frame(UAs$ID,ua_parse(UAs$user_agent))
#breath count
#br_count <- tapply(tasks_raw$breath,tasks_raw$breath$ID,nrow)
#br_n <- data.frame(br=br_count,ID=names(br_count))
#merge
#ua_br <- merge(uas,br_n,by.x='UAs.ID',by.y='ID')
#ua_br <- ua_br[order(ua_br$br, ua_br$userAgent),]
#write.csv(ua_br,'UAs_br_count.csv')

##################################################################### Clean ----
# Survey ------
cleanSurvey <- function(survey){
  TMS_r <- c('Not at all','A little','Moderately','Quite a bit','Very much')
  SSO_r <- c('Never or very rarely true','Rarely true','Sometimes true','Often true','Very often or always true')
  
  recode <- function(var_name,var_opts, dataset){ # var_name<-"TMS";var_opts<-TMS_r;dataset<-survey
    var_opts <- tolower(var_opts)
    var_opts <- gsub(" ","_",var_opts)
    for(i in 1:length(var_opts)){ #i<-1
      dataset[,grep(var_name,names(dataset))] <- apply(dataset[,grep(var_name,names(dataset))],2,function(x){ gsub(var_opts[i],i,x) })
    }
    dataset[,grep(var_name,names(dataset))] <- apply(dataset[,grep(var_name,names(dataset))],2,as.numeric)
    return(dataset)
  }
  
  survey <- recode("TMS",TMS_r,survey)
  survey <- recode("SSO",SSO_r,survey)
  #survey <- survey[,-which(names(survey)%in%c("user_agent","platform"))]
  survey$TMS <- apply(survey[grep("TMS",names(survey))],1,mean)
  survey$SSO <- apply(survey[grep("SSO",names(survey))],1,mean)
  survey$nationality[survey$nationality==F] <- 'Taiwan'
  survey$platform[survey$platform==F] <- 'sona'
  return(survey)
}

survey <- cleanSurvey(tasks_raw$survey)
# Breath -----
cleanBreath <- function(breath){ #breath<-tasks_raw$breath
  breath[,c("br","conf")] <- apply(breath[,c("br","conf")],2,function(x){ tolower(gsub("Arrow","",x)) })
  breath$conf <- gsub("left",1,breath$conf)
  breath$conf <- gsub("right",2,breath$conf)
  
  ## reset deletes everything up to last down press
  downs <- c(0,which(grepl("down",breath$br)|!duplicated(breath$ID,fromLast=T)))
  resets <- which(grepl("reset",breath$br))
  reset_seq <- unlist(lapply(resets,function(x){
    dlx <- downs[downs<x]
    (dlx[length(dlx)]+1):x
  })) #check: breath$reset <- F; breath[reset_seq,"reset"] <- T
  breath <- breath[-c(reset_seq,which(breath$conf=='no_conf')),]
  seq_ends <- diff(c(0,which(grepl("down",breath$br)|!duplicated(breath$ID,fromLast=T))))
  breath$count <- unlist(lapply(seq_ends,seq))
  breath$sequence <- unlist(tapply(grepl("down",breath$br),breath$ID,function(x){
    x[length(x)]<-T
    rep(1:sum(x),diff(c(0,which(x))))
  }))
  
  breath$target <- ifelse(breath$count==9,'down','up')
  breath$correct <- breath$br==breath$target
  names(breath) <- c('response','resp_rt','confidence','conf_rt','ID','count','sequence','target','correct')
  breath <- breath[,c('ID','count','sequence','target','response','resp_rt','correct','confidence','conf_rt')]
  return(breath)
}

breath <- cleanBreath(tasks_raw$breath)
# Dots, Span, Gabor; Combine -----
createTaskList <- function(tasks_raw){
  #extract for ease
  dots <- tasks_raw$dots
  gabor <- tasks_raw$gabor
  span <- tasks_raw$span
  
  #intensity
  dots$intensity <- ((dots$n_dots - 0.1) * 100) / (5.743 - 0.1)
  gabor$contrast <- as.numeric(gabor$max)-as.numeric(gabor$min)
  gabor$intensity <- ((gabor$contrast - 2) * 100) / (255 - 2)
  span$span <- as.numeric(span$span)
  span_max <- 30 #max(span$span)==28
  span$intensity <- ((((span_max+1)-span$span) - 1) * 100) / (span_max - 1)
  
  #separate 2AFC and YN tasks
  gabor_YN <- gabor[gabor$trial_n<=70,]
  gabor_2AFC <- gabor[gabor$trial_n>70,]
  gabor_2AFC$trial_n <- gabor_2AFC$trial_n-70
  span_YN <- span[span$trial_n<=70,]
  span_2AFC <- span[span$trial_n>70,]
  span_2AFC$trial_n <- span_2AFC$trial_n-70
  
  #list
  tasks <- list(dots=dots,gabor_YN=gabor_YN,gabor_2AFC=gabor_2AFC,
                span_YN=span_YN,span_2AFC=span_2AFC,breath=breath)
  #lapply(tasks,function(x){unique(x$ID)})
  #select valid ids
  lapply(tasks,function(x){ 
    x <- x[x$ID%in%valid_ids,]
    x$ID <- as.factor(x$ID);x})
}

#tasks <- createTaskList(tasks_raw)
#saveRDS(tasks,'rds_files/tasks.rds')
tasks <- readRDS('rds_files/tasks.rds')
############################################################# SDT functions ----
# NR_SX: Stimulus-response-confidence counts ----
NR_SX <- function(task,s1_s2=unique(task$response)){ # task<-split_data[[1,6]];s1_s2=s_rs[[6]]
  S1<-s1_s2[1];S2<-s1_s2[2]
  task$ID <- droplevels(task$ID)
  task <- task[task$confidence!=-1,] #remove practice trials
  task$confidence <- as.factor(task$confidence)
  #counts of confidece ratings for each stimulus-response pairing
  S1_R1 <- table(task[task$target==S1 & task$response==S1,c("ID","confidence")]) #note ID must be a factor, confidence too?
  S1_R2 <- table(task[task$target==S1 & task$response==S2,c("ID","confidence")])
  S2_R1 <- table(task[task$target==S2 & task$response==S1,c("ID","confidence")])
  S2_R2 <- table(task[task$target==S2 & task$response==S2,c("ID","confidence")])
  S1_R1 <- S1_R1[,ncol(S1_R1):1] #reverse
  S2_R1 <- S2_R1[,ncol(S2_R1):1]
  
  #bind as seaprate counts for each stimulus
  nr_s1 <- cbind(S1_R1,S1_R2) #bind
  nr_s2 <- cbind(S2_R1,S2_R2)
  colnames(nr_s1) <- paste0(rep(c("S1_R1_C","S1_R2_C"),each=ncol(S1_R1)), colnames(nr_s1))
  colnames(nr_s2) <- paste0(rep(c("S2_R1_C","S2_R2_C"),each=ncol(S1_R1)), colnames(nr_s2))
  nr_s1 <- data.frame(nr_s1) #maybe take this out?
  nr_s2 <- data.frame(nr_s2)
  
  #list for use in SDT calculations
  nr_sx <- list(nr_s1=nr_s1,nr_s2=nr_s2)
  return(nr_sx)
}

# Individual-level SDT ----

## TYPE-1 ##
type1 <- function(nr_sx.5,c_levs){
  #counts
  t1 <- data.frame(sapply(1:4,function(x){
    apply(nr_sx.5[(c_levs*(x-1)+1):(c_levs*x)],1,sum)
  }))+.5
  names(t1) <- c("CR","FA","miss","hit")
  #rates
  rates <- cbind(HR=t1$hit/(t1$hit+t1$miss), FAR=t1$FA/(t1$FA+t1$CR))
  rates <- qnorm(rates)
  #stats
  d1 <- rates[,1]-rates[,2]
  c1 <- -.5*apply(rates,1,sum)
  #bind
  data.frame(d1,c1)
}


## TYPE-2 ##
type2 <- function(nr_sx.5,c_levs){
  #sums
  C <- apply(nr_sx.5[,c(1:c_levs,(c_levs*3+1):(c_levs*4))],1,sum) #correct
  I <- apply(nr_sx.5[,c((c_levs+1):(c_levs*3))],1,sum) #incorrect
  #rates
  t2 <- nr_sx.5[,F]
  for(x in 2:c_levs){
    t2[,paste0("Hc",x)] <- apply(nr_sx.5[,c(1:(c_levs+1-x),(c_levs*3+x):(c_levs*4))],1,sum)/C
    t2[,paste0("Fc",x)] <- apply(nr_sx.5[,c((c_levs+x):(c_levs*2),(c_levs*2+1):(c_levs*3+1-x))],1,sum)/I
  }
  t2_z <- apply(t2,2,qnorm)
  # formula switch depending on if conf_l>2 for regression model
  if(c_levs>2){
    d2 <- apply(t2_z,1,function(x){
      t2_lm <- lm(x[grep("H",names(x))]~x[grep("F",names(x))]) #;plot(x[grep("FA",names(x))],x[grep("H",names(x))]); abline(t2_lm)
      return(sqrt(2/(1+(t2_lm$coefficients[2]^2)))*t2_lm$coefficients[1]) #see p. 62-63: Hautus, Macmillan & Creelman (2021) Detection Theory: A User's Guide 3rd ed.
    })
    c2 <- -.5*(t2_z[,paste0("Hc",ceiling(c_levs/2))]-t2_z[,paste0("Fc",ceiling(c_levs/2))])
  } else {
    d2 <- t2_z[,"Hc2"]+t2_z[,"Fc2"]
    c2 <- -.5*(t2_z[,"Hc2"]-t2_z[,"Fc2"])
  }
  
  data.frame(d2,c2)
}


## Meta-d' ##
Hmetad <- function(nr_sx,c_levs,est_dp){
  SDT <- nr_sx[,F] #empty dataframe preserving rows
  
  mcmc_params <- list(response_conditional=0,estimate_dprime=est_dp,nchains=3,
                      nadapt=1000,nburnin=1000,nsamples=10000,nthin=1,
                      dic=0,rhat=0,parallel=1,monitorparams=0,saveallsamples=0,
                      estimate_mratio=0)
  if(est_dp){
    SDT[c("fit_d1","fit_c1","meta_d","M_ratio","fit_d1_rs","fit_c1_rs",'M_ratio_rS1','M_ratio_rS2')]<-NA
  } else {
    SDT[c("meta_d","M_ratio",'M_ratio_rS1','M_ratio_rS2')]<-NA
  }
  for(i in 1:nrow(nr_sx)){ #i<-1
    ## RESPONSE GENERAL
    message(paste0('Running participant ',i,' of ',nrow(nr_sx)))
    mcmc_params$response_conditional<-0
    fit <- fit_meta_d_mcmc(nr_sx[i,1:(c_levs*2)],nr_sx[i,(c_levs*2+1):(c_levs*4)],mcmc_params)
    if(est_dp){   
      SDT[i,c("fit_d1","fit_c1","meta_d","M_ratio")] <- fit[c("d1","c1","meta_d","M_ratio")]
    } else {
      SDT[i,c("meta_d","M_ratio")] <- fit[c("meta_d","M_ratio")]
    }
    
    ## RESPONSE SPECIFIC
    mcmc_params$response_conditional<-1
    fit <- fit_meta_d_mcmc(nr_sx[i,1:(c_levs*2)],nr_sx[i,(c_levs*2+1):(c_levs*4)],mcmc_params)    #nR_S1=nR_S1[s,] ; nR_S2=nR_S2[s,]
    if(est_dp){   
      SDT[i,c("fit_d1_rs","fit_c1_rs",'M_ratio_rS1','M_ratio_rS2')] <- fit[c("d1","c1",'M_ratio_rS1','M_ratio_rS2')]
    } else {
      SDT[i,c('M_ratio_rS1','M_ratio_rS2')] <- unlist(fit[c('M_ratio_rS1','M_ratio_rS2')])
    }
  }
  
  return(SDT)
}

MLmetad <- function(nr_sx,c_levs,est_dp){
  SDT <- nr_sx[,F] #empty dataframe preserving rows
  for(i in 1:nrow(nr_sx)){ #i<-1
    fit <- fit_meta_d_MLE(nr_sx[i,1:(c_levs*2)],nr_sx[i,(c_levs*2+1):(c_levs*4)])
    SDT[i,c("ML_meta_d","ML_M_ratio")] <- fit[c("meta_da","M_ratio")]
  }
  return(SDT)
}

## Call SDT functions ##
SDT <- function(nr_sx,est_dp=0){ # nr_sx<-lapply(nrsxs[,1],function(x){x[1:3,]});nr_sx<-nrsxs[,1]
  nr_sx <- cbind(nr_sx$nr_s1,nr_sx$nr_s2) #bind
  nr_sx.5 <- nr_sx+.5 #adjusted version
  c_levs <- ncol(nr_sx)/4 #n confience levels
  t1 <- type1(nr_sx,c_levs)
  t2 <- type2(nr_sx.5,c_levs)
  md <- Hmetad(nr_sx,c_levs,est_dp)
  data.frame(t1,t2,md)
}

# Reliability analysis -----

#extract odd-even split halves
splitData <- function(dataset){ # dataset <- c_tasks1
  if("breath" %in% names(dataset)){ #do breath sequences rather than breaths
    names(dataset$breath)[names(dataset$breath)=='sequence'] <- 'trial_n'
  }
  split_data <- sapply(dataset,function(x){ list( #x<-dataset[[6]]
    even=x[!x$trial_n%%2,],
    odd=x[!!x$trial_n%%2,]
  )})
  split_nrsx <- apply(split_data,1,function(x){ mapply(NR_SX,x,s_rs[names(dataset)]) })
  do.call(rbind,split_nrsx)
}

## Hierarchical reliability ##
relibility <- function(SX,split_dat,ntask=3){ # SX<-c(0,0);split_dat<-c_split[,6,drop=F];x<-1
  mcmc_params <- list(response_conditional=2,rc_alternating=c(SX,SX),estimate_dprime=0,
                      nchains=ntask,nadapt=ntask*1000,nburnin=ntask*1000,nsamples=ntask*10000,
                      nthin=ntask,dic=0,rhat=1,parallel=1,monitorparams=0,saveallsamples=0)
  if(length(colnames(split_dat))>1){
    lapply(colnames(split_dat),function(x){
      fit_meta_d_mcmc_groupCorr(split_dat[c(1,3),x],split_dat[c(2,4),x],mcmc_params)
    })
  } else {
    fit_meta_d_mcmc_groupCorr(split_dat[c(1,3),],split_dat[c(2,4),],mcmc_params)
  }
}

########################################################## Individual-level ----
# Inidividual-level SDT analysis ----------------
#order of S1 and S2 responses
s_rs <- list(dots=c('left','right'),
             gabor_YN=c("absent","present"),
             gabor_2AFC=c("present","absent"), # gabor 2AFC time 1==present
             span_YN=c("absent","present"),
             span_2AFC=c("present","absent"), # span 2AFC left==present
             breath=c("up","down"))
#calculate nR
nrsxs <- mapply(NR_SX,tasks,s_rs)
saveRDS(nrsxs,'rds_files/nrsxs.rds')
#nrsxs <- readRDS('rds_files/nrsxs.rds')

#SDTs <- apply(nrsxs,2,SDT)
#saveRDS(SDTs,'rds_files/sdt.rds')
SDTs <- readRDS('rds_files/sdt.rds')

# Summarise Responses/Performance -----
respSum <- function(dataset){#dataset<-gfac[[6]]; # dataset<-tasks[[6]]
  t <- dataset #short name
  t <- t[t$confidence!=-1,] #removes practice trials
  t$ID <- droplevels(t$ID)
  t_sum <- data.frame(ID = unique(t$ID))[,F]
  rownames(t_sum) <- unique(t$ID) 
  perf_name <- colnames(t)[colnames(t)%in%c('n_dots','contrast','span','count')]
  
  if("count" %in% names(t)){  #breath task
    br_raw <- tasks_raw$breath[tasks_raw$breath$ID%in%valid_ids,]
    t_sum[,c('n_Reset','n_Breath')] <- data.frame(sapply(c('reset','ArrowUp|ArrowDown'), function(x){ 
      tapply(grepl(x,br_raw$br),br_raw$ID,sum) 
    }))
    trial_n <- t_sum$n_Breath
  } else { 
    t_sum[,"intensity_M"] <- tapply(t$intensity,t$ID,mean)
    t_sum[,"intensity_SD"] <- tapply(t$intensity,t$ID,sd)
    t_sum[,"stim_min"] <- tapply(t[,perf_name],t$ID,min)
    t_sum[,"stim_mean"] <- tapply(t[,perf_name],t$ID,mean)
    trial_n <- sum(t$ID==t$ID[1])
  }
  
  t_sum[,"stim_max"] <- tapply(t[,perf_name],t$ID,max)
  t_sum[,"accuracy_pctg"] <- (tapply(t$correct,t$ID,sum)/trial_n)*100
  t_sum[,"response_pctg"] <- (tapply(t$response==t$response[1],t$ID,sum)/trial_n)*100
  t_sum[,"confidence_SD"] <- tapply(t$confidence,t$ID,sd)
  
  return(t_sum)
}

resps <- lapply(tasks,respSum) #resps <- do.call("cbind", resps)

################################################################ Exclusions ----
# Plot Descriptives ----
plotStat <- function(dataset,stat=NA,save=F,name=stat,filename=name,bins=50,txt_sz=1){ #dataset=SDTs
  if(save==T){ 
    jpeg(file=paste0("graphs/",filename,".jpeg"),900,600,quality=100) 
    txt_sz=2
  }
  print(class(dataset)[1])
  if(class(dataset)[1]=='list'){
    if(!stat %in% names(dataset[[6]])){
      dataset<-dataset[1:5]
    }
    dataset <- sapply(dataset,function(x){x[,stat]})
    if(save==T){
      par(mfrow=c(2,3), oma=c(1,1,5,1),mar=c(2,2,3,2))
    } else {
      par(mfrow=c(3,2), oma=c(0,0,2,0),mar=rep(2,4))
    }
  }
  sapply(colnames(dataset),function(x){
    data <- round(dataset[,x],2)
    name <- paste(toupper(substr(x,1,1)), substr(x,2,nchar(x)),sep="")
    title <- paste0(gsub("_"," ",name)," {",min(data),", ",max(data),"}")
    hist(dataset[,x],bins,main=title,cex.axis=txt_sz,cex.main=txt_sz,col='lightskyblue')}) #
  name <- paste0(name," (N=",NROW(dataset),")")
  mtext(name, line=txt_sz-.5, side=3, outer=TRUE, cex=txt_sz,font=2)
  cormat <- round(cor(dataset),2)
  print(stat)
  print(cormat)
  if(save==T){ 
    dev.off()
    write.csv(cormat,paste0("graphs/",filename,".csv"))
  }
}
plotStat(SDTs,'M_ratio')
pp_plot <- plotStat(SDTs,'M_ratio',save=T,name="Per-protocol M-ratios",filename='Hist_mratios_pp')

#SDTz <- apply(sapply(SDTs,function(x){x[,'M_ratio']}),2,function(x){gmean<-mean(x); sds<-sd(x); (x-gmean)/sds})
#plotStat(SDTz,'M_ratio') #in SD units
#invisible(lapply(names(resps[[1]]),function(x){plotStat(resps,x)})) #performance estimates
#invisible(lapply(names(SDTs[[1]]),function(x){plotStat(SDTs,x)})) #all SDT vars
#plotStat(resps[[6]],'Breath Task Results') #save=T) #all breath vars

# Plot Exclusions ----
excl <- function(dataset,var,lb,ub,vec=F,SD=F){ # dataset<-SDTs;var='M_ratio';lb<- -3;ub<-3;SD=T
  var_data <- sapply(dataset,function(x){x[,var]})
  apply(var_data>5,2,which)
  if(missing(lb)){lb<- min(var_data)}
  if(missing(ub)){ub<- max(var_data)}
  if(SD==T){
    var_dataz <- apply(var_data,2,function(x){gmean<-mean(x); sds<-sd(x); (x-gmean)/sds})
    if(ub==max(var_data)){ ub=max(var_dataz) }
    if(lb==min(var_data)){ lb=min(var_dataz) }
    var_data <- var_dataz
  }
  exclusions <- apply(var_data>= lb & var_data<=ub,1,sum)<NCOL(var_data)
  print(sum(exclusions))
  if(vec==T){ return(exclusions) #TRUE= doesn't meet criteria
  } else { 
    dataset <- lapply(dataset,function(x){x[!exclusions,]})
    return(dataset)
  }
}

#Guggenmos example
guggenmos <- excl(SDTs,'M_ratio',0,1.6,vec=T)

### d'<5 ###
mind1 <- round2(qnorm(.6)-qnorm(.4),1) # d' cutoff
#plot
excl_d <- excl(SDTs,'d1',mind1,vec=T) # vector of exclusions
c_SDT_dp <- lapply(SDTs,function(x){x[excl_d,]}) 
plotStat(c_SDT_dp,'M_ratio',save=T,name="M-ratios (d'<.5)",filename='Hist_mratios_excl_d.5')

# Run Exclusions -----------
#dp
min_d <- round2(qnorm(.6)-qnorm(.4),1) # d' cutoff
dp_excl <- apply(sapply(SDTs,function(x){ x$d1>=min_d }),1,sum)==6 #vector of participants to keep
sum(!dp_excl) #number exlcuded on d'
minmax_mrat_dp_excl <- round2(apply(sapply(SDTs,function(x){x[dp_excl,"M_ratio"]}),2,range))
#breath
br_crit <- c(12-5,20+5)*15 #+/-5bmp buffer on 12-20bmp over 15 mins (Sapra et al., 2022)
cl_br <- tasks$breath[tasks$breath$ID%in%rownames(SDTs$breath)[dp_excl],] #remove dp exclusions first
n_br <- tapply(cl_br,cl_br$ID,nrow) #count row numbers
br_excl <- n_br>br_crit[1] & n_br<br_crit[2]
sum(!br_excl) #number excluded on breath count
final_ids <- names(n_br[br_excl]) # get final valid IDs
length(final_ids) #final sample size
saveRDS(final_ids,'rds_files/final_ids.rds')
final_ids <- readRDS('rds_files/final_ids.rds')
# clean everything
c_resp <- lapply(resps,function(x){x[rownames(x) %in% final_ids,]})
c_SDTs <- lapply(SDTs,function(x){x[rownames(x) %in% final_ids,]})
c_tasks <- lapply(tasks,function(x){
  x <- x[x$ID %in% final_ids,]
  x$ID <- droplevels(x$ID);x
})
c_nrsxs <- mapply(NR_SX,c_tasks,s_rs)
saveRDS(c_nrsxs,'rds_files/c_nrsxs.rds')
c_nrsxs <- readRDS('rds_files/c_nrsxs.rds')
c_split <- splitData(c_tasks)
#plot final dataset
plotStat(c_SDTs,'M_ratio')
final_hist <- plotStat(c_SDTs,'M_ratio',save=T,name="Histograms of M-ratios after exclusions",filename='Hist_Mratios_excl_br')

#sapply(c_resp,function(x){x[x$response_pctg<10|x$response_pctg>90,]})
#sapply(c_resp,function(x){x[x$accuracy_pctg<=55,]})
#c_SDTs$dots[rownames(c_SDTs$dots)==401231,]

############################################################### Staircasing ----
# Staircasing plots ----

#flexible plotting function
plotTrials <- function(dataset,stat="intensity",colour=1,analysis=mean,title,ylims,plotter,cut_t=10,S2,leg_loc='topright',labs=F,tasklist=clean){ # dataset<-dots;cut_t<-10;analysis<-mean;stat='response';ylims=c(0,1)
  #see https://www.science.org/doi/10.1126/science.1191883#supplementary-materials
  #replace missing args
  if(missing(title)){
    title <- paste0(match.call()$dataset)
    if(deparse(substitute(analysis)=='sd')){ title <- paste0(stat," ",deparse(substitute(analysis))) }
  }
  title <- paste0(title," (+/-1SE)")
  if(missing(plotter)){
    if(match.call()$dataset=='tasks$dots'||match.call()$dataset=='clean$dots'||match.call()$dataset=='dots'||(is.numeric(dataset)&&dataset==1)){ 
      plotter = plot
      S1='left'
    } else { plotter = lines }
  }
  #flexible to lapply numbers
  if(is.numeric(dataset)||match.call()$dataset=='x'){
    #if(dataset==5){ legend(leg_loc,col=1:5,lwd=2,legend=names(tasks)[1:5],cex=.7) }
    S2=s_rs[[dataset]][[2]]
    dataset <- tasklist[[dataset]]
  }
  
  #prep data
  dataset <- dataset[dataset$ID%in%valid_ids,]
  n <- length(unique(dataset$ID))
  if(min(dataset$trial_n)==71){ dataset$trial_n <- dataset$trial_n-70 }
  dataset$trial_n <- cut_t*(dataset$trial_n%/%cut_t + as.logical(dataset$trial_n%%cut_t))
  #stat dependent changes
  name <- stat
  if(stat=='confidence'){ dataset <- dataset[dataset$trial_n>20,] 
  } else if(stat=="intensity_ou"){ #intensity original units
    name <- colnames(dataset)[colnames(dataset)%in%c("n_dots","contrast","span")]
  }
  if(stat=="response"){
    data <- dataset[,name]==S2
  } else { data <- dataset[,name] }
  
  #analyse
  pmeans <- tapply(data,list(dataset$ID, dataset$trial_n),analysis)
  means <- apply(pmeans,2,mean)
  ses <- apply(pmeans,2,sd)/sqrt(nrow(pmeans))
  trial_num <- as.numeric(names(means))
  ylab <- paste0(stat," (",deparse(substitute(analysis)),")")
  if(missing(ylims)){ ylims <- round(c(min(means-ses),max(means+ses)),3) }
  
  #plot
  plotter(trial_num,means,type="l",ylim=ylims,col=colour,lwd=4,ann=labs,ylab=ylab,xlab='Trial Number',main=title,cex.axis=2.5) #
  arrows(trial_num, means-ses, trial_num, means+ses, length=0.05, angle=90, code=3,col=colour)
}


# call plotting function on final vars
staircasePlot <- function(dataset){ #dataset<-clean
  txt_sz <- 1.7
  if(length(dev.list()!=0)){dev.off()}
  jpeg("graphs/staircasing_trials.jpeg", height=250, width=500,units='mm',quality=200,res=300)
  par(mfrow=c(2,4), oma=c(1,3,11,0),xpd=NA,mar=c(2,4,2,2))
  #Prop Correct
  invisible(lapply(1:(length(dataset)-1),function(x){plotTrials(x,'correct',x,mean,'Prop. Correct Mean',c(.65,.85),labs=F,tasklist=dataset)}))
  mtext('Mean', line=1.5, side=3,cex=txt_sz,font=2)
  mtext('Prop. Correct', line=4, side=2,cex=txt_sz,font=2)
  invisible(lapply(1:(length(dataset)-1),function(x){plotTrials(x,'correct',x,sd,'Prop. Correct SD',c(.34,.48),leg_loc='bottomright',tasklist=dataset)}))
  mtext('SD', line=1.5, side=3,cex=txt_sz,font=2)
  #Confidence
  invisible(lapply(1:(length(dataset)-1),function(x){plotTrials(x,'confidence',x,mean,'Confidence Mean',c(2.7,3.7),tasklist=dataset)}))
  mtext('Confidence', line=3.5, side=2,cex=txt_sz,font=2)
  mtext('Mean', line=1.5, side=3,cex=txt_sz,font=2)
  invisible(lapply(1:(length(dataset)-1),function(x){plotTrials(x,'confidence',x,sd,'Confidence SD',c(.7,1.2),tasklist=dataset)}))
  mtext('SD', line=1.5, side=3,cex=txt_sz,font=2)
  #Response
  invisible(lapply(1:(length(dataset)-1),function(x){plotTrials(x,'response',x,mean,'Prop. S2 Response Mean',c(.35,.65),tasklist=dataset)}))
  mtext('Prop. RS2', line=4, side=2,cex=txt_sz,font=2)
  invisible(lapply(1:(length(dataset)-1),function(x){plotTrials(x,'response',x,sd,'Prop. S2 Response SD',c(.44,.51),tasklist=dataset)}))
  #Intensity
  invisible(lapply(1:(length(dataset)-1),function(x){plotTrials(x,'intensity',x,mean,'Intensity % Mean',c(30,80),tasklist=dataset)}))
  mtext('Intensity %', line=3.5, side=2,cex=txt_sz,font=2)
  invisible(lapply(1:(length(dataset)-1),function(x){plotTrials(x,'intensity',x,sd,'Intensity % SD',c(0,15),tasklist=dataset)}))
  #Title and legend
  mtext('Staircasing', line=7, side=3, outer=TRUE, cex=2,font=2)
  par(fig=c(0,1,0,1),oma=c(1,3,0,0),mar=c(0,0,5,0),new=T)
  plot(0,0,type='l',bty='n',xaxt='n',yaxt='n')
  legend(x='top',col=1:5,lwd=6,legend=tnames[1:5],cex=2,horiz=T,xpd=T)#legend=names(tasks)[1:5]
  #segments(x0=-.03,y0=-1.1,x1=-.03,y1=.9,lwd=2)
  dev.off()
}

staircasePlot(c_tasks)


#instensity in original units
intensityPlots <- function(dataset,stat=mean,ylim_gab=c(80,170),ylim_span=c(7.5,12)){ #plot intensity in original units ;dataset<-clean
  plotTrials(dataset$dots,'intensity_ou',1,stat,'Dot exp() mean',labs=T,plotter=plot)
  plotTrials(dataset$gabor_YN,'intensity_ou',2,stat,'Gabor Contrast Mean',ylim_gab,labs=T,plotter=plot)
  plotTrials(dataset$gabor_2AFC,'intensity_ou',3,stat,plotter=lines) #8-10.5
  legend("topright",col=2:3,lwd=2,legend=names(tasks)[2:3],cex=.7)
  plotTrials(dataset$span_YN,'intensity_ou',4,stat,'Span Length Mean',ylim_span,labs=T,plotter=plot)
  plotTrials(dataset$span_2AFC,'intensity_ou',5,stat,plotter=lines) #8-10.5
  legend("topleft",col=4:5,lwd=2,legend=names(tasks)[4:5],cex=.7)
}
#intensityPlots(c_tasks,sd,c(2,40),c(.8,1.9))
#intensityPlots(c_tasks,mean)

#table(c_tasks$gabor_YN$target,c_tasks$gabor_YN$correct)

# Trial-cuts on heirarhcical reliabilities ----

#Test effects of trial cut on full M-ratio heirarchical reliabilities
Hrels <- function(dataset,cut_seq=seq(20,40,10),cols=F,itr=3){
  h_rels <- lapply(cut_seq,function(t){ #t=20
    tcut <- lapply(dataset[1:5],function(x){x[x$trial_n>t,]}) #cut earlier trials
    nr_cut <- splitData(tcut) # odd-even split half nrsx counts
    if(!cols){ relibility(0,nr_cut,itr) #cols are the task number?? can't remember....
    } else { relibility(0,nr_cut[,cols,drop=F],itr) }
  })
  names(h_rels) <- cut_seq # names from trail lower limit
  h_rels <- lapply(h_rels,function(x){ names(x)<-names(dataset)[1:5]; x }) #rename cols
  return(h_rels)
}

#h_rels <- Hrels(c_tasks)
#saveRDS(h_rels,'rds_files/heirarchical_reliabilities.rds')
h_rels <- readRDS('rds_files/heirarchical_reliabilities.rds')

# Plot heirarchical reliabilities
#level 1 is starting trial number, level 2 is reliabilities
plotRelsCut <- function(reliabilities,tseq=c(20,30,40),name='rels_trial_cuts'){ # reliabilities <- h_relNEW
  trels <- sapply(reliabilities,function(x){sapply(x,function(y){ y$rho$list })})
  conv <- sapply(reliabilities,function(x){sapply(x,function(y){ y$mcmc$rhat$full$mpsrf })})
  dimnames(trels) <- list(tnames[1:5],tseq)
  print(round2(trels))
  jpeg(file=paste0("graphs/",name,".jpeg"), height=80, width=150,units='mm',quality=200,res=200)
  par(mfrow=c(1,2),mar=rep(1.5,4),oma=c(1,1,0,1)) # ,mgp=c(0,0,0)), oma=c(1,1,2,1)
  plot(tseq,trels[1,],cex.lab=1.1,cex.axis=1.1,type="b",pch=19,lwd=2,xaxt='n',ylim=c(0,1),ylab='Reliability',xlab='Starting Trial',main='Reliability')
  axis(1, at=tseq,cex=2)
  lapply(2:5,function(x){lines(tseq,trels[x,],type="b",pch=19,lwd=2,col=x)})
  #legend('topleft',legend=rownames(trels),col=1:5,lwd=5,horiz=F,cex=.7,text.font=2,y.intersp=1)
  plot(tseq,conv[1,],cex.lab=1.1,cex.axis=1.1,type="b",pch=19,lwd=2,xaxt='n',ylim=c(1,10),ylab='Rhat',xlab='Starting Trial',main='Rhat')
  axis(1, at=tseq,cex=2)
  abline(h=1.1, xpd=F, lwd=1.5,lty='dotted')
  lapply(2:5,function(x){lines(tseq,conv[x,],type="b",pch=19,lwd=2,col=x)})
  legend('topleft',legend=rownames(trels),col=1:5,pch=16,horiz=F,cex=.7,text.font=2,y.intersp=1)
  #mtext("Reliability over starting trials", side=3, line=0, outer=TRUE,cex=1.5,font=2)
  dev.off()
}
plotRelsCut(h_rels,c(0,10,20),'H_rels')

# Trial-cuts on individual-level reliabilities ----

# get SDT analysis on split halves for various trial cuts
splitSDT <- function(dataset,cut_seq=seq(20,40,10)){  # dataset<-c_tasks[6]
  if("breath" %in% names(dataset)){ #do breath sequences rather than breaths
    names(dataset$breath)[names(dataset$breath)=='sequence'] <- 'trial_n' #column standardisation
  }
  cut_rels <- lapply(cut_seq,function(x){ # x<-30
    cutt <- lapply(dataset,function(y){y[y$trial_n>x,]})
    splitt <- splitData(cutt) #odd even split halves
    e <- apply(splitt[1:2,],2,SDT) #get SDT on split cols
    o <- apply(splitt[3:4,],2,SDT)
    list(even=e,odd=o)
  })
  names(cut_rels) <- cut_seq
  return(cut_rels)
}

#indiv_rels <- splitSDT(tasks[1:5]) #do breath later as it won't hold up to trial cuts
#saveRDS(indiv_rels,'rds_files/indiv_rels.rds')
indiv_rels <- readRDS('rds_files/indiv_rels.rds')

# Correlate indiv split-halves
indivRels <- function(split_SDT){ # split_SDT <- indiv_rels
  SDT_vars <- c("d1","c1","d2","c2","meta_d","M_ratio","M_ratio_rS1","M_ratio_rS2")
  tvr <- lapply(SDT_vars,function(x){ # x<-SDT_vars[[1]]
    var_rels <- sapply(split_SDT,function(y){ # y<-split_SDT[[1]]
      even <- lapply(y$even,function(z){z[rownames(z)%in%final_ids,x]}) #lists for mapply
      odd <- lapply(y$odd,function(z){z[rownames(z)%in%final_ids,x]})
      rels <- mapply(corr.test,even,odd)
      do.call(rbind,rels['r',])
    })
    rownames(var_rels) <- tnames[1:nrow(var_rels)]
    return(var_rels)
  })
  names(tvr) <- SDT_vars
  
  return(tvr)
}

rels_p <- indivRels(indiv_rels) #individual reliabilities
rels_r <- lapply(rels_p,round2) #round
rels_b <- cbind(rels_r$d1,rels_r$meta_d,rels_r$M_ratio) # exctract cols
write.csv(rels_b,'csv_files/rels_mat.csv')

# Plit individual reliabilities
plotIndivRels <- function(got_rels){
  tcuts <- as.numeric(colnames(got_rels$d1))-20
  plot(tcuts,got_rels$d1[1,],type="l",pch=1,lwd=2,ylim=c(-.4,.8),xaxt='n',lty=3)
  axis(1, at=tcuts,cex=2)
  lapply(2:5,function(x){ lines(tcuts,got_rels$d1[x,],lwd=2,col=x,lty=3) })
  lapply(1:5,function(x){ lines(tcuts,got_rels$meta_d[x,],lwd=2,col=x,lty=5) })
  lapply(1:5,function(x){ lines(tcuts,got_rels$M_ratio[x,],lwd=2,col=x,lty=1) })
}
plotIndivRels(rels_p)
  
# Add Breath task to reliabilities ----
#add breath task
#split_SDT <- lapply(indiv_rels$'20',lapply,function(x){ x[rownames(x)%in%final_ids,] })
#split_SDT$even$breath <- SDT(c_split[1:2,6])
#split_SDT$odd$breath <- SDT(c_split[3:4,6])
#saveRDS(split_SDT,'rds_files/split_SDT.rds')
split_SDT <- readRDS('rds_files/split_SDT.rds')
rel_cors <- mapply(corr.test,split_SDT$even,split_SDT$odd)
i_rels <- sapply(rel_cors['r',],diag)
############################################################### Group-level ----
# Group-Level Hmeta-d ------
c_nrsxs <- readRDS('rds_files/c_nrsxs.rds')
#setup
mcmc_params <- list(response_conditional=0,estimate_dprime=0,rhat=1,parallel=1,
                    nchains=3,nadapt=1000,nburnin=1000,nsamples=10000,nthin=1)
# Response-general
   #group_nrc <- apply(c_nrsxs,2,function(x){ fit_meta_d_mcmc_group(x[[1]],x[[2]],mcmc_params) })
   #saveRDS(group_nrc,'rds_files/group_nrc.rds')
group_nrc <- readRDS('rds_files/group_nrc.rds')

# Response-Specific
  mcmc_params$response_conditional <- 1
  #group_rc <- apply(c_nrsxs,2,function(x){ fit_meta_d_mcmc_group(x[[1]],x[[2]],mcmc_params) })
  #saveRDS(group_rc,'rds_files/group_rc.rds')
group_rc <- readRDS('rds_files/group_rc.rds')

# Correct HDIs: originally forgot to exp() the group-level results.
HDIs_nrc <- lapply(group_nrc,function(x){HPDinterval(as.mcmc(exp(x$samples)))})
HDIs_rc <- lapply(group_rc,function(x){HPDinterval(as.mcmc(exp(x$samples)))})
HDIs <- Map(rbind,HDIs_nrc,HDIs_rc)

#correlate c and responses
c_resp_cor <- round2(mapply(function(x,y){corr.test(x$c1,y$response_pctg)$r},c_SDTs,c_resp))

# Group-Level Bayes Factors ------
groupBFs <- function(){
  nsubj <- length(group_nrc[[1]]$observed_d1)
  rg_est <- exp(sapply(group_nrc,function(x){x$group_means})['mu_logMratio',])
  rs_est <- exp(sapply(group_rc,function(x){x$group_means})[c('mu_logMratio_rS1','mu_logMratio_rS2'),])
  #analysis of whichever was higher
    #m_dir <- apply(rs_est,2,function(x){ifelse(x[1]>x[2],1,2)})
    #rs_samples <- lapply(group_rc,function(x){exp(x$samples)})
    #diff_samp <- mapply(function(x,y){x[,y]-x[,-y]},rs_samples,m_dir)
  
  diff_samp <- sapply(group_rc,function(x){exp(x$samples[,2])-exp(x$samples[,1])})

  g_BF <- round2(cbind(
    Mean = apply(diff_samp,2,mean),
    SE = apply(diff_samp,2,sd),
    H1 = rg_est/2, BF=NA,RRl=NA,RRu=NA))
  rownames(g_BF) <- tnames

  for(i in 1:nrow(g_BF)){ #i<-1
    BF <- bfrr(sample_mean=g_BF[i,'Mean'],sample_se=g_BF[i,'SE'],sample_df=nsubj-1,
                model="normal", mean=0, sd=g_BF[i,'H1'], tail=2, criterion=3,
                rr_interval=list(mean=c(0,0),sd=c(0,5),precision=.01))
    g_BF[i,c('BF','RRl','RRu')] <- c(round2(BF$BF),BF$RR$sd)
  }
  
  g_BF[g_BF[,'BF']>1000,'BF'] <- format(g_BF[g_BF[,'BF']>1000,'BF'], scientific=T,digits=2)
  g_BF[g_BF[,'RRu']==5,'RRu'] <- ">5"
  g_BF[,'RRl'] <- paste0('[',g_BF[,'RRl'],', ',g_BF[,'RRu'],']')
  g_BF <- g_BF[,colnames(g_BF)!='RRu']
  write.csv(g_BF,'csv_files/group_BFs.csv')
  return(g_BF)
}
group_BFs <- groupBFs()

# Group-Level Plots ------
ci <- function(x,bound='lower'){
  m <- mean(x)
  cib <- (sd(x)/sqrt(length(x)))*1.96
  if(bound=='lower'){ m-cib } else { m+cib }
}

# Plot group-level bias and M-ratios
groupPlot <- function(){
  group_means <- t(rbind(
    'Type-1 c' = sapply(c_SDTs,function(x){mean(x$c1)}),
    'M-ratio' = sapply(group_nrc,function(x){exp(x$group_means['mu_logMratio'])}),
    'M-ratio S1' = sapply(group_rc,function(x){exp(x$group_means['mu_logMratio_rS1'])}),
    'M-ratio S2' = sapply(group_rc,function(x){exp(x$group_means['mu_logMratio_rS2'])})
  ))
  
  group_err <- lapply(c('lower','upper'),function(x){ #x<-'lower'
    t(rbind(sapply(c_SDTs,function(y){ci(y$c1,x)}),
            sapply(HDIs,function(y){y[,x]}) ))
  })

  par(mar=rep(2,4))
  bias_plot <- barplot(group_means,beside=T,col=c('darkgrey',2:6),main='Group-level Bias and M-ratios',ylim=c(min(group_err[[1]]),max(group_err[[2]])))
  segments(bias_plot[1,2]-1,1,bias_plot[6,4]+1,1,lty=2)
  arrows(bias_plot, group_err[[1]], bias_plot, group_err[[2]],angle=90,code=3,length=.05)
  legend('topleft',legend=tnames,fill=c('darkgrey',2:6))#,horiz=T,xpd=T,cex=.7)#,inset=c(-.09,-.09))
}

jpeg("graphs/group_mratios.jpeg", height=150, width=200,units='mm',quality=200,res=300)
groupPlot()
dev.off()


# Plot type-1 performance
comparePerf <- function(task_arr=1:6){
  group_means <- t(rbind(
    'Accuracy' = sapply(c_resp[task_arr],function(x){mean(x$accuracy_pctg/100)}),
    "d'"       = sapply(c_SDTs[task_arr],function(x){mean(x$d1)}),
    "meta-d'"  = sapply(c_SDTs[task_arr],function(x){mean(x$meta_d)})
    #,"M-ratio"  = sapply(group_nrc[task_arr],function(x){exp(x$group_means['mu_logMratio'])})
  ))
  
  group_err <- lapply(c('lower','upper'),function(x){ #x<-'lower'
    t(rbind(sapply(c_resp[task_arr],function(y){ci(y$accuracy_pctg/100,x)}),
            sapply(c_SDTs[task_arr],function(y){ci(y$d1,x)}),
            sapply(c_SDTs[task_arr],function(y){ci(y$meta_d,x)})
            #,sapply(HDIs[task_arr],function(y){y['mu_logMratio',x]})
    ))
  })
  par(mar=rep(2,4))
  bias_plot <- barplot(group_means,beside=T,col=sub(1,'darkgrey',task_arr),main='Group-level Performance',ylim=c(0,max(group_err[[2]])))
  arrows(bias_plot, group_err[[1]], bias_plot, group_err[[2]],angle=90,code=3,length=.05)
  legend('topleft',legend=tnames[task_arr],fill=sub(1,'darkgrey',task_arr))#,horiz=T,xpd=T,cex=.7)#,inset=c(-.09,-.09))
  #abline(h=1, lty=2)
  segments(bias_plot[1,1]-1,.71,bias_plot[6,1]+1,.71,lty=2)
  segments(bias_plot[1,2]-1,1.1,bias_plot[6,3]+1,1.1,lty=2)
  #abline(h=1.1, lty=2)
  if(6 %in% task_arr){
    round2(rbind(
      mean=apply(group_means[1:5,1:3],2,mean),
      SD=apply(group_means[1:5,1:3],2,sd), 
      breath=group_means[6,1:3]
    ))
  }
}

jpeg('graphs/performance.jpeg', height=150, width=200,units='mm',quality=200,res=300)
comparePerf(1:6)
dev.off()

################################################# Hierarchical Correlations ----
# Hierarchical Correlations --------------
  # Takes several hours to run and may crash computer (overload RAM?)
  # Can just have c_nrsxs in your env and run and save each model individually. Consider reducing samples
# Setup
c_nrsxs <- readRDS('rds_files/c_nrsxs.rds')
Ntask<-6
mcmc_params <- list(response_conditional=0,rc_alternating=rep(c(0,2),3),estimate_dprime=0,
  nchains=3, nadapt=Ntask*1000, nburnin=Ntask*1000, nsamples=Ntask*10000,
  nthin=Ntask, dic=0, rhat=1, parallel=1, monitorparams=0, saveallsamples=0) #monitorparams=c('mu_logMratio', 'sigma_logMratio', 'rho'),

# Full Response Model #
  #c_Hcor_S0 <- fit_meta_d_mcmc_groupCorr(c_nrsxs[1,],c_nrsxs[2,],mcmc_params)
  #saveRDS(c_Hcor_S0,'rds_files/c_Hcor_S0.rds',compress=F);beep()
  c_Hcor_S0 <- readRDS('rds_files/c_Hcor_S0.rds')
  
# YN tasks in R.S. 1 #
  mcmc_params$response_conditional <- 2
  mcmc_params$rc_alternating <- rep(c(0,1),3)
  #c_Hcor_S1 <- fit_meta_d_mcmc_groupCorr(c_nrsxs[1,],c_nrsxs[2,],mcmc_params)
  #saveRDS(c_Hcor_S1,'rds_files/c_Hcor_S1.rds',compress=F);beep()
  c_Hcor_S1 <- readRDS('rds_files/c_Hcor_S1.rds')
  
# YN tasks in R.S. 2 #
  mcmc_params$response_conditional <- 2
  mcmc_params$rc_alternating <- rep(c(0,2),3)
  #c_Hcor_S2 <- fit_meta_d_mcmc_groupCorr(c_nrsxs[1,],c_nrsxs[2,],mcmc_params)
  #saveRDS(c_Hcor_S2,'rds_files/c_Hcor_S2.rds',compress=F);beep()
  c_Hcor_S2 <- readRDS('rds_files/c_Hcor_S2.rds')

# RS Highest M-ratio #
  mcmc_params$response_conditional <- 2
  mcmc_params$rc_alternating <- c(0,2,1,1,0,0)
  #c_Hcor_S3 <- fit_meta_d_mcmc_groupCorr(c_nrsxs[1,],c_nrsxs[2,],mcmc_params)
  #saveRDS(c_Hcor_S3,'rds_files/c_Hcor_S3.rds',compress=F);beep()
  c_Hcor_S3 <- readRDS('rds_files/c_Hcor_S3.rds')

# Per-Protocol data #
  # full <- fit_meta_d_mcmc_groupCorr(nrsx[1,],nrsx[2,],mcmc_params)
  
# 2AFC only (model test) #
  #mcmc_params$response_conditional <- 0
  #cor_2AFC <- fit_meta_d_mcmc_groupCorr(c_nrsxs[1,c(1,3,5)],c_nrsxs[2,c(1,3,5)],mcmc_params)
  #saveRDS(cor_2AFC,'rds_files/cor_2AFC.rds',compress=F);beep()
  #cor_2AFC <- readRDS('rds_files/cor_2AFC.rds')

# Matrix of sample plots
plotMat <- function(Hcor){
  n_task <- NROW(Hcor$group_level)
  rho_idx <- t(combn(1:n_task, 2))
  rho_s <- Hcor$mcmc$samples$rho
  rho_m <- Hcor$rho$mat
  rho_hdi <- Hcor$mcmc$HDI$rho
  
  par(mfrow=c(Ntask-1,Ntask-1), oma=c(0,3,3,0),mar=rep(1,4))
  i<-0
  for(x in (Ntask-1):1){
    if(x<(Ntask-1)){
      for(y in 1:((Ntask-1)-x)){ 
        plot(0, type="n",axes=F,ann=F)
      }
    }
    for(y in 1:x){ i<-i+1
    #if(x==(n_task-1)){title <- names(tasks)[y+1];ann=T} else { ann=F }
    hist(rho_s[,i],freq=F,breaks=50, ann=F,cex.axis=1) #,main=title)#main=colnames(rho_s)[i])#ann=F)#
    #to plot probabilities rather than density: https://stackoverflow.com/a/17433446/7705626
    abline(v=0,col=2, lwd=3, lty='dotted')
    abline(v=rho_m[rho_idx[i,1],rho_idx[i,2]],col=1,lwd=2)
    abline(v=rho_hdi[i,1],col=4, lwd=3, lty='dashed',lend=2)
    abline(v=rho_hdi[i,2],col=4, lwd=3, lty='dashed',lend=2)
    if(y==1){
      mtext(tnames[n_task-x+1], line=-.7, side=3, outer=TRUE, at=grconvertX(0.5,'npc','nic'),cex=.7,font=2)
      mtext(tnames[n_task-x], line=1, side=2, outer=TRUE, at=grconvertY(0.5,'npc','nic'),cex=.7,font=2) #mtext(names(tasks)[n_task-x], line=2, side=2,cex=.6) #
    }
    }
  }
  mtext('Posterior Densities on Correlations', line=1, side=3, outer=TRUE, cex=1,font=2)
}
#plotMat(c_Hcor_S0)

# Hierarchical Pooling ----

hPool <- function(H_cor,file_name){ #H_cor<-c_Hcor_S1;file_name<-'abc'
  #correlation between domains vs within
  #correlation between 2AFC and YN tasks vs between 2AFC tasks
  
  # datasets
  rhos <- r2z(H_cor$rho$mat)
  HDIs <- r2z(H_cor$mcmc$HDI$rho)
  #SEs <- (HDIs[,2]-HDIs[,1])/4
  # idxs
  AFC_idx <- which(names(tasks)%in%c('dots','gabor_2AFC','span_2AFC'))
  YN_idx <- which(names(tasks)%in%c('gabor_YN','span_YN'))
  V_idx <- which(names(tasks)%in%c('dots','gabor_YN','gabor_2AFC'))
  M_idx <- which(names(tasks)%in%c('span_YN','span_2AFC'))
  G_idx <- which(names(tasks)%in%c('gabor_YN','gabor_2AFC'))
  B_idx <- which(names(tasks)%in%'breath')
  
  # row and col idx pairs
  effect_idxs <- list(#within task type
                      'AFC'=list(AFC_idx,AFC_idx),
                      'YN'=list(YN_idx,YN_idx),
                      #between task type
                      'YN_2AFC'=list(AFC_idx,YN_idx),
                      #within domain
                      'gabor'=list(G_idx,G_idx), #within task variant
                      'visual'=list(V_idx,V_idx),
                      'memory'=list(M_idx,M_idx), #within task variant
                      #between domain
                      'between_domains'=list(V_idx,M_idx),
                      'vis_b'=list(V_idx,B_idx), #not used
                      'mem_b'=list(M_idx,B_idx)) #not used
  
  #selcted cors
  selected_cors <- lapply(effect_idxs, function(x){ #x<-rc_idxs[['FC']]
    sub_mat <- rhos[x[[1]],x[[2]]]
    if(identical(x[[1]],x[[2]])){
      sub_mat <- sub_mat[lower.tri(sub_mat)]
    }
    return(sub_mat)
    #ses <- SEs[which((idx[,1]%in%x[[1]] & idx[,2]%in%x[[2]]) | (idx[,2]%in%x[[1]] & idx[,1]%in%x[[2]]))]
    #cbind(rho=c(sub_mat),se=ses)
  })
  
  #concatenate within-domain effects
  selected_cors <- c(selected_cors,list(within_domains = do.call('rbind',selected_cors[c('visual','memory')]),
                                        task_variants = do.call('rbind',selected_cors[c('gabor','memory')])))
  
  #mean differences
  categories <- list(YN_2AFC=c('YN','AFC'),
                     YN2AFC_AFC=c('YN_2AFC','AFC'),
                     #task_variants=c('task_variants','within_domains'),
                     domains=c('between_domains','within_domains'))
  
  stuff <- t(sapply(categories,function(x){ #x<-categories[[1]]
    lower <- selected_cors[[x[1]]]
    higher <- selected_cors[[x[2]]]
    #ses <- c(lower[,'se'],higher[,'se'])
    l_m <- round2(mean(lower)) #[,'rho']
    h_m <- round2(mean(higher))
    m <- r2z(h_m-l_m)
    #SE <- sqrt(sum(ses^2))
    z_SE <- round2(1/sqrt(nrow(H_cor$Mratio)-3))
    SE <- sqrt((z_SE^2)*2)
    H1 <- (h_m)/2
    BF_inputs <- round2(c(m=m,SE=SE,H1=H1))
    BF <- bfrr(sample_mean=BF_inputs['m'],sample_se=SE,sample_df=nrow(H_cor$Mratio)-1,
               model="normal", mean=0, sd=BF_inputs['H1'], tail=1, criterion=3,
               rr_interval=list(mean=c(0,0),sd=c(0,1)),precision=.01)
    c(lower=x[1],higher=x[2],l_m=l_m,h_m=h_m,BF_inputs,BF=BF$BF,RR=BF$RR$sd)
  }))
  stuff[,3:8] <- round2(as.numeric(stuff[,3:8]))
  #stuff[,'m'] <- paste0(stuff[,'m'],' (',stuff[,'SE'],')')
  stuff[stuff[,'RR1']==0.01,'RR1'] <- '0'
  stuff[stuff[,'RR2']==1,'RR2'] <- '>1'
  stuff[,'RR1'] <- paste0('[',stuff[,'RR1'],' ,',stuff[,'RR2'],']')
  stuff <- stuff[,-which(colnames(stuff)%in%c('SE','RR2'))]
  return(stuff)
}

cat_s0 <- hPool(c_Hcor_S0,"S0_category_test.csv")
cat_s1 <- hPool(c_Hcor_S1,'s1_category_test.csv')
cat_s2 <- hPool(c_Hcor_S2,'s2_category_test.csv')
cat_s3 <- hPool(c_Hcor_S3,'s3_category_test.csv')

write.csv(cat_s0,'csv_files/s0_category_test.csv')
write.csv(cat_s1,'csv_files/s1_category_test.csv')
write.csv(cat_s2,'csv_files/s2_category_test.csv')
write.csv(cat_s3,'csv_files/s3_category_test.csv')
# Hierarchical Relibailities -----
#odd even splits
#pp_split <- splitData(tasks) #per-protocol

#Standard S0 Reliabilities
  #c_rels_S0 <- h_rels$'20'
  #c_rels_S0$breath <- relibility(1,c_split[,6,drop=F],3)
  #saveRDS(c_rels_S0,'rds_files/c_rels_S0.rds')
c_rels_S0 <- readRDS('rds_files/c_rels_S0.rds')

#RS1 reliabilities
  #c_rels_S1 <- relibility(1,c_split,3)
  #saveRDS(c_rels_S1,'rds_files/c_rels_S1.rds')
c_rels_S1 <- readRDS('rds_files/c_rels_S1.rds')

#RS2 reliabilities
  #c_rels_S2 <- relibility(2,c_split,3)
  #saveRDS(c_rels_S2,'rds_files/c_rels_S2.rds')
c_rels_S2 <- readRDS('rds_files/c_rels_S2.rds')

#Examine convergence
rel_rhat <- rbind(s0=sapply(c_rels_S0,function(x){x$mcmc$rhat$full$mpsrf}),
                  s1=sapply(c_rels_S1,function(x){x$mcmc$rhat$full$mpsrf}),
                  s2=sapply(c_rels_S2,function(x){x$mcmc$rhat$full$mpsrf}))

# Correct Relibailities -----
#oops - names would be good
names(c_rels_S0) <- names(tasks)
names(c_rels_S1) <- names(tasks)
names(c_rels_S2) <- names(tasks)
# select tasks with highest M-ratio (do first)
c_rels_S3 <- c_rels_S0
c_rels_S3[c('gabor_2AFC','span_YN')] <- c_rels_S1[c('gabor_2AFC','span_YN')]
c_rels_S3['gabor_YN'] <- c_rels_S2['gabor_YN']
# select 2AFC tasks
c_rels_S1[AFC] <- c_rels_S0[AFC]
c_rels_S2[AFC] <- c_rels_S0[AFC]


# get hrels in correct format
relList <- function(rels){
  rel_l <- t(rbind(r=sapply(rels,function(x){x$rho$list}),
                   sapply(rels,function(x){x$mcmc$HDI$rho[1,]})))
  round2(sb(rel_l))
}

rels_h <- list(M_ratio = relList(c_rels_S0),
               M_ratio_rS1 = relList(c_rels_S1),
               M_ratio_rS2 = relList(c_rels_S2),
               h_mr = relList(c_rels_S3))

# Hierarchical BFs ---------
Hcor <- function(cors,rels){ # cors<-c_Hcor_S3; rels<-rels_h[[3]]
  ntask <- ncol(cors$d1)
  # combined reliabilities
  H1s <- ( rels[idx[,1],'r'] * rels[idx[,2],'r'] ) **.5 /2
  # actual correlations
  ztran <- round2(r2z(cbind(r=cors$rho$list, cors$mcmc$HDI$rho, H1=H1s))) #ztransfrom stuff for BFs
  rhos <- cbind(ztran,SE=(ztran[,'high']-ztran[,'low'])/4)#(qnorm(.975)*2)) #new SE transform 
  #BFs
  BFs <- cbind(rhos[,F],BF=NA,RRl=NA,RRu=NA)
  for(i in 1:nrow(rhos)){ # i<-1
    BF <- bfrr(sample_mean=rhos[i,'r'],sample_se=rhos[i,'SE'],sample_df=nrow(cors$Mratio)-1,
               model="normal", mean=0, sd=H1s[i], tail=1, criterion=3,
               rr_interval=list(mean=c(0,0),sd=c(0,1)),precision=.01)
    BFs[i,] <- c(BF$BF,BF$RR$sd)
  }
  output <- round2(cbind(r=cors$rho$list, 
                         cors$mcmc$HDI$rho, 
                         h1=rhos[,'H1'], 
                         BFs))
  #sort formatting issues with my Hmeta-d function
  colnames(output)[colnames(output)%in%c('low','high')] <- c("lower","upper")
  return(output)
}

#run
bf_h <- list(M_ratio = Hcor(c_Hcor_S0,rels_h[[1]]),
             M_ratio_rS1 = Hcor(c_Hcor_S1,rels_h[[2]]),
             M_ratio_rS2 = Hcor(c_Hcor_S2,rels_h[[3]]),
             highest_mratio = Hcor(c_Hcor_S3,rels_h[[4]]))


############################################# Non-hierarchical Correlations ----
# Non-Hierarchical Correlations -----

# List variables to correlate
corList <- function(c_SDTs){
  # datasets of variables to correlate
  cor_vars <- lapply(var_names,function(x){ sapply(c_SDTs,function(y){y[,x]}) })
  names(cor_vars) <- var_names
  
  # Replace tasks for different analyses
  cor_vars$highest_mratio <- cor_vars$M_ratio
  cor_vars$highest_mratio[,c('gabor_2AFC','span_YN')] <- cor_vars$M_ratio_rS1[,c('gabor_2AFC','span_YN')]
  cor_vars$highest_mratio[,'gabor_YN'] <- cor_vars$M_ratio_rS2[,'gabor_YN']
  cor_vars$M_ratio_rS1[,AFC] <- cor_vars$M_ratio[,AFC] # 'No' for YN tasks only
  cor_vars$M_ratio_rS2[,AFC] <- cor_vars$M_ratio[,AFC] # 'Yes' for YN tasks only
  
  return(cor_vars)
}
# get main correlations
data_nh <- corList(c_SDTs)
list_nh <- lapply(data_nh,corr.test)
z_SE <- round2(1/sqrt(list_nh$M_ratio$n-3))
cor_nh <- lapply(list_nh,function(x){ cbind(x$ci[,1:3], se=z_SE) }) #x$se[lower.tri(x$se)]

# Non-Hierarchical Relibailities -----

# reliabilities
reliabilitiesNh <- function(split_SDT){
  nh_even <- corList(split_SDT$even) #list vars to correlate
  nh_odd <- corList(split_SDT$odd)
  nh_cors <- mapply(corr.test,nh_even,nh_odd)['ci',] # correlate odd and even
  allv <- expand.grid(tnames,tnames) # all possible combinations of vars
  nh_rels <- lapply(nh_cors,function(x){ x[which(allv[,1]==allv[,2]),1:3]}) #reduce to reliabilities
  nh_rels <- lapply(nh_rels,sb) #Spearman-Brown Adjust
  lapply(nh_rels,round2) #round
}

#get reliabilities list
rels_nh <- reliabilitiesNh(split_SDT)

# Non-Hierarchical BFs -----
#Bayes Factors
bfsNh <- function(rels_nh,cor_nh){
  #reliabilities to H1 #note these are sb adjusted above
  h1_nh <- sapply(rels_nh,function(x){ (x[idx[,1],'r'] * x[idx[,2],'r']) **.5 /2 }) #combined reliabilities
  #r to z transform
  h1_nh <- r2z(h1_nh) 
  rz <- lapply(cor_nh,r2z)
  # list
  eff_nh <- lapply(1:length(cor_nh),function(x){cbind(r=rz[[x]][,'r'], #z transormed correlation
                                                      se=cor_nh[[x]][,'se'], #standard se
                                                      h1=h1_nh[,x])}) #combined reliabilities
  names(eff_nh) <- names(cor_nh) #nicer names
  #get bayes factors
  lapply(names(eff_nh),function(x){ # x<-eff_nh[[1]]
    eff_tab <- eff_nh[[x]]
    BFs <- c()
    for(i in 1:nrow(eff_tab)){ # i<-1
      if(is.nan(eff_tab[i,'h1']) | eff_tab[i,'h1']<0){ # skip if H1 is NA or <0
        BFs <- rbind(BFs,rep(NA,3))
        next 
        }
      #Bfs
      BF <- bfrr(sample_mean=eff_tab[i,'r'], sample_se=eff_tab[i,'se'],
                 sample_df=nrow(c_SDTs$dots)-1, model="normal", mean=0, 
                 sd=eff_tab[i,'h1'], tail=1, criterion=3, 
                 rr_interval=list(mean=c(0,0),sd=c(-1,1)),precision=.01)
      BFs <- rbind(BFs,c(BF=BF$BF,RR=BF$RR$sd))
    }
    colnames(BFs) <- c('BF','RRl','RRu')
    cbind(cor_nh[[x]],round2(BFs)) #final tables
  })
  
}

#get bayes factors
bf_nh <- bfsNh(rels_nh,cor_nh)

################################################################ Formatting ----
# Concatenate columns -----
matChars <- function(mat,cortype='r',err='CI',ntask=6){ # mat=bf_h
  lapply(mat,function(x){ # x<-mat[[1]]
    x <- round2(data.frame(x))
    p <- paste0(cortype,'=',x$r,'\r',err,'=[',x$lower,', ',x$upper,']')
    if('BF'%in%colnames(mat[[1]])){
      x[is.na(x)] <- ''
      p <- paste0(p,'\rBF=',x$BF,'\rRR=[',x$RRl,', ',x$RRu,']')
      p <- gsub('\\rBF=\\rRR=\\[,\\s]','\r-',p)
    }
    return(p)
  })
}

# Print as matrices -----
matify <- function(h_chars,nh_chars,r_chars){ #h_chars=h_list[[1]];nh_chars=nh_list[[1]];r_chars=r_list
  ntask <- length(r_chars[[1]])
  #Make mat
  lapply(1:length(h_chars),function(x){ # x<-1
    mat <- matrix(1,ntask,ntask)
    for(i in 1:nrow(idx)){ #i<-6
      mat[idx[i,1],idx[i,2]] <- nh_chars[[x]][i] #top
      if(i<=ntask){mat[i,i] <- r_chars[[x]][i]} #diag
      mat[idx[i,2],idx[i,1]] <- h_chars[[x]][i] #bottom
    }
    dimnames(mat) <- list(tnames,tnames)
    write_excel_csv(data.frame(mat),paste0('csv_files/',names(h_chars)[x],'_mat.csv'))
    return(mat)
  })
}

#matrix has H on bottom, NH on top, both inside the 
correlationMatrices <- function(bf_h,rels_h,bf_nh,rels_nh){
  bf_h_chars <- matChars(bf_h,'\u03C1','HDI')
  rels_h_chars <- matChars(rels_h,'\u03C1','HDI')
  bf_nh_chars <- matChars(bf_nh,'r','CI')
  rels_nh_chars <- matChars(rels_nh,'r','CI')
  r_list <- lapply(1:4,function(x){paste0(rels_nh_chars[[x]],'\n\n',rels_h_chars[[x]])}) #mapply(paste0,nh_mats[[2]],h_mats[[2]])
  mats <- matify(bf_h_chars,bf_nh_chars,r_list)
}

correlationMatrices(bf_h,rels_h,bf_nh,rels_nh)

lapply(names(bf_h),function(x){write.csv(bf_h[[x]],paste0('csv_files/',x,'_list.csv'))})
# Reliabilities barplot ----
# done here so to add in breath
#c_rels_S0 <- readRDS('rds_files/c_rels_S0.rds')

reliabilityPlot <- function(i_rels,h_rels){
  rels_i <- i_rels[c('d1','meta_d','M_ratio'),]
  rels_h <- c(sapply(h_rels$'20',function(x){x$rho$list}),breath=c_rels_S0$breath$rho$list)
  rels_a <- t(rbind(rels_i,rels_h))
  colnames(rels_a) <- c("d'","meta-d'","M-ratio","Hmeta-d")
  jpeg('graphs/reliabilities.jpeg',height=200,width=200,units='mm',res=100,quality=200)
  par(mar=rep(3,4))
  barplot(rels_a,beside=T,col=c(1:6),main='Reliabilities',cex.axis=2,cex.names=2,cex.main=2)
  legend('topleft',tnames,col=1:6,pch=15,text.font=2)
  dev.off()
}

reliabilityPlot(i_rels,h_rels)
####################################################################### Fin ----
# End timer & save -----
save.image(file='gfactor.RData')
end_time <- Sys.time()-start_time
print(end_time)
