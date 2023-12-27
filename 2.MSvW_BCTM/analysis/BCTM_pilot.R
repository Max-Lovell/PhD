# Analysis code for MSvW BCT-M Pilot study by Max Lovell #
# Setup ----
#make working dir same as script file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Note have BCTM env file from run of main study's code first
load('BCTM.RData')
# or, can source, recommend commenting out anything that can be loaded as rds instead
start <- Sys.time()
source("BCTM.R")
start - Sys.time()
source("fit_meta_d_mcmc.R")
# Helpers -------------------------
pnames <- c("2SO"="SOW 1&2","n_breath"="N Breath","n_reset"="N Reset","adj_dp2rs"="Adj. Type-2 d' (rs)",
            #same as main study
            "TMS"="TMS-D", "SOW_F"="SOW 1st", "SOW_S"="SOW 2nd",
            "PHQ"="PHQ-8", "GAD"="GAD-7", "RRS"="RRS",
            "RRS_D"="RRS Depression", "RRS_B"="RRS Brooding", "RRS_R"="RRS Rumination",
            "WBSI"="WBSI", "WBSI_sup"="WBSI Supression", "WBSI_int"="WBSI Intrusion",
            "acc"="Accuracy", "dp"="d'", "c"="c", "adj_dp"="Adj. d'",
            "dp2"="Type-2 d'", "c2"="Type-2 c", "adj_dp2"="Adj. type-2 d'",
            "adj_meta_d"="Adj. meta-d'", "H_meta_d"="Hmeta-d", "H_M_ratio"="HM-ratio",
            "ln_M_ratio"="log(M-ratio)","dp2rs"="Type-2 d' (rs)")

###################################################################### Data ----
# Task Data ------------
readPilot <- function(){
  task_files <- list.files("pilot_data/task_data/")
  task <- c()
  for(file_n in task_files){ #file_n<-task_files[1]
      #note warnings when reading data files: had issues creating a CSV file in JavaScript with \r\n and URIs being read literally
        #probably due to sanitizing JSON string for secure sending with JQuery AJAX to PHP file so only data and commas were used.
    file <- tryCatch({
        suppressWarnings(read.table(paste0("pilot_data/task_data/",file_n),header=FALSE,sep=",")) 
      }, error = function(e) { 
        data.frame("no_data","-1") 
        })
    if(ncol(file)==1){ file <- data.frame("no_data","-1") } #don's use NA as resp NAs are resets and need as.numeric conversion
    file <- gsub("\\[|\\]","",file) #file construction error
    if(length(file)==0){print(f);next} #print names of incorrect files
    file <- data.frame(file)
    
    p_data <- data.frame(ID=sub("\\D+","",file_n),
                         time=sub(".*(pre|post).*","\\1",file_n),
                         resp=file[c(T,F),],
                         rt=file[c(F,T),])
    task <- rbind(task,p_data)
  }
  task$resp[task$resp=="NA"] <- "reset"
  task$resp <- tolower(gsub("Arrow","",task$resp))
  task$rt <- as.numeric(task$rt)
  return(task)
}

pilot_raw <- readPilot()

cleanPilot <- function(task){ # task <- br_raw
  # delete repeated confidence
  confs <- grepl("right|left",task$resp) & duplicated(task$ID) # confidences, excluding when first row of each id
  has_conf <- which(confs) # rows with confidence ratings
  conf_reps <- has_conf[diff(c(0,has_conf))==1] # rows with successive confidence presses
  task <- task[-(conf_reps-1),] # delete repeated confidence ratings
  # move confidence to corresponding breath
  confs <- grepl("right|left",task$resp) & duplicated(task$ID) # confidences, excluding when first row of each id
  task[which(confs)-1,c('conf','conf_rt')] <- task[confs,c('resp','rt')]
  task <- task[!confs,]
  # breath count
  ends <- grepl("down|reset",task$resp) | !duplicated(task$ID,fromLast=T) # confidences, excluding when first row of each id
  counts <-  diff(c(0,which(ends)))
  task$count <- unlist(lapply(counts,seq))
  # resets
  task[task$resp=='reset','count'] <- NA
  # breath sequence
  downs <- task$resp=='down'| !duplicated(task$ID,fromLast=T)
  task$sequence <- unlist(tapply(downs,task$ID,function(x){
    rep(1:sum(x),diff(c(0,which(x)))) # number of presses between TRUEs (starting at 0), with sequence number repeated the relevant number of times
  }))
  #correct and incorrect responses
  task$stim <- ifelse(task$count==9,'down','up')
  task$conf <- ifelse(task$conf=='left',1,2)
  task$correct <- task$resp==task$stim
  task$ID <- factor(task$ID)
  return(task)
}

pilot_cl <- cleanPilot(pilot_raw)

write.csv(pilot_cl,"csv_files/pilot_raw_BCT_data.csv")

# SDT -----
pilotSDT <- function(dataset,resps=c('up','down')){ # dataset<-br_cl
  #breath accuracy
  n_reset <- tapply(dataset$resp=='reset',dataset$ID,sum)
  dataset <- dataset[dataset$resp!='reset',]
  n_breath <- as.numeric(table(dataset$ID))
  acc <- tapply(dataset$stim=='down' & dataset$resp=='down', dataset$ID, sum)/
            tapply(dataset$resp=='down', dataset$ID, sum)

  # SDT
  nr_sx <- NRSX(dataset,resps)
  nconf <- ncol(nr_sx)/4
  #adjust for 0 counts
  nr_sx_adj <- nr_sx+.5 #from Hmeta-d is adjusted by (1/(n_rat*2)) instead - should standardise this?
  t1 <- type1(nr_sx_adj,nconf)
  
  ## new t2 ##
  dataset <- dataset[!is.na(dataset$conf),]
  ###########
  nr_sx <- NRSX(dataset,resps)
  nr_sx_adj <- nr_sx+.5 #from Hmeta-d is adjusted by (1/(n_rat*2)) instead - should standardise this?
  t2 <- type2(nr_sx_adj,nconf)
  t2_rs <- type2RS(nr_sx_adj,nconf)
  metad <- Hmetad(nr_sx)
  metad_rs <- HmetadRS(nr_sx)
  
  #combine
  output <- data.frame(n_reset,n_breath,acc,
                       t1,t2,t2_rs,
                       metad,metad_rs) #dataframe for easier column assignment below

  #Add ID and time
  output$ID <- rownames(nr_sx)
  output$time <- tapply(dataset$time,dataset$ID,function(x){x[1]})
  return(output)
}

pilot_SDT <- pilotSDT(pilot_cl)
#saveRDS(pilot_SDT,'rds_files/pilot_BCT_performance.rds')
#write.csv(pilot_SDT,'csv_files/pilot_BCT_performance.csv')
pilot_SDT <- readRDS('rds_files/pilot_BCT_performance.rds')

SDT_vars <- colnames(pilot_SDT)[!colnames(pilot_SDT)%in%c("ID","time",'n_breath','n_reset')]

# Survey data ------------------
pilotSurvey <- function(){
  survey_files <- list.files("pilot_data/survey_data")
  file_rgxs <- paste("^breath",c("pre",paste0("[w|m]",1:10),"post"),"",sep="_")
  for(file_rgx in file_rgxs){ #loop through the files for pre, post, and each day separately #file_rgx<-"^breath_[w|m]1_"
    day_files <- survey_files[grep(file_rgx,survey_files)] #read in first dataset to add the others to
    day_data <- read.csv(paste0("pilot_data/survey_data/",day_files[1]))
    day_data <- day_data[3:nrow(day_data),]
    
    for(file_name in 2:length(day_files)){ #loop through the rest and add the rows
      file_name <- day_files[file_name]
      day_file <- read.csv(paste0("pilot_data/survey_data/",file_name))
      day_file <- day_file[3:nrow(day_file),]
      day_data <- rbind(day_data,day_file) #add file to rest of data for that day
    }
    suffix <- gsub('\\^|(breath)|_|(\\[w\\|m\\])','',file_rgx)
    names(day_data)[!grepl("Email",names(day_data))] <- paste(names(day_data)[!grepl("Email",names(day_data))],suffix,sep="_")
    
    #bit of data cleaning, fix variable names, etc
    if(file_rgx=="^breath_pre_"){
      #remove unfinished responses
      day_data <- day_data[day_data$Progress_pre==100 & day_data$Duration..in.seconds._pre>0 & day_data$Finished_pre==1 & day_data$DistributionChannel_pre!="preview",]
      day_data$Gender_pre[day_data$Gender_pre == 3] <- "male"
      day_data$Gender_pre[day_data$Gender_pre == 4] <- "female"
      day_data$Gender_pre[day_data$Gender_pre == 5] <- "other"
      day_data$Gender_pre[day_data$Gender_pre == 6] <- "opt_out"
      #day_data[day_data$Email=="[redacted]@redacted","Email"] <- "[redacted]@redacted" #personally contacted about not recieving emails by this person, who entered their email wrong
      BCT_pre <- pilot_SDT[pilot_SDT$time=="pre",]
      names(BCT_pre) <- paste(names(BCT_pre),suffix,sep="_")
      survey_data <- merge(day_data,BCT_pre, by.x="Random_ID_pre", by.y="ID_pre", all=T) #add task data
    } else {
      if(file_rgx=="^breath_post_"){#fixing some badly named vars in qualtrics:
        # remove unfinished responses
        day_data <- day_data[day_data$Progress_post==100 & day_data$Duration..in.seconds._post>0 & day_data$Finished_post==1 & day_data$DistributionChannel_post!="preview",]
        TMS_names <- grepl("^Q1_[0-9]{1,2}_post$",names(day_data))
        names(day_data)[TMS_names] <- gsub("Q1","TMS",names(day_data)[TMS_names])
        GAD_names <- grepl("^(Q1_)[0-9]{1,2}\\.1{1}_post$",names(day_data))
        names(day_data)[GAD_names] <- gsub("Q1","GAD",names(day_data)[GAD_names])
        names(day_data)[GAD_names] <- gsub("\\.1","",names(day_data)[GAD_names])
        names(day_data) <- gsub("PHQ.8","PHQ",names(day_data))
        names(day_data) <- gsub("Q16","RRS",names(day_data))
        names(day_data) <- gsub("Q17","WBSI",names(day_data))
        BCT_post <- pilot_SDT[pilot_SDT$time=="post",]
        names(BCT_post) <- paste(names(BCT_post),suffix,sep="_")
        day_data <- merge(day_data,BCT_post, by.x="Random_ID_post", by.y="ID_post", all=T) #add breath task data
      } else {
        if(file_rgx=="^breath_[w|m]1_"){
          phq_exp <- grepl("^exp_GAD_([8,9]|[0-9]{2})_1$",names(day_data))
          names(day_data)[phq_exp] <- paste0("exp_PHQ_",seq(1:8),"_1")
        }
        names(day_data) <- gsub("Q1_1|focus_1","focus",names(day_data))
        names(day_data) <- gsub("Q2","activities",names(day_data))
        #clean before merge
        day_data <- day_data[day_data[,grep('Progress',colnames(day_data))]==100,]
        day_data <- day_data[day_data[,grep('Finished',colnames(day_data))]==1,]
        day_data <- day_data[day_data[,grep('Duration..in.seconds.',colnames(day_data))]>0,]
        day_data <- day_data[day_data[,grep('DistributionChannel',colnames(day_data))]!="preview",]
      }
      names(day_data)[names(day_data) == 'RecipientEmail'] <- 'Email'
      #add to the other days' data by email
      survey_data <- merge(survey_data,day_data,by="Email",all=T)
    }
  }
  
  #Remove duplicate and test emails
  #Remove duplicate emails: people who completed the same survey more than once.
  survey_data$Email <- tolower(survey_data$Email)
  dups <- survey_data$Email[duplicated(survey_data$Email)]
  print(paste0('n with duplicated emails: ',length(unique(dups))))
  survey_data <- survey_data[!survey_data$Email %in% dups,]
  #test_emails <- c("redacted@redacted","redacted@redacted")
  #survey_data <- survey_data[!(survey_data$Email %in% test_emails),]
  rownames(survey_data) <- NULL #reset row numbers
  
  #drop unneeded variables
  drops <- c("Status","IPAddress","Progress","Finished","RecordedDate","ResponseId",
             "Name","RecipientEmail","ExternalReference", "Location","DistributionChannel","UserLanguage",
             "Email.validation","Prolific.ID","SC0","PROLIFIC_PID","always1",
             #add the following back in later?:
             "Date","Duration","Consent","Resolution","Click","timer",".Submit",
             "Q32_post","Q33_post","Q34_post","Q34_4_TEXT_post","Q22_post","Q22_3_TEXT_post","Q35_post" 
             #post-test feedback qs: Q32 = time_total, Q33 = course_content, Q34 = useful, Q34_4_TEXT = useful_text, 
             #Q22 = authentic, Q22_3_TEXT = authentic_text, Q35 = comments
  )
  survey_data <- survey_data[!grepl(paste(drops, collapse = "|"),names(survey_data))]
  
  #SORT TYPE ISSUES
  non_numeric <- grepl("^(TMS|GAD|PHQ|SOW|RRS|WBSI)_[0-9]{1,2}_(pre|post)$|^exp_",names(survey_data))
  survey_data[,non_numeric] <- sapply(survey_data[,non_numeric],as.numeric)
  survey_data <- survey_data[survey_data$condition_pre!='',] #not sure what happened here
  survey_data$condition_pre <- factor(survey_data$condition_pre)
  #survey_data <- survey_data[!is.na(survey_data$condition_pre),]
  return(survey_data)
}

pilot_data <- pilotSurvey()
#saveRDS(pilot_data,'rds_files/pilot_anon_raw_data.rds')
#write.csv(pilot_data,'csv_files/pilot_anon_raw_data.csv')
pilot_data <- readRDS('rds_files/pilot_anon_raw_data.rds')
nrow(pilot_data)
# Exclusions ---------------------------
#sdt_cols<-paste0(rep(SDT_vars,each=2),c('_pre','_post'))
# missing post
missing_post <- is.na(pilot_data$Random_ID_post)
pilot_data <- pilot_data[!missing_post,] 
# not missing some post
nmissing_some_post <- complete.cases(pilot_data[,which(names(pilot_data)=="TMS_1_post"):which(names(pilot_data)=="WBSI_15_post")])
pilot_data <- pilot_data[nmissing_some_post,] 
# not missing some expectancies
nmissing_exp <- complete.cases(pilot_data[,grep("exp",names(pilot_data))]) | pilot_data$condition_pre=="control"
pilot_data <- pilot_data[nmissing_exp,] 
# accuracy
bad_acc <- which(pilot_data$acc_pre==0|is.na(pilot_data$acc_pre)|is.nan(pilot_data$acc_pre)|
                 pilot_data$acc_post==0|is.na(pilot_data$acc_post)|is.nan(pilot_data$acc_post))
pilot_data <- pilot_data[-bad_acc,]
# d-prime
bad_dp <- which(pilot_data$dp_pre==0|is.na(pilot_data$dp_pre)|is.nan(pilot_data$dp_pre)|
                   pilot_data$dp_post==0|is.na(pilot_data$dp_post)|is.nan(pilot_data$dp_post))
#pilot_data[bad_dp,sdt_cols] 
pilot_data <- pilot_data[-bad_dp,] 
#breath rate
br_range <- c(12-5,20+5)*10 #+/-5bpm buffer around 12-20bpm (Sapra et al., 2022) for 10 mins
br_rate <- which(pilot_data$n_breath_pre<br_range[1]|pilot_data$n_breath_pre>br_range[2]|
                   pilot_data$n_breath_post<br_range[1]|pilot_data$n_breath_post>br_range[2])
pilot_data <- pilot_data[-br_rate,]
nrow(pilot_data)
table(pilot_data$condition_pre)
exclusions <- rbind(same_twice=23,
                    missing_cond=1,
                    missing_some_post=sum(!nmissing_some_post),
                    missing_exp=sum(!nmissing_exp),
                    missing_post=sum(missing_post),
                    bad_acc=length(bad_acc),
                    bad_dp=length(bad_dp),
                    br_rate=length(br_rate))

# Clean pilot data ------------------------------
# recode vars
pilot_data <- scale_recode(pilot_data,"TMS",c(1:7)) #  TMS: https://sci-hub.yncjkj.com/10.1891/0889-8391.23.3.185.
pilot_data <- scale_recode(pilot_data,"GAD",c(1:7))
pilot_data <- scale_recode(pilot_data,"PHQ",c(1:8)) # PHQ-8: https://asset-pdf.scinapse.io/prod/1985329916/1985329916.pdf

# convert to numeric
for(i in 1:length(scales)){
  sc_names <- c(scales[[i]]$pre,scales[[i]]$post)
  pilot_data[,sc_names] <- sapply(pilot_data[,sc_names], as.numeric)
}

# diff vars
pilot_data <- meanDiffs(pilot_data)

saveRDS(pilot_data,"rds_files/pilot_dataset.rds")
write.csv(pilot_data,"csv_files/pilot_dataset.csv")
pilot_data <- readRDS("rds_files/pilot_dataset.rds")

#pilot_data[pilot_data$H_M_ratio_pre>20000,sdt_cols]
################################################################## Analysis ----
# Windsorising functions -----------
trimMean <- function(x,q=.1){
  x <- x[!is.na(x)&!is.infinite(x)]
  trim <- x[x>=quantile(x,q) & x<=quantile(x,1-q)]
  mean(trim)
}

trimSE <- function(x,q=.1){
  x <- x[!is.na(x)&!is.infinite(x)]
  n <- length(x)
  # trimmed mean
  trim <- x[x>=quantile(x,q) & x<=quantile(x,1-q)]
  nr <- length(trim)
  # windsorised min/max replacement for SE
  windsor_low <- rep(min(trim),sum(x<quantile(x,q)))
  windsor_high <- rep(max(trim),sum(x>quantile(x,1-q)))
  windsor <- c(windsor_low,trim,windsor_high)
  S2w <- var(windsor)
  sqrt(S2w/nr*(n-1)/(nr-1))
}

#trim_ms <- apply(pilot_data[,paste0(SDT_vars,"_diff")],2,tapply,pilot_data$condition_pre,trimMean)
#trim_ses <- apply(pilot_data[,paste0(SDT_vars,"_diff")],2,tapply,pilot_data$condition_pre,trimSE)

# Descriptives ----------
shortlist <- c("TMS","SOW_F","SOW_S","PHQ","GAD","RRS","WBSI","acc","dp","dp2rs","H_M_ratio")

groupDesc <- function(dataset){ # dataset<-pilot_data
  cells <- function(vars,mf,sef){ #vars=names(scales);mf=mean;sef=se
    #vars <- vars[vars%in%shortlist] #shortlist
    col_names <- paste0(rep(levels(dataset$condition_pre),each=2),c('_pre','_post'))
    means <- matrix(NA,length(vars),length(col_names),dimnames=list(vars,col_names))
    ses <- means
    for(time in c('_pre','_post')){
      for(v in vars){
        means[v,paste0(names(abc),time)] <- tapply(as.numeric(dataset[,paste0(v,time)]),dataset$condition_pre,mf)
        ses[v,paste0(names(abc),time)] <- tapply(as.numeric(dataset[,paste0(v,time)]),dataset$condition_pre,sef)
      }
    }
    list(mean=round2(means),se=round2(ses))
  }
  scale_desc <- cells(names(scales),mean,se)
  sdt_desc <- cells(SDT_vars,trimMean,trimSE)
  Map(rbind,scale_desc,sdt_desc)
}

pilot_desc <- groupDesc(pilot_data)

prettyDesc <- function(mses){ #mses<-pilot_desc
  #pretty print
  tab <- matrix(paste0(mses$mean,' (',mses$se,')'),nrow(mses[[1]]),ncol(mses[[1]]),dimnames=dimnames(mses[[1]]))
  tab <- tab[shortlist,]
  rownames(tab) <- pnames[rownames(tab)]
  colnames(tab) <- tools::toTitleCase(sub("_"," ",colnames(tab)))
  return(tab)
}

pilot_desc_pp <- prettyDesc(pilot_desc)
write.csv(pilot_desc_pp,'csv_files/pilot_group_means.csv')

# Expectancies ---------
expects <- expectancies(pilot_data)
write.csv(expects,"csv_files/pilot_expectancies.csv")
# Analysis -----------
shortlist_group <-c("TMS","SOW_F","SOW_S","PHQ","GAD","RRS","WBSI","adj_dp2","adj_dp2rs","adj_meta_d","adj_meta_d_rs","H_M_ratio")

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
                'adj_dp2rs' = adj(data_c,'dp2rs',c('c2','dp','c')), #adjusted r.s. type-2 d'
                'adj_meta_d'= adj(data_c,'H_meta_d',c('dp','c')), #adjusted Hmeta-d
                'adj_meta_d_rs' = adj(data_c,'H_meta_d_rS2',c('c2','dp','c'))) #adjusted r.s. type-2 d'
    
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
  comps <- comps[!rownames(comps)%in%c('n_reset','c','c1','c2','c2rs','ln_M_ratio','H_meta_d_rS1','H_M_ratio_rS1','ln_M_ratio_rS1','ln_M_ratio_rS2'),] # c is too small to run - log(M) has NAs
  # Bayes Factors
  for(comp in c('mvw','wvwl','mvwl')){ #comp<-'mvw'
    for(i in rownames(comps)){ #i<-rownames(comps)[1] 
      #print(paste0(comp,': ',i))
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

res_pilot <- BFs(pilot_data)

# Correlations ------------
#all_vars <- c(names(scales),SDT_vars)
cor_vars <- c("TMS","GAD","PHQ","SOW_F","SOW_S","RRS","WBSI",
              "n_breath","n_reset","acc","dp","dp2","dp2rs","H_M_ratio")
pilot_cors <- round2(cor(pilot_data[,paste0(cor_vars,'_pre')]))
pilot_cors <- pilot_cors[,c('TMS_pre','acc_pre','dp_pre','dp2_pre','dp2rs_pre','H_M_ratio_pre')]
rownames(pilot_cors) <- pnames[gsub("_pre","",rownames(pilot_cors))]
colnames(pilot_cors) <- pnames[gsub("_pre","",colnames(pilot_cors))]
write.csv(pilot_cors,"csv_files/pilot_correlations.csv")
############################################################## Presentation ----
# Format results ------
shortlist_group <-c("TMS","SOW_F","SOW_S","PHQ","GAD","RRS","WBSI","acc","adj_dp","adj_dp2rs","adj_meta_d","H_M_ratio")
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

res_pilot_f <- formatResults(res_pilot)
write.csv(res_pilot_f,'csv_files/pilot_results.csv')

# Bar graphs -----------
dir.create('graphs', showWarnings = FALSE)
shorts <- unique(c(shortlist,shortlist_group))
colours <- turbo(length(shorts))
names(colours) <- shorts
pconds <- c(control='Waitlist',world='World',mental='Mental States')

# plot diffs
jpeg('graphs/pilot_group_means.jpeg', height=150, width=200,units='mm',quality=200,res=300)
plotDiff(pilot_data)
dev.off()

jpeg('graphs/pilot_interactions.jpeg', height=150, width=200,units='mm',quality=200,res=300)
plotRes(res_pilot,'Time*Condition Interactions (Per-Protocol)')
dev.off()

################################################################## Save env ----
save.image(file='BCTM_pilot.RData')
