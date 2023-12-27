#SETUP#######################---\\\\\\\\
# Packages----
list.of.packages <- c("tidyverse", "dplyr","magrittr","eeptools","mice","MASS",
                      "aod","psych","papaja","emmeans","bfrr","bayesplay","readxl",
                      "foreign") 
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
  #install.packages("devtools")
  #devtools::install_github("debruine/bfrr")
  #devtools::install_github("crsh/papaja")

# Read in data----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
raw_data <- read.csv("MSvW_raw_data.csv")
  #str(raw_data, list.len=ncol(raw_data))

# Reorder and remove raw variables ---------------------------
#calculate age
days0x <- paste0(0,raw_data$Q3_pre)
days0 <- ifelse(raw_data$Q3_pre < 10,days0x,raw_data$Q3_pre)
months0x <- paste0(0,raw_data$Q4_pre)
months0 <- ifelse(raw_data$Q4_pre < 10,months0x,raw_data$Q4_pre)
dobx <- ifelse(!is.na(raw_data$Q5_pre),paste0(raw_data$Q5_pre,"-",months0,"-",days0),NA)
dob <- as.Date(dobx)
exp_date <- as.Date(substr(raw_data$StartDate_pre, 1,10))
ages <- age_calc(na.omit(dob), enddate =  na.omit(exp_date), units = "years", precise = F)
raw_data$age[!is.na(dob)] <- ages

#rename factor levels
raw_data$Condition <- as.factor(raw_data$Condition)
levels(raw_data$Condition)[levels(raw_data$Condition)=="mental states"|
                           levels(raw_data$Condition)=="mind"] <- "mental"
levels(raw_data$Condition)[levels(raw_data$Condition)==""] <- NA

raw_data$Q6_pre <- as.factor(raw_data$Q6_pre)
levels(raw_data$Q6_pre)[levels(raw_data$Q6_pre)=="cis female"|
                        levels(raw_data$Q6_pre)=="f"|
                          levels(raw_data$Q6_pre)=="female "|
                          levels(raw_data$Q6_pre)=="woman"|
                          levels(raw_data$Q6_pre)=="fem"] <- "female"
levels(raw_data$Q6_pre)[levels(raw_data$Q6_pre)=="male "|
                          levels(raw_data$Q6_pre)=="m"] <- "male"

#vars for missing data on entire pre or post-test
raw_data$URN_pre <- ifelse(is.na(raw_data$URN_pre),0,1)
raw_data$URN_pre <- as.factor(raw_data$URN_pre)
raw_data$URN_post <- ifelse(is.na(raw_data$URN_post),0,1)
raw_data$URN_post <- as.factor(raw_data$URN_post)


#remove and reorder variables
raw_data <- raw_data[,!colSums(is.na(raw_data)) == nrow(raw_data)]

raw_data <- raw_data[, 
  c(which(colnames(raw_data) == "Condition"),
  which(colnames(raw_data) == "age"),
  which(colnames(raw_data) == "Q6_pre"),
  which(colnames(raw_data) == "X1a"):which(colnames(raw_data) == "X1b"),
  which(colnames(raw_data) == "URN_pre"),
  which(colnames(raw_data) == "URN_post"),
  which(colnames(raw_data) == "Q7_pre"):which(colnames(raw_data) == "Q14_18_pre"),
  which(colnames(raw_data) == "Q6_post"):which(colnames(raw_data) == "Q13_18_post"))]

#rename variables
names(raw_data)[which(names(raw_data) == "Condition"):which(names(raw_data) == "URN_post")] <- 
  c("condition","age","gender","exp_anx","exp_dep","missing_pre","missing_post")

raw_data <- raw_data %>%
  rename_with(~gsub("_a_", "_", .)) %>%
  rename_with(~gsub("Q7|Q6","AWH",.), grep("(^Q7_.*pre$)|(^Q6_.*post$)",names(.))) %>%
  rename_with(~gsub("Q8|Q7","AMH",.), grep("(^Q8_.*pre$)|(^Q7_.*post$)",names(.))) %>%
  rename_with(~gsub("Q9|Q8","FFMQ",.), grep("(^Q9_.*_pre$)|(^Q8_.*_post$)",names(.))) %>%
  rename_with(~gsub("Q10|Q9","PSS",.), grep("(^Q10_.*_pre$)|(^Q9_.*_post$)",names(.))) %>%
  rename_with(~gsub("Q11|Q10","PHQ",.), grep("(^Q11_.*_pre$)|(^Q10_.*_post$)",names(.))) %>%
  rename_with(~gsub("Q12|Q11","RRS",.), grep("(^Q12_.*_pre$)|(^Q11_.*_post$)",names(.))) %>%
  rename_with(~gsub("Q13|Q12","WBSI",.), grep("(^Q13_.*_pre$)|(^Q12_.*_post$)",names(.))) %>%
  rename_with(~gsub("Q14|Q13","PWB",.), grep("(^Q14_.*_pre$)|(^Q13_.*_post$)",names(.)))


# count missing data---------------------------
#select pre-test and post-test variables
pre_range <- c(which(colnames(raw_data) == "AWH_pre"):which(colnames(raw_data) == "PWB_18_pre"))
post_range <- c(which(colnames(raw_data) == "AWH_post"):which(colnames(raw_data) == "PWB_18_post"))

m_cond <- sum(is.na(raw_data$condition))
m_some_pre <- table(raw_data[raw_data$missing_pre==1 & !complete.cases(raw_data[pre_range]),"condition"])
m_all_pre <- table(raw_data[raw_data$missing_pre==0,"condition"])
m_some_post <- table(raw_data[raw_data$missing_post==1 & !complete.cases(raw_data[post_range]),"condition"])
m_all_post <- table(raw_data[raw_data$missing_post==0,"condition"])

missing <- rbind(m_some_pre,m_all_pre,m_some_post,m_all_post)
row.names(missing) <- c("Some Pre-test","Entire Pre-test","Some Post-test","Entire Post-test")

missing_sums <- colSums(is.na(raw_data))
#graphs of missing data pattern. commented out for code speed.
#jpeg("pre_listwise_missing_pattern.jpg", width = 15000, height = 15000)
#md.pattern(raw_data)
#dev.off()

#write.csv(raw_data, file="all_pre_listwise.csv")
# Round2 ----
round2 <- function(x, digits = 2) {  # Function to always round 0.5 up
  posneg <- sign(x)
  z <- abs(x) * 10^digits
  z <- z + 0.5
  z <- trunc(z)
  z <- z / 10^digits
  z * posneg
}
#INITIAL DATA CHANGES#######################
# Remove missing pre-test data (listwise deletion) ---------------------------
#remove missing data we don't plan to impute (i.e. if not whole just post-test missing)
all_data <- raw_data 
incomplete_index_post <- raw_data #create duplicate dataset
incomplete_index_post[raw_data$missing_post==0,post_range] <- incomplete_index_post[raw_data$missing_post==0,pre_range] #Dummy LOCF

raw_data <- raw_data[
                    !is.na(raw_data$condition)
                    & raw_data$missing_pre == 1 #people not missing pre-test
                    & complete.cases(raw_data[pre_range])
                    #& complete.cases(incomplete_index_post[post_range]) #REMOVES PARTICIPANTS MISSING SOME BUT NOT ALL POST TEST DATA
                    ,!(names(raw_data) %in% "missing_pre")
                    ]

raw_data$gender <- droplevels(raw_data$gender)

sum(raw_data$missing_post==0)/length(raw_data$missing_post==0)*100 #~29% of entire post-test sessions are missing
# Scale Re-numbering ---------------------------
scale_recode <- function(dataset, scale_name, q_nums, scale_max=0) {
  q_nums <- as.character(q_nums)
  q_names <- c(paste0(scale_name,"_",q_nums,"_pre"), paste0(scale_name,"_",q_nums,"_post"))
  if(scale_max == 0){ dataset[q_names] <- dataset[q_names] - 1
  } else {dataset[q_names] <- (scale_max+1) - dataset[q_names]} #scale_max > 0 for reverse code
  return(dataset)
}

raw_data <- scale_recode(raw_data,"PSS",c(1:4,9))
raw_data <- scale_recode(raw_data,"PHQ",c(1:4))
raw_data <- scale_recode(raw_data,"FFMQ",c(4,5,7,8,11,12,14,17,19,22:24), 5)
raw_data <- scale_recode(raw_data,"PSS",c(4,5,7,8), 4)
raw_data <- scale_recode(raw_data,"PWB", c(2,4,6,10,11,14,16), 6)
#write.csv(raw_data, file="raw_data.csv")
#DATA AGGREGATING#######################
# Scale List ---------------------------  
agg_data <- raw_data

scales_raw <- list(AMH=c(""),AWH=c(""),
            FFMQ_OB=c(6,10,15,20),FFMQ_NR=c(3,9,13,18,21),FFMQ_DS=c(1,2,5,11,16),
            FFMQ_NJ=c(4,7,14,19,24),FFMQ_AA=c(8,12,17,22,23),
            PHQ_anx=c(1,2),PHQ_dep=c(3,4),
            PWB_env=c(1,2,4,7,11,12,13,16,18),PWB_acc=c(3,5,6,8,9,10,14,15,17),
            RRS_D=c(1:4,6,8,9,14,17:19,22),RRS_B=c(5,10,13,15,16),RRS_R=c(7,11,12,20,21),
            WBSI_sup=c(13,1,8,10,11,14),WBSI_int=c(2,3,4,5,6,7,9,12,15),
            PSS=c(1:10), RRS=c(1:22), WBSI=c(1:15), PHQ=c(1:4), FFMQ=c(1:24))

expandScales <- function(scale_name, q_nums){
  scale_name <- sub("_.*", "", scale_name)
    scale_name <- sub("OB|NR|DS|NJ|AA.*", "FFMQ", scale_name)
    q_nums <- as.character(q_nums)
  pre_scale <- paste(scale_name,q_nums,"pre",sep="_")
    post_scale <- paste(scale_name,q_nums,"post",sep="_")
  pre_scale <- gsub('__','_',pre_scale)
    post_scale <- gsub('__','_',post_scale)
  scale <- list(pre = pre_scale, post = post_scale)
  return(scale)
}

scales <- Map(expandScales, names(scales_raw), scales_raw)
# Mean Scores & scale names list ---------------------------
agg_data$exp_mean <- rowMeans(agg_data[c("exp_anx","exp_dep")])

scaleMeans <- function(dataset){
  for(i in 1:length(scales)){
    scale <- names(scales[i])
    if(grepl("^A.H",scale)){next}
    for(j in 1:2){
      measure <- dataset[,scales[[i]][[j]]]
      varname <- paste(scale,names(scales[[i]][j]),sep='_')
      dataset <- dplyr::mutate(dataset,!!varname := rowMeans(measure))
      dataset[,varname] <- as.numeric(dataset[,varname])
    }
  }
  return(dataset)
}

agg_data <- scaleMeans(agg_data)

#####ANALYSIS####################
#H1 estimates#######################
# Scale Measures-----
scales_dat <- read_excel("previous_studies/MSvW_prior_effects.xlsx", sheet = "scales_data")
scales_dat <- scales_dat[,!(names(scales_dat) %in% c("alt_control","sd"))]
scales_dat[scales_dat==""]<-NA
scales_dat <- scales_dat %>% fill(everything(), .direction = "down")
unique(scales_dat$study)

wide_scales <- scales_dat %>% pivot_wider(
  id_cols = c(study, measure, scale_range, num_q),
  names_from = c(condition, time), 
  values_from = m) 

short_scales <- wide_scales[!(wide_scales$measure == "Scale of Positive and Negative Experience - affect balance" |
                                wide_scales$measure == "Brief Irritability Test" |
                                wide_scales$measure == "Brief Resilience Scale")
                            ,] #SPANE subtracted one subscale from another. 
#BITe and BRS already averaged and coding this out was difficult, attempt below:

#adjusted_scales <- short_scales %>%
#  mutate(scale_size = as.numeric(substr(scale_range, 3,3)) - as.numeric(substr(scale_range, 1,1)),
#         across(Experimental_Pre:Control_Post, ~(((.x/num_q)/scale_size)*4)),
#         #attempt to change analysis for only some measures
#         #across(Experimental_Pre:Control_Post, ~ifelse(measure == "Brief Irritability Test"|| measure == "Brief Resilience Scale",
#         #(((.x/scale_size)*4)),(((.x/num_q)/scale_size)*4)))
#         adj_diff = abs((Experimental_Post-Experimental_Pre)-(Control_Post-Control_Pre))) %>%
#  select(-scale_range,-scale_size,-c(Experimental_Pre:Control_Post))
#
#adjusted_scales[adjusted_scales$study == "Cavanagh et al 2018" | adjusted_scales$study == "Cavanagh et al 2013",]$adj_diff %>% sd
#scales_prior <- round(mean(adjusted_scales$adj_diff),2)

# Meta-d'----

# schmidt study----

#schmidt from raw data
#convert for the raw data files from the OSF directory
  #txt_files <- list.files("schmidt_data/Data/R",".txt$")
  #file_names <- gsub(".txt$","",txt_files)
  #schmitR <- vector(mode = "list", length = length(txt_files))
  #names(schmitR) <- file_names
  #for(file in 1:length(txt_files)){
  #  txt_file <- read.delim(paste0("schmidt_data/Data/R/",txt_files[file]), header=TRUE, sep=" ")
  #  schmitR[[file]] <- txt_file #assign(txt_files[file],txt_file)
  #  #write.csv(txt_file, file=paste0(file_names[file],".csv"))
  #}
#
#summary data, pre-post separately
  #MFN <- read.spss("schmidt_data/Data/SPSS/MFN_.sav", to.data.frame=TRUE)
  #write.csv(MFN, file="MFN.csv")
#change data over time
  #DeltaBases <- read.spss("schmidt_data/Data/SPSS/DeltaBases.sav", to.data.frame=TRUE)
  #write.csv(DeltaBases, file="DeltaBases.csv")

DeltaBases <- read.csv("previous_studies/DeltaBases.csv")
deltas <- DeltaBases[DeltaBases$Task_ == "gab",c("Grupo_","deltaMeta_Dprime","mdRatio_Delta21")]

#meta-d'
schmidt_metad_m <- mean(deltas[deltas$Grupo_=="in","deltaMeta_Dprime"])-mean(deltas[deltas$Grupo_=="ex","deltaMeta_Dprime"])
metad_in_se <- sd(deltas[deltas$Grupo_=="in","deltaMeta_Dprime"])/sqrt(sum(deltas$Grupo_=="in"))
metad_ex_se <- sd(deltas[deltas$Grupo_=="ex","deltaMeta_Dprime"])/sqrt(sum(deltas$Grupo_=="ex"))
schmidt_metad_se <- sqrt(metad_in_se^2 + metad_ex_se^2)
schmidt_metad <- data.frame(schmidt_metad_m,schmidt_metad_se)
names(schmidt_metad) <- c("m_diff","se_diff")

#meta-d'/d'
schmidt_mratio_m <- mean(deltas[deltas$Grupo_=="in","mdRatio_Delta21"])-mean(deltas[deltas$Grupo_=="ex","mdRatio_Delta21"])
mratio_in_se <- sd(deltas[deltas$Grupo_=="in","mdRatio_Delta21"])/sqrt(sum(deltas$Grupo_=="in"))
mratio_ex_se <- sd(deltas[deltas$Grupo_=="ex","mdRatio_Delta21"])/sqrt(sum(deltas$Grupo_=="ex"))
schmidt_mratio_se <- sqrt(mratio_in_se^2 + mratio_ex_se^2)
schmidt_mratio <- data.frame(schmidt_mratio_m,schmidt_mratio_se)
names(schmidt_mratio) <- c("m_diff","se_diff")

#log(meta-d'/d')
schmidtp <- read.csv("previous_studies/MFN.csv")
schmidtp <- schmidtp[schmidtp$Task_ == "gab",c("Grupo_","Session__","Id","mdRatio")]
schmidtp$mdRatio <- log(schmidtp$mdRatio)
schmidtp <- pivot_wider(schmidtp,id_cols=c("Id","Grupo_"),names_from="Session__",values_from="mdRatio")
schmidtp$diff <- schmidtp$post-schmidtp$`pre `
schmidt_ex <- schmidtp[schmidtp$Grupo_=="ex","diff"]
schmidt_in <- schmidtp[schmidtp$Grupo_=="in","diff"]
schmidt_log_m <- mean(schmidt_in$diff) - mean(schmidt_ex$diff)
log_ex_se <- sd(schmidt_ex$diff)/sqrt(length(schmidt_ex$diff))
log_in_se <- sd(schmidt_in$diff)/sqrt(length(schmidt_in$diff))
schmidt_log_se <- sqrt(log_ex_se^2+log_in_se^2)
schmidt_log <- data.frame(schmidt_log_m,schmidt_log_se)
names(schmidt_log) <- c("m_diff","se_diff")

# schmidt adjusted meta-d' means----
#data extraction
schmidt_MFN <- read.csv("previous_studies/MFN.csv")
schmidt_MFN <- schmidt_MFN[schmidt_MFN$Task_ == "gab",c("Id","Session__", "Task_", "Grupo_", "dprime","metaDprime")]
schmidt_MFN <- pivot_wider(schmidt_MFN,id_cols=c("Id","Grupo_"),names_from="Session__",values_from=c("metaDprime","dprime"))
schmidt_MFN$metad_diff <- schmidt_MFN$metaDprime_post-schmidt_MFN$`metaDprime_pre `
schmidt_MFN$d_diff <- schmidt_MFN$dprime_post-schmidt_MFN$`dprime_pre `
schmidt_MFN$Grupo_ <- as.factor(schmidt_MFN$Grupo_)
schmidt_MFN <- schmidt_MFN[,c("Grupo_","metad_diff","d_diff")]

#with emmeans package
schmidt_lm <- lm(metad_diff ~ Grupo_ + d_diff, data=schmidt_MFN)
schmidt_em <- emmeans(schmidt_lm,"Grupo_")
#check by by hand as emmeans may not be what is expected
schmidt_sum <- summary(schmidt_lm)
in_adj  <- schmidt_sum$coefficients[1,1] + schmidt_sum$coefficients[2,1] + 
          (schmidt_sum$coefficients[3,1]*mean(schmidt_MFN$d_diff))
ex_adj  <- schmidt_sum$coefficients[1,1] + 
          (schmidt_sum$coefficients[3,1]*mean(schmidt_MFN$d_diff))
#original means
in_orig <- mean(schmidt_MFN$metad_diff[schmidt_MFN$Grupo_=="in"])
ex_orig <- mean(schmidt_MFN$metad_diff[schmidt_MFN$Grupo_=="ex"])

schmidt_adj <- in_adj-ex_adj

# Carpenter study -----
#data extracted in MATLAB code
carpenter <- read.csv("previous_studies/carpenter_diffs.csv")
carpenter$m_diff <- carpenter$exp_mdiff-carpenter$control_mdiff
carpenter$se_diff <- sqrt((carpenter$control_se ^2)+(carpenter$exp_se ^2))
  #write.csv(carpenter,file="carpenter.csv")

#schmidt_2 <- read.csv("previous_studies/MFN.csv")
#schmidt_2$dprime
#carpenter_raw$s10_d %>% mean

#adjusted meta-d' means
carpenter_raw <- read.csv("previous_studies/carpenter_raw.csv")
carpenter_raw$diff_d <- carpenter_raw$s10_d-carpenter_raw$s1_d     
carpenter_raw$diff_metad <- carpenter_raw$s10_metad-carpenter_raw$s1_metad

carpenter_lm <- lm(diff_metad ~ condition + diff_d, data=carpenter_raw)
carpenter_sum <- summary(carpenter_lm)
exp_adj  <- carpenter_sum$coefficients[1,1] + carpenter_sum$coefficients[2,1] + 
            (carpenter_sum$coefficients[3,1]*mean(carpenter_raw$diff_d))
cont_adj  <- carpenter_sum$coefficients[1,1] + 
            (carpenter_sum$coefficients[3,1]*mean(carpenter_raw$diff_d))
carpenter_adj <- exp_adj-cont_adj

#carpenter mratio and log by hand to check
carpenter_raw2 <- read.csv("previous_studies/carpenter_raw.csv")
carpenter_raw2$metad_diff <- carpenter_raw2$s10_metad - carpenter_raw2$s1_metad
carpenter_raw2$mratio_pre <- carpenter_raw2$s1_metad / carpenter_raw2$s1_d
carpenter_raw2$mratio_post <- carpenter_raw2$s10_metad / carpenter_raw2$s10_d 
carpenter_raw2$mratio_diff <- carpenter_raw2$mratio_post - carpenter_raw2$mratio_pre
carpenter_raw2$logmratio_pre <- log(carpenter_raw2$mratio_pre)
carpenter_raw2$logmratio_post <- log(carpenter_raw2$mratio_post)
carpenter_raw2$logmratio_diff <- carpenter_raw2$logmratio_post - carpenter_raw2$logmratio_pre
mean(carpenter_raw2$metad_diff[carpenter_raw2$condition == "exp"])-
  mean(carpenter_raw2$metad_diff[carpenter_raw2$condition == "control"])
mean(carpenter_raw2$mratio_diff[carpenter_raw2$condition == "exp"])-
  mean(carpenter_raw2$mratio_diff[carpenter_raw2$condition == "control"])
mean(carpenter_raw2$logmratio_diff[carpenter_raw2$condition == "exp"])-
  mean(carpenter_raw2$logmratio_diff[carpenter_raw2$condition == "control"])

hist(carpenter_raw2$metad_diff, breaks=20)
hist(carpenter_raw2$mratio_diff, breaks=20)
hist(carpenter_raw2$logmratio_diff, breaks=20)

plot(density(carpenter_raw2$metad_diff))
plot(density(carpenter_raw2$mratio_diff))
plot(density(carpenter_raw2$logmratio_diff))


# combined Schmidt and Carpenter----
#sample sizes.
  #Carpenter control:12+20=32 exp:17+12=29 total=61
  #Schmidt control:13 exp:14 total=27
m_both <- round(mean(c(carpenter$m_diff[1],schmidt_metad$m_diff)),1)
schmidt_metad_se_adj <- schmidt_metad$se_diff/(sqrt(27)/sqrt(61))
se_both <- round(sqrt(mean(schmidt_metad_se_adj^2,carpenter$se_diff[1]^2)),2)
m_both_adj <- round(mean(c(carpenter_adj,schmidt_adj)),1)

m_both_ratio <- round(mean(c(carpenter$m_diff[2],schmidt_mratio$m_diff)),1)
schmidt_mratio_se_adj <- schmidt_mratio$se_diff/(sqrt(27)/sqrt(61))
se_both_ratio <- round(sqrt(mean(schmidt_mratio_se_adj^2,carpenter$se_diff[2]^2)),2)

m_both_log <- round(mean(c(carpenter$m_diff[3],schmidt_log$m_diff)),1)
schmidt_log_se_adj <- schmidt_log$se_diff/(sqrt(27)/sqrt(61))
se_both_log <- round(sqrt(mean(schmidt_log_se_adj^2,carpenter$se_diff[3]^2)),2)



#table
both_tab <- carpenter[c("row_names","m_diff","se_diff")]
names(both_tab) <- c("row_names","carpenter_m","carpenter_se")
both_tab$schmidt_m <- c(schmidt_metad$m_diff, schmidt_mratio$m_diff, schmidt_log$m_diff)
both_tab$schmidt_se <- c(schmidt_metad$se_diff, schmidt_mratio$se_diff, schmidt_log$se_diff)
both_tab$schmidt_se_adj <- c(schmidt_metad_se_adj,schmidt_mratio_se_adj, schmidt_log_se_adj)
both_tab$combined_m <- c(m_both,m_both_ratio,m_both_log)
both_tab$combined_se <- c(se_both,se_both_ratio,se_both_log)

#write.csv(both_tab,"metad_priors.csv")
#MULTIPLE IMPUTATION#######################
# Auxiliary variables & missingness predictors ----
 pre_vars <- names(scales[1:which(names(scales) == "PSS")])

  m_count_cond <- xtabs(~condition+missing_post, data = agg_data)
  m_count_gender <- xtabs(~gender+missing_post, data = agg_data)
  #xtabs(~missing_post+age, data = agg_data)
#check if auxiliary and pre-test variables are predictive of missingness: http://stats.idre.ucla.edu/r/dae/logit-regression/
aux_logit <- glm(missing_post ~ condition+age+gender, data = agg_data, family = "binomial")
  summary(aux_logit)
  aux_ORs2 <- exp(cbind(OR = coef(aux_logit), confint(aux_logit)))
    #wald <- wald.test(b = coef(mylogit), Sigma = vcov(mylogit), Terms = 2:3)
#only condition is significant, analyse with just condition
  cond_logit <- glm(missing_post ~ condition, data = agg_data, family = "binomial")
  aux_ORs <- exp(cbind(OR = coef(cond_logit), confint(cond_logit)))
  log_tab <- with(agg_data, data.frame(condition = factor(levels(condition))))
  log_tab$cond_pred <- predict(cond_logit, newdata = log_tab, type = "response")
  wald.test(b = coef(cond_logit), Sigma = vcov(cond_logit), L = cbind(0,1,-1))
#check pre-test variables on missingness
  eq <- as.formula(paste("missing_post~",paste(pre_vars,collapse="+")))
  pre_logit <- glm(eq, data = agg_data, family = "binomial")
  summary(pre_logit)
  
# Scatter plots----

scatterPlot <- function(var){
  #jpeg(paste0(var,"_scatters.jpeg"), width = 2000, height = 2000)
  post <- paste0(var,"_post")
    #scatter plots
  par(mfrow=c(3,7))
  for(scale in 1:length(scales)){
    pre <- paste0(names(scales[scale]),"_pre")
    scat_data <- agg_data[!is.na(agg_data[,post]),c(post,pre)]
    plot(scat_data[,post], scat_data[,pre], main=pre, xlab=post, ylab=pre)
    fit <- lm(get(post) ~ get(pre), data = scat_data)
    abline(fit, col="red")
  }
  #dev.off()
}
#lapply(names(scales),scatterPlot)

# Select both scale and item level predictors for each scale----
# Setup and run multiple imputation----
#Hybrid Non-AIC
  subscales <- 3:which(names(scales) == "WBSI_int")
  subscales_pre <- paste0(names(scales[subscales]),"_pre")
  full_scales_pre <- paste0(names(scales[-subscales]),"_pre")
  
scalePreds <- function(var){
  #var <- names(scales)[2]
  item_level <- scales[var][[1]][[1]]
  scale_name <- sub("_.*", "", var)
  subscale_level <- subscales_pre[!grepl(var, subscales_pre)]
  scale_level <- full_scales_pre[grepl("AMH_pre|AWH_pre|PSS_pre", full_scales_pre) & 
                                   !grepl(var, full_scales_pre)]
  #hybrid
    #subscale_level <- subscales_pre[grepl(scale_name, subscales_pre) & !grepl(var, subscales_pre)]
    #scale_level <- full_scales_pre[!grepl(scale_name, full_scales_pre)]
  scale_preds <- list()
  scale_preds <- c(item_level,subscale_level,scale_level,"condition")
  return(scale_preds)
}

sc_len <- which(names(scales) == "PSS")
scale_preds <- lapply(names(scales[1:sc_len]), scalePreds)
  names(scale_preds) <- names(scales[1:sc_len])

#change the matrix of predictor-outcome analyses
pred_mat <- quickpred(agg_data)
pred_mat[,] <- 0
#original version
for(i in seq_along(scale_preds)){
  predictors <- scale_preds[[i]]
  outcomes <- scales[[i]][[2]]
  pred_mat[outcomes,predictors] <- 1
}

dataset[which(names(dataset)==""):which(names(dataset)=="")]

exp_preds <- c("AMH_pre","AWH_pre","PSS_pre",subscales_pre[!grepl("PHQ", subscales_pre)], "condition")
  exp_anx_pred <- c("PHQ_1_pre","PHQ_2_pre","PHQ_dep_pre",exp_preds)
    pred_mat["exp_anx",exp_anx_pred] <- 1
  exp_dep_pred <- c("PHQ_3_pre","PHQ_4_pre","PHQ_anx_pre",exp_preds)
    pred_mat["exp_dep",exp_dep_pred] <- 1
  exp_mean_pred <- c("PHQ_1_pre","PHQ_2_pre","PHQ_3_pre","PHQ_4_pre",exp_preds)
    pred_mat["exp_mean",exp_mean_pred] <- 1

#list of methods
meths <- rep('',length(names(agg_data)))
  names(meths) <- names(agg_data)
  post_vars <- c("exp_anx","exp_dep","exp_mean",names(raw_data[,grepl("[0-9|A.H]_post",names(raw_data))]))
  meths[post_vars] <- 'norm'

#run imputations (this takes a while with a good computer, lower m and maxit to test)
imp <- mice(agg_data,m=100,maxit=30,meth=meths,pred=pred_mat,print=F,seed=123)
#original data: imp$data ; imputed data: imp$imp[outcomes] ; first complete dataset: complete(imp,1)
#write.csv(imp$imp[outcomes],file="imp.csv")
  #CHECK MAX VALS: #as.character(imp_new$loggedEvents[, "out"]) #x <- c() #for(i in imp_new$imp[post_vars]){for(j in i){x <- c(x,max(j))}} #max(x)
comps <- mice::complete(imp,"all")
  #write.csv(comps$`1`,file="comp.csv")
# Scale means and difference scores----
scaleDiffs <- function(dataset){
  for(i in 1:length(scales)){
    scale <- names(scales[i])
    scale_pre <- paste0(scale,"_pre")
    scale_post <- paste0(scale,"_post")
    varname <- paste0(scale,'_diff')
    dataset[[varname]] <- with(dataset, eval(parse(text=scale_post))-eval(parse(text=scale_pre)))
  }
  return(dataset)
}

comps_agg <- map(comps, scaleMeans)
comps_agg <- map(comps_agg, scaleDiffs)

#DESCRIPTIVES#######################
# Demographics----
descTable <- function(dataset){
  #setup
    descriptives <- as.data.frame(matrix(ncol=1,nrow=7))
    colnames(descriptives) <- c("n")
  #categorical
    counts <- c(summary(dataset$condition),summary(dataset$gender))
    descriptives[1:length(counts),] <- counts
    rownames(descriptives)[1:length(counts)] <- names(counts)
  #age
    i <- c(length(counts)+1,length(counts)+2)
    rownames(descriptives)[i] <- c("Age: Mean (SD)","Age Range")
    age <- round(summary(dataset$age))
    descriptives[i,] <- c(paste0(age["Mean"]," (",round(sd(dataset$age),2),")"),
                          paste(age["Min."],age["Max."],sep="-"))
  return(descriptives)
}
descriptives <- descTable(raw_data)

# Cronbach's alpha----
cronAlpha <- function(dataset){
  scales <- scales[!grepl("^A.H",names(scales))]
  cronbach <- as.data.frame(matrix(ncol=2,nrow=length(scales)))
  colnames(cronbach) <- c("Pre","Post")
  rownames(cronbach) <- names(scales)
  for(scalei in 1:length(scales)){
    scale <- scales[[scalei]]
    cronbach[scalei,1] <- psych::alpha(dataset[,scale$pre])$total$raw_alpha
    cronbach[scalei,2] <- psych::alpha(dataset[,scale$post])$total$raw_alpha
  }
  return(cronbach)
}
cron_alpha <- map(comps_agg, cronAlpha)
cron_reduce <- Reduce("+", cron_alpha)/length(cron_alpha)
cron_reduce <- round2(cron_reduce,2)

# Conditional cell means----
cellMeans <- function(dataset){
  n_stats <- length(c("mean","SD","SE"))
  n_times <- length(c("pre","post"))
  n_cond <- length(levels(agg_data$condition))
  cell_means <- as.data.frame(matrix(ncol=n_cond*n_times*n_stats,nrow=length(scales))) 
  cols <- expand.grid(c("m","sd","se"),c("pre","post"),levels(agg_data$condition))
  colnames(cell_means) <- paste(cols$Var1,cols$Var2,cols$Var3,sep="_")
  rownames(cell_means) <- c(names(scales))
  for(i in 1:length(scales)){
    for(group in levels(dataset$condition)){
      cond_log <- dataset$condition==group
      for(test_time in c("pre","post")){
        scale <- names(scales[i])
        time <- paste0(scale,"_",test_time)
        time_group_col <- grep(paste0(test_time,"_",group),names(cell_means))
        cell_means[i,time_group_col] <- c(mean(dataset[cond_log,time],na.rm=T), 
                                          sd(dataset[cond_log,time],na.rm=T),
                                          sd(dataset[cond_log,time],na.rm=T)/sqrt(sum(cond_log))
                                          )
      }
    }
  }
  cell_means <- round(cell_means,2)
  return(cell_means)
}


cell_means <- map(comps_agg, cellMeans)
means_reduce <- Reduce("+", cell_means)/length(cell_means)
means_reduce <- round2(means_reduce,2)
rownames(means_reduce) <- c(rownames(means_reduce)[1:7],"PHQ Anxiety","PHQ Depression",
                           "PWB Environmental Mastery", "PWB Self-Acceptance",
                           "RRS Depression","RRS Brooding","RRS Reflection",
                           "WBSI Supression","WBSI Intrusion","PSS-10", "RRS","WBSI","PHQ-4","FFMQ")

means_final <- means_reduce
means_final$m_pre_control <-paste0(means_final$m_pre_control," (",means_final$sd_pre_control,")")
means_final$m_post_control <-paste0(means_final$m_post_control," (",means_final$sd_post_control,")")
means_final$m_pre_mental <-paste0(means_final$m_pre_mental," (",means_final$sd_pre_mental,")")
means_final$m_post_mental <-paste0(means_final$m_post_mental," (",means_final$sd_post_mental,")")
means_final$m_pre_world <-paste0(means_final$m_pre_world," (",means_final$sd_pre_world,")")
means_final$m_post_world <-paste0(means_final$m_post_world," (",means_final$sd_post_world,")")
means_final <- means_final[,!grepl("^s",names(means_final))]


# Bargraphs----
barGraph <- function(means_table){
  mean_cols <- grepl("^m_",names(means_table))
  mean_tab <- means_table[mean_cols]
  tibby <- data.frame(Control = as.numeric(mean_tab[,grepl("control$",names(mean_tab))]),
                      Mental = as.numeric(mean_tab[,grepl("mental$",names(mean_tab))]),
                      World = as.numeric(mean_tab[,grepl("world$",names(mean_tab))]),
                      row.names = c("Pre","Post"))
  tibby <- as.matrix(tibby)
  
  se_cols <- grepl("^se_",names(means_table))
  CI_lo <- as.numeric(means_table[mean_cols]-(means_table[se_cols]*1.96))
  CI_up <- as.numeric(means_table[mean_cols]+(means_table[se_cols]*1.96))
  
  tab_plot <- barplot(tibby, main = row.names(means_table),
                      xlab = "Condition", ylab = "Score +/-CI95%", ylim = c(0, max(CI_up)+1),
                      beside=T, legend.text=rownames(tibby),
                      args.legend=list(x="topright", inset=c(-0.1, -0.15)))
  co <- seq(ncol(means_table[mean_cols]-1))
  arrows(tab_plot[co],CI_lo[co],tab_plot[co],CI_up[co], angle=90, code=3, length=.08)

  return(tab_plot)
}
#INFERENTIALS#######################
# Expectations ---------------------------
expect <-function(dataset){
  #data setup
  dataset <- dataset[!is.na(dataset$exp_mean) & 
                    (dataset$condition=='mental'|
                    dataset$condition=='world'),] #Note: data loss: these were collected with emails
  mental <- dataset$condition=='mental'
  world <- dataset$condition=='world'
  #dataset$condition <- droplevels(dataset$condition)
  exp_vars <- names(dataset)[grep("^exp_",names(dataset))]
  #table setup
  exp_tab <- as.data.frame(matrix(ncol=6,nrow=length(exp_vars))) 
  colnames(exp_tab) <- c("m_diff","SE","B","RR","t","p")
  rownames(exp_tab) <- c(exp_vars)
  #prior calculation
  exp_reg <- summary(lm(PHQ_post ~ PHQ_pre + exp_mean, data=dataset))$coefficients["exp_mean",c("Estimate","Std. Error")]
  up_conf <- exp_reg[[1]]+exp_reg[[2]]*1.96
  dataset$PHQ_diff <- dataset$PHQ_post -  dataset$PHQ_pre
  phq_diff <- mean(dataset[mental,"PHQ_diff"],na.rm=T) - mean(dataset[world,"PHQ_diff"],na.rm=T)
  original_h1 <- (abs(round2(phq_diff,1))/round2(up_conf))
  #room to move heuristic
  r2m <- round((5-mean(dataset[world,"exp_mean"]))/2,1) -.1
  #priors
  h0_mod <- prior(family="point", point=0)
  h1_mod <- prior(family="normal", mean=r2m, sd=r2m/2)
  
  #analysis
  for(i in 1:length(exp_vars)){
    var <- exp_vars[i]
    #mean diff, SEs
    exp_tab[i,1] <- round2(mean(dataset[mental,var]) - mean(dataset[world,var]))
    exp_tab[i,2] <- round2(sqrt(
      sd(dataset[mental,var])/sqrt(nrow(dataset[mental,])) ^2 +
        sd(dataset[world,var])/sqrt(nrow(dataset[world,])) ^2))
    #bayes
    data_mod <- bayesplay::likelihood("normal", mean=exp_tab[i,1], sd=exp_tab[i,2])
    exp_tab[i,3] <- round2(integral(data_mod * h1_mod)/integral(data_mod * h0_mod))
    #Robustness Regions
    exp_rr <-  bfrr(exp_tab[i,1],exp_tab[i,2],nrow(dataset[mental,])+nrow(dataset[world,])-2,
                    model="normal", mean=r2m, sd=r2m/2, tail=2, criterion=3, 
                    rr_interval=list(mean=c(exp_tab[i,1],exp_tab[i,1]),sd=c(0, 11)), precision = 0.01)
    exp_tab[i,4] <- paste0("[",exp_rr$RR$sd[1],"-",exp_rr$RR$sd[2],"]")
    #t-test
    var_t <- t.test(get(var) ~ condition, data=dataset,var.equal=T)[c("statistic","p.value")]
    exp_tab[i,5] <- round2(var_t[[1]])
    exp_tab[i,6] <- round2(var_t[[2]],3)
    print(data.frame(mental = paste(round2(mean(dataset[mental,var])),round2(sd(dataset[mental,var]))),
                     world =  paste(round2(mean(dataset[world,var])),round2(sd(dataset[world,var])))))
  }

  cor.test(dataset$exp_anx, dataset$exp_dep, method = "pearson", conf.level = 0.95)
  nrow(dataset[mental,])
  nrow(dataset[world,])
  return(exp_tab)
}
expect_table <- expect(agg_data)

# Main analysis----
BFs <- function(dataset){
  #SETUP
    #table setup
      Bf_table <- as.data.frame(matrix(ncol=14,nrow=length(scales))) 
      cols <- expand.grid(c("m_diff","SE","B","RR","RR_up","t","p"),c("active","control"))
      colnames(Bf_table) <- paste(cols$Var1,cols$Var2,sep="_")
      rownames(Bf_table) <- c(names(scales))
      contrasts = list(Main = c(0, -1, 1), Control = c(-1, 1, 0))
    #condition rows
      mental <- dataset$condition == "mental"
      world <- dataset$condition == "world"
      control <- dataset$condition == "control"
    #Bayes models
      h0_mod <- prior(family="point", point=0)
  
  #ANALYSIS FOR EACH VARIABLE
  for(i in 1:length(scales)){
    #SETUP
      scale <- paste0(names(scales[i]),"_diff")
      #SEs
        SE_mental <- sd(dataset[mental,scale])/sqrt(nrow(dataset[mental,]))
        SE_world <- sd(dataset[world,scale])/sqrt(nrow(dataset[world,]))
        SE_control <- sd(dataset[control,scale])/sqrt(nrow(dataset[control,]))
      #bayes setup
        ifelse(grepl("^A.H",scale), H1<-2.5, H1<-.2)
        h1_mod <- prior(family="normal", mean=0, sd=H1, range=c(0, Inf))
      #frequentist setup
        contrast_lm <- lm(get(scale) ~ condition, dataset)
        contrast_mm <- summary(emmeans(contrast_lm, "condition", contr=contrasts)$contrasts)
      #loop setup
        analyses <- list(active=c("mental","world",c(1:7)),control=c("world","control",c(8:14)))

    #TABLE
    for(j in 1:length(analyses)){
      #extract analysis vars
        analysis <- analyses[[j]][1:2]
        col <- as.numeric(analyses[[j]][3:9])
      #mean and SD
        Bf_table[i,col[1]] <- mean(dataset[get(analysis[1]),scale]) - mean(dataset[get(analysis[2]),scale])
        Bf_table[i,col[2]] <- sqrt(SE_mental^2+SE_world^2)
        if(grepl("^A.H",scale)){Bf_table[i,col[1]] <- Bf_table[i,col[1]]*10
                                Bf_table[i,col[2]] <- Bf_table[i,col[2]]*10
                                RR_SDs <- c(0, 25)}
      #Bayes Factor
        active_mean <- Bf_table[i,col[1]]
        if(grepl("^PSS|PHQ|RRS|WBSI",scale)){active_mean <- active_mean*-1}
        bayes_active <- bayesplay::likelihood("normal", mean=active_mean, sd=Bf_table[i,col[2]])
        Bf_table[i,col[3]] <- integral(bayes_active * h1_mod)/integral(bayes_active * h0_mod)
            #Bf_table[i,col[3]] <- Bf(Bf_table[i,col[2]],active_mean,40, "normal", "normal", modeoftheory = 0, 
            #                         scaleoftheory = .2, dftheory = 102, tail = 1)
      #Robustness Regions
        RR_SDs <- c(0, 2)
        if(grepl("^A.H",scale)){RR_SDs <- c(0, 10)}
        active_rr <-  bfrr(active_mean,Bf_table[i,col[2]],nrow(dataset[get(analysis[1]),])+nrow(dataset[get(analysis[2]),])-2,
                           model="normal", mean=0, sd=H1, tail=1, criterion=3, 
                           rr_interval=list(mean=c(active_mean,active_mean),sd=RR_SDs), precision = 0.01)
        Bf_table[i,col[4]] <- active_rr$RR$sd[1]
        Bf_table[i,col[5]] <- active_rr$RR$sd[2]
      #T-test
        Bf_table[i,col[6]] <- contrast_mm$t.ratio[j]
        Bf_table[i,col[7]] <- contrast_mm$p.value[j]
    }
  }
  return(Bf_table)
}
# Final table edits----
BF_tabs <- map(comps_agg, BFs)
  main_reduce <- Reduce("+", BF_tabs)/length(BF_tabs)
  
  main_reduce[,!grepl("^p_",names(main_reduce))] <- round2(main_reduce[,!grepl("^p_",names(main_reduce))])
  main_reduce[,grepl("^p_",names(main_reduce))] <- round2(main_reduce[,grepl("^p_",names(main_reduce))],3)
  main_reduce$RR_active <- ifelse(main_reduce$RR_active==0.01,0,main_reduce$RR_active)
  main_reduce$RR_control <- ifelse(main_reduce$RR_control==0.01,0,main_reduce$RR_control)
  main_reduce$RR_active <- paste0('[',main_reduce$RR_active,', ',main_reduce$RR_up_active,']')
    main_reduce$RR_control <- paste0('[',main_reduce$RR_control,', ',main_reduce$RR_up_control,']')
  main_reduce$m_diff_active <- paste0(main_reduce$m_diff_active,' (',main_reduce$SE_active,')')
    main_reduce$m_diff_control <- paste0(main_reduce$m_diff_control,' (',main_reduce$SE_control,')')
  main_reduce <- main_reduce[,!grepl("^RR_up|SE_",names(main_reduce))]
  
  rownames(main_reduce) <- c(rownames(main_reduce)[1:7],"PHQ Anxiety","PHQ Depression",
                             "PWB Environmental Mastery", "PWB Self-Acceptance",
                             "RRS Depression","RRS Brooding","RRS Reflection",
                             "WBSI Supression","WBSI Intrusion","PSS-10", rownames(main_reduce)[18:21])
  
  main_reduce

#SAMPLE SIZE ESTIMATION#######################
  
  BF_estimate <- function(dataset,n_q,conds,n_total,H,H1m,H1SD){
    #Returns the Bayes Factor for a given dataset, 
      #n_q = number of questions to sample of FFMQ
      #conds = conditions to compare (MSvW or WvControl)
      #n_total = suggested total sample size to estimate on
      #H = mean of model of data
      #H1m = mean of H1 model, H1SD = SD of H1 model
  
    #calulate SEs: SD of 100 means from a random selection of FFMQ questions
    samples <- replicate(100,sample(c(1:24),n_q)) 
    mean_diffs <- c()
    for(i in 1:ncol(samples)){
      qs <- samples[,i]
      pre_qs <- paste("FFMQ",qs,"pre", sep="_")
      post_qs <- paste("FFMQ",qs,"post", sep="_")
      pre <- rowMeans(dataset[pre_qs])
      post <- rowMeans(dataset[post_qs])
      diffs <- data.frame(diff = post-pre, cond = dataset$condition)
      mean_diff <- mean(diffs[diffs$cond == conds[1],"diff"]) - mean(diffs[diffs$cond == conds[2],"diff"])
      mean_diffs <- c(mean_diffs,mean_diff)
    }
    diffs_sd <- round(sd(mean_diffs),2)
    
    #Calculate BFs from SE
    n_group = n_total/3 #calculation is on a per-group basis
    n_cond <- c(sum(dataset$condition==conds[1]), sum(dataset$condition==conds[2]))
    harmonic_mean <- round(1/mean(1/n_cond))
    SE_estimate <- diffs_sd*sqrt(harmonic_mean/n_group)
    h0_mod <- bayesplay::prior(family="point", point=0)
    if(H1m==0){
      h1_mod <- bayesplay::prior(family="normal", mean=H1m, sd=H1SD, range=c(0, Inf))
    } else {
      h1_mod <- bayesplay::prior(family="normal", mean=H1m, sd=H1SD)
    }
    bayes_mod <- bayesplay::likelihood("normal", mean=H, sd=SE_estimate)
    bayes_factor <- integral(bayes_mod * h1_mod)/integral(bayes_mod * h0_mod)
    return(bayes_factor)
  }
  
  N_estimate <- function(dataset,n_qs=7,H1m=0,H1SD=.2,sizes=seq(120,350,10)){
    #creates a table (from a single dataset) of sample sizes and bayes factors for different sample sizes
    active <- c("mental","world")
    control <- c("world","control")
    bns <- as.data.frame(matrix(ncol=4,nrow=1))
    colnames(bns) <- c("comparison","n_q","n_p","bf")
    i=0
    for(cond_pair in list(active,control)){
      for(n_q in n_qs){
        for(ns in sizes){
          b <- BF_estimate(dataset,n_q=n_q,conds=cond_pair,n_total=ns,H=0,H1m=H1m,H1SD=H1SD)
          #bns[i,1] <- j
          i = i+1
          bns[i,1] <- cond_pair[1]
          bns[i,2] <- n_q
          bns[i,3] <- ns
          bns[i,4] <- b
        }
      }
    }
    return(bns)
  }
  
  N_table <- function(N_estimate_output){
    #reduces the output from the function above into a single meaned table
      #and then extracts only the relevant info
    a_list <- vector(mode = "list", length = length(N_estimate_output))
    i=0
    for(dataset in N_estimate_output){
      i=i+1
      a_list[i] <- dataset[4]
    }
    NEst_reduce <- N_estimate_output[[1]][1:3]
    NEst_reduce$bfs <- Reduce("+", a_list)/length(a_list)
    NEst_reduce$comparison <- paste0(NEst_reduce$comparison, NEst_reduce$n_q)
    n_need <- NEst_reduce[NEst_reduce$bf<.33,c(1,3,4)]
    a <- vector(length = 2); j <- 0
    for(i in unique(n_need$comparison)){
      j = j+1
      abc <- n_need[n_need$comparison==i,"n_p"]
      a[j] <- paste(i,min(abc))
    }
    print(min(n_need[,"n_p"]))
    #print(n_need[,"n_p"] %>% density %>% plot)
    #print(n_need[,"n_p"] %>% hist(breaks =6))
    return(a)
  }
  
  #Call functions to estimate sample size for TMS and GAD-7
  TMS_ns <- lapply(comps_agg,N_estimate)
  TMS_n <- N_table(TMS_ns)
  #PHQ-8
  PHQ_ns <- lapply(comps_agg,N_estimate,n_qs=8)
  PHQ_n <- N_table(PHQ_ns)
  #Expectancies
  exp_ns <- lapply(comps_agg,function(d){N_estimate(d, n_qs=c(4:7), H1m=0.2,H1SD=0.1)})
  exp_n <- N_table(exp_ns)
  #meta-d' BHN(0,.01)
  metad_est <- function(sizes, m_diff, se_diff){
    bns <- as.data.frame(matrix(ncol=2,nrow=1))
    colnames(bns) <- c("conds_n","bf")
    h0_mod <- bayesplay::prior(family="point", point=0)
    #dprime$m_diff/2
    h1_mod <- bayesplay::prior(family="normal", mean=0, sd=m_diff, range=c(0, Inf))
    #h1_mod <- prior(family = "uniform", min = 0, max = schmidt_metad$m_diff)
    schmidt_control <- 13
    schmidt_exp <- 14
    carpenter_control <- 32
    carpenter_exp <- 29
    harmonic_mean <- round(1/mean(1/c(schmidt_control,schmidt_exp,carpenter_control,carpenter_exp)))
    i=0
    for(ns in sizes){ 
      i = i+1
      bns[i,1] <- ns
      n_group = ns/3 #calculation is on a per-group basis
      SE_estimate <- se_diff*sqrt(harmonic_mean/n_group)
      bayes_mod <- bayesplay::likelihood("normal", mean=0, sd=SE_estimate)
      bayes_factor <- integral(bayes_mod * h1_mod)/integral(bayes_mod * h0_mod)
      bns[i,2] <- bayes_factor
    }
    bns_out <- bns#[bns$bf < .33]
    return(bns_out)
  }

  #mean of schmidt and carpenter meta-d'
  both_metad_n <- metad_est(seq(50,120,10),both_tab$combined_m[1],both_tab$combined_se[1])
  both_metad_nc <- both_metad_n[both_metad_n$bf<.33,][1,1]
  #mean of schmidt and carpenter for meta-d'/d'
  both_mratio_n <- metad_est(seq(120,300,10),both_tab$combined_m[2],both_tab$combined_se[2])
  both_mratio_nc <- both_mratio_n[both_mratio_n$bf<.33,][1,1]
  #log(meta-d'/d')
  both_log_n <- metad_est(seq(300,350,10),both_tab$combined_m[3],both_tab$combined_se[3])
  both_log_nc <- both_log_n[both_log_n$bf<.33,][1,1]
  #adjusted means
  both_adj_n <- metad_est(seq(120,200,10),m_both_adj,both_tab$combined_se[1])
  both_adj_nc <- both_adj_n[both_adj_n$bf<.33,][1,1]
  
  both_tab$est_n <- c(both_metad_nc,both_mratio_nc,both_log_nc)
  write.csv(both_tab,"est_n.csv")
  
#CALL RESULTS#######################
  nrow(all_data) - nrow(comps_agg$`1`)
  #missing data
  missing
  #write.csv(missing, file="output/MISSING_DATA.csv")
  #descriptives
  descriptives
  #write.csv(descriptives, file="output/DESCRIPTIVES.csv")
  #cronbach's alpha
  cron_reduce
  #write.csv(cron_reduce, file="output/CRONBACHS_ALPHA.csv")
  #condition*time cell means
  means_final
  #write.csv(means_final, file="output/CELL_MEANS.csv")
  #expectations
  expect_table
  #write.csv(expect_table,file="output/EXPECTANCIES.CSV")
  #bargraphs
  invisible(lapply(1:nrow(means_reduce),function(x){barGraph(means_reduce[x,])}))
  #jpeg("bargraphs.jpeg", width = 2000, height = 2000)
  #par(mfrow=c(4,5))
  #dev.off()
  #results table
  main_reduce
  #write.csv(main_reduce, file="output/RESULTS.csv")
  #minimum estimated sample size
  TMS_n
  metad_n
  #amount of expectancy questions needed
  exp_n
  
  