#SETUP#######################
# Load packages and data ---------------------------
  #NOTE: if there are issue with code, makes sure dplyr functions haven't been masked  
  #invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))
list.of.packages <- c("car", "pastecs", "stats", "data.table", "emmeans", "Rmisc", "tidyverse", "dplyr","magrittr","psych", "papaja") 
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  lapply(list.of.packages, require, character.only = TRUE)
  #remotes::install_github("crsh/papaja", dependencies = TRUE)
  #library("papaja")
  
  
mindful_data <- read.csv("mindful_data.csv") %>% 
  mutate_if(is.character, factor) %>% 
  dplyr::mutate(condition = fct_relevel(condition, "Waitlist", "World", "Mental"))

# Missing data ---------------------------

#survey website allowed people to continue without answering every question
  pre_range <- c(which(colnames(mindful_data) == "AWH_pre"):which(colnames(mindful_data) == "PWB_18_pre"))
    missing_pre <- mindful_data[!complete.cases(mindful_data[pre_range]),]%>%
                         group_by(condition) %>% summarise(count = length(condition)) #save count for later
  mindful_data <- mindful_data[complete.cases(mindful_data[pre_range]), ] #delete pre-test missing
#delete post-test missing questions, but not those who missed the whole second session
  post_range <- c(which(colnames(mindful_data) == "AWH_post"):which(colnames(mindful_data) == "PWB_18_post"))
      missing_post <- mindful_data[!complete.cases(mindful_data[!is.na(mindful_data$URN_post),post_range]),] %>%
                        group_by(condition) %>% summarise(count = length(condition)) #save count for later
    index_missing <- is.na(mindful_data$URN_post) #index of those that missed second session
    mindful_data_index <- mindful_data #create duplicate dataset
    mindful_data_index[index_missing,post_range] <- mindful_data_index[index_missing,pre_range] #Dummy ITT data
    mindful_data <- mindful_data[complete.cases(mindful_data_index[post_range]), ] #those missing some but not all data
    
    m_all_post <- mindful_data %>% filter(URN_post %>% is.na) %>% group_by(condition) %>% summarise(count = length(condition))
#First few rows were calculated when converging the different datasets, although this code is in a separate file as this included emails.
missing_data <-tribble(
                 ~Variable,       ~Count,                                                    ~Condition,
                 "Condition",      62,                                                         "-",        
                 "E-mail",         1,                                                          "Mental",   
                 "All Pre-test",   1,                                                          "Mental",   
                 "All Pre-test",   3,                                                          "World",    
                 "All Pre-test",   3,                                                          "Waitlist",
                 "Some Pre-test",  filter(missing_pre, condition=="Mental") %>% pull(count),   "Mental",
                 "Some Pre-test",  filter(missing_pre, condition=="Waitlist") %>% pull(count), "Waitlist",
                 "Some Post-test", filter(missing_post, condition=="Mental") %>% pull(count),  "Mental",
                 "Some Post-test", filter(missing_post, condition=="World") %>% pull(count),   "World",
                 "All Post-test",  filter(m_all_post, condition=="Mental") %>% pull(count),    "Mental",
                 "All Post-test",  filter(m_all_post, condition=="World") %>% pull(count),     "World",
                 "All Post-test",  filter(m_all_post, condition=="Waitlist") %>% pull(count),  "Waitlist"
                )

mindful_data %>% group_by(condition) %>% summarise(count = length(condition))

missing_data %>% write.csv(file="missing_data.csv")    
# Scale Re-coding ---------------------------
PSS_recode <- c(1:4,9) %>% as.character() %>%
    {c(paste0("PSS_",.,"_pre"), paste0("PSS_",.,"_post"))}
  mindful_data[PSS_recode] <- mindful_data[PSS_recode] - 1
PHQ_recode <- c(1:4) %>% as.character() %>%
    {c(paste0("PHQ_",.,"_pre"), paste0("PHQ_",.,"_post"))}
  mindful_data[PHQ_recode] <- mindful_data[PHQ_recode] - 1

#Reverse Coding
FFMQ_reverse <- c(4,5,7,8,11,12,14,17,19,22:24) %>% as.character() %>%
    {c(paste0("FFMQ_",.,"_pre"), paste0("FFMQ_",.,"_post"))}
  mindful_data[FFMQ_reverse] <- 6 - mindful_data[FFMQ_reverse]
PSS_reverse <- c(4,5,7,8) %>% as.character() %>%
    {c(paste0("PSS_",.,"_pre"), paste0("PSS_",.,"_post"))}
  mindful_data[PSS_reverse] <- 5 - mindful_data[PSS_reverse]
PWB_reverse <- c(2,4,6,10,11,14,16) %>% as.character() %>%
    {c(paste0("PWB_",.,"_pre"), paste0("PWB_",.,"_post"))}
  mindful_data[PWB_reverse] <- 7 - mindful_data[PWB_reverse]

# Scale Objects ---------------------------      

#Awareness last hour
  AH_pre <- c("AMH_pre","AWH_pre")
  AH_post <- c("AMH_post","AWH_post")
  
#FFMQ
DS <- c(1,2,5,11,16) %>% as.character()
  DS_pre <- {paste0("FFMQ_",DS,"_pre")}
  DS_post <- {paste0("FFMQ_",DS,"_post")}
NR <- c(3,9,13,18,21) %>% as.character()
  NR_pre <- {paste0("FFMQ_",NR,"_pre")}
  NR_post <- {paste0("FFMQ_",NR,"_post")}
NJ <- c(4,7,14,19,24) %>% as.character()
  NJ_pre <- {paste0("FFMQ_",NJ,"_pre")}
  NJ_post <- {paste0("FFMQ_",NJ,"_post")}
OB <- c(6,10,15,20) %>% as.character()
  OB_pre <- {paste0("FFMQ_",OB,"_pre")}
  OB_post <- {paste0("FFMQ_",OB,"_post")}
AA <- c(8,12,17,22,23) %>% as.character()
  AA_pre <- {paste0("FFMQ_",AA,"_pre")}
  AA_post <- {paste0("FFMQ_",AA,"_post")}

#Mental Health scales
PHQ_anx <- c(1,2) %>% as.character()
  PHQ_anx_pre <- {paste0("PHQ_",PHQ_anx,"_pre")}
  PHQ_anx_post <- {paste0("PHQ_",PHQ_anx,"_post")}
PHQ_dep <- c(3,4) %>% as.character()
  PHQ_dep_pre <- {paste0("PHQ_",PHQ_dep,"_pre")}
  PHQ_dep_post <- {paste0("PHQ_",PHQ_dep,"_post")}
PWB_env <- c(1,2,4,7,11,12,13,16,18) %>% as.character()
  PWB_env_pre <- {paste0("PWB_",PWB_env,"_pre")}
  PWB_env_post <- {paste0("PWB_",PWB_env,"_post")}
PWB_acc <- c(3,5,6,8,9,10,14,15,17) %>% as.character()
  PWB_acc_pre <- {paste0("PWB_",PWB_acc,"_pre")}
  PWB_acc_post <- {paste0("PWB_",PWB_acc,"_post")}
  
PSS <- c(1:10) %>% as.character()
  PSS_pre <- {paste0("PSS_",PSS,"_pre")}
  PSS_post <- {paste0("PSS_",PSS,"_post")}
RRS <- c(1:22) %>% as.character()
  RRS_pre <- {paste0("RRS_",RRS,"_pre")}
  RRS_post <- {paste0("RRS_",RRS,"_post")}
WBSI <- c(1:15) %>% as.character()
  WBSI_pre <- {paste0("WBSI_",WBSI,"_pre")}
  WBSI_post <- {paste0("WBSI_",WBSI,"_post")}

#used for expectations
PHQ_all <- c(1:4) %>% as.character()
  PHQ_all_pre <- {paste0("PHQ_",PHQ_all,"_pre")}
  PHQ_all_post <- {paste0("PHQ_",PHQ_all,"_post")}
  
#subscale splits for appendix
RRS_D <- c(1:4,6,8,9,14,17:19,22) %>% as.character
  RRS_D_pre <- {paste0("RRS_",RRS_D,"_pre")}
  RRS_D_post <- {paste0("RRS_",RRS_D,"_post")}
RRS_B <- c(5,10,13,15,16) %>% as.character
  RRS_B_pre <- {paste0("RRS_",RRS_B,"_pre")}
  RRS_B_post <- {paste0("RRS_",RRS_B,"_post")}
RRS_R <- c(7,11,12,20,21) %>% as.character
  RRS_R_pre <- {paste0("RRS_",RRS_R,"_pre")}
  RRS_R_post <- {paste0("RRS_",RRS_R,"_post")}
WBSI_sup <- c(13,1,8,10,11,14) %>% as.character()
  WBSI_sup_pre <- {paste0("WBSI_",WBSI_sup,"_pre")}
  WBSI_sup_post <- {paste0("WBSI_",WBSI_sup,"_post")}
WBSI_int <- c(2,3,4,5,6,7,9,12,15) %>% as.character()
  WBSI_int_pre <- {paste0("WBSI_",WBSI_int,"_pre")}
  WBSI_int_post <- {paste0("WBSI_",WBSI_int,"_post")}

#list of all scales, and pre and post separately
  scales <- c("AMH","AWH",
              "OB","NR","DS","NJ","AA",
              "PHQ_anx","PHQ_dep",
              "PWB_env","PWB_acc",
              "PSS","RRS","WBSI",
              "RRS_D","RRS_B","RRS_R",
              "WBSI_sup","WBSI_int",
              "PHQ_all")

    for(scale in scales){
      if(scale == scales[1]){
        scales_pre <- paste0(scale,"_pre")
        scales_post <- paste0(scale,"_post")
        scales_time <- c(scales_pre,scales_post)
      }
      else{
        scales_pre <- c(scales_pre,paste0(scale,"_pre"))
        scales_post <- c(scales_post,paste0(scale,"_post"))
        scales_time <- c(scales_time,paste0(scale,"_pre"),paste0(scale,"_post"))
      }
    }
# Mean Scores --------------------------- 
  for(scale in scales_time){
    if(grepl("^A.H",scale)){next}
    measure <- mindful_data[get(scale)]
    varname <- paste0(scale)
    mindful_data <- dplyr::mutate(mindful_data,!!varname := rowMeans(measure), na.rm=TRUE)
  } 
  
mindful_data <- mindful_data[,!(names(mindful_data) %in% c("na.rm"))] %>% 
  dplyr::mutate(expect_mean = rowMeans(mindful_data[c("expect_anx","expect_dep")], na.rm=TRUE))

#Imputation-----
mindful_means <- mindful_data[which(colnames(mindful_data)=="condition"):which(colnames(mindful_data)=="expect_mean")]
#22% of post-test values are missing
sum(is.na(mindful_means["OB_post"]))/length(mindful_means[,"OB_post"])*100


#LOCF
  ##mindful_data[is.na(mindful_data$URN_post),scales_post] <- mindful_data[is.na(mindful_data$URN_post),scales_pre]
  
#Difference scores
  for(scale in scales){
    mindful_data[paste0(scale,"_diff")] <- mindful_data[paste0(scale,"_post")] - mindful_data[paste0(scale,"_pre")]
  }
  
# Filtered Datasets ---------------------------   
  waitlist <- mindful_data[mindful_data$condition == "Waitlist",]
  world <- mindful_data[mindful_data$condition == "World",]
  mental <- mindful_data[mindful_data$condition == "Mental",]

  exp_data <- mindful_data %$%
    {!(is.na(expect_anx) | is.na(expect_dep)) & (condition == "World" | condition == "Mental")} %>% 
    {mindful_data[.,grepl("^expect_|^condition$|^PHQ_(anx|dep|all)", names(mindful_data))]} %>% 
    droplevels
#DATA EXPLORATION#######################
# Descriptives ---------------------------
#Demographics
mindful_data %>% {tibble(Variable = c(levels(.$condition),levels(.$gender),"Age","Range"),
                         Description = c(sapply(Variable[1:3], function(x) sum(.$condition == x)),
                                       sapply(Variable[4:5], function(x) sum(.$gender == x)),
                                       printnum(mean(.$age)),paste0(min(.$age),"-",max(.$age))))} %>% 
    write.csv("demographics.csv")

# Missing data
  missing_pre_cond <- 69 #missing all pre-test and condition
tibble(Reason = c("Missing-at-Random", "Dropouts", "Intention-to-Treat"),
       Count  = c(missing_pre_cond + missing_pre + missing_post, sum(is.na(mindful_data$URN_post)), 
                  NROW(mindful_data)))

# Descriptives
lapply(levels(mindful_data$condition), function(x){
  tibble("{x}_Pre" := sapply(scales_pre, function(y){paste0(printnum(mean(mindful_data[which(mindful_data$condition==x),y])),
                                                      " (",printnum(sd(mindful_data[which(mindful_data$condition==x),y])),")")}),
         "{x}_Post" := sapply(scales_post, function(z){paste0(printnum(mean(mindful_data[which(mindful_data$condition==x),z])),
                                                        " (",printnum(sd(mindful_data[which(mindful_data$condition==x),z])),")")})
         )}) %>% 
  bind_cols %>% dplyr::mutate(Scale = scales, .before = Waitlist_Pre) %>% 
  write.csv("descriptives.csv")

# Exectation Descriptives
lapply(levels(exp_data$condition), function(x){
  tibble("{x}" := sapply(names(exp_data[c("expect_anx","expect_dep")]), function(y){paste0(printnum(mean(exp_data[which(exp_data$condition==x),y])),
                                                            " (",printnum(sd(exp_data[which(exp_data$condition==x),y])),")")})
  )}) %>% bind_cols %>% dplyr::mutate(Scale = names(exp_data[c("expect_anx","expect_dep")]), .before = World) %>%
  write.csv("exp_desc.csv")

# Cronbach's Alpha
lapply(scales, function(x){
  if(grepl("^A.H",x)){return()}
  tibble(Scale = x,
         Pre = psych::alpha(mindful_data[get(paste0(x,"_pre"))])$total$raw_alpha,
         Post = psych::alpha(mindful_data[get(paste0(x,"_post"))])$total$raw_alpha
  )}) %>% bind_rows %>% printnum %>%
  write.csv("cronbachs.csv")

lapply(scales, function(x){
  if(grepl("^A.H",x)){return()}
  psych::alpha(mindful_data[get(paste0(x,"_pre"))])
  psych::alpha(mindful_data[get(paste0(x,"_post"))])
  })
# Bargraphs ---------------------------   
melt_se <- function(dataset){
  scale_melt <- gather(dplyr::select(dataset, scales_time, condition, URN), scale, score, -URN, -condition)
  scale_melt$time <- gl(2, NROW(dataset), labels=c("Pre-test","Post-test"))
  scale_melt$measure <- gl(NROW(scales), NROW(dataset)*2, labels=scales)
  scale_SE <- summarySE(scale_melt, measurevar="score", groupvars=c("measure","time","condition"))
  return(scale_SE)
}
mindful_melt <- melt_se(mindful_data) %>% 
                mutate(measure = fct_recode(measure, #rename graphed variables
                                            "Awareness: Mental States" = "AMH",
                                            "Awareness: World" = "AWH",
                                            "FFMQ: Observe" = "OB",
                                            "PHQ: Anxiety" = "PHQ_anx",
                                            "PHQ: Depression" =  "PHQ_dep",
                                            "PWB: Environmental Mastery" = "PWB_env",
                                            "RRS: Depression-related" = "RRS_D",                     
                                            "RRS: Brooding" = "RRS_B",                      
                                            "RRS: Reflection" = "RRS_R",
                                            "WBSI: Intrusive Thoughts" = "WBSI_int"
                                            ))
bargraphs <- function(dataset,free_fixed){
  ggplot2::ggplot(dataset, 
    aes(condition, score, fill=time))+
    geom_bar(position=position_dodge(),width=.7,stat="identity")+
    geom_errorbar(aes(ymin=score-ci, ymax=score+ci), width=.2, position=position_dodge(.7))+
    facet_wrap(~measure, scales = free_fixed)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values=c( "#0072B2", "#D55E00"))+
    guides(fill=guide_legend(title="Time"))+ 
    labs(x = "Condition", y = "Score (+/- C.I. 95%)")
}

bargraphs(mindful_melt[mindful_melt$measure == "Awareness: Mental States"| 
                       mindful_melt$measure == "Awareness: World",], "fixed")
ggsave("Graph_Aware.png",  width = 20, height = 10, units = "cm")

bargraphs(mindful_melt[mindful_melt$measure == "FFMQ: Observe",], "fixed")
ggsave("Graph_FFMQ.png",  width = 15, height = 10, units = "cm")

bargraphs(mindful_melt[mindful_melt$measure == "PHQ: Anxiety"|
                       mindful_melt$measure == "PWB: Environmental Mastery",], "free")
ggsave("Graph_MH.png",  width = 20, height = 10, units = "cm")

bargraphs(mindful_melt[mindful_melt$measure == "RRS: Depression-related"|
                       mindful_melt$measure == "RRS: Brooding"|
                       mindful_melt$measure == "RRS: Reflection"|
                       mindful_melt$measure == "WBSI: Intrusive Thoughts"
                       ,], "free")
ggsave("Graph_MH_Appendix.png",  width = 20, height = 20, units = "cm")
#INFERENTIAL#######################
# Bayes Factor function ---------------------------
Bf<-function(sd, obtained, dfdata = 1, likelihood = c("normal", "t"), 
             modeloftheory= c("normal","t","cauchy", "uniform") ,lower =0, 
             upper=1, modeoftheory = 0, scaleoftheory = 1, dftheory = 1, tail = 2)
{
  if(likelihood=="normal"){
    dfdata=10^10
  }
  if(modeloftheory=="normal"){
    dftheory = 10^10
  } else if(modeloftheory=="cauchy"){
    dftheory = 1
  }
  area <- 0
  normarea <- 0
  if(modeloftheory=="uniform"){
    theta <- lower
    range <- upper - lower
    incr <- range / 2000
    for (A in -1000:1000){
      theta <- theta + incr
      dist_theta <- 1 / range
      height <- dist_theta * dt((obtained-theta)/sd, df=dfdata)
      area <- area + height * incr
    }
    LikelihoodTheory <- area
  }else{
    theta <- modeoftheory - 10 * scaleoftheory
    incr <- scaleoftheory/200
    for (A in -2000:2000){
      theta <- theta + incr
      dist_theta <- dt((theta-modeoftheory)/scaleoftheory, df=dftheory)
      if(identical(tail, 1)){
        if (theta <= modeoftheory){
          dist_theta <- 0
        } else {
          dist_theta <- dist_theta * 2
        }
      }
      height <- dist_theta * dt((obtained-theta)/sd, df = dfdata)
      area <- area + height * incr
      normarea <- normarea + dist_theta*incr
    }
    LikelihoodTheory <- area/normarea
  }
  Likelihoodnull <- dt(obtained/sd, df = dfdata)
  BayesFactor <- LikelihoodTheory/Likelihoodnull
  BayesFactor
}

# Expectations ---------------------------
#Note: data loss: these were collected with emails
#Datasets
#T-test
exp_t <- lapply(exp_data[c("expect_anx","expect_dep","expect_mean")], function(x){
    t.test(x ~ condition, exp_data, var.equal=TRUE) %$% 
      c(printnum(statistic[[1]]), printp(p.value), printnum(conf.int)) %>%
      set_names(c("t","p","CI-","CI+"))
  }) %>% bind_rows %>% unite("CI", "CI-","CI+", sep=", ")

##observed differences 
#obs_diff_anx <- {mean(exp_data[exp_data$condition == "Mental","PHQ_anx_diff"]) -
#        mean(exp_data[exp_data$condition == "World","PHQ_anx_diff"])} %>% {.*-1} %>% round(2)
#obs_diff_dep <- {mean(exp_data[exp_data$condition == "Mental","PHQ_dep_diff"]) -
#        mean(exp_data[exp_data$condition == "World","PHQ_dep_diff"])} %>% {.*-1}%>% round(2)
##original H1 calculation
#exp_prior_anx <- lm(PHQ_anx_post ~ PHQ_anx_pre + expect_anx, data=exp_data) %>% summary() %$%
#  {coefficients["expect_anx","Estimate"] + (coefficients["expect_anx","Std. Error"]*1.96)} %>% round(2) %>% 
#  {obs_diff_anx/.} %>% round(2) %>% printnum
#exp_prior_dep <- lm(PHQ_dep_post ~ PHQ_dep_pre + expect_dep, data=exp_data) %>% summary() %$%
#  {coefficients["expect_dep","Estimate"] + (coefficients["expect_dep","Std. Error"]*1.96)} %>% round(2) %>% 
#  {obs_diff_dep/.} %>% round(2) %>% printnum

##PRIOR
#exp_prior_all <- lm(PHQ_all_post ~ PHQ_all_pre + expect_mean, data=exp_data) %>% summary() %$%
#  {coefficients["expect_mean","Estimate"] + (coefficients["expect_mean","Std. Error"]*1.96)} %>% round(2) %>% 
#  {.2/.} %>% round(2) %>% printnum
#
##room-to-move heuristic
#r2m_anx <- mean(exp_data[exp_data$condition == "World","expect_anx"]) %>% 
#            round(2) %>% {5-.} %>% {./2}
#r2m_dep <- mean(exp_data[exp_data$condition == "World","expect_dep"]) %>% 
#            round(2) %>% {5-.} %>% {./2}

#Table and analyses
exp_data %>%
  group_by(condition) %>%
  dplyr::summarise(across(starts_with("expect"), 
                   list(mean = mean, n = length, se = ~sd(.x)/sqrt(length(.x))
                   ))) %>% printnum %>%
  {tibble(Expectancy = c("Anxiety", "Depression", "Both"),
          Diff = printnum(abs(as.numeric(unlist(.[2,grepl("_mean$",names(.))])) -
                                as.numeric(unlist(.[1,grepl("_mean$",names(.))])))),
          SE = printnum(sqrt((as.numeric(unlist(.[2,grepl("_se",names(.))]))^2) +
                               (as.numeric(unlist(.[1,grepl("_se",names(.))]))^2))),
          Diff_SE = paste0(Diff," (",SE,")"),
          Bf = printnum(Bf(as.numeric(SE), as.numeric(Diff), 39, "normal", 
                           "normal", modeoftheory = 1.5, scaleoftheory = .75, tail = 2,
                           dftheory = 100))
  )} %>% select(-Diff,-SE) %>% add_column(exp_t) %>%
  write.csv("expectations.csv")
#"uniform",lower =0, upper=3,dftheory = 100, tail = 1
#"normal", modeoftheory = 3, scaleoftheory = 1.5, dftheory = 100, tail = 2
# Main Analyses ---------------------------
contrasts = list(Main    = c(0, -1, 1),
                 Control = c(-1, 1, 0))

lapply(paste0(scales,"_diff"), function(x){
  contrast_lm <- lm(get(x) ~ condition, mindful_data) %>%
        emmeans(., "condition", contr = contrasts) %$% contrasts %>% 
        c(summary(.), apa_print(.))
  
  ifelse(grepl("^A.H",x), H1<-2.5, H1<-.2)
  
  tibble(Scales          = sub("_diff", "", x),
         Active_mean     = printnum(mean(mental[,x]) - mean(world[,x])),
         Active_SE       = printnum(sqrt((sd(mental[,x])/sqrt(NROW(mental))) ^2+
                                (sd(world[,x])/sqrt(NROW(world))) ^2)),
         Active          = paste0(Active_mean," (",Active_SE,")"),
         Active_BF       = Bf(as.numeric(Active_SE), 
                              ifelse(x%in%c("PSS_diff","PHQ_anx_diff","PHQ_dep_diff","PHQ_all_diff",
                                            "RRS_diff","RRS_D_diff","RRS_B_diff","RRS_R_diff",
                                            "WBSI_diff","WBSI_sup_diff","WBSI_int_diff"),as.numeric(Active_mean)*-1,as.numeric(Active_mean)),
                              NROW(mental) + NROW(world) - 2, "normal", "normal", modeoftheory = 0, 
                              scaleoftheory = H1, dftheory = 102, tail = 1),
        "Active_t(119)" := contrast_lm$t.ratio[1],
         Active_p        = contrast_lm$p.value[1],

         Control_mean    = printnum(mean(world[,x]) - mean(waitlist[,x])),
         Control_SE      = printnum(sqrt((sd(world[,x])/sqrt(NROW(world))) ^2+
                                (sd(waitlist[,x])/sqrt(NROW(waitlist))) ^2)),
         Control         = paste0(Control_mean," (",Control_SE,")"),
         Control_BF      = Bf(as.numeric(Control_SE), 
                              ifelse(x%in%c("PSS_diff","PHQ_anx_diff","PHQ_dep_diff","PHQ_all_diff",
                                            "RRS_diff","RRS_D_diff","RRS_B_diff","RRS_R_diff",
                                            "WBSI_diff","WBSI_sup_diff","WBSI_int_diff"),as.numeric(Control_mean)*-1,as.numeric(Control_mean)), 
                              NROW(world) + NROW(waitlist) - 2, "normal","normal", modeoftheory = 0,
                              scaleoftheory = H1, dftheory = 102, tail = 1),
        "Control_t(119)" = contrast_lm$t.ratio[2],
         Control_p       = contrast_lm$p.value[2]
  )}) %>% bind_rows %>% select(-ends_with("_mean"), -ends_with("_SE")) %>%
  dplyr::mutate(across(ends_with("_p"), printp), across(!ends_with("_p"), printnum)) %>%
  write.csv("analysis.csv")

#Waitlist-World for Awareness and PWB Env
lapply(paste0(c("AMH","AWH","PWB_env"),"_diff"), function(x){
  ifelse(grepl("^A.H",x), H1<-2.5, H1<-.2)
  tibble(Scales          = sub("_diff", "", x),
         Control_mean    = printnum(mean(waitlist[,x]) - mean(world[,x])),
         Control_SE      = printnum(sqrt((sd(world[,x])/sqrt(NROW(world))) ^2+
                                           (sd(waitlist[,x])/sqrt(NROW(waitlist))) ^2)),
         Control         = paste0(Control_mean," (",Control_SE,")"),
         Control_BF      = Bf(as.numeric(Control_SE), as.numeric(Control_mean), 
                              NROW(world) + NROW(waitlist) - 2, "normal","normal", modeoftheory = 0,
                              scaleoftheory = H1, dftheory = 102, tail = 1)
  )}) %>% bind_rows %>% select(-ends_with("_mean"), -ends_with("_SE")) %>%
  dplyr::mutate(across(!ends_with("_p"), printnum))

# Robustness Regions ----
lapply(paste0(scales,"_diff"), function(x){
  tibble(Scales           = sub("_diff", "", x),
         Active_SE        = sqrt((sd(mental[,x])/sqrt(NROW(mental))) ^2+ (sd(world[,x])/sqrt(NROW(world))) ^2),
         Active_mean_ur   = mean(mental[,x]) - mean(world[,x]),
         Active_mean      = ifelse(x%in%c("PSS_diff","PHQ_anx_diff","PHQ_dep_diff","PHQ_all_diff",
                                          "RRS_diff","RRS_D_diff","RRS_B_diff","RRS_R_diff",
                                          "WBSI_diff","WBSI_sup_diff","WBSI_int_diff"),Active_mean_ur*-1,Active_mean_ur),
         Control_SE       = sqrt((sd(world[,x])/sqrt(NROW(world))) ^2+ (sd(waitlist[,x])/sqrt(NROW(waitlist))) ^2),
         Control_mean_ur  = mean(world[,x]) - mean(waitlist[,x]),
         Control_mean     = ifelse(x%in%c("PSS_diff","PHQ_anx_diff","PHQ_dep_diff","PHQ_all_diff",
                                          "RRS_diff","RRS_D_diff","RRS_B_diff","RRS_R_diff",
                                          "WBSI_diff","WBSI_sup_diff","WBSI_int_diff"),Control_mean_ur*-1,Control_mean_ur)
  )}) %>% bind_rows %>% select(-ends_with("_ur")) %>% printnum

#Scales   Active                    Control
#AMH      0.45, -0.14 [.99,∞]       0.55, -0.66  [.70,∞]
#AWH      0.45,  0.62 [0,6.32]      0.57, -1.74  [.34,∞]
#OB       0.13,  0.23 [.17,.24]     0.13, -0.01  [0,.35]
#NR       0.15,  0.10 [0,.82]       0.14,  0.00  [0,.4]
#DS       0.12,  0.07 [0,.59]       0.13, -0.08  [0,.23]
#NJ       0.16, -0.05 [0,.34]       0.15,  0.12  [0,.94]
#AA       0.16,  0.25 [0,.22]       0.15,  0.04  [0,.53]
#PSS      0.10,  0.03 [0,.37]       0.11,  0.05  [0,.47]
#PHQ_anx  0.15,  0.39 [.07,2.89]    0.19,  0.01  [0,.56]
#PHQ_dep  0.14,  0.16 [0,1.38]      0.16,  0.08  [0,.74]
#RRS      0.09,  0.12 [0,1.18]      0.09,  0.06  [0,.48]
#WBSI     0.13,  0.15 [0,1.31]      0.11,  0.05  [0,.47]
#PWB_env  0.13,  0.24 [.14,.34]     0.12, -0.12  [.17,∞]
#PWB_acc  0.10,  0.10 [0,.82]       0.09, -0.02  [0,.21]

#RRS_D    0.10, 0.13  [0,1.24]      0.10, 0.06  [0,.51]
#RRS_B    0.11, 0.19  [0,2.80]      0.10, 0.00  [0,.29]
#RRS_R    0.13, 0.01  [0,.39]       0.13, 0.13  [0,1.07]
#WBSI_sup 0.15, 0.12  [0,.94]       0.13, 0.05  [0,.53]
#WBSI_int 0.13, 0.17  [0,1.65]      0.13, 0.05  [0,.53]

Bf(0.13, 0.05, scaleoftheory = .53,
   39,"normal", "normal", modeoftheory = 0,
   dftheory = 100, tail = 1) %>% printnum

#Expectations  
#Anxiety .35,.14     U: [1.42,>10]   N: [.80,<10]
#Depression	.32,.28  U: [2.22,>10]   N: [1.08,<10]
#Combined .30,.21    U: [2.24,>10]   N: [.88,<10]


mean_h1 = .8
sd_h1 = mean_h1/2
Bf(.35,.14 ,
   dfdata = 39, likelihood = "normal", 
   "normal", modeoftheory = mean_h1, 
   scaleoftheory = sd_h1, 
   tail = 2, dftheory = 100) %>% printnum
#"uniform", lower =0, upper=3,  tail = 1,
  #"normal", modeoftheory = .2, scaleoftheory = .1, tail = 2,

# Sample Size Estimation ----

#To calculate expected posterior SD for TMS and OB
  #Get the SD of 100 differences between mental-world means for n random FFMQ questions
SE_mindful <- function(n_questions, analysis){
    replicate(100, sample(c(1:24),n_questions)) %>% as.tibble %>%
    lapply(function(x){
      cond <- mindful_data$condition
      pre <- mindful_data[paste("FFMQ",x,"pre", sep="_")] %>% rowMeans
      post <- mindful_data[paste("FFMQ",x,"post", sep="_")] %>% rowMeans
      diff <- bind_cols(pre,post) %>% transmute(diff = ...2 - ...1) %>% replace_na(list(diff = 0)) %>%
        mutate(cond = cond)
      if(analysis == "active"){
      mean(diff$diff[diff$cond == "Mental"]) - mean(diff$diff[diff$cond == "World"])
      } else if(analysis == "control"){
      mean(diff$diff[diff$cond == "World"]) - mean(diff$diff[diff$cond == "Waitlist"])
      }
    }) %>% unlist %>% sd %>% round(2)
}

#Estimate sample size needed, see: https://www.youtube.com/watch?v=10Lsm_o_GRg
n_size <- function(SE, N_total, H = .2, analysis = "active"){
  N_group = N_total/3 #calculation is on a per-group basis
  
 if(analysis == "active"){
    df_theory <- NROW(mental) + NROW(world) - 2
    active_n <- c(NROW(mental), NROW(world))
    harmonic_mean <- 1/mean(1/active_n)
    harmonic_mean <- round(harmonic_mean)
  } else if(analysis == "control"){
    df_theory <- NROW(world) + NROW(waitlist) - 2
    control_n <- c(NROW(world), NROW(waitlist))
    harmonic_mean <- 1/mean(1/control_n)
    harmonic_mean <- round(harmonic_mean)
  }

  SE_estimate <- SE*sqrt(harmonic_mean/N_group)

  Bf(SE_estimate, H, 100,"normal", "normal", modeoftheory = 0, 
     scaleoftheory = .2, dftheory = df_theory, tail = 1) %>% 
    round(digits = 2) %>% paste("BF =",.)
}

#TMS decentering and Full FFMQ OB (minus 1 irrelevant Q), and the GAD all have 7 questions
set.seed(1)
active_7 <- SE_mindful(7, "active")
set.seed(1)
control_7 <- SE_mindful(7, "control")

#for estimating how many expectancy qs to put in
set.seed(1)
exp_4 <- SE_exp(4, "active")
set.seed(1)
exp_5 <- SE_exp(5, "active")
set.seed(1)
exp_6 <- SE_exp(6, "active")
set.seed(1)
exp_7 <- SE_exp(7, "active")

#PHQ-9
set.seed(1)
active_9 <- SE_mindful(9, "active")
set.seed(1)
control_9 <- SE_mindful(9, "control")

n_size(control_7,1, H=0, "control") #55 original, round to 57 so divisible by 3
n_size(active_7,24, H=0, "active") #92
n_size(control_9,36, H=0, "control") #35
n_size(active_9,66, H=0, "active") #64

#estimation for expectations
n_size_exp <- function(SE, N_total, H = .2, analysis = "active"){
  N_group = N_total/3 #calculation is on a per-group basis
  
  if(analysis == "active"){
    df_theory <- NROW(mental) + NROW(world) - 2
    active_n <- c(NROW(mental), NROW(world))
    harmonic_mean <- 1/mean(1/active_n)
    harmonic_mean <- round(harmonic_mean)
  } else if(analysis == "control"){
    df_theory <- NROW(world) + NROW(waitlist) - 2
    control_n <- c(NROW(world), NROW(waitlist))
    harmonic_mean <- 1/mean(1/control_n)
    harmonic_mean <- round(harmonic_mean)
  }
  
  SE_estimate <- SE*sqrt(harmonic_mean/N_group)
  
  Bf(SE_estimate, H, 100,"normal", "normal", modeoftheory = .2, scaleoftheory = .1, tail = 2
     #"uniform", lower =0, upper=.4, tail = 1
     #"normal", modeoftheory = .2, scaleoftheory = .1, tail = 2
     dftheory = df_theory) %>% 
    round(digits = 2) %>% paste("BF =",.)
}

SE_exp <- function(n_questions, analysis){
  replicate(100, sample(c(1:24),n_questions)) %>% as.tibble %>%
    lapply(function(x){
      cond <- mindful_data$condition
      pre <- mindful_data[paste("FFMQ",x,"pre", sep="_")] %>% rowMeans
      post <- mindful_data[paste("FFMQ",x,"post", sep="_")] %>% rowMeans
      diff <- bind_cols(pre,post) %>% transmute(diff = ...2 - ...1) %>% replace_na(list(diff = 0)) %>%
        mutate(cond = cond)
      if(analysis == "active"){
        mean(diff$diff[diff$cond == "Mental"]) - mean(diff$diff[diff$cond == "World"])
      } else if(analysis == "control"){
        mean(diff$diff[diff$cond == "World"]) - mean(diff$diff[diff$cond == "Waitlist"])
      }
    }) %>% unlist %>% sd %>% round(2)
}



n_size(exp_4,161, H=0, "active") #U: 170, N: 161
n_size(exp_5,124, H=0, "active") #both same
n_size(exp_6,92, H=0, "active")
n_size(exp_7,92, H=0, "active")




#Ignore#######################
# OB question means -----
lapply(levels(mindful_data$condition), function(x){
  tibble("{x}" := sapply(OB, function(y){mean(mindful_data[which(mindful_data$condition==x),paste0("FFMQ_",y,"_post")] - 
                           mindful_data[which(mindful_data$condition==x),paste0("FFMQ_",y,"_pre")], na.rm=T)}))}) %>%
  bind_cols %>% printnum %>% dplyr::mutate(Q = paste0("#", OB), .before = Waitlist) %>%
  write.csv("observe.csv")

# WBSI separate sup and int -----
mindful_data_WBSI <- mindful_data

WBSI_sup <- c(13,1,8,10,11,14) %>% as.character()
WBSI_sup_pre <- {paste0("WBSI_",WBSI_sup,"_pre")}
WBSI_sup_post <- {paste0("WBSI_",WBSI_sup,"_post")}
WBSI_int <- c(2,3,4,5,6,7,9,12,15) %>% as.character()
WBSI_int_pre <- {paste0("WBSI_",WBSI_int,"_pre")}
WBSI_int_post <- {paste0("WBSI_",WBSI_int,"_post")}

WBSI_scales <- c("WBSI_sup_pre", "WBSI_sup_post", "WBSI_int_pre", "WBSI_int_post")
#means pre and post
for(scale in WBSI_scales){
  measure <- mindful_data_WBSI[get(scale)]
  varname <- paste0(scale)
  mindful_data_WBSI <- dplyr::mutate(mindful_data_WBSI,!!varname := rowMeans(measure), na.rm=TRUE)
} 

#Intention-to-treat
mindful_data_WBSI[is.na(mindful_data_WBSI$URN_post),"WBSI_sup_post"] <- mindful_data_WBSI[is.na(mindful_data_WBSI$URN_post),"WBSI_sup_pre"]
mindful_data_WBSI[is.na(mindful_data_WBSI$URN_post),"WBSI_int_post"] <- mindful_data_WBSI[is.na(mindful_data_WBSI$URN_post),"WBSI_int_pre"]

#Difference scores
mindful_data_WBSI$WBSI_sup_diff <- mindful_data_WBSI$WBSI_sup_post - mindful_data_WBSI$WBSI_sup_pre
mindful_data_WBSI$WBSI_int_diff <- mindful_data_WBSI$WBSI_int_post - mindful_data_WBSI$WBSI_int_pre


waitlist_WBSI <- mindful_data_WBSI[mindful_data_WBSI$condition == "Waitlist",]
world_WBSI <- mindful_data_WBSI[mindful_data_WBSI$condition == "World",]
mental_WBSI <- mindful_data_WBSI[mindful_data_WBSI$condition == "Mental",]


contrasts = list(Main    = c(0, -1, 1),
                 Control = c(-1, 1, 0))

lapply(c("WBSI_sup_diff","WBSI_int_diff"), function(x){
  print(x)
  contrast_lm <- lm(get(x) ~ condition, mindful_data_WBSI) %>%
    emmeans(., "condition", contr = contrasts) %$% contrasts %>% 
    c(summary(.), apa_print(.))
  H1<-.2
  tibble(Scales          = sub("_diff", "", x),
         Active_mean     = mean(mental_WBSI[,x]) - mean(world_WBSI[,x]),
         Active_SE       = sqrt((sd(mental_WBSI[,x])/sqrt(NROW(mental_WBSI))) ^2+
                                  (sd(world_WBSI[,x])/sqrt(NROW(world_WBSI))) ^2),
         Active          = paste0(printnum(Active_mean)," (",printnum(Active_SE),")"),
         Active_BF       = Bf(round(Active_SE,2), round(Active_mean*-1,2), NROW(mental_WBSI) + NROW(world_WBSI) - 2, 
                              "normal", "normal", modeoftheory = 0, scaleoftheory = H1, dftheory = 102, tail = 1),
         "Active_t(119)" := contrast_lm$t.ratio[1],
         Active_p        = contrast_lm$p.value[1],
         
         Control_mean    = mean(world_WBSI[,x]) - mean(waitlist_WBSI[,x]),
         Control_SE      = sqrt((sd(world_WBSI[,x])/sqrt(NROW(world_WBSI))) ^2+
                                  (sd(waitlist_WBSI[,x])/sqrt(NROW(waitlist_WBSI))) ^2),
         Control          = paste0(printnum(Control_mean)," (",printnum(Control_SE),")"),
         Control_BF      = Bf(round(Control_SE,2), round(Control_mean*-1,2), NROW(world_WBSI) + NROW(waitlist_WBSI) - 2, 
                              "normal", "normal", modeoftheory = 0, scaleoftheory = H1, dftheory = 102, tail = 1),
         "Control_t(119)" = contrast_lm$t.ratio[2],
         Control_p       = contrast_lm$p.value[2]
  )}) %>% bind_rows %>% select(-ends_with("_mean"), -ends_with("_SE")) %>%
  dplyr::mutate(across(ends_with("_p"), printp), across(!ends_with("_p"), printnum))%>%
  write.csv("WBSI_split.csv")

# SE for Scales ----
#Posterior SD for the GAD-7 and PHQ-9 is the obtained SE diff for PHQ_anx and PHQ_dep respectively
#PWB is retained for the main study.
SE_mhealth <- function(scale_var){
  diff <- paste0(scale_var,"_diff") #variable
  SE_world <- sd(world[,diff])/sqrt(NROW(world))  
  SE_mental <- sd(mental[,diff])/sqrt(NROW(mental)) 
  round(sqrt(SE_mental^2+SE_world^2), digits = 2)
}

# Round2 ----
round2 <- function(x, digits = 0) {  # Function to always round 0.5 up
  posneg <- sign(x)
  z <- abs(x) * 10^digits
  z <- z + 0.5
  z <- trunc(z)
  z <- z / 10^digits
  z * posneg
}

