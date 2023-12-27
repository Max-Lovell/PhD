#Packages----
list.of.packages <- c("dplyr","magrittr","eeptools") 
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

#Clean raw dataset----
raw_data <- read.csv("combined_datasets_raw.csv")
  #str(raw_data, list.len=ncol(raw_data))

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

#remove and reorder variables
  #raw_data[,!colSums(is.na(raw_data)) == nrow(raw_data)]
col_order <- c(which(colnames(raw_data) == "email"),which(colnames(raw_data) == "Condition"),
               which(colnames(raw_data) == "age"),which(colnames(raw_data) == "Q6_pre"),
               which(colnames(raw_data) == "X1a"):which(colnames(raw_data) == "X1b"),
               which(colnames(raw_data) == "URN_pre"),
               which(colnames(raw_data) == "Q7_pre"):which(colnames(raw_data) == "Q14_18_pre"),
               which(colnames(raw_data) == "URN_post"),               
               which(colnames(raw_data) == "Q6_post"):which(colnames(raw_data) == "Q13_18_post"))
raw_data <- raw_data[,col_order]

#rename columns
demo <- c(which(colnames(raw_data) == "email"):which(colnames(raw_data) == "X1b"))
names(raw_data)[demo] <- c("email.ran_num","condition","age","gender","exp_anx","exp_dep")

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
  
#count missing data
m_cond <- sum(is.na(raw_data$condition))
m_email <- table(raw_data[is.na(raw_data$email.ran_num),"condition"])
m_pre <- table(raw_data[is.na(raw_data$URN_pre),"condition"])
m_post <- table(raw_data[is.na(raw_data$URN_post),"condition"])
pre_range <- c(which(colnames(raw_data) == "AWH_pre"):which(colnames(raw_data) == "PWB_18_pre"))
m_pre_item <- table(raw_data[!is.na(raw_data$URN_pre) & !complete.cases(raw_data[pre_range]),"condition"])
post_range <- c(which(colnames(raw_data) == "AWH_post"):which(colnames(raw_data) == "PWB_18_post"))
m_post_item <- table(raw_data[!is.na(raw_data$URN_post) & !complete.cases(raw_data[post_range]),"condition"])


merge(m_email,m_pre,m_pre_item,m_post,m_post_item, all.x)

merge(m_email,m_pre)



missing_data <-tribble(
  ~Variable,       ~Count,                                                    ~Condition,
  "Condition",      m_cond,                                                         "-",        
  "E-mail",         m_email,                                                          "Mental",   
  "All Pre-test",   filter(missing_pre, condition=="Mental") %>% pull(count),     "Mental",   
  "All Pre-test",   filter(missing_pre, condition=="World") %>% pull(count),     "World",    
  "All Pre-test",   filter(missing_pre, condition=="Waitlist") %>% pull(count),     "Waitlist",
  "Some Pre-test",  filter(missing_pre, condition=="Mental") %>% pull(count),   "Mental",
  "Some Pre-test",  filter(missing_pre, condition=="Waitlist") %>% pull(count), "Waitlist",
  "Some Post-test", filter(missing_post, condition=="Mental") %>% pull(count),  "Mental",
  "Some Post-test", filter(missing_post, condition=="World") %>% pull(count),   "World",
  "All Post-test",  filter(m_all_post, condition=="Mental") %>% pull(count),    "Mental",
  "All Post-test",  filter(m_all_post, condition=="World") %>% pull(count),     "World",
  "All Post-test",  filter(m_all_post, condition=="Waitlist") %>% pull(count),  "Waitlist"
)

