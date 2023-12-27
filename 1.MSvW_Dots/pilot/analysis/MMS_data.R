# Packages ----
list.of.packages <- c("eeptools", "magrittr", "tidyverse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

# Pre-test ----
#CSV data exported fromm onlinesurveys.ac.uk. tick: "include unique response number", 
  #"date and time response started"/"submitted", "code responses","combine scale/rank values".
pre_test <- list.files(recursive = TRUE, pattern = "^first_201(.*)csv$") %>% 
  lapply(., read.csv) %>% bind_rows()

names(pre_test) <- paste0(names(pre_test),"_pre")

pre_test <- pre_test %>%
  mutate_all(list(tolower)) %>%
  dplyr::rename(email = Q2_pre) %>%
  mutate(email = trimws(email)) %>%
  filter(email != "a" ) %>%
  mutate(email = dplyr::recode(email, #rename matchable emails
                        'redacted@gmaill.com'='redacted@gmail.com'))



# Post-test ----
#Note: second_2019 (about N=35) is missing the WBSI, which is usually found in Q12. 
  #manually changed all Q12 to Q13 and added in new Q12 (values ITT LOCF imputation from session 1)
post_test <- list.files(recursive = TRUE, pattern = "^second_201(.*)csv$") %>%
  lapply(., read.csv) %>% bind_rows()

names(post_test) <- paste0(names(post_test),"_post")

post_test <- post_test %>%
  mutate_all(list(tolower)) %>%
  dplyr::rename(email = Q1_post) %>%
  mutate(email = trimws(email)) %>%
  filter(!duplicated(email) & !email %in% c("1@a.com", "[redacted]@sussex.ac.uk"))

#Combine datasets
both_sessions <- merge(pre_test, post_test, by="email", all= T)

# Emails ----
email_data <- list.files(recursive = TRUE, pattern = "*_emails.csv$") %>%
  lapply(., read.csv) %>% bind_rows()

#Email Data
email_data <- email_data %>%
  mutate_all(list(tolower)) %>%
  mutate(Participant.email = trimws(Participant.email))

#keep separate for counting missing data                
email_data2 <- email_data %>%
  filter(!duplicated(Participant.email)) %>%
  dplyr::rename(email = Participant.email) %>%
  mutate(email = dplyr::recode(email,
                        'redacted@gmaill.com'='redacted@gmail.com',
                        ))

combined_datasets_raw <- merge(both_sessions,email_data2,by="email",all=T)

combined_datasets_raw["email"] <- ifelse(combined_datasets_raw$email=="","",
                                       sample(1:1000, nrow(combined_datasets_raw)))
write.csv(combined_datasets_raw,file="combined_datasets_raw.csv")
# Clean data ---------------------------
#Recode factors for counting
combined_datasets_raw <- combined_datasets_raw %>%
  mutate(Condition = dplyr::recode(Condition,'mental states'='Mental','mind'='Mental','mental'='Mental',
                                   'world'='World','control'='Waitlist'),
         Q6_pre = dplyr::recode(Q6_pre,'female'='Female','female '='Female','cis female'='Female','f'='Female',
                                'fem'='Female','woman'='Female','male'='Male','m'='Male','male '='Male'))

#separate to count missing data
combined_datasets <- combined_datasets_raw %>%
  filter(!is.na(Participant.email) & Participant.email != "") %>%
  filter(!is.na(Condition), !Condition == "", !is.na(URN_pre))  %>%
  mutate(StartDate_pre = as.Date(substr(StartDate_pre, 1,10)),
         dob      = as.Date(paste0(Q5_pre,"-", Q4_pre,"-", Q3_pre)),
         age = age_calc(dob, enddate =  StartDate_pre, units = "years", precise = F))  %>%
  select_if(function(x){!all(is.na(x))})  %>% #remove empty variables
  dplyr::select(-c(email, Q1_1_pre, Q1_2_pre, Q3_pre,	Q4_pre,	Q5_pre, Q12_pre,	Q13_pre, StartDate_pre, 
                   CompletionDate_pre, Q11_post,	Q12_post, StartDate_post, CompletionDate_post, Q2_post, 
                   Q3_post, Q4_post, Q5_post, X, X1c, X2a, X2b,	X3a, X3b, X4a, X4b, X5a, X5b,	X6a, X6b,	
                   X7a, X7b, X8a, X8b,	X9a, X9b, X10a,	X10b, dob)) %>%
  rename_with(~gsub("_a_", "_", .)) %>%
  rename_with(~gsub("Q7|Q6","AWH",.), grep("(^Q7_.*pre$)|(^Q6_.*post$)",names(.))) %>%
  rename_with(~gsub("Q8|Q7","AMH",.), grep("(^Q8_.*pre$)|(^Q7_.*post$)",names(.))) %>%
  rename_with(~gsub("Q9|Q8","FFMQ",.), grep("(^Q9_.*_pre$)|(^Q8_.*_post$)",names(.))) %>%
  rename_with(~gsub("Q10|Q9","PSS",.), grep("(^Q10_.*_pre$)|(^Q9_.*_post$)",names(.))) %>%
  rename_with(~gsub("Q11|Q10","PHQ",.), grep("(^Q11_.*_pre$)|(^Q10_.*_post$)",names(.))) %>%
  rename_with(~gsub("Q12|Q11","RRS",.), grep("(^Q12_.*_pre$)|(^Q11_.*_post$)",names(.))) %>%
  rename_with(~gsub("Q13|Q12","WBSI",.), grep("(^Q13_.*_pre$)|(^Q12_.*_post$)",names(.))) %>%
  rename_with(~gsub("Q14|Q13","PWB",.), grep("(^Q14_.*_pre$)|(^Q13_.*_post$)",names(.))) %>%
  dplyr::rename(URN = URN_pre, gender = Q6_pre, condition = Condition, expect_anx  = X1a, expect_dep = X1b) %>%
  mutate(condition = dplyr::recode(condition,'mental states'='Mental','mind'='Mental','mental'='Mental',
                            'world'='World','control'='Waitlist'),
         gender = dplyr::recode(gender,'female'='Female','female '='Female','cis female'='Female','f'='Female',
                         'fem'='Female','woman'='Female','male'='Male','m'='Male','male '='Male'))

write.csv(combined_datasets, file = "mindful_data.csv")




#Missing data----
m_email <- email_data[email_data$Participant.email == "",]
m_cond <- sum(is.na(combined_datasets_raw$Condition) | combined_datasets_raw$Condition == "")
m_pre <- combined_datasets_raw[is.na(combined_datasets_raw$URN_pre),] %>%
  group_by(Condition) %>% summarise(count = length(Condition))

tribble(
  ~Variable, ~Count, ~Condition,
  "Condition", m_cond, "-",
  "E-mail", nrow(m_email), "Mental",
  "All Pre-test", filter(m_pre, Condition=="Mental") %>% pull(count), "Mental",
  "All Pre-test", filter(m_pre, Condition=="World") %>% pull(count), "World",
  "All Pre-test", filter(m_pre, Condition=="Waitlist") %>% pull(count), "Waitlist"
)

#initial intervention allocation
combined_datasets_raw %>% filter(!is.na(Condition), !Condition == "") %>% group_by(Condition) %>% summarise
