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
  filter(!duplicated(email) & !email %in% c("1@a.com", "redacted@redacted.co.uk"))

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

combined_datasets_raw <- combined_datasets_raw %>% select(-c("Q12_pre", "Q13_pre",
                                                             "Q11_post", "Q12_post"))

write.csv(combined_datasets_raw,file="combined_datasets_raw.csv")


