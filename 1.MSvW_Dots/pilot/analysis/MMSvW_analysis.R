#Packages----
list.of.packages <- c("car", "pastecs", "stats", "data.table", "emmeans", "Rmisc", 
                      "tidyverse", "dplyr","magrittr","psych","papaja","mice") 
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
#Clean raw dataset----
#initial intervention allocation
combined_datasets_raw %>% filter(!is.na(Condition), !Condition == "") %>% group_by(Condition) %>% summarise

mindful_data <- read.csv("mindful_data.csv") %>% 
  mutate_if(is.character, factor) %>% 
  dplyr::mutate(condition = fct_relevel(condition, "Waitlist", "World", "Mental"))



#select pre-test and post-test variables
pre_range <- c(which(colnames(mindful_data) == "AWH_pre"):which(colnames(mindful_data) == "PWB_18_pre"))
post_range <- c(which(colnames(mindful_data) == "AWH_post"):which(colnames(mindful_data) == "PWB_18_post"))

missing_pre <- mindful_data[!complete.cases(mindful_data[pre_range]),]%>%
  group_by(condition) %>% summarise(count = length(condition)) #save count for later

