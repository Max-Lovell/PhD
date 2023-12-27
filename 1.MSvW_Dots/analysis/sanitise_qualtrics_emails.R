sanitise <- function(folder){ #folder<-'survey_data'
  files <- list.files(folder,full.names=T)
  emails <- c()
  for(f in files){ # f<-files[1]
    #read data
    csv <- read.csv(f)
    #extract email
    email <- csv[-c(1,2),colnames(csv) %in% c('Email', 'Email validation', 'RecipientEmail')]
    #convert to lowercase
    email_l <- sapply(email,tolower)
    #save lowercase back into original location
    csv[-c(1,2),colnames(csv) %in% c('Email', 'Email validation', 'RecipientEmail')] <- email_l
    #save converted data
    write.csv(csv,f)
    #save emails in list
    emails <- c(emails, email_l)
  }
  
  #get unique emails
  emails <- unname(unlist(emails))
  emails <- emails[emails!='']
  u_emails <- unique(emails)
  
  #generate random Ids
  random_ids <- c()
  while(length(random_ids)<length(u_emails)){
    random_ids <- unique(replicate(length(u_emails),paste(sample(c(LETTERS,0:9), 20, replace=T),collapse='')))
  }
  ref <- cbind(email=u_emails,id=random_ids)
  
  write.csv(ref,paste0(sub('/','-',folder),'-email_reference.csv'))
  
  #swamp emails for random ids
  for(f in files){ # f<-files[100]
    csv <- read.csv(f) #print(which(files==f))
    #dataset columns
    e_cols <- colnames(csv) %in% c('Email', 'Email validation', 'RecipientEmail')
    email <- unlist(csv[-c(1,2),e_cols]) #extract each email to replace
    email <- email[email != ""]
    #loop through actual emails and replace with ids
    for(e in email){ #e<-email[1]
      csv[csv==e] <- ref[ref[,'email']==e,'id']
    }
    write.csv(csv,f)
  }
  return(ref)
}

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dots_ref <- sanitise('survey_data')
