setwd('emtpb/')
library(readr)
library(dplyr)
library(tidyr)
library(tibble)

EMTscores <- read_csv("metadata/EMTscores.csv")
cancertypes <- unique(c('PANCAN', unlist(EMTscores['TCGA Desc'])))
run <- "run2"
resp = read_csv("metadata/matrix_resp.csv") %>% dplyr::select(-c("COSMIC ID","TCGA Desc"))
df <- tibble(drug = character(length = ncol(resp)),
             values_na = lapply(1:ncol(resp), list),
             values = lapply(1:ncol(resp), list),
             cancertype = rep(cancertypes[which_run],ncol(resp)),
             .rows = ncol(resp)) <- paste0("metadata/",run,"/results")
dir.create(results_path)
overwrite <- FALSE
# costum
#which_run <- 4
# costum


# iterate over cancer
for( which_run in 1:length(cancertypes)){
  save_path <- paste0(results_path,"/",cancertypes[which_run],"_performances.csv")
  if(overwrite & file.exists(save_path)){
    message(paste0("File for cancertype ",cancertypes[which_run]," exists, skipping.."))
  }else{
    df <- tibble(drug = character(length = ncol(resp)),
                 values_na = lapply(1:ncol(resp), list),
                 values = lapply(1:ncol(resp), list),
                 cancertype = rep(cancertypes[which_run],ncol(resp)),
                 .rows = ncol(resp))
    # iteratre over run
    for( which in 1:ncol(resp)){
      which <- 1
      preds_dir <- paste0("metadata/",run,"/predictions/",cancertypes[which_run],"/")
      files <- list.files(path = preds_dir, pattern = paste0("._",as.character(which-1),"_."), full.names = TRUE)
      
      model_false <- read_csv(files[1])
      model_true <- read_csv(files[2])
      if((ncol(model_false)<2) | (ncol(model_true)<2)){
        l <- c('false'=NA,"true"=NA)
        l_na <- c('false'=NA,"true"=NA)
      }else{
        model_false <- model_false %>%
          group_by(repeatfold) %>%
          summarize(cor_na = cor(preds, truth)) %>% 
          mutate(cor = replace_na(cor_na, 0))
        model_true <- model_true %>%
          group_by(repeatfold) %>%
          summarize(cor_na = cor(preds, truth)) %>% 
          mutate(cor = replace_na(cor_na, 0))
        
        l = list()
        l[['false']] = model_false$cor
        l[['true']] = model_true$cor
       
        
        l_na = list()
        l_na[['false']] = model_false$cor_na
        l_na[['true']] = model_true$cor_na
      }
      
      df[which, "values"][[1]] <- list(l)
      df[which, "values_na"][[1]] <- list(l_na)
      df[which,"drug"] <- colnames(resp)[which]
      
      write_csv(df, path = save_path)
    }
  }
}
