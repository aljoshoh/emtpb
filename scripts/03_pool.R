#setwd('emtpb/')
library(readr)
library(dplyr)
library(tidyr)
library(tibble)

experiment <- read_csv("metadata/paper/benchmark_paper_exp1.csv", na = character())
overwrite <- FALSE


# iterate over cancer
pb <- txtProgressBar(min=1, max=nrow(experiment), style = 3, width = 20) 
for( which_run in 1:nrow(experiment)){
  
  ###
  # costum
  #which_run <- 25
  # costum
  
  score <- experiment$score[which_run]
  run <- paste0("run",as.character(experiment$run[which_run]))
  model <- experiment$model_type[which_run]
  response_type <- experiment$response_type[which_run]
  cancertype_index <- experiment$cancer_type[which_run+1]
  
  EMTscores <- suppressMessages(read_csv(paste0("metadata/EMTscores",score,".csv")))
  cancertypes <- unique(c('PANCAN', unlist(EMTscores['TCGA Desc'])))
  resp = suppressMessages(read_csv("metadata/matrix_resp.csv") %>% dplyr::select(-c("COSMIC ID","TCGA Desc")))
  df <- tibble(drug = character(length = ncol(resp)),
               values_na = lapply(1:ncol(resp), list),
               values = lapply(1:ncol(resp), list),
               cancertype = rep(cancertypes[cancertype_index],ncol(resp)),
               .rows = ncol(resp))
  
  results_path <- paste0("metadata/",run,"/results")
  dir.create(results_path)
  ###
  
  
  setTxtProgressBar(pb, which_run)
  message(which_run)
  #CHANGING RESULTS PATH
  #save_path <- paste0(results_path,"/",cancertypes[cancertype_index],"_performances.rds")
  save_path <- paste0(results_path,"/",
                      "performances_",
                      cancertypes[cancertype_index],"_",
                      model,"_",
                      score,"_",
                      response_type,
                      ".rds")
  if(overwrite & file.exists(save_path)){
    message(paste0("File for cancertype ",cancertypes[cancertype_index]," exists, skipping.."))
  }else{
    df <- tibble(drug = character(length = ncol(resp)),
                 values_na = lapply(1:ncol(resp), list),
                 values = lapply(1:ncol(resp), list),
                 cancertype = rep(cancertypes[cancertype_index],ncol(resp)),
                 .rows = ncol(resp))
    # iteratre over run
    #pb <- txtProgressBar(min=1, max=ncol(resp), style = 3, width = 20) 
    for( which in 1:ncol(resp)){
      #which <- 1
      #setTxtProgressBar(pb, which)
      preds_dir <- paste0("metadata/",run,"/predictions/",cancertypes[cancertype_index],"/")
      files <- list.files(path = preds_dir, pattern = paste0("._",as.character(which-1),"_."), full.names = TRUE)
      rm(model_false); rm(model_true)
      model_false <- read_csv(files[1], 
                              col_types = cols(
                                `COSMIC ID` = col_double(),
                                preds = col_double(),
                                truth = col_double(),
                                `repeat` = col_double(),
                                folds = col_double(),
                                repeatfold = col_character()
                              )
                              )
      model_true <- tryCatch(read_csv(files[2], 
                             col_types = cols(
                               `COSMIC ID` = col_double(),
                               preds = col_double(),
                               truth = col_double(),
                               `repeat` = col_double(),
                               folds = col_double(),
                               repeatfold = col_character()
                             )
                             ), error = function(e) data.frame(NA) )
      if((ncol(model_false)<2) | (ncol(model_true)<2)){
        l <- c('false'=NA,"true"=NA)
        l_na <- c('false'=NA,"true"=NA)
      }else{
        model_false <- model_false %>%
          group_by(repeatfold) %>%
          summarize(cor_na = suppressWarnings(cor(preds, truth)), .groups = 'drop') %>% 
          mutate(cor = replace_na(cor_na, 0))
        model_true <- model_true %>%
          group_by(repeatfold) %>%
          summarize(cor_na = suppressWarnings(cor(preds, truth)), .groups = 'drop') %>% 
          mutate(cor = replace_na(cor_na, 0))
        
        l = list()
        l[['false']] = model_false$cor
        l[['true']] = model_true$cor
       
        
        l_na = list()
        l_na[['false']] = model_false$cor_na
        l_na[['true']] = model_true$cor_na
      }
      
      if(all(is.na(model_true))){
        l <- unique(model_false$repeatfold)
        l <- tryCatch(as.numeric(strsplit(gsub(",,",",",gsub(",,",",",gsub("\\)",",",gsub("\\(",",",gsub("\\]",",",gsub("\\[",",",gsub(" ",",",l))))))),",")[[1]][-1]), error = function(e) l)
      }
      
      df[which, "values"][[1]] <- list(l)
      df[which, "values_na"][[1]] <- list(l_na)
      df[which,"drug"] <- colnames(resp)[which]
    }
    saveRDS(df, file = save_path)
    
    if(max(unlist(lapply(df$values, length))) == 4){ # for grf if 'values' column has CIs
      df <- unnest_wider(df, col = values)
      saveRDS(df, file = save_path)
    }
  }
}


# for grf models:

# 


