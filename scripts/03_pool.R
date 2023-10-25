setwd('emtpb/')
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(yaml)

config <- yaml.load_file("scripts/config.yaml")

# local paths
path <- config$r_local_path
path_calcs <- config$r_local_path_calcs

# server paths
#path <- config$r_server_path
#path_calcs <- config$r_server_path_calcs

t_test <- function(x, vadj) {
  x <- na.omit(x)
  n <- length(x)
  stat <- mean(x) / sqrt(var(x)/n*vadj)
  2 * pt(abs(stat), n - 1, lower.tail = FALSE)
}
# -> vadj=7.25 for 5 times 5-fold CV (details below). 
# Implementation reduces to standard t-test for vadj=1. 
# for repeated cv T times K-fold, it is vadj = 1+T*K/(K-1)

# get drug names
tmp1 <- read_csv(paste0(path,"data/PANCANCER_IC_Wed Feb 8 15_46_36 2023_GDSC1.csv"))
tmp2 <- read_csv(paste0(path,"data/PANCANCER_IC_Wed Feb 8 15_46_31 2023_GDSC2.csv"))
tmp_drugs <- full_join(tmp1, tmp2) %>%
  mutate(`drug` = paste0(`Drug ID`,"-",`Dataset Version`)) %>%
  dplyr::select(c("Drug Name","drug","Max Conc")) %>%
  distinct()

experiment_full <- read_csv(paste0(path_calcs,"metadata/paper/benchmark_paper_exp3.csv"), na = character())
experiment_second <- read_csv(paste0(path_calcs,"metadata/paper/benchmark_paper_exp6.csv"), na = character())
experiment_full <- bind_rows(experiment_full, experiment_second)

experiment <- experiment_full %>% dplyr::select(-c(...1,drug)) %>% distinct
write_csv(experiment, file = (paste0(path,"metadata/paper/benchmark_paper_exp3+6_03_pool.csv")))
overwrite <- FALSE

# iterate over cancer
pb <- txtProgressBar(min=1, max=nrow(experiment), style = 3, width = 20) 
for( which_run in 1:nrow(experiment)){ #1:nrow(experiment) # c(25,81,473)
  setTxtProgressBar(pb, which_run)
  #message(which_run)
  
  ###
  # costum
  #which_run <- 25
  # costum
  
  score <- experiment$score[which_run]
  run <- paste0("run",as.character(experiment$run[which_run]))
  model <- experiment$model_type[which_run]
  response_type <- experiment$response_type[which_run]
  cancertype_index <- experiment$cancer_type[which_run] + 1
  
  EMTscores <- suppressMessages(read_csv(paste0(path,"metadata/EMTscores","",".csv"))) # BENCHMARKING ORDER FOR GDSC CANCER TYPES
  cancertypes <- unique(c('PANCAN', unlist(EMTscores['TCGA Desc'])))
  EMTscores <- suppressMessages(read_csv(paste0(path,"metadata/EMTscores",score,".csv")))
  resp = suppressMessages(read_csv(paste0(path,"metadata/matrix_resp",response_type,".csv")) %>% dplyr::select(-c("COSMIC ID","TCGA Desc")))
  resp_help = suppressMessages(read_csv(paste0(path,"metadata/matrix_resp",response_type,".csv"))) # %>% dplyr::select(-c("COSMIC ID","TCGA Desc")))
  mut = suppressMessages(read_csv(paste0(path,"metadata/matrix_mut_",cancertypes[cancertype_index],".csv")))
  df <- tibble(drug = character(length = ncol(resp)),
               values_na = lapply(1:ncol(resp), list),
               values = lapply(1:ncol(resp), list),
               cancertype = rep(cancertypes[cancertype_index],ncol(resp)),
               .rows = ncol(resp))
  
  results_path <- paste0(path_calcs,"metadata/summaries")
  dir.create(results_path)
  ###
  
  #CHANGING RESULTS PATH
  #save_path <- paste0(results_path,"/",cancertypes[cancertype_index],"_performances.rds")
  save_path <- paste0(results_path,"/",
                      which_run,"_",
                      "performances_",
                      cancertypes[cancertype_index],"_",
                      model,"_",
                      score,"_",
                      response_type,
                      "_v3.rds")
  if(!overwrite & file.exists(save_path)){
    #message(paste0("File for cancertype ",cancertypes[cancertype_index]," exists, skipping.."))
    message("skipping..")
  }else{
    #df <- tibble(drug = character(length = ncol(resp)),
    #             values_na = lapply(1:ncol(resp), list),
    #             values = lapply(1:ncol(resp), list),
    #             effect = lapply(1:ncol(resp), list),
    #             cancertype = rep(cancertypes[cancertype_index],ncol(resp)),
    #             .rows = ncol(resp))
    # iteratre over run
    #pb <- txtProgressBar(min=1, max=ncol(resp), style = 3, width = 20) 
    for( which in 1:ncol(resp)){
      #which <- 1
      #setTxtProgressBar(pb, which)
      preds_dir <- paste0(path_calcs,"metadata/",run,"/predictions/",cancertypes[cancertype_index],"/")
      files <- list.files(path = preds_dir, pattern = paste0("._",as.character(which-1),"_."), full.names = TRUE)
      rm(model_false); rm(model_true)
      if(length(files)!=0){ # if running files do not exist, because of missing cancertype in secrier
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
      }else{
        model_false <- as.data.frame(NA)
      }
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
      if(model == "eln"){
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
        # no effect for eln
        p <- NA
      }
      
      if(model %in% c("grf","grfo")){
        if(ncol(model_false)<2){
          l <- c('false'=NA,"true"=NA)
          l_na <- c('false'=NA,"true"=NA)
          p <- NA
        }else{
          model_false_h <- model_false %>%
            group_by(repeatfold) %>%
            summarize(cor_na = suppressWarnings(cor(preds, truth)), .groups = 'drop') %>% 
            mutate(cor = replace_na(cor_na, 0))
          #model_true <- model_true %>%
          #  group_by(repeatfold) %>%
          #  summarize(cor_na = suppressWarnings(cor(preds, truth)), .groups = 'drop') %>% 
          #  mutate(cor = replace_na(cor_na, 0))
          
          l = list()
          l[['false']] = model_false_h$cor
          l[['true']] = model_true$cor
          
          l_na = list()
          l_na[['false']] = model_false_h$cor_na
          l_na[['true']] = model_true$cor_na
          
          # get effect (p)
          model_false_hh <- model_false %>% group_by(repeatfold) %>% summarize(p_raw = unique(effects))
          tmp <- lapply(model_false_hh$p_raw, function(x) eval(parse(text = tryCatch(gsub("\\]",")",gsub("\\[","list(",(gsub("\\)list","),list",gsub("\\(list\\(",",list(list(",gsub("\\(\\[","list(",gsub("\\]\\)",")",gsub("array\\(","\\(",gsub(" ","",gsub("\n","",x)))))))))), error = function(e) x))))
          tmp <- map(tmp, ~ map(., matrix, nrow = 1, byrow = TRUE))
          model_false_hh$p <- tmp
          #model_false_hh$p <- lapply(model_false_hh$p_raw, function(x) tryCatch(as.numeric(strsplit(gsub(",,",",",gsub(",,",",",gsub("\\)",",",gsub("\\(",",",gsub("\\]",",",gsub("\\[",",",gsub(" ",",",x))))))),",")[[1]][-1]), error = function(e) x))
          p <- list(model_false_hh$p)#, model_false_hh$p_raw)
        }
      }
      
      
      df[which, "values"][[1]] <- list(l)
      df[which, "values_na"][[1]] <- list(l_na)
      df[which,"drug"] <- colnames(resp)[which]
      df[which,"effect"] <- list(p)
    }
    
    
    ### some analysis on the whole df (cors, fraction_extrapolated)
    performances <- df
    performances$delta <- unlist(lapply(1:nrow(performances),function(i) tryCatch(mean(performances[i,]$values[[1]]$false-performances[i,]$values[[1]]$true), error = function(e) NA)))
    performances$cor_emt <- unlist(lapply(1:nrow(performances),function(i) tryCatch(mean(performances[i,]$values[[1]]$false), error = function(e) NA)))
    performances$cor_mut <- unlist(lapply(1:nrow(performances),function(i) tryCatch(mean(performances[i,]$values[[1]]$true), error = function(e) NA)))
    performances$cor_emt_sd <- unlist(lapply(1:nrow(performances),function(i) tryCatch(sd(performances[i,]$values[[1]]$false), error = function(e) NA)))
    performances$cor_mut_sd <- unlist(lapply(1:nrow(performances),function(i) tryCatch(sd(performances[i,]$values[[1]]$true), error = function(e) NA)))
    performances$pvalue <- unlist(lapply(1:nrow(performances),function(i) tryCatch(t_test(performances[i,]$values[[1]]$false-performances[i,]$values[[1]]$true, vadj = 7.25), error = function(e) NA)))
    performances$TCGA <- cancertypes[cancertype_index]
    performances$Drug <- performances$drug
    performances$label <- ""
    performances$drugname <- unlist(lapply(1:nrow(performances), function(x){tmp_drugs$`Drug Name`[tmp_drugs$drug == performances[x,]$drug]}))
    performances <- left_join(x = performances, y = tmp_drugs %>% dplyr::select(-c("Max Conc")))
    
    fraction_extrapolated <- c()
    if(cancertypes[cancertype_index] == "PANCAN"){
      j <- full_join(EMTscores, mut)
    }else{
      j <- full_join(EMTscores[EMTscores$`TCGA Desc` == cancertypes[cancertype_index],], mut)
    }
    for( i in 1:nrow(performances)){
      #print(i) # do not print the 700 drugs for each of the 300k runs....
      if(cancertypes[cancertype_index] == "PANCAN"){
        df <- full_join(
          resp_help[,c(1,i+2),drop = F], j, by = "COSMIC ID") %>% 
          distinct() %>%
          dplyr::select(-c(`COSMIC ID`,`TCGA Desc`)) %>%
          na.omit()
      }else{
        df <- full_join(
          resp_help[resp_help$`TCGA Desc` == cancertypes[cancertype_index],c(1,i+2),drop = F], j, by = "COSMIC ID") %>% 
          distinct() %>%
          dplyr::select(-c(`COSMIC ID`,`TCGA Desc`)) %>%
          na.omit()
      }
      maxc <- tmp_drugs$`Max Conc`[tmp_drugs$drug == colnames(df)[1]]
      fraction_extrapolated[which(performances$drug == colnames(df)[1])] <- length(which((exp(unlist(df[1])) > 2*maxc)))/nrow(df)
    }; performances$fraction_extrapolated <- fraction_extrapolated
    ### some analysis on the whole df (cors, fraction_extrapolated)
    ###############################################################
    
    saveRDS(performances, file = save_path)
    
    #if(max(unlist(lapply(df$values, length))) == 4){ # for grf if 'values' column has CIs
    #  df <- unnest_wider(df, col = values)
    #  saveRDS(df, file = save_path)
    #}
  }
}




































if(F){
      ###
      # OLD version
      ####################################################
      # iterate over cancer
      #
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
        
        results_path <- paste0(path_calcs,"metadata/",run,"/results")
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
      ####################################################
}