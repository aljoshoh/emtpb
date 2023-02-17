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
which_run <- 25
# costum

# iterate over cancer
#for( which_run in 1:length(cancertypes)){
save_path <- paste0(results_path,"/",cancertypes[which_run],"_performances.rds")
performances <- read_rds(save_path)
performances$delta <- unlist(lapply(1:nrow(performances),function(i) tryCatch(mean(performances[i,]$values[[1]]$false-performances[100,]$values[[1]]$true), error = function(e) NA)))

