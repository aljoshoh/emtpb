# IS THIS DEPRECATED ?

setwd('emtpb/')
library(readr)
library(dplyr)
library(tidyr)
library(tibble)

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
tmp1 <- read_csv("data/PANCANCER_IC_Wed Feb 8 15_46_36 2023_GDSC1.csv")
tmp2 <- read_csv("data/PANCANCER_IC_Wed Feb 8 15_46_31 2023_GDSC2.csv")
tmp_drugs <- full_join(tmp1, tmp2) %>%
  mutate(`drug` = paste0(`Drug ID`,"-",`Dataset Version`)) %>%
  dplyr::select(c("Drug Name","drug","Max Conc")) %>%
  distinct()



cancertypes <- unique(c('PANCAN', unlist(EMTscores['TCGA Desc'])))
run <- "run5"
which_score <- ""#"gsva" # "" 
results_path <- paste0("metadata/",run,"/results")
dir.create(results_path)

# costum
which_run <- 25
# costum

EMTscores <- read_csv(paste0("metadata/EMTscores",which_score,".csv"))
resp = read_csv("metadata/matrix_resp.csv") #%>% dplyr::select(-c("COSMIC ID","TCGA Desc"))

# iterate over cancer
for( which_run in 1:length(cancertypes)){
  mut = read_csv(paste0("metadata/matrix_mut_",cancertypes[which_run],".csv"))
  save_path <- paste0(results_path,"/",cancertypes[which_run],"_performances.rds")
  performances <- read_rds(save_path)
  performances$delta <- unlist(lapply(1:nrow(performances),function(i) tryCatch(mean(performances[i,]$values[[1]]$false-performances[i,]$values[[1]]$true), error = function(e) NA)))
  performances$cor_emt <- unlist(lapply(1:nrow(performances),function(i) tryCatch(mean(performances[i,]$values[[1]]$false), error = function(e) NA)))
  performances$cor_mut <- unlist(lapply(1:nrow(performances),function(i) tryCatch(mean(performances[i,]$values[[1]]$true), error = function(e) NA)))
  performances$cor_emt_sd <- unlist(lapply(1:nrow(performances),function(i) tryCatch(sd(performances[i,]$values[[1]]$false), error = function(e) NA)))
  performances$cor_mut_sd <- unlist(lapply(1:nrow(performances),function(i) tryCatch(sd(performances[i,]$values[[1]]$true), error = function(e) NA)))
  performances$pvalue <- unlist(lapply(1:nrow(performances),function(i) tryCatch(t_test(performances[i,]$values[[1]]$false-performances[i,]$values[[1]]$true, vadj = 7.25), error = function(e) NA)))
  performances$TCGA <- cancertypes[which_run]
  performances$Drug <- performances$drug
  performances$label <- ""
  performances$drugname<- unlist(lapply(1:nrow(performances), function(x){tmp_drugs$`Drug Name`[tmp_drugs$drug == performances[x,]$drug]}))
  performances <- left_join(x = performances, y = tmp_drugs %>% dplyr::select(-c("Max Conc")))
  
  fraction_extrapolated <- c()
  if(cancertypes[which_run] == "PANCAN"){
    j <- full_join(EMTscores, mut)
  }else{
    j <- full_join(EMTscores[EMTscores$`TCGA Desc` == cancertypes[which_run],], mut)
  }
  for( i in 1:nrow(performances)){
    if(cancertypes[which_run] == "PANCAN"){
      df <- full_join(
        resp[,c(1,i+2),drop = F], j, by = "COSMIC ID") %>% 
        distinct() %>%
        dplyr::select(-c(`COSMIC ID`,`TCGA Desc`)) %>%
        na.omit()
    }else{
      df <- full_join(
          resp[resp$`TCGA Desc` == cancertypes[which_run],c(1,i+2),drop = F], j, by = "COSMIC ID") %>% 
        distinct() %>%
        dplyr::select(-c(`COSMIC ID`,`TCGA Desc`)) %>%
        na.omit()
    }
    maxc <- tmp_drugs$`Max Conc`[tmp_drugs$drug == colnames(df)[1]]
    fraction_extrapolated[which(performances$drug == colnames(df)[1])] <- length(which((exp(unlist(df[1])) > 2*maxc)))/nrow(df)
  }; performances$fraction_extrapolated <- fraction_extrapolated
  
  save_path <- paste0(results_path,"/",cancertypes[which_run],"_performances_v2.rds")
  saveRDS(performances, file = save_path)
}
