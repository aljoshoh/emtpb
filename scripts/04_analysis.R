# local paths
path <- "/Users/alexander.ohnmacht/research/marisa/emtpb/"
path_calcs <- "/Volumes/pheb/lustre/groups/cbm01/code/alexander.ohnmacht/emtpb/"

# server paths
#path <- "/lustre/groups/cbm01/code/alexander.ohnmacht/emtpb/"
#path_calcs <- path

library(dplyr)
# Import summaries
files <- list.files(path = paste0(path_calcs,"metadata/summaries"),full.names = TRUE, pattern = "_performances") #  pattern = "eln", 
#files_grf <- list.files(path = paste0(path_calcs,"metadata/summaries"), pattern = "grf", full.names = TRUE)



names <- c("drug","drugname","TCGA","pvalue","fdr","effectsize","delta","cor_emt","cor_mut","cor_emt_sd","cor_mut_sd","method","resp_type","score","fraction_extrapolated","kept_extrapolated")
ll <- list()
for( i in 1:length(files) ){ ## c(168,429)
  print(i)
  extrapolated_threshold <- 0.70 #0.65
  tmp <- readRDS(files[i]) #SKCM # tmp <- readRDS(files_eln[1]) #PANCAN

  if(grepl("eln",files[i], fixed = TRUE)&(!grepl("auc",files[i], fixed = TRUE))){
    tmp <- tmp %>% 
      dplyr::mutate(kept_extrapolated = fraction_extrapolated < extrapolated_threshold) %>%
      group_by(kept_extrapolated) %>%
      dplyr::mutate(fdr = p.adjust(pvalue, method = "BH")) %>%
      ungroup %>%
      dplyr::mutate(fdr = ifelse(kept_extrapolated, fdr, NA)) %>%
      dplyr::mutate(effectsize = NA) %>%
      dplyr::mutate(method = "eln") %>%
      dplyr::mutate(resp_type = "ic50")
  }
  
  if(grepl("eln",files[i], fixed = TRUE)&(grepl("auc",files[i], fixed = TRUE))){
    tmp <- tmp %>% 
      dplyr::mutate(kept_extrapolated = fraction_extrapolated < extrapolated_threshold) %>%
      #group_by(kept_extrapolated) %>%
      dplyr::mutate(fdr = p.adjust(pvalue, method = "BH")) %>%
      #ungroup %>%
      #dplyr::mutate(fdr = ifelse(kept_extrapolated, fdr, NA)) %>%
      dplyr::mutate(effectsize = NA) %>%
      dplyr::mutate(method = "eln") %>%
      dplyr::mutate(resp_type = "auc")
  }
  
  if(grepl("grf",files[i], fixed = TRUE)&(!grepl("auc",files[i], fixed = TRUE))){
    tmp$effect_p <- lapply(tmp$effect, function(y) unlist(lapply(y, function(x) unique(x)[2])))
    tmp$pvalue <- unlist(lapply(tmp$effect_p, function(x) as.numeric(mean(x)))) # harmonicmeanp::hmp.stat
    tmp$effectsize <- unlist(lapply(tmp$effect, function(y) mean(unlist(lapply(y, function(x) unique(x)[5])))))
    tmp <- tmp %>% 
      dplyr::mutate(kept_extrapolated = fraction_extrapolated < extrapolated_threshold) %>%
      group_by(kept_extrapolated) %>%
      dplyr::mutate(fdr = p.adjust(pvalue, method = "BH")) %>%
      ungroup %>%
      dplyr::mutate(fdr = ifelse(kept_extrapolated, fdr, NA)) %>%
      #dplyr::mutate(effectsize = cor_emt - cor_mut) %>%
      dplyr::mutate(method = "grf") %>%
      dplyr::mutate(resp_type = "ic50")
  }
  
  if(grepl("grf",files[i], fixed = TRUE)&(grepl("auc",files[i], fixed = TRUE))){
    tmp$effect_p <- lapply(tmp$effect, function(y) unlist(lapply(y, function(x) unique(x)[2])))
    tmp$pvalue <- unlist(lapply(tmp$effect_p, function(x) as.numeric(mean(x)))) # harmonicmeanp::hmp.stat
    tmp$effectsize <- unlist(lapply(tmp$effect, function(y) mean(unlist(lapply(y, function(x) unique(x)[5])))))
    tmp <- tmp %>% 
      dplyr::mutate(kept_extrapolated = fraction_extrapolated < extrapolated_threshold) %>%
      #group_by(kept_extrapolated) %>%
      dplyr::mutate(fdr = p.adjust(pvalue, method = "BH")) %>%
      #ungroup %>%
      #dplyr::mutate(fdr = ifelse(kept_extrapolated, fdr, NA)) %>%
      #dplyr::mutate(effectsize = cor_emt - cor_mut) %>%
      dplyr::mutate(method = "grf") %>%
      dplyr::mutate(resp_type = "auc")
  }
  
  if(grepl("tan_088",files[i], fixed = TRUE)){
    tmp$score <- "tan"
  }else{
    if(grepl("secrier_065",files[i], fixed = TRUE)){
      tmp$score <- "secrier"
    }else{
      if(grepl("gsva",files[i], fixed = TRUE)){
        tmp$score <- "gsva"
      }else{
        tmp$score <- "mak"
      }
    }
  }
  
  tmp <- tmp[,names]
  ll[[i]] <- tmp
}
ll_df <- do.call(rbind, ll)

saveRDS(object = ll_df, paste0(path_calcs,"metadata/summaries/PERFORMANCES_v1.rds"))
saveRDS(object = ll_df, paste0(path,"metadata/summaries/PERFORMANCES_v1.rds"))









