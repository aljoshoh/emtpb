# local paths
path <- "/Users/alexander.ohnmacht/research/marisa/emtpb/"
path_calcs <- "/Volumes/pheb/lustre/groups/cbm01/code/alexander.ohnmacht/emtpb/"

# server paths
#path <- "/lustre/groups/cbm01/code/alexander.ohnmacht/emtpb/"
#path_calcs <- path

library(dplyr)
myfnc1 <- function(arg1){
  if(all(is.na(arg1))){
    temp <- NA
  }else{
    temp <- ((lapply(arg1, function(x) unlist(x[[2]]))))
  }
  return(temp)
}
myfnc2 <- function(arg1,arg2){
  if(all(is.na(arg1))){
    temp <- NA
  }else{
    temp <- ((lapply(arg1, function(x) unlist(x[[3]][[arg2]]))))
  }
  return(temp)
}

# Import summaries
files <- list.files(path = paste0(path_calcs,"metadata/summaries"),full.names = TRUE, pattern = "_performances.*v3") #  pattern = "eln", 
#files_grf <- list.files(path = paste0(path_calcs,"metadata/summaries"), pattern = "grf", full.names = TRUE)



names <- c("drug","drugname","TCGA","pvalue","fdr","effectsize","effectsize_mean","lower","upper","inference","delta","cor_emt","cor_mut","cor_emt_sd","cor_mut_sd","method","resp_type","score","fraction_extrapolated","kept_extrapolated")
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
      dplyr::mutate(inference = NA) %>%
      dplyr::mutate(effectsize_mean = NA) %>%
      dplyr::mutate(effectsize = NA) %>%
      dplyr::mutate(upper = NA) %>%
      dplyr::mutate(lower = NA) %>%
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
      dplyr::mutate(inference = NA) %>%
      dplyr::mutate(effectsize_mean = NA) %>%
      dplyr::mutate(effectsize = NA) %>%
      dplyr::mutate(upper = NA) %>%
      dplyr::mutate(lower = NA) %>%
      dplyr::mutate(method = "eln") %>%
      dplyr::mutate(resp_type = "auc")
  }
  
  if(grepl("grf",files[i], fixed = TRUE)&(!grepl("auc",files[i], fixed = TRUE))){
    tmp$effect_p <- lapply(tmp$effect, function(y) lapply(y, function(x) unlist(x[[1]])))
    tmp$pvalue <- lapply(tmp$effect_p, function(y) mean(unlist(lapply(y, function(x) unique(x)))))
    tmp$inference <- tmp$effect
    tmp$effectsize_mean <- lapply(tmp$effect, function(y) mean(unlist(myfnc1(y))))
    tmp$effectsize <- lapply(tmp$effect, function(y)  myfnc1(y))
    tmp$lower <- lapply(tmp$effect, function(y)  mean(unlist(myfnc2(y,1))))
    tmp$upper <- lapply(tmp$effect, function(y)  mean(unlist(myfnc2(y,2))))
    #tmp$effect_p <- lapply(tmp$effect, function(y) unlist(lapply(y, function(x) unique(x)[2])))
    #tmp$pvalue <- unlist(lapply(tmp$effect_p, function(x) as.numeric(mean(x)))) # harmonicmeanp::hmp.stat
    #tmp$inference <- tmp$effect#(lapply(tmp$effect, function(y) (unlist(lapply(y, function(x) unique(x))))))
    #tmp$effectsize_mean <- unlist(lapply(tmp$effect, function(y) mean(unlist(lapply(y, function(x) unique(x)[3])))))
    #tmp$effectsize <- (lapply(tmp$effect, function(y) (unlist(lapply(y, function(x) unique(x)[3])))))
    #tmp$lower <- (lapply(tmp$effect, function(y) (unlist(lapply(y, function(x) unique(x)[4])))))
    #tmp$upper <- (lapply(tmp$effect, function(y) (unlist(lapply(y, function(x) unique(x)[5])))))
    tmp <- tmp %>% 
      dplyr::mutate(kept_extrapolated = fraction_extrapolated < extrapolated_threshold) %>%
      group_by(kept_extrapolated) %>%
      dplyr::mutate(fdr = p.adjust(pvalue, method = "BH")) %>%
      ungroup %>%
      dplyr::mutate(fdr = ifelse(kept_extrapolated, fdr, NA)) %>%
      #dplyr::mutate(effectsize = cor_emt - cor_mut) %>%
      dplyr::mutate(method = ifelse(grepl("grfo",files[i], fixed = TRUE),"grfo","grf")) %>%
      dplyr::mutate(resp_type = "ic50")
  }
  
  if(grepl("grf",files[i], fixed = TRUE)&(grepl("auc",files[i], fixed = TRUE))){
    tmp$effect_p <- lapply(tmp$effect, function(y) lapply(y, function(x) unlist(x[[1]])))
    tmp$pvalue <- lapply(tmp$effect_p, function(y) mean(unlist(lapply(y, function(x) unique(x)))))
    tmp$inference <- tmp$effect
    tmp$effectsize_mean <- lapply(tmp$effect, function(y) mean(unlist(myfnc1(y))))
    tmp$effectsize <- lapply(tmp$effect, function(y)  myfnc1(y))
    tmp$lower <- lapply(tmp$effect, function(y)  mean(unlist(myfnc2(y,1))))
    tmp$upper <- lapply(tmp$effect, function(y)  mean(unlist(myfnc2(y,2))))
    #tmp$effect_p <- lapply(tmp$effect, function(y) unlist(lapply(y, function(x) unique(x)[2])))
    #tmp$pvalue <- unlist(lapply(tmp$effect_p, function(x) as.numeric(mean(x)))) # harmonicmeanp::hmp.stat
    #tmp$inference <- tmp$effect#(lapply(tmp$effect, function(y) (unlist(lapply(y, function(x) unique(x))))))
    #tmp$effectsize_mean <- unlist(lapply(tmp$effect, function(y) mean(unlist(lapply(y, function(x) unique(x)[3])))))
    #tmp$effectsize <- (lapply(tmp$effect, function(y) (unlist(lapply(y, function(x) unique(x)[3])))))
    #tmp$lower <- (lapply(tmp$effect, function(y) (unlist(lapply(y, function(x) unique(x)[4])))))
    #tmp$upper <- (lapply(tmp$effect, function(y) (unlist(lapply(y, function(x) unique(x)[5])))))
    tmp <- tmp %>% 
      dplyr::mutate(kept_extrapolated = fraction_extrapolated < extrapolated_threshold) %>%
      #group_by(kept_extrapolated) %>%
      dplyr::mutate(fdr = p.adjust(pvalue, method = "BH")) %>%
      #ungroup %>%
      #dplyr::mutate(fdr = ifelse(kept_extrapolated, fdr, NA)) %>%
      #dplyr::mutate(effectsize = cor_emt - cor_mut) %>%
      dplyr::mutate(method = ifelse(grepl("grfo",files[i], fixed = TRUE),"grfo","grf")) %>%
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
ll_df$index <- 1:nrow(ll_df)


saveRDS(object = ll_df, paste0(path_calcs,"metadata/summaries/PERFORMANCES_v3.rds"))
saveRDS(object = ll_df, paste0(path,"metadata/summaries/PERFORMANCES_v3.rds"))









if(F){ # test combine CI methods....
  test <- ll_df[(ll_df$drug == "1559-GDSC2")&(ll_df$TCGA == "SKCM")&(ll_df$method == "grf")&(ll_df$resp_type == "ic50")&(ll_df$score == "mak"),]
  u <- test$upper[[1]]
  l <- test$lower[[1]]
  e <- test$effectsize[[1]]
  
  #### way too small variance
  Variance <- ((u - l) / (2 * qnorm(0.975)))^2
  InvVariance <- 1 / Variance
  weighted_mean <- sum(e * InvVariance) / sum(InvVariance)
  overall_variance <- 1 / sum(InvVariance)
  
  z_score <- qnorm(1 - (1 - 0.95) / 2)
  margin_of_error <- z_score * sqrt(overall_variance)
  
  l_o <- weighted_mean - margin_of_error
  u_o <- weighted_mean + margin_of_error
  e_o <- weighted_mean
  ###########################
  
  #### percentile method
  overall_lower <- quantile(l, probs = 0.025)
  overall_upper <- quantile(u, probs = 0.975)
  
  l_o <- overall_lower
  u_o <- overall_upper
  e_o <- mean(e)
  ###########################
  
  #### percentile method
  # Combine confidence intervals using weighted average
  weights <- rep(40, 25)  # Vector of sample sizes
  combined_lower <- sum(l * weights) / sum(weights)
  combined_upper <- sum(u * weights) / sum(weights)
  
  # Adjust confidence interval for cross-validation
  num_folds <- length(l/5)
  adjustment_factor <- sqrt(num_folds)
  adjusted_interval <- c((combined_lower + combined_upper) / 2,
                         (combined_upper - combined_lower) / (2 * adjustment_factor))
  
  l_o <- adjusted_interval[1] - adjusted_interval[2]
  u_o <- adjusted_interval[1] + adjusted_interval[2]
  e_o <- mean(e)
  ###########################
  
  
  df <- data.frame(label = c("overall",1:25), mean = c(e_o,e), lower = c(l_o,l), upper = c(u_o,u))
  
  # reverses the factor level ordering for labels after coord_flip()
  df$label <- factor(df$label, levels=rev(df$label))
  
  library(ggplot2)
  fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
    geom_pointrange() + 
    geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
    coord_flip() +  # flip coordinates (puts labels on y axis)
    xlab("Label") + ylab("Mean (95% CI)") +
    theme_bw()  # use a white background
  print(fp)
}

