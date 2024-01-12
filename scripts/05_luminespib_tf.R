setwd('emtpb/')
library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggbeeswarm)
library(tibble)
library(data.table)
library(yaml)
config <- yaml.load_file("scripts/config.yaml")
path <- config$r_local_path

# specifics
library(enrichR)
library(limma)
#rna <- read_csv(paste0(path,"metadata/matrix_exp.csv")) %>% dplyr::filter(`COSMIC ID` %in% (EMTscores$`COSMIC ID`[EMTscores$`TCGA Desc` == "SKCM"]))
#EMT_gs <- readRDS(paste0(path,"metadata/EMT_gs.rds"))
rna <- read_csv(paste0(path,"metadata/matrix_exp.csv"))

# thresholds
fdr_threshold <- 0.2
howmanyhits <- 3

# read results
df <- readRDS(paste0(path,"metadata/summaries/PERFORMANCES_v3.rds")) # path_calcs
df$pvalue <- unlist(lapply(df$pvalue, function(x){x}))
df$effectsize_mean <- unlist(lapply(df$effectsize_mean, function(x){x}))
df$lower <- unlist(lapply(df$lower, function(x){x}))
df$upper <- unlist(lapply(df$upper, function(x){x}))
df <- df %>%
  dplyr::select(-c(effectsize, inference))
df_3 <- df %>%
  #dplyr::filter(resp_type == "ic50") %>%
  dplyr::filter(!is.na(fdr)) %>%
  #dplyr::filter(score == "mak") %>%
  dplyr::filter(TCGA != "PANCAN") %>%
  dplyr::filter(method %in% c("grfo","eln")) %>%
  dplyr::filter(TCGA %in% as.character(unlist(unique(na.omit(df[(df$method == "eln")&(df$fdr < fdr_threshold),"TCGA"]))))) 
df_3 <- df_3 %>%
  dplyr::select(c(method, drugname, drug,TCGA,score,pvalue,fdr,resp_type,cor_emt,effectsize_mean,delta, lower, upper)) %>%
  tidyr::pivot_wider(names_from = method, values_from = c(pvalue,fdr,cor_emt,effectsize_mean, lower, upper)) %>%
  #dplyr::filter(drug == "1559-GDSC2") %>%
  mutate(id = paste0(drugname," (",TCGA,", ",resp_type,", ",score,")")) %>%
  mutate(strat = paste0(drug,"-",score,"-",resp_type,"-",TCGA)) %>% 
  mutate(filt = paste0(drugname,"-",TCGA))
label_data <- df_3[(df_3$fdr_eln < fdr_threshold)&(df_3$filt %in% names(which(table(df_3[df_3$fdr_eln < fdr_threshold,"filt"] ) >=howmanyhits))),] %>%
  filter(!is.na(fdr_eln)) %>% 
  dplyr::select(c(drugname,drug,TCGA,resp_type)) %>% distinct
  

# loop over results
lenriched <- list()
lenriched_full <- list()
for(j in 1:nrow(label_data)){
  print(round(j/nrow(label_data)*100))
  ct <- label_data$TCGA[j]
  dr <- label_data$drug[j]
  resp_type <- label_data$resp_type[j]
  response_type <- ifelse(label_data$resp_type[j] == "ic50","","_auc")
  score <- ifelse(label_data$score[j] == "MAK","",paste0("_",label_data$score[j]))
  resp <- suppressMessages(read_csv(paste0(path,"metadata/matrix_resp",response_type,".csv"))) #%>% dplyr::select(-c("COSMIC ID","TCGA Desc"))
  EMTscores <- suppressMessages(read_csv(paste0(path,"metadata/EMTscores",score,".csv"))) #%>% dplyr::select(-c("COSMIC ID","TCGA Desc"))
  gex_gdsc_t <- rna %>% dplyr::filter(`COSMIC ID` %in% (EMTscores$`COSMIC ID`[EMTscores$`TCGA Desc` == ct]))
  skcm_emt <- EMTscores[EMTscores$`TCGA Desc` == ct,]
  skcm_emt <- left_join(skcm_emt, resp[,c("COSMIC ID",dr),drop = F])
  tmp <- left_join(skcm_emt, gex_gdsc_t, by = c("COSMIC ID" = "COSMIC ID"))
  tmp <- na.omit(tmp)
  exp_skcm <- tmp[,-c(1:5)]
  luminespib_resp_skcm <- unlist(tmp[dr])
  
  my_design <- model.matrix(~luminespib_resp_skcm)
  my_fit <- lmFit(exp_skcm%>%as.matrix%>%t, my_design)
  my_ebayes <- eBayes(my_fit)
  my_results <- suppressMessages(topTable(my_ebayes, n=nrow(my_ebayes))%>%rownames_to_column("gene"))
  
  limma_luminespib_pos <- my_results$gene[(my_results$adj.P.Val < 0.1)&(my_results$logFC > 0)]
  limma_luminespib_neg <- my_results$gene[(my_results$adj.P.Val < 0.1)&(my_results$logFC < 0)]
  
  if((ct == "SKCM") & (dr == "1559-GDSC2") & (resp_type == "auc")){
    limma_luminespib_auc <- limma_luminespib_pos
  }
  if((ct == "SKCM") & (dr == "1559-GDSC2") & (resp_type == "ic50")){
    limma_luminespib_ic50 <- limma_luminespib_pos
  }
  
  # enrichment
  #terms <- c("ARCHS4_Tissues","ChEA_2022","RNAseq_Automatic_GEO_Signatures_Human_Up","RNAseq_Automatic_GEO_Signatures_Human_Down","LINCS_L1000_CRISPR_KO_Consensus_Sigs")
  #terms <- "ChEA_2022"
  terms <- c("ChEA_2022")
  enriched_pos <- enrichr(limma_luminespib_pos, terms); enriched_pos <- do.call(rbind, enriched_pos); enriched_pos$sign <- "positive"
  enriched_neg <- enrichr(limma_luminespib_neg, terms); enriched_neg <- do.call(rbind, enriched_neg); enriched_neg$sign <- "negative"
  enriched <- rbind(enriched_pos, enriched_neg)
  enriched$cancertype <- ct
  enriched$drug <- dr
  enriched$drugname <-  label_data$drugname[j]
  enriched$resp_type <-  label_data$resp_type[j]
  
  lenriched[[j]] <- enriched
  
  # enrichment
  #terms <- c("ARCHS4_Tissues","ChEA_2022","RNAseq_Automatic_GEO_Signatures_Human_Up","RNAseq_Automatic_GEO_Signatures_Human_Down","LINCS_L1000_CRISPR_KO_Consensus_Sigs")
  #terms <- "ChEA_2022"
  terms <- c("GO_Biological_Process_2023")
  enriched_pos <- enrichr(limma_luminespib_pos, terms); enriched_pos <- do.call(rbind, enriched_pos); enriched_pos$sign <- "positive"
  enriched_neg <- enrichr(limma_luminespib_neg, terms); enriched_neg <- do.call(rbind, enriched_neg); enriched_neg$sign <- "negative"
  enriched <- rbind(enriched_pos, enriched_neg)
  enriched$cancertype <- ct
  enriched$drug <- dr
  enriched$drugname <-  label_data$drugname[j]
  enriched$resp_type <-  label_data$resp_type[j]
  
  lenriched_full[[j]] <- enriched 
}

# save go results
enr_full <- rbindlist(lenriched_full)
saveRDS(enr_full, file = paste0(path,"metadata/paper/tf_target_enrichments_v3.rds"))

# save tf results
enr <- rbindlist(lenriched)
saveRDS(enr, file = paste0(path,"metadata/paper/tf_target_enrichments_v2.rds"))
saveRDS(limma_luminespib_auc, file = paste0(path,"metadata/paper/tf_diff_genes_luminespib_auc.rds"))
saveRDS(limma_luminespib_ic50, file = paste0(path,"metadata/paper/tf_diff_genes_luminespib_ic50.rds"))





# DEPRECATED
if(F){
  # genes MITF targets
  #strsplit(enr_luminespib$Genes[enr_luminespib$Term == "MITF 21258399 ChIP-Seq MELANOMA Human"][1],";")
  
  colors_genes <- c(limma_luminespib_neg, limma_luminespib_pos)
  names(colors_genes) = colors_genes
  colors_genes[rep(T, length(colors_genes))] <- "black"
  colors_genes[names(colors_genes) %in% strsplit(enriched$Genes[enriched$Term == "MITF 21258399 ChIP-Seq MELANOMA Human" & (enriched$sign == "positive")],";")[[1]]] <- "deepskyblue"
  colors_genes[names(colors_genes) %in% strsplit(enriched$Genes[enriched$Term == "FIBROBLAST" & (enriched$sign == "negative")],";")[[1]]] <- "darkorange"
  colors_genes[names(colors_genes) %in% c("CDH1","CDH2")] <- "deepskyblue4"
  saveRDS(colors_genes, file = "metadata/paper/colors_genes_enrichr_luminespib_melanoma.rds")
  
  
  
  
  
  # luad rock inhibitor
  gex_gdsc_t <- read_csv(paste0(path,"metadata/matrix_exp.csv")) %>% dplyr::filter(`COSMIC ID` %in% (EMTscores$`COSMIC ID`[EMTscores$`TCGA Desc` == "LUAD"]))
  luad_emt <- EMTscores[EMTscores$`TCGA Desc` == "LUAD",]
  luad_emt <- left_join(luad_emt, resp[,c("COSMIC ID","1192-GDSC1"),drop = F])
  tmp <- left_join(luad_emt, gex_gdsc_t, by = c("COSMIC ID" = "COSMIC ID"))
  tmp <- na.omit(tmp)
  exp_luad <- tmp[,-c(1:5)]
  rock_resp_luad <- tmp$`1192-GDSC1`
  
  my_design <- model.matrix(~rock_resp_luad)
  library(limma)
  my_fit <- lmFit(exp_luad%>%as.matrix%>%t, my_design)
  my_ebayes <- eBayes(my_fit)
  my_results <- topTable(my_ebayes, n=nrow(my_ebayes))%>%rownames_to_column("gene")
  
  
  limma_rock_pos <- my_results$gene[(my_results$adj.P.Val < 0.1)&(my_results$logFC > 0)]
  limma_rock_neg <- my_results$gene[(my_results$adj.P.Val < 0.1)&(my_results$logFC < 0)]
  
  library(enrichR)
  terms <- c("ARCHS4_Tissues","ChEA_2022","RNAseq_Automatic_GEO_Signatures_Human_Up","RNAseq_Automatic_GEO_Signatures_Human_Down","LINCS_L1000_CRISPR_KO_Consensus_Sigs")
  enriched_pos <- enrichr(limma_rock_pos, terms); enriched_pos <- do.call(rbind, enriched_pos); enriched_pos$sign <- "positive"
  enriched_neg <- enrichr(limma_rock_neg, terms); enriched_neg <- do.call(rbind, enriched_neg); enriched_neg$sign <- "negative"
  enriched <- rbind(enriched_pos, enriched_neg)
  
  
  
  
  
  # exploration of gene lists in LUAD (rock inhibitor)
  if(F){
    dbs <- listEnrichrDbs()
    terms <- c("RNAseq_Automatic_GEO_Signatures_Human_Up","RNAseq_Automatic_GEO_Signatures_Human_Down","LINCS_L1000_CRISPR_KO_Consensus_Sigs")
    enriched_pos <- enrichr(limma_rock_pos, terms); enriched_pos <- do.call(rbind, enriched_pos); enriched_pos$sign <- "positive"
    enriched_neg <- enrichr(limma_rock_neg, terms); enriched_neg <- do.call(rbind, enriched_neg); enriched_neg$sign <- "negative"
    enriched <- rbind(enriched_pos, enriched_neg)
  }
  #
}