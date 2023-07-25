library(ggplot2)
library(ggrepel)
# local paths
path <- "/Users/alexander.ohnmacht/research/marisa/emtpb/"
path_calcs <- "/Volumes/pheb/lustre/groups/cbm01/code/alexander.ohnmacht/emtpb/"

# analysis
df <- readRDS(paste0(path,"metadata/summaries/PERFORMANCES_v1.rds")) # path_calcs
drug_meta <- read_csv("data/screened_compounds_rel_8.4.csv") %>% 
  mutate(drugname = DRUG_NAME, target = TARGET, pathway = TARGET_PATHWAY) %>%
  select(c(drugname, pathway)) %>% # target
  distinct()
drug_meta <- drug_meta[!drug_meta$drugname %>% duplicated,]
df <- left_join(df, drug_meta, by = "drugname")



# 1. density plot of two methods
ggplot(data = df)+
  geom_density(aes(x = cor_emt, color = method))+
  theme_minimal()

# 2. volcano of all stuff
df_2 <- df %>% 
  dplyr::filter(resp_type == "ic50") %>%
  #dplyr::filter(score == "mak") %>%
  dplyr::filter(TCGA != "PANCAN") %>%
  dplyr::select(c(method, drugname, drug,TCGA,score,pvalue,fdr,resp_type,cor_emt)) %>%
  tidyr::pivot_wider(names_from = method, values_from = c(pvalue,fdr,cor_emt))
ggplot(data = df_2)+
  geom_point(aes(x = -log10(pvalue_eln), y = -log10(pvalue_grf)), color = "black", alpha = 0.2)+
  geom_point(data = df_2[(df_2$fdr_eln < 0.2) | (df_2$fdr_grf < 0.2),], aes(x = -log10(pvalue_eln), y = -log10(pvalue_grf), color = TCGA), alpha = 0.8)+
  geom_text_repel(data = df_2[(df_2$fdr_eln < 0.2) | (df_2$fdr_grf < 0.2),],
                  aes(x = -log10(pvalue_eln), y = -log10(pvalue_grf), label = drugname), color = "black", size = 4, show_guide = F, min.segment.length = 0.02)+
  theme_minimal()+
  xlab("-log10(p-value) standard")+
  ylab("-log10(p-value) causal")+
  geom_abline(intercept = 0, linetype = "dashed")+
  theme(aspect.ratio=1)

ggplot(data = df_2)+
  geom_point(aes(x = cor_emt_eln, y = cor_emt_grf), color = "black", alpha = 0.2)+
  geom_point(data = df_2[(df_2$fdr_eln < 0.2),], aes(x = cor_emt_eln, y = cor_emt_grf, color = TCGA), alpha = 0.8)+
  geom_text_repel(data = df_2[(df_2$fdr_eln < 0.2),],
                  aes(x = cor_emt_eln, y = cor_emt_grf, label = drugname), color = "black", size = 4, show_guide = F, min.segment.length = 0.02)+
  theme_minimal()+
  xlab("Pearson's r standard")+
  ylab("Pearson's r causal")+
  geom_abline(intercept = 0, linetype = "dashed")+
  theme(aspect.ratio=1)


# 3. barplot for one drug and cancer type
df_3 <- df %>%
  dplyr::filter(resp_type == "auc") %>%
  dplyr::filter(!is.na(fdr)) %>%
  #dplyr::filter(score == "mak") %>%
  dplyr::filter(TCGA == "SKCM") %>%
  dplyr::select(c(method, drugname, drug,TCGA,score,pvalue,fdr,resp_type,cor_emt,effectsize,delta)) %>%
  tidyr::pivot_wider(names_from = method, values_from = c(pvalue,fdr,cor_emt,effectsize)) %>%
  #dplyr::filter(drug == "1559-GDSC2") %>%
  mutate(id = paste0(drugname,"-",score,"-",resp_type))
df_3$effectsize_grf[df_3$score == "mak"] <- -df_3$effectsize_grf[df_3$score == "mak"]
fdr_threshold <- 0.2
x <- "delta"
h0 <- sort(unlist(na.omit(df_3[,"pvalue_eln"])))[max(which(sort(unlist(na.omit(df_3[,"fdr_eln"])))<fdr_threshold))]
p <- ggplot()+
  geom_point(data = df_3, aes(x = !!sym(x), y = -log10(pvalue_eln), color = score), alpha = 0.2)+
  geom_point(data = df_3[(df_3$delta > 0.3),], aes(x = !!sym(x), y = -log10(pvalue_eln), color = score), alpha = 0.8)+
  geom_text_repel(data = df_3[(df_3$delta > 0.3),],
                  aes(x = !!sym(x), y = -log10(pvalue_eln), label = id), color = "black", size = 4, show_guide = F, min.segment.length = 0.02)+
  theme_minimal()+
  xlab("treatment effect estimate")+
  ylab("unadjusted -log10(p)")+
  geom_hline(yintercept=-log10(h0), linetype='dotted', color = "grey66", size = 0.9)+
  annotate("text", x = min(df_3[,x]%>%na.omit), y = -log10(h0), label = paste0("FDR<",as.character(round(fdr_threshold*100)),"%"), hjust = 0, vjust = 0); p





# EMT correlations
EMTscores_mak <- suppressMessages(read_csv(paste0(path,"metadata/EMTscores","",".csv"))) #%>% mutate(EMT_score = -EMT_score) #%>% mutate(method = "mak")
EMTscores_gsva <- suppressMessages(read_csv(paste0(path,"metadata/EMTscores","_gsva",".csv"))) #%>% mutate(method = "gsva")
EMTscores_tan <- suppressMessages(read_csv(paste0(path,"metadata/EMTscores","_tan_088",".csv"))) #%>% mutate(method = "tan")
EMTscores_secrier <- suppressMessages(read_csv(paste0(path,"metadata/EMTscores","_secrier_065",".csv"))) #%>% mutate(method = "secrier")
emt <- full_join(full_join(full_join(EMTscores_mak, EMTscores_gsva, by = "COSMIC ID"), EMTscores_tan, by = "COSMIC ID"), EMTscores_secrier, by = "COSMIC ID")
emtp <- select_if(emt, is.numeric)
colnames(emtp) <- c("cosmic","Mak","GSVA","Tan","Secrier")
cor(emtp[emt$`TCGA Desc.x` == "SKCM",-1], use = "complete.obs")





en_df <- list();for(i in 1:nrow(hm_fig1_en)){
  en <- enrichment(tab = table(drug_targets[drug_targets %in% names(which(table(drug_targets)>=1))], 
                               hm_fig1_en[i,drug_targets %in% names(which(table(drug_targets)>=1))]%>%unlist), 
                   col_label = "biomarker", row_label = "target",
                   digits = 3) #### digits !
  en <- en[en$biomarker != "none",]
  en$fdr <- p.adjust(abs(en$value), method = "BH")
  en_df[[i]] <- if(nrow(en) == 0){NULL}else{en$type <- row.names(hm_fig1_en)[i];en}
}; en_df <- do.call(rbind, en_df); en_df$FDR <- p.adjust(abs(en_df$value), method = "BH")



