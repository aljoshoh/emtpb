library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)
library(RColorBrewer)
library(plotly)
library(GGally)

# local paths
path <- "/Users/alexander.ohnmacht/research/marisa/emtpb/"
path_calcs <- "/Volumes/pheb/lustre/groups/cbm01/code/alexander.ohnmacht/emtpb/"

# analysis
df <- readRDS(paste0(path,"metadata/summaries/PERFORMANCES_v1.rds")) # path_calcs
drug_meta <- read_csv(paste0(path,"data/screened_compounds_rel_8.4.csv")) %>% 
  mutate(drugname = DRUG_NAME, target = TARGET, pathway = TARGET_PATHWAY) %>%
  select(c(drugname, pathway)) %>% # target
  distinct()
drug_meta <- drug_meta[!drug_meta$drugname %>% duplicated,]
df <- left_join(df, drug_meta, by = "drugname")


# estimates for luminespib 
test <- resp_help[resp_help$`TCGA Desc` == "SKCM","1559-GDSC2",drop= F]
exp(c(1,0.5,1.5)*sd(test$`1559-GDSC2`  %>% na.omit)+mean(test$`1559-GDSC2` %>% na.omit))

# figure S1 (performances (fdr) dependent on all 'drugname'&'TCGA'&'resp_type','emtscore')
# only showing cancer types and drugs with at least one FDR<0.2
# gaps are not enough samples and extrapolated IC50 (Methods)
df_4 <- df %>%
  #dplyr::filter(resp_type == "ic50") %>%
  #dplyr::filter(!is.na(fdr)) %>%
  #dplyr::filter(score == "mak") %>%
  dplyr::filter(method == "eln") %>%
  dplyr::select(c(fdr, pvalue, delta, score,drug, drugname, TCGA, resp_type)) %>%
  #dplyr::filter(drugname %in% unique(df$drugname)[1:100]) %>%
  mutate(id = paste0(score," (",toupper(resp_type),")")) %>%
  mutate(drugname2 = paste0(drugname," (",sub(".*-","",drug),")"))
df_4 <- df_4 %>%
  dplyr::filter(TCGA %in% unique((df_4 %>% group_by(TCGA) %>% dplyr::filter(fdr < 0.2))$TCGA)) %>% # filter for TCGA with at least one < 0.2
  dplyr::filter(drug %in% unique((df_4 %>% group_by(drug) %>% dplyr::filter(fdr < 0.2))$drug)) %>% # filter for drugname with at least one < 0.2 
  mutate(txt = ifelse(fdr < 0.2,"-","")) %>%
  mutate(txt = ifelse(is.na(txt), "", txt))
h <- ggplot(df_4, aes(x = id, y = drugname2, fill = -log10(pvalue))) +
  geom_tile() +
  facet_grid(. ~ TCGA, scales = "free", space = "free")+
  scale_fill_continuous(na.value = 'grey88')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_text(aes(label=txt), color = "white"); h



# Fig. 2 volcano for all cancer types and drugs (BRCA, SKCM, GBM, LUAD, SCLC, COREAD)
fdr_threshold <- 0.2
df_3 <- df %>%
  #dplyr::filter(resp_type == "ic50") %>%
  dplyr::filter(!is.na(fdr)) %>%
#dplyr::filter(score == "mak") %>%
  dplyr::filter(TCGA != "PANCAN") %>%
  dplyr::filter(TCGA %in% as.character(unlist(unique(na.omit(df[(df$method == "eln")&(df$fdr < fdr_threshold),"TCGA"]))))) 
df_3 <- df_3 %>%
  dplyr::select(c(method, drugname, drug,TCGA,score,pvalue,fdr,resp_type,cor_emt,effectsize,delta)) %>%
  tidyr::pivot_wider(names_from = method, values_from = c(pvalue,fdr,cor_emt,effectsize)) %>%
  #dplyr::filter(drug == "1559-GDSC2") %>%
  mutate(id = paste0(drugname," (",TCGA,",",resp_type,", ",score,")")) %>%
  mutate(strat = paste0(drug,"-",score,"-",resp_type,"-",TCGA)) %>% 
  mutate(filt = paste0(drugname,"-",TCGA))
df_3 <- df_3 %>% group_by(strat) %>% summarize(effectsize_grf = -na.omit(effectsize_grf), 
                                                  TCGA = TCGA, 
                                                  score = score, 
                                                  drugname = drugname,
                                                  resp_type = resp_type,
                                                  pvalue_eln = na.omit(pvalue_eln),
                                                  pvalue_grf = na.omit(pvalue_grf),
                                                  delta = na.omit(delta),
                                                  fdr_eln = na.omit(fdr_eln),
                                                  filt = filt,
                                                  id = id) %>% 
  distinct %>% 
  ungroup()
#df_3$effectsize_grf[df_3$score == "mak"] <- -df_3$effectsize_grf[df_3$score == "mak"] ### FIX THIS !
x <- "delta"
y <- "pvalue_eln"
h0 <- max((df_3[,c(y,"fdr_eln")] %>% arrange(fdr_eln) %>% filter(fdr_eln < fdr_threshold))[,y])
label_data <- df_3[(df_3$fdr_eln < fdr_threshold)&(df_3$filt %in% names(which(table(df_3[df_3$fdr_eln < fdr_threshold,"filt"] ) >=3))),]
#h0 <- sort(unlist(na.omit(df_3[,"pvalue_eln"])))[max(which(sort(unlist(na.omit(df_3[,"fdr_eln"])))<=fdr_threshold))]
p <- ggplot()+
  geom_point(data = df_3, aes(x = !!sym(x), y = -log10(!!sym(y)), color = TCGA), alpha = 0.2)+
  geom_point(data = df_3[(df_3$fdr_eln < fdr_threshold),], aes(x = !!sym(x), y = -log10(!!sym(y)), color = TCGA), alpha = 1)+
  geom_text_repel(data = label_data,
                  aes(x = !!sym(x), y = -log10(!!sym(y)), label = id), color = "black", show_guide = F, min.segment.length = 0.02,
                  hjust = 1, direction = "y", nudge_x = 0.1- as.numeric(unlist(label_data[,x])), segment.size = 0.1, size = 3
                  )+
  theme_minimal()+
  xlab("delta Pearson correlation")+
  ylab("unadjusted -log10(p)")+
  scale_fill_brewer(name = "cancer type",palette="Pastel2", direction=-1) + scale_color_brewer(name = "cancer type",palette="Pastel2", direction=-1)+
  geom_hline(yintercept=-log10(h0), linetype='dotted', color = "grey66", size = 0.9)+
  annotate("text", x = min(df_3[,x]%>%na.omit), y = -log10(h0), label = paste0("FDR<",as.character(round(fdr_threshold*100)),"%"), hjust = 0, vjust = 0); p



# supplement volcano for cancer types and drugs (PANCAN)
fdr_threshold <- 0.05
df_3 <- df %>%
  #dplyr::filter(resp_type == "ic50") %>%
  dplyr::filter(!is.na(fdr)) %>%
  #dplyr::filter(score == "mak") %>%
  dplyr::filter(TCGA %in% "PANCAN") 
df_3 <- df_3 %>%
  dplyr::select(c(method, drugname, drug,TCGA,score,pvalue,fdr,resp_type,cor_emt,effectsize,delta)) %>%
  tidyr::pivot_wider(names_from = method, values_from = c(pvalue,fdr,cor_emt,effectsize)) %>%
  #dplyr::filter(drug == "1559-GDSC2") %>%
  mutate(id = paste0(drugname," (",TCGA,",",resp_type,", ",score,")")) %>%
  mutate(strat = paste0(drug,"-",score,"-",resp_type,"-",TCGA)) %>% 
  mutate(filt = paste0(drugname,"-",TCGA))
df_3 <- df_3 %>% group_by(strat) %>% summarize(effectsize_grf = -na.omit(effectsize_grf), 
                                               TCGA = TCGA, 
                                               score = score, 
                                               drugname = drugname,
                                               resp_type = resp_type,
                                               pvalue_eln = na.omit(pvalue_eln),
                                               pvalue_grf = na.omit(pvalue_grf),
                                               delta = na.omit(delta),
                                               fdr_eln = na.omit(fdr_eln),
                                               filt = filt,
                                               id = id) %>% 
  distinct %>% 
  ungroup()
#df_3$effectsize_grf[df_3$score == "mak"] <- -df_3$effectsize_grf[df_3$score == "mak"] ### FIX THIS !
x <- "delta"
y <- "pvalue_eln"
h0 <- max((df_3[,c(y,"fdr_eln")] %>% arrange(fdr_eln) %>% filter(fdr_eln < fdr_threshold))[,y])
label_data <- df_3[(df_3$fdr_eln < fdr_threshold)&(df_3$filt %in% names(which(table(df_3[df_3$fdr_eln < fdr_threshold,"filt"] ) >=3))),]
#h0 <- sort(unlist(na.omit(df_3[,"pvalue_eln"])))[max(which(sort(unlist(na.omit(df_3[,"fdr_eln"])))<=fdr_threshold))]
p <- ggplot()+
  geom_point(data = df_3, aes(x = !!sym(x), y = -log10(!!sym(y))), alpha = 0.2, color = "black")+
  geom_point(data = df_3[(df_3$fdr_eln < fdr_threshold),], aes(x = !!sym(x), y = -log10(!!sym(y))), color = "black", alpha = 0.8)+
  geom_text_repel(data = label_data,
                  aes(x = !!sym(x), y = -log10(!!sym(y)), label = id), color = "black", show_guide = F, min.segment.length = 0.02,
                  hjust = 1, direction = "y", nudge_x = 0.- as.numeric(unlist(label_data[,x])), segment.size = 0.1, size = 3
                  )+
  theme_minimal()+
  xlab("delta Pearson correlation")+
  ylab("unadjusted -log10(p)")+
  geom_hline(yintercept=-log10(h0), linetype='dotted', color = "grey66", size = 0.9)+
  annotate("text", x = min(df_3[,x]%>%na.omit), y = -log10(h0), label = paste0("FDR<",as.character(round(fdr_threshold*100)),"%"), hjust = 0.5, vjust = 0); p




# Figure 1 volcano for cancer types and drugs (SKCM)
pastel_colors <- brewer.pal(6, "Pastel2")
df_3 <- df %>%
  #dplyr::filter(resp_type == "ic50") %>%
  dplyr::filter(!is.na(fdr)) %>%
  #dplyr::filter(score == "mak") %>%
  dplyr::filter(TCGA == "SKCM") 
df_3 <- df_3 %>%
  dplyr::select(c(method, drugname, drug,TCGA,score,pvalue,fdr,resp_type,cor_emt,effectsize,delta)) %>%
  tidyr::pivot_wider(names_from = method, values_from = c(pvalue,fdr,cor_emt,effectsize)) %>%
  #dplyr::filter(drug == "1559-GDSC2") %>%
  mutate(id = paste0(drugname," (",TCGA,",",resp_type,", ",score,")")) %>%
  mutate(strat = paste0(drug,"-",score,"-",resp_type,"-",TCGA)) %>% 
  mutate(filt = paste0(drugname,"-",TCGA))
df_3 <- df_3 %>% group_by(c(strat)) %>% summarize(effectsize_grf = -na.omit(effectsize_grf), 
                                               TCGA = TCGA, 
                                               score = score, 
                                               drugname = drugname,
                                               resp_type = resp_type,
                                               pvalue_eln = na.omit(pvalue_eln),
                                               pvalue_grf = na.omit(pvalue_grf),
                                               delta = na.omit(delta),
                                               fdr_eln = na.omit(fdr_eln),
                                               filt = filt,
                                               id = id) %>% 
  distinct %>% 
  ungroup()
#df_3$effectsize_grf[df_3$score == "mak"] <- -df_3$effectsize_grf[df_3$score == "mak"] ### FIX THIS !
fdr_threshold <- 0.2
x <- "delta"
y <- "pvalue_eln"
h0 <- max((df_3[,c(y,"fdr_eln")] %>% arrange(fdr_eln) %>% filter(fdr_eln < fdr_threshold))[,y])
label_data <- df_3[(df_3$fdr_eln < fdr_threshold)&(df_3$filt %in% names(which(table(df_3[df_3$fdr_eln < 0.2,"filt"] ) >=3))),]
#h0 <- sort(unlist(na.omit(df_3[,"pvalue_eln"])))[max(which(sort(unlist(na.omit(df_3[,"fdr_eln"])))<=fdr_threshold))]
p <- ggplot()+
  geom_point(data = df_3, aes(x = !!sym(x), y = -log10(!!sym(y))), alpha = 0.2, color = pastel_colors[1])+
  geom_point(data = df_3[(df_3$fdr_eln < fdr_threshold),], aes(x = !!sym(x), y = -log10(!!sym(y))), alpha = 0.8, color = pastel_colors[1])+
  geom_text_repel(data = label_data,
                  aes(x = !!sym(x), y = -log10(!!sym(y)), label = id), color = "black", show_guide = F, min.segment.length = 0.02,
                  hjust = 1, direction = "y", nudge_x = 0.1- as.numeric(unlist(label_data[,x])), segment.size = 0.1, size = 3
  )+
  theme_minimal()+
  xlab("delta Pearson correlation")+
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
emtp$TCGA <- factor(emt$`TCGA Desc.x`)
emtp <- emtp[emtp$TCGA %in% c("BRCA","SKCM","LUAD","SCLC","GBM","COREAD"),]
emtp$TCGA <- factor(emtp$TCGA)
colnames(emtp) <- c("cosmic","mak","gsva","tan","secrier","TCGA")

# Create a grid of scatter plots
my_dens <- function(data, mapping, ...) {
  ggplot(data = data, mapping=mapping) +
    geom_density(..., alpha = 1) 
}
ggp <- ggpairs(emtp, columns = 2:ncol(emtp), aes(colour=TCGA),
        diag = list(continuous = my_dens)
        ) + scale_fill_manual(values = rev(pastel_colors))+scale_color_manual(values = rev(pastel_colors))+
  theme_minimal()
for(i in 1:ggp$nrow) {
  for(j in 1:ggp$ncol){
    if((j == 4)|((j==5)&(i == 4)))
      ggp[i,j] <- ggp[i,j] + 
      scale_fill_manual(values=rev(pastel_colors)[-3]) +
      scale_color_manual(values=rev(pastel_colors)[-3])  
  }
}












if(F){
# DEPRECATED


# 2. volcano of all stuff
df_2 <- df %>% 
  dplyr::filter(resp_type == "ic50") %>%
  #dplyr::filter(score == "secrier") %>%
  #dplyr::filter(TCGA == "HNSC") %>%
  dplyr::select(c(method, drugname, drug,TCGA,score,pvalue,fdr,resp_type,cor_emt)) %>%
  tidyr::pivot_wider(names_from = method, values_from = c(pvalue,fdr,cor_emt))
ggplot(data = df_2)+
  geom_point(aes(x = -log10(pvalue_eln), y = -log10(pvalue_grf)), color = "black", alpha = 0.2)+
  geom_point(data = df_2[(df_2$fdr_eln < 0.2),], aes(x = -log10(pvalue_eln), y = -log10(pvalue_grf), color = TCGA), alpha = 0.8)+
  geom_text_repel(data = df_2[(df_2$fdr_eln < 0.2),],
                  aes(x = -log10(pvalue_eln), y = -log10(pvalue_grf), label = drugname), color = "black", size = 4, show_guide = F, min.segment.length = 0.02)+
  theme_minimal()+
  xlab("-log10(p-value) standard")+
  ylab("-log10(p-value) causal")+
  geom_abline(intercept = 0, linetype = "dashed")+
  theme(aspect.ratio=1)


# 1. density plot of two methods
ggplot(data = df)+
  geom_density(aes(x = cor_emt, color = method))+
  theme_minimal()


# failing enrichments
en_df <- list();for(i in 1:nrow(hm_fig1_en)){
  en <- enrichment(tab = table(drug_targets[drug_targets %in% names(which(table(drug_targets)>=1))], 
                               hm_fig1_en[i,drug_targets %in% names(which(table(drug_targets)>=1))]%>%unlist), 
                   col_label = "biomarker", row_label = "target",
                   digits = 3) #### digits !
  en <- en[en$biomarker != "none",]
  en$fdr <- p.adjust(abs(en$value), method = "BH")
  en_df[[i]] <- if(nrow(en) == 0){NULL}else{en$type <- row.names(hm_fig1_en)[i];en}
}; en_df <- do.call(rbind, en_df); en_df$FDR <- p.adjust(abs(en_df$value), method = "BH")

}

