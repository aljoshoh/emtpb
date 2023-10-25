setwd('emtpb/')
library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)
library(RColorBrewer)
library(plotly)
library(GGally)
library(yaml)
config <- yaml.load_file("scripts/config.yaml")

# local paths
path <- config$r_local_path
path_calcs <- config$r_local_path_calcs

# analysis
df <- readRDS(paste0(path,"metadata/summaries/PERFORMANCES_v3.rds")) # path_calcs
df$pvalue <- unlist(lapply(df$pvalue, function(x){x}))
df$effectsize_mean <- unlist(lapply(df$effectsize_mean, function(x){x}))
df$lower <- unlist(lapply(df$lower, function(x){x}))
df$upper <- unlist(lapply(df$upper, function(x){x}))
df <- df %>%
  dplyr::select(-c(effectsize, inference))


drug_meta <- read_csv(paste0(path,"data/screened_compounds_rel_8.4.csv")) %>% 
  mutate(drugname = DRUG_NAME, target = TARGET, pathway = TARGET_PATHWAY) %>%
  select(c(drugname, pathway)) %>% # target
  distinct()
drug_meta <- drug_meta[!drug_meta$drugname %>% duplicated,]
df <- left_join(df, drug_meta, by = "drugname")


# estimates for luminespib 
resp_help = suppressMessages(read_csv(paste0(path,"metadata/matrix_resp","",".csv"))) # %>% dplyr::select(-c("COSMIC ID","TCGA Desc")))
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
  dplyr::filter(method %in% c("grfo","eln")) %>%
  dplyr::filter(TCGA %in% as.character(unlist(unique(na.omit(df[(df$method == "eln")&(df$fdr < fdr_threshold),"TCGA"]))))) 
df_3 <- df_3 %>%
  dplyr::select(c(method, drugname, drug,TCGA,score,pvalue,fdr,resp_type,cor_emt,effectsize_mean,delta, lower, upper)) %>%
  tidyr::pivot_wider(names_from = method, values_from = c(pvalue,fdr,cor_emt,effectsize_mean, lower, upper)) %>%
  #dplyr::filter(drug == "1559-GDSC2") %>%
  mutate(id = paste0(drugname," (",TCGA,", ",resp_type,", ",score,")")) %>%
  mutate(strat = paste0(drug,"-",score,"-",resp_type,"-",TCGA)) %>% 
  mutate(filt = paste0(drugname,"-",TCGA))
df_3 <- df_3 %>% group_by(strat) %>% summarize(effectsize_grfo = na.omit(effectsize_mean_grfo), 
                                                  TCGA = TCGA, 
                                                  score = score, 
                                                  drugname = drugname,
                                                  resp_type = resp_type,
                                                  pvalue_eln = na.omit(pvalue_eln),
                                                  pvalue_grf = na.omit(pvalue_grfo),
                                                  delta = na.omit(delta),
                                                  fdr_eln = na.omit(fdr_eln),
                                                  filt = filt,
                                                  id = id,
                                                  lower = na.omit(lower_grfo),
                                                  upper = na.omit(upper_grfo)) %>% 
  distinct %>% 
  ungroup()
#df_3$effectsize_grf[df_3$score == "mak"] <- -df_3$effectsize_grf[df_3$score == "mak"] ### FIX THIS !
x <- "delta"
y <- "pvalue_eln"
howmanyhits <- 3
h0 <- max((df_3[,c(y,"fdr_eln")] %>% arrange(fdr_eln) %>% filter(fdr_eln < fdr_threshold))[,y])
label_data <- df_3[(df_3$fdr_eln < fdr_threshold)&(df_3$filt %in% names(which(table(df_3[df_3$fdr_eln < fdr_threshold,"filt"] ) >=howmanyhits))),]
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
  dplyr::filter(method %in% c("grfo","eln")) %>%
  dplyr::filter(TCGA %in% "PANCAN") 
df_3 <- df_3 %>%
  dplyr::select(c(method, drugname, drug,TCGA,score,pvalue,fdr,resp_type,cor_emt,effectsize_mean,delta, lower, upper)) %>%
  tidyr::pivot_wider(names_from = method, values_from = c(pvalue,fdr,cor_emt,effectsize_mean, lower, upper)) %>%
  #dplyr::filter(drug == "1559-GDSC2") %>%
  mutate(id = paste0(drugname," (",TCGA,",",resp_type,", ",score,")")) %>%
  mutate(strat = paste0(drug,"-",score,"-",resp_type,"-",TCGA)) %>% 
  mutate(filt = paste0(drugname,"-",TCGA))
df_3 <- df_3 %>% group_by(strat) %>% summarize(effectsize_grf = na.omit(effectsize_mean_grfo), 
                                               TCGA = TCGA, 
                                               score = score, 
                                               drugname = drugname,
                                               resp_type = resp_type,
                                               pvalue_eln = na.omit(pvalue_eln),
                                               pvalue_grf = na.omit(pvalue_grfo),
                                               delta = na.omit(delta),
                                               fdr_eln = na.omit(fdr_eln),
                                               filt = filt,
                                               id = id,
                                               lower = na.omit(lower_grfo),
                                               upper = na.omit(upper_grfo)) %>% 
  distinct %>% 
  ungroup()
#df_3$effectsize_grf[df_3$score == "mak"] <- -df_3$effectsize_grf[df_3$score == "mak"] ### FIX THIS !
x <- "delta"
y <- "pvalue_eln"
howmanyhits <- 3
h0 <- max((df_3[,c(y,"fdr_eln")] %>% arrange(fdr_eln) %>% filter(fdr_eln < fdr_threshold))[,y])
label_data <- df_3[(df_3$fdr_eln < fdr_threshold)&(df_3$filt %in% names(which(table(df_3[df_3$fdr_eln < fdr_threshold,"filt"] ) >=howmanyhits))),]
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
  dplyr::filter(method %in% c("grfo","eln")) %>%
  dplyr::filter(TCGA == "SKCM") 
df_3 <- df_3 %>%
  dplyr::select(c(method, drugname, drug,TCGA,score,pvalue,fdr,resp_type,cor_emt,effectsize_mean,delta, lower, upper)) %>%
  tidyr::pivot_wider(names_from = method, values_from = c(pvalue,fdr,cor_emt,effectsize_mean, lower, upper)) %>%
  #dplyr::filter(drug == "1559-GDSC2") %>%
  mutate(id = paste0(drugname," (",TCGA,",",resp_type,", ",score,")")) %>%
  mutate(strat = paste0(drug,"-",score,"-",resp_type,"-",TCGA)) %>% 
  mutate(filt = paste0(drugname,"-",TCGA))
df_3 <- df_3 %>% group_by(c(strat)) %>% summarize(effectsize_grf = na.omit(effectsize_mean_grfo), 
                                               TCGA = TCGA, 
                                               score = score, 
                                               drugname = drugname,
                                               resp_type = resp_type,
                                               pvalue_eln = na.omit(pvalue_eln),
                                               pvalue_grf = na.omit(pvalue_grfo),
                                               delta = na.omit(delta),
                                               fdr_eln = na.omit(fdr_eln),
                                               filt = filt,
                                               id = id,
                                               lower = na.omit(lower_grfo),
                                               upper = na.omit(upper_grfo)) %>% 
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
                  hjust = 1, direction = "y", nudge_x = 0.18- as.numeric(unlist(label_data[,x])), segment.size = 0.1, size = 3
  )+
  theme_minimal()+
  xlab("delta Pearson correlation")+
  ylab("unadjusted -log10(p)")+
  geom_hline(yintercept=-log10(h0), linetype='dotted', color = "grey66", size = 0.9)+
  annotate("text", x = min(df_3[,x]%>%na.omit), y = -log10(h0), label = paste0("FDR<",as.character(round(fdr_threshold*100)),"%"), hjust = 0, vjust = 0); p

x <- "id"
y <- "effectsize_grf"
df_4 <- (df_3[df_3$resp_type == "auc",])
df_4 <- df_4[(df_4$`c(strat)` %in% label_data$`c(strat)`),]

p <- ggplot(data = df_4)+
  geom_errorbar(aes(x = reorder(!!sym(x), !!sym(y)), ymin = lower, ymax=upper), color="black", alpha = 1) +
  geom_line(aes(x = reorder(!!sym(x), !!sym(y)), y = !!sym(y), group = 1), color = "grey")+ 
  geom_point(aes(x = reorder(!!sym(x), !!sym(y)), y = !!sym(y)), color = "black", alpha = 1) +
  #geom_ribbon(aes(x = 1:length(!!sym(x)), ymin = lower, ymax=upper), fill = "lightgrey", alpha = 0.4)+
  #geom_point(data = df_3, aes(x = !!sym(x), y = -log10(!!sym(y))), alpha = 0.2, color = pastel_colors[1])+
  #geom_point(data = df_3[(df_3$fdr_eln < fdr_threshold),], aes(x = !!sym(x), y = -log10(!!sym(y))), alpha = 0.8, color = pastel_colors[1])+
  #geom_text_repel(data = df_4, aes(x = reorder(!!sym(x), !!sym(y)), y = !!sym(y), label = id), color = "black", size = 3)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), plot.margin = margin(t = 0, r = 0, b = 0, l = 50, unit = "pt"))+
  xlab("")+
  ylab("Standardised estimated causal effect");p

df_5 <- df_3
df_5 <- df_5[(df_5$`c(strat)` %in% label_data$`c(strat)`),]
resp_help = suppressMessages(read_csv(paste0(path,"metadata/matrix_resp","",".csv"))) # %>% dplyr::select(-c("COSMIC ID","TCGA Desc")))
resp_help_auc = suppressMessages(read_csv(paste0(path,"metadata/matrix_resp","_auc",".csv"))) # %>% dplyr::select(-c("COSMIC ID","TCGA Desc")))
for(i in 1:nrow(df_5)){
  d <- paste0(strsplit(df_5$`c(strat)`[i],"-")[[1]][1:2], collapse = "-")
  if(df_5$resp_type[i] == "ic50"){
    test <- resp_help[resp_help$`TCGA Desc` == df_5$TCGA[i],,drop= F]
    #tt <- exp(c(df_5$effectsize_grf[i],df_5$lower[i],df_5$upper[i])*sd(unlist(test[,d])  %>% na.omit)+mean(unlist(test[,d]) %>% na.omit)) -> transform back EMT
    tt <- exp(abs(c(df_5$effectsize_grf[i],df_5$lower[i],df_5$upper[i])))
  }
  if(df_5$resp_type[i] == "auc"){
    test <- resp_help_auc[resp_help_auc$`TCGA Desc` == df_5$TCGA[i],,drop= F]
    #tt <- c(df_5$effectsize_grf[i],df_5$lower[i],df_5$upper[i])*sd(unlist(test[,d]) %>% na.omit)+mean(unlist(test[,d]) %>% na.omit) -> transform back EMT
    tt <- (abs(c(df_5$effectsize_grf[i],df_5$lower[i],df_5$upper[i])))
  }
  df_5$effectsize_grf_e[i] <- tt[1]
  df_5$lower_e[i] <- tt[2]
  df_5$upper_e[i] <- tt[3]
}

# plot confidence band for luminespib AUC
df_6 <- df_5[df_5$id == "Luminespib (SKCM,auc, mak)",]
error.measurement <- abs(df_6$effectsize_grf-df_6$upper)
ribbondata<- data.frame(x=c(0, 2),
                        ymin=c(  0 - error.measurement,
                                 2*df_6$effectsize_grf - error.measurement),
                        ymax=c( 0  + error.measurement,
                                2*df_6$effectsize_grf + error.measurement)
)
ggplot()+
    geom_point(data = df_6, aes(0, effectsize_grf), shape = 4, size = 3) +  # Scatter plot points
  geom_abline(slope = df_6$effectsize_grf, intercept = df_6$effectsize_grf-df_6$lower,linetype = "dashed") +
  geom_abline(slope = df_6$effectsize_grf, intercept = df_6$effectsize_grf-df_6$upper,linetype = "dashed") +
  geom_ribbon(data=ribbondata, aes(x=x,ymin=ymin,ymax=ymax), inherit.aes = FALSE, alpha = 0.3) +
  geom_abline(intercept = 0, slope = df_6$effectsize_grf, color = "black") +  # Fixed slope line
  theme_minimal()+
  xlim(c(-0.5,2))+
  ylim(c(-0.25,0.05))+
  ylab("delta AUC")+
  xlab("sigma EMT change")+
  scale_x_continuous(breaks=c(0,1,2))
sp
  
  
  
  

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

