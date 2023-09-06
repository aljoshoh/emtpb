library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)
library(RColorBrewer)
library(plotly)
library(GGally)
library(ggbeeswarm)
library(ggpubr)
library(tidyr)
library(rstatix)
library(tibble)
library(data.table)
library(limma)
library(enrichR)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(shades)
library(viridis)

EMTscores <- suppressMessages(read_csv(paste0(path,"metadata/EMTscores","",".csv"))) # BENCHMARKING ORDER FOR GDSC CANCER TYPES
cancertypes <- unique(c('PANCAN', unlist(EMTscores['TCGA Desc'])))

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


# EMT scores
EMTscores_mak <- suppressMessages(read_csv(paste0(path,"metadata/EMTscores","",".csv"))) #%>% mutate(EMT_score = -EMT_score) #%>% mutate(method = "mak")
EMTscores_gsva <- suppressMessages(read_csv(paste0(path,"metadata/EMTscores","_gsva",".csv"))) #%>% mutate(method = "gsva")
EMTscores_tan <- suppressMessages(read_csv(paste0(path,"metadata/EMTscores","_tan_088",".csv"))) #%>% mutate(method = "tan")
EMTscores_secrier <- suppressMessages(read_csv(paste0(path,"metadata/EMTscores","_secrier_065",".csv"))) #%>% mutate(method = "secrier")
emt <- full_join(full_join(full_join(EMTscores_mak, EMTscores_gsva, by = "COSMIC ID"), EMTscores_tan, by = "COSMIC ID"), EMTscores_secrier, by = "COSMIC ID")
emtp <- select_if(emt, is.numeric)
emtp$TCGA <- factor(emt$`TCGA Desc.x`)
emtp <- emtp[emtp$TCGA %in% c("BRCA","SKCM","LUAD","SCLC","GBM","COREAD"),]
emtp$TCGA <- factor(emtp$TCGA)
colnames(emtp) <- c("COSMIC ID","mak","gsva","tan","secrier","TCGA Desc")

# responses
resp_ic50 = suppressMessages(read_csv(paste0(path,"metadata/matrix_resp","",".csv"))) # %>% dplyr::select(-c("COSMIC ID","TCGA Desc")))
resp_auc = suppressMessages(read_csv(paste0(path,"metadata/matrix_resp","_auc",".csv"))) # %>% dplyr::select(-c("COSMIC ID","TCGA Desc")))

# get experiments
experiment <- read_csv(paste0(path,"metadata/paper/benchmark_paper_exp3+6_03_pool.csv"), na = character())
#experiment_full <- read_csv(paste0(path,"metadata/paper/benchmark_paper_exp3.csv"), na = character())
#experiment_second <- read_csv(paste0(path,"metadata/paper/benchmark_paper_exp6.csv"), na = character())
#experiment_full <- bind_rows(experiment_full, experiment_second)

# get expression
gex <- read_csv(paste0(path,"metadata/matrix_exp.csv"))
mut_skcm <-read_csv(paste0(path,"metadata/matrix_mut_SKCM.csv"))
### init






# prepare validation plots
df_3 <- df %>%
  dplyr::filter(!is.na(fdr)) %>%
  dplyr::filter(score == "mak") %>%
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


df_4 <- df_3[df_3$drugname == "Luminespib",]
df_5 <- df_4
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




#
# auc
#


# plot confidence band for luminespib AUC
df_6 <- df_5[df_5$`c(strat)` == "1559-GDSC2-mak-auc-SKCM",]
error.measurement <- abs(df_6$effectsize_grf-df_6$upper)
ribbondata<- data.frame(x=c(0, 1),
                        ymin=c(  0 - error.measurement,
                                 1*df_6$effectsize_grf - error.measurement),
                        ymax=c( 0  + error.measurement,
                                1*df_6$effectsize_grf + error.measurement))

# import validation
V <- read_excel(path = paste0(path,"metadata/goeksu/exp1/AUCs.xlsx"))
V$tgf <- factor(unlist(lapply(V$CELL_DRUG, function(x) substr(x,1,1))))
V$emt <- c("mes.","mes.","mes.","mes.","epi.","epi.","epi.","epi.")
V$name <- c("RPMI-7951","RPMI-7951","A375","A375","SK-MEL-5","SK-MEL-5","IGR-37","IGR-37")
V <- V[!V$name %in% "A375",]
V <- V %>% 
  mutate(label = factor(paste(emt,tgf)))
  
# summary
ggplot(V, aes(x = name, y = AUC, fill = label, color = label, shape = tgf)) + 
  geom_point(position = position_dodge(width = 0.4), size = 2.5) +
  scale_color_manual(values = c("grey30","darkorange","grey30","navyblue","red")) +
  scale_fill_manual(values = alpha(c("grey30","darkorange","grey30","navyblue","red"), 0.3)) +
  scale_shape_manual(values = c(21, 24)) +
  theme_minimal() + 
  ggtitle("Response to luminespib") +
  guides(fill = guide_legend(override.aes = list(size = 3)))+
  guides(shape = guide_legend(override.aes = list(size = 3)))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

Vred <- V %>%
  group_by(name) %>%
  mutate(delta = AUC - lag(AUC, order_by = tgf)) %>%
  na.omit

# paper fig
Vredci <- full_join(Vred %>% 
                  as.data.frame %>%
                  mutate(upper = NA) %>%
                  mutate(lower = NA),
                df_6 %>% 
                  dplyr::select(c("effectsize_grf","lower","upper")) %>% 
                  mutate(label = "n.a.") %>%
                  mutate(name = "n.a.") %>%
                  mutate(emt = "n.a.") %>%
                  mutate(tgf = "n.a.") %>%
                  mutate(upper = ribbondata$ymax[ribbondata$x == 0]) %>%
                  mutate(lower = ribbondata$ymin[ribbondata$x == 0]) %>%
                  mutate(delta = 0)
                ) %>%
  mutate(name = factor(name,levels = c("RPMI-7951","IGR-37","SK-MEL-5","n.a.")))

ggplot(data = Vredci)+
  geom_errorbar(aes(x = name, y = -delta, ymin = lower, ymax=upper), color="black", width = 0.4) +
  geom_point(aes(x = name, y = -delta, color = emt, shape = tgf), size = 2.5) +
  theme_minimal()+
  ylim(c(-0.05,0.05))+
  ylab("delta AUC")+
  xlab("Cell Line")+
  scale_color_manual(values = c("mes."="navyblue","epi."="darkorange","n.a."="black"))+
  scale_shape_manual(values = c(21, 24)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))







#
# ic50
#

# plot confidence band for luminespib AUC
df_6 <- df_5[df_5$`c(strat)` == "1559-GDSC2-mak-ic50-SKCM",]
error.measurement <- abs(df_6$effectsize_grf-df_6$upper)
ribbondata<- data.frame(x=c(0, 1),
                        ymin=c(  0 - error.measurement,
                                 1*df_6$effectsize_grf - error.measurement),
                        ymax=c( 0  + error.measurement,
                                1*df_6$effectsize_grf + error.measurement))

# import validation
V <- read_excel(path = paste0(path,"metadata/goeksu/exp1/IC50s.xlsx"))
V$tgf <- factor(unlist(lapply(V$CELL_DRUG, function(x) substr(x,1,1))))
V$emt <- c("epi.","epi.","epi.","mes.","epi.","mes.","mes.","mes.")
V$name <- c("IGR-37","SK-MEL-5","SK-MEL-5","A375","IGR-37","RPMI-7951","RPMI-7951","A375")
V <- V %>% 
  mutate(label = factor(paste(emt,tgf))) %>%
  mutate(IC50 = `IC50(nM)`)
V <- V[!V$name %in% "A375",]

# summary
ggplot(V, aes(x = name, y = IC50, fill = label, color = label, group = label, shape = tgf)) + 
  geom_point(position = position_dodge(width = 0.4), size = 2.5) +
  scale_color_manual(values = c("grey30","darkorange","grey30","navyblue","red")) +
  scale_fill_manual(values = alpha(c("grey30","darkorange","grey30","navyblue","red"), 0.3)) +
  scale_shape_manual(values = c(21, 24)) +
  theme_minimal() + 
  ggtitle("Response to luminespib") +
  guides(fill = guide_legend(override.aes = list(size = 3)))+
  guides(shape = guide_legend(override.aes = list(size = 3)))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

Vred <- V %>%
  group_by(name) %>%
  mutate(delta = log(IC50 / lag(IC50, order_by = tgf))) %>% # calc log(c_t/c_u) for comparison
  na.omit

# paper fig
Vredci <- full_join(Vred %>% 
                      as.data.frame %>%
                      mutate(upper = NA) %>%
                      mutate(lower = NA),
                    df_6 %>% 
                      dplyr::select(c("effectsize_grf","lower","upper")) %>% 
                      mutate(label = "n.a.") %>%
                      mutate(name = "n.a.") %>%
                      mutate(emt = "n.a.") %>%
                      mutate(tgf = "n.a.") %>%
                      mutate(upper = ribbondata$ymax[ribbondata$x == 0]) %>%
                      mutate(lower = ribbondata$ymin[ribbondata$x == 0]) %>%
                      mutate(delta = 0)
) %>%
  mutate(name = factor(name,levels = c("RPMI-7951","IGR-37","SK-MEL-5","n.a.")))

ggplot(data = Vredci)+
  geom_errorbar(aes(x = name, y = -delta, ymin = lower, ymax=upper), color="black", width = 0.4) +
  geom_point(aes(x = name, y = -delta, color = emt, shape = tgf), size = 2.5) +
  theme_minimal()+
  #ylim(c(-0.05,0.05))+
  ylab("delta log(IC50)")+
  xlab("Cell Line")+
  scale_color_manual(values = c("mes."="navyblue","epi."="darkorange","n.a."="black"))+
  scale_shape_manual(values = c(21, 24)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
































