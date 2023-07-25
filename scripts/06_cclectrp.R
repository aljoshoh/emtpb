#########################################################################################
## Data Wrangling                                                                      ##
#########################################################################################
library(tidyverse)
library(dplyr)
library(ggplot2)
source("R/functions.R")

dir.create("data/depmap")
# sample_info.csv
# https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/35020903/sample_info.csv?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAIYCQYOYV5JSSROOA/20230306/eu-west-1/s3/aws4_request&X-Amz-Date=20230306T165227Z&X-Amz-Expires=10&X-Amz-SignedHeaders=host&X-Amz-Signature=4e2d2e93bddad7043510a0c439ee7454a95a1ed7715e52326024e24cf130c0d8

# OmicsExpressionProteinCodingGenesTPMLogp1.csv
# https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/38357462/OmicsExpressionProteinCodingGenesTPMLogp1.csv?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAIYCQYOYV5JSSROOA/20230306/eu-west-1/s3/aws4_request&X-Amz-Date=20230306T165923Z&X-Amz-Expires=10&X-Amz-SignedHeaders=host&X-Amz-Signature=

# CTRPv2.0_2015_ctd2_ExpandedDataset.zip
# https://ctd2-data.nci.nih.gov/Public/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/CTRPv2.0_2015_ctd2_ExpandedDataset.zip
unzip("data/depmap/CTRPv2.0_2015_ctd2_ExpandedDataset.zip")

sample_info <- read.csv("data/depmap/sample_info.csv")
expression <- read.csv("data/depmap/OmicsExpressionProteinCodingGenesTPMLogp1.csv")
colnames(expression) <- gsub("\\.\\..*","",colnames(expression))

drug_response <- read.delim("data/depmap/v20.data.curves_post_qc.txt")
meta_experiment <- read.delim("data/depmap/v20.meta.per_experiment.txt")
meta_cell <- read.delim("data/depmap/v20.meta.per_cell_line.txt")
meta_compound <- read.delim("data/depmap/v20.meta.per_compound.txt")
meta_assay <- read.delim("data/depmap/v20.meta.per_assay_plate.txt")
column_description <- read.delim("data/depmap/v20._COLUMNS.txt")

###check if drugs with good EMT models are in CCLE CTRP
### good drugs in GDSC were: Luminespib, Bleomycin, PARP_9482, GSK1059615, Cytarabine, CUDC_101, 
### GSK269962A, Dactolisib, Staurosporine, Topotecan, Teniposide, Nutlin_3a (PAAD), Daporinad, TW_37, 

check_drugs <- c("luminespib", "bleomycin", "PARP_9482", "GSK1059615", "cytarabine", "CUDC_101", 
                 "GSK269962A", "dactolisib", "staurosporine", "topotecan", "teniposide", "nutlin_3", "daporinad", "TW-37")
check_drugs[which(check_drugs %in% meta_compound$cpd_name)]

drugs_in_ccle <- c()
for(i in 1:length(check_drugs)){
  drugs_in_ccle <- c(drugs_in_ccle, meta_compound$cpd_name[grep(check_drugs[i], meta_compound$cpd_name, ignore.case = T)])
}
drugs_in_ccle


#################################################################################
## prepare data for EMT score calculation                                      ##
#################################################################################

#add master cell line id to expression dataframe and drug respose

sample_info <- left_join(meta_cell[,1:2], sample_info, by = c("ccl_name" = "stripped_cell_line_name"))
expression <- inner_join(sample_info[,c(1,3)], expression, by = c("DepMap_ID" = "X"))
rownames(expression) <- expression$master_ccl_id
expression <- expression[,-c(1,2)]
expression <- merge(meta_cell[,c(1,4)], expression, by.x = 1, by.y = 0, all = T)
marker_genes <- c("CDH1", "CDH2", "FN1", "VIM")
#all(marker_genes %in% colnames(expression))

EMTscores <- data.frame(cell_id=character(), tissue=character(), EMT_score=numeric())
EMT_gs <- list()
for(tissue in 1:length(unique(expression$ccle_primary_site))){ # was 12 ?
  curr_tissue <- unique(expression$ccle_primary_site)[tissue]
  print(curr_tissue)
  if(nrow(expression[expression$ccle_primary_site == curr_tissue,]) <= 2){next}
  tmp <- expression[expression$ccle_primary_site == curr_tissue,]
  tmp <- tmp %>% 
    dplyr::rename(`COSMIC ID` = master_ccl_id) %>%
    dplyr::rename(`TCGA Desc` = ccle_primary_site)
  EMT_output <- EMTscore_full(gex = tmp, marker_genes = marker_genes)
  EMTscores <- rbind(EMTscores, EMT_output$EMT_score)
  EMT_gs[[paste(curr_tissue)]] <- EMT_output$EMT_genes
}

#View(expression[expression$ccle_primary_site == curr_tissue,])
saveRDS(EMT_gs, file = "metadata/EMT_gs_ccle.rds")
colnames(EMTscores) = c("COSMIC ID","TCGA Desc","EMT_score")
write_csv(EMTscores, file = "metadata/EMTscores_ccle.csv")


ggplot(EMTscores, aes(x = EMT_score)) + geom_histogram() + facet_wrap(vars(`TCGA Desc`)) +
  theme_classic() + xlab("EMT score")




############################################################################################
## find overlapping cell lines between GDSC and CCLE                                      ##
## correlation of EMT scores?                                                             ##
############################################################################################
# Import resp (drug data) and EMTscores
EMTscores_ccle <- read_csv("metadata/EMTscores_ccle.csv")
#EMTscores_ccle$EMT_score <- -EMTscores_ccle$EMT_score
EMTscores <- read_csv("metadata/EMTscores.csv")
#EMTscores$EMT_score <- -EMTscores$EMT_score
resp = read_csv(paste0(path,"metadata/matrix_resp.csv")) #%>% dplyr::select(-c("COSMIC ID","TCGA Desc"))
GDSC_dr_wide2 <- full_join(resp, EMTscores)
#GDSC_dr_wide2 <- read.csv("~/Dropbox/Dokumente/Master Human Biology/Master Thesis/GDSC/GDSC_drugresponse_mutations.csv")
#GDSC_dr_wide2$COSMIC_ID <- as.character(GDSC_dr_wide2$COSMIC_ID)





tmp <- right_join(sample_info[,c(1,7)], EMTscores_ccle, by = c("master_ccl_id" = "COSMIC_ID"))
tmp$COSMICID <- as.character(tmp$COSMICID)
tmp <- tmp %>%
  mutate(in_GDSC = case_when(COSMICID %in% GDSC_dr_wide2$`COSMIC ID` ~ "yes",
                             T ~ "no"))
nrow(tmp[tmp$in_GDSC == "yes",]) # 644 cell lines in GDSC
nrow(tmp) # 1104 cell lines in total

GDSC_dr_wide2$`COSMIC ID` <- as.character(GDSC_dr_wide2$`COSMIC ID`)
tmp <- tmp %>%
  rename("COSMIC ID" = "COSMICID")
test <- right_join(GDSC_dr_wide2[,c("COSMIC ID", "EMT_score")], tmp, by = "COSMIC ID")
tmp <- right_join(GDSC_dr_wide2[,c("COSMIC ID", "EMT_score")], tmp, by = c("COSMIC ID"))
colnames(tmp) # EMT_score.x (GDSC), EMT_score.y (CCLE)
tmp <- tmp %>%
  rename(EMT_GDSC = EMT_score.x) %>%
  rename(EMT_CCLE = EMT_score.y)

ggplot(tmp, aes(x = EMT_GDSC, y = EMT_CCLE)) + geom_point() +
  stat_smooth(method="lm", color = "gray50", fill = "gray65") + theme_classic() + xlab("EMT score from CCLE data") + 
  ylab("EMT score from GDSC data") + 
  annotate("text", x = 3, y = -3.5, size = 8,
           label = paste("Pearson's R:\n", round(cor(tmp$EMT_GDSC[!is.na(tmp$EMT_GDSC)], 
                                                     na.omit(tmp$EMT_GDSC)), 5))) +
  theme(plot.title = element_text(size = 20), 
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15))





