setwd('emtpb/')
library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)

source("R/functions.R")
dir.create("data/")
dir.create("metadata/")

# download drug metadata release 8.4
emtpb_download(path = "data/screened_compounds_rel_8.4.csv",
               url = 'https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.4/screened_compounds_rel_8.4.csv',
               unzp = F)

# download mutational data gdsc1/gdsc2 pancancer
emtpb_download(path = "data/PANCANCER_Genetic_features_Sat Feb  4 15_33_58 2023.csv",
               url = 'https://www.cancerrxgene.org/downloads/download/genetic_feature?mutation=both&tissue=PANCANCER&screening_set=GDSC1',
               unzp = F)
emtpb_download(path = "data/PANCANCER_Genetic_features_Sat Feb  4 15_42_40 2023.csv",
               url = 'https://www.cancerrxgene.org/downloads/download/genetic_feature?mutation=both&tissue=PANCANCER&screening_set=GDSC2',
               unzp = F)

# prepare matrix for pancancer
tmp1 <- read_csv("data/PANCANCER_Genetic_features_Sat Feb  4 15_33_58 2023.csv")
tmp2 <- read_csv("data/PANCANCER_Genetic_features_Sat Feb  4 15_42_40 2023.csv")
tmp <- full_join(tmp1, tmp2)
tmp <- tmp[,c("COSMIC ID","Genetic Feature","IS Mutated","TCGA Desc")] %>% 
  pivot_wider(names_from = `Genetic Feature`, values_from = `IS Mutated`) %>%
  mutate(`COSMIC ID` = as.character(`COSMIC ID`))
matrix_mut_pancan <- tmp



# download mutational data gdsc1 cancer-specific
tmp1 <- read_csv("data/PANCANCER_Genetic_features_Sat Feb  4 15_42_40 2023.csv")
tmp1 <- tmp1$`TCGA Desc` %>% factor %>% levels %>% setdiff(c("UNCLASSIFIED","OTHER"))
tmp2 <- read_csv("data/PANCANCER_Genetic_features_Sat Feb  4 15_33_58 2023.csv")
tmp2 <- tmp2$`TCGA Desc` %>% factor %>% levels %>% setdiff(c("UNCLASSIFIED","OTHER"))
cancertypes <- unique(c(tmp1, tmp2))

list_mut_ct <- list()
for(cancertype in cancertypes){
  path1 <- paste0("data/PANCANCER_Genetic_features_Sat Feb  4 15_33_58 2023_",cancertype,"_GDSC1.csv")
  path2 <- paste0("data/PANCANCER_Genetic_features_Sat Feb  4 15_42_40 2023_",cancertype,"_GDSC2.csv")
  emtpb_download(path = path1,
                 url = paste0('https://www.cancerrxgene.org/downloads/download/genetic_feature?mutation=both&tissue=',cancertype,'&screening_set=GDSC1'),
                 unzp = F
  )
  emtpb_download(path = path2,
                 url = paste0('https://www.cancerrxgene.org/downloads/download/genetic_feature?mutation=both&tissue=',cancertype,'&screening_set=GDSC2'),
                 unzp = F
                 )
  tmp1 <- read_csv(path1)
  tmp2 <- read_csv(path2)
  tmp <- full_join(tmp1, tmp2)
  tmp <- tmp[,c("COSMIC ID","Genetic Feature","IS Mutated","TCGA Desc")] %>% pivot_wider(names_from = `Genetic Feature`, values_from = `IS Mutated`)
  list_mut_ct[[cancertype]] <- tmp
}


# download gdsc1 and gdsc2
emtpb_download(path = "data/PANCANCER_IC_Wed Feb 8 15_46_36 2023_GDSC1.csv",
               url = 'https://www.cancerrxgene.org/downloads/download/ic?tissue=PANCANCER&screening_set=GDSC1',
               unzp = F)
emtpb_download(path = "data/PANCANCER_IC_Wed Feb 8 15_46_31 2023_GDSC2.csv",
               url = 'https://www.cancerrxgene.org/downloads/download/ic?tissue=PANCANCER&screening_set=GDSC2',
               unzp = F)
tmp1 <- read_csv("data/PANCANCER_IC_Wed Feb 8 15_46_36 2023_GDSC1.csv")
tmp2 <- read_csv("data/PANCANCER_IC_Wed Feb 8 15_46_31 2023_GDSC2.csv")
tmp <- full_join(tmp1, tmp2)
tmp_save <- tmp
tmp <- tmp[,c("Cosmic ID","IC50","TCGA Classification","Drug Name","Drug ID","Dataset Version")] %>% 
  mutate(`DRUG ID` = paste0(`Drug ID`,"-",`Dataset Version`)) %>%
  mutate(`COSMIC ID` = `Cosmic ID`) %>%
  mutate(`TCGA Desc` = `TCGA Classification`) %>%
  dplyr::select(c("COSMIC ID","IC50","TCGA Desc","DRUG ID")) %>%
  pivot_wider(names_from = `DRUG ID`, values_from = `IC50`, ) %>%
  mutate(`COSMIC ID` = as.character(`COSMIC ID`))
matrix_resp <- tmp

tmp <- tmp_save[,c("Cosmic ID","AUC","TCGA Classification","Drug Name","Drug ID","Dataset Version")] %>% 
  mutate(`DRUG ID` = paste0(`Drug ID`,"-",`Dataset Version`)) %>%
  mutate(`COSMIC ID` = `Cosmic ID`) %>%
  mutate(`TCGA Desc` = `TCGA Classification`) %>%
  dplyr::select(c("COSMIC ID","AUC","TCGA Desc","DRUG ID")) %>%
  pivot_wider(names_from = `DRUG ID`, values_from = `AUC`, ) %>%
  mutate(`COSMIC ID` = as.character(`COSMIC ID`))
matrix_resp_auc <- tmp


# download processed gene expression data from the GDSC
emtpb_download(path = "data/Cell_line_RMA_proc_basalExp.txt.zip",
               url = 'https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/preprocessed/Cell_line_RMA_proc_basalExp.txt.zip',
               unzp = T
)
tmp <- read_delim("data/Cell_line_RMA_proc_basalExp.txt", delim = "\t")
matrix_exp <- tmp %>% 
  dplyr::select(-c("GENE_title")) %>% 
  dplyr::filter(!is.na(GENE_SYMBOLS)) %>% 
  as.data.frame %>% 
  column_to_rownames("GENE_SYMBOLS") %>% 
  t %>% 
  as.data.frame %>% 
  rownames_to_column("COSMIC ID") %>% 
  as_tibble() %>% 
  mutate(`COSMIC ID` = gsub("DATA.", "", `COSMIC ID`))
matrix_exp <- full_join(matrix_mut_pancan[,c("COSMIC ID","TCGA Desc")], matrix_exp)


# add EMT score (1: marisa) ####################
EMTscores <- data.frame(`COSMIC ID`=character(), `TCGA Desc`=character(), EMT_score=numeric())
EMT_gs <- list()
tcga_used <- c()
cancertypes <- cancertypes[cancertypes %in% names(which(table(matrix_exp$`TCGA Desc`) > 5))] # cancer types with more than 5 samples
for(tcga in 1:length(cancertypes)){
  print(tcga)
  if(cancertypes[tcga] != "UNCLASSIFIED" & !is.na(cancertypes[tcga])){
    curr_tcga <- cancertypes[tcga]} else{next}
  tcga_used[length(tcga_used) + 1] <- curr_tcga
  gex_tmp <- as.data.frame(matrix_exp[matrix_exp$`TCGA Desc` == curr_tcga &!is.na(matrix_exp$`TCGA Desc`),])
  EMT_output <- EMTscore_full(gex = gex_tmp)
  EMTscores <- rbind(EMTscores, EMT_output$EMT_score)
  EMT_gs[[tcga]] <- EMT_output$EMT_genes
}
names(EMT_gs) <- cancertypes
saveRDS(EMT_gs, file = "metadata/EMT_gs.rds")
colnames(EMTscores) = c("COSMIC ID","TCGA Desc","EMT_score")
write_csv(EMTscores, file = "metadata/EMTscores.csv")
##############################################


# add EMT score (2: GSVA) ####################
library(msigdbr)
library(GSVA)
EMTscores <- data.frame(`COSMIC ID`=character(), `TCGA Desc`=character(), EMT_score=numeric())
EMT_gs <- list()
tcga_used <- c()
cancertypes <- cancertypes[cancertypes %in% names(which(table(matrix_exp$`TCGA Desc`) > 5))] # cancer types with more than 5 samples

for(tcga in 1:length(cancertypes)){
  print(tcga)
  if(cancertypes[tcga] != "UNCLASSIFIED" & !is.na(cancertypes[tcga])){
    curr_tcga <- cancertypes[tcga]} else{next}
  tcga_used[length(tcga_used) + 1] <- curr_tcga
  gex_tmp <- as.data.frame(matrix_exp[matrix_exp$`TCGA Desc` == curr_tcga &!is.na(matrix_exp$`TCGA Desc`),])
  
  Hallmark_genesets <- msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::select(gs_name, gene_symbol)
  hallmarks_list <- unlist(split(Hallmark_genesets, f = Hallmark_genesets$gs_name), recursive = F)
  gsva_hallmarks <- as.data.frame(gsva(t(as.matrix(gex_tmp %>% dplyr::select(-c("TCGA Desc")) %>% column_to_rownames("COSMIC ID"))%>%na.omit), hallmarks_list, verbose=FALSE))
  EMT_output <- as.data.frame(t(gsva_hallmarks))[,"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.gene_symbol",drop = F]
  EMT_output <- EMT_output %>%
    mutate(`TCGA Desc` = curr_tcga) %>%
    mutate(EMT_score = HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.gene_symbol) %>%
    rownames_to_column("COSMIC ID") %>% 
    dplyr::select(-c(HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.gene_symbol))

  EMTscores <- rbind(EMTscores, EMT_output)
  EMT_gs[[tcga]] <- hallmarks_list$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.gene_symbol
}
names(EMT_gs) <- cancertypes
saveRDS(EMT_gs, file = "metadata/EMT_gs_gsva.rds")
write_csv(EMTscores, file = "metadata/EMTscores_gsva.csv")
##############################################


# add EMT score (3: tan et al 2014) #################### -> using matrix_resp for tcga matching
if(F){
  anno <- read_csv("data/Model.csv") %>% mutate(CellLine = ModelID)
  scores <- read_excel("data/Table_S4C_EMT_curated.xlsx")
  scores$`COSMIC ID` <- as.numeric(unlist(lapply(scores$`CosmicID Label`, function(x) strsplit(x," ")[[1]][1])))
  scores <- scores %>% 
    mutate(EMT_score = as.numeric(gsub(" ","", `Generic.EMT .\nCL.Ksscore`))) %>%
    mutate(`COSMIC ID` = as.character(`COSMIC ID`))
  scores <- left_join(scores, (matrix_resp%>%mutate(gdsc_resp = "yes"))[,c("COSMIC ID","TCGA Desc","gdsc_resp")])
  
  # PROBLEM only 692/789 cell lines overlap (88%).... 
  if(T){
    EMTscores <- scores[,c("COSMIC ID","TCGA Desc","EMT_score")]
    print(nrow(na.omit(EMTscores))/nrow(EMTscores))
  }
  #EMT_gs <- list()
  #cancertypes <- cancertypes[cancertypes %in% names(which(table(matrix_exp$`TCGA Desc`) > 5))] # cancer types with more than 5 samples
  #names(EMT_gs) <- cancertypes
  #saveRDS(EMT_gs, file = "metadata/EMT_gs_gsva.rds")
  write_csv(EMTscores, file = "metadata/EMTscores_tan_088.csv")
}
##############################################


# add EMT score (4: secrier) ####################  -> using matrix_resp for tcga matching
if(F){
  emtpb_download(path = "data/secrier_supplementary_tables.xlsx",
                 url = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-36439-7/MediaObjects/41467_2023_36439_MOESM8_ESM.xlsx',
                 unzp = F)
  #emtpb_download(path = "data/Model.csv",
  #               url = 'https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/38466923/Model.csv?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAIYCQYOYV5JSSROOA/20230227/eu-west-1/s3/aws4_request&X-Amz-Date=20230227T153029Z&X-Amz-Expires=10&X-Amz-SignedHeaders=host&X-Amz-Signature=fbccd46f6f002d8246596cc74cdc35c70e8337b0baf4d4c3d27209f2f209765a',
  #               unzp = F)
  
  # get ccle data
  anno <- read_csv("data/Model.csv") %>% mutate(CellLine = ModelID)
  scores <- read_excel("data/secrier_supplementary_tables.xlsx", sheet = 5)
  scores <- scores[,8:ncol(scores)] %>% na.omit
  
  scores$cell_line_name_ccle <- unlist(lapply(scores$CCLE_ID, function(x) strsplit(x,"_")[[1]][2]))
  
  scores <- left_join(scores, anno, by = c("cell_line_name_ccle"="StrippedCellLineName")) %>%
    mutate(`COSMIC ID` = as.character(COSMICID)) %>%
    mutate(EMT_score = emt_score)
  scores <- left_join(scores, (matrix_resp%>%mutate(gdsc_resp = "yes"))[,c("COSMIC ID","TCGA Desc","gdsc_resp")])
  scores_ccle <- scores[,c("COSMIC ID","TCGA Desc","EMT_score")]
  
  # get gdsc data -> using their inquired supplementary data
  load("data/EMT_ic50__secrier.RData") # inquired supplementary file # contains emt_ic50
  emt_ic50 <- emt_ic50[,c("gdsc_ID","emt_score")]
  tmp1 <- read_csv("data/PANCANCER_Genetic_features_Sat Feb  4 15_33_58 2023.csv")
  tmp2 <- read_csv("data/PANCANCER_Genetic_features_Sat Feb  4 15_42_40 2023.csv")
  tmp <- full_join(tmp1, tmp2)
  anno_gdsc <- distinct(tmp[,c("COSMIC ID","Cell Line Name","TCGA Desc")]) %>% rename(gdsc_ID = `Cell Line Name`)
  emt_ic50 <- left_join(emt_ic50, anno_gdsc) %>%
    rename(EMT_score = emt_score)
  scores_gdsc <- emt_ic50[,c("COSMIC ID","TCGA Desc","EMT_score")]
  
  
  
  # PROBLEM only 477/734 (65%) cell lines overlap.... table(!is.na(scores_ccle$`COSMIC ID`))
  if(T){
    EMTscores <- rbind(scores_ccle, scores_gdsc)
    EMTscores <- EMTscores[!(duplicated(EMTscores$`COSMIC ID`) & !is.na(EMTscores$`COSMIC ID`)),]
    print(nrow(na.omit(EMTscores))/nrow(EMTscores))
  }
  #EMT_gs <- list()
  #cancertypes <- cancertypes[cancertypes %in% names(which(table(matrix_exp$`TCGA Desc`) > 5))] # cancer types with more than 5 samples
  #names(EMT_gs) <- cancertypes
  #saveRDS(EMT_gs, file = "metadata/EMT_gs_gsva.rds")
  write_csv(EMTscores, file = "metadata/EMTscores_secrier_065.csv")
}
##############################################


# save
write_csv(matrix_exp, file = "metadata/matrix_exp.csv")
write_csv(matrix_mut_pancan, file = "metadata/matrix_mut_PANCAN.csv")
write_csv(matrix_resp, file = "metadata/matrix_resp.csv")
write_csv(matrix_resp_auc, file = "metadata/matrix_resp_auc.csv")
for(cancertype in cancertypes){
  write_csv(list_mut_ct[[cancertype]], file = paste0("metadata/matrix_mut_",cancertype,".csv"))
}


#
#
#



















# DEPRECATED
# add EMT score (4: secrier) #################### -> using matrix_resp for tcga matching
if(F){
  emtpb_download(path = "data/secrier_supplementary_tables.xlsx",
                 url = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-36439-7/MediaObjects/41467_2023_36439_MOESM8_ESM.xlsx',
                 unzp = F)
  #emtpb_download(path = "data/Model.csv",
  #               url = 'https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/38466923/Model.csv?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAIYCQYOYV5JSSROOA/20230227/eu-west-1/s3/aws4_request&X-Amz-Date=20230227T153029Z&X-Amz-Expires=10&X-Amz-SignedHeaders=host&X-Amz-Signature=fbccd46f6f002d8246596cc74cdc35c70e8337b0baf4d4c3d27209f2f209765a',
  #               unzp = F)
  anno <- read_csv("data/Model.csv") %>% mutate(CellLine = ModelID)
  scores <- read_excel("data/secrier_supplementary_tables.xlsx", sheet = 5)
  scores <- scores[,8:ncol(scores)] %>% na.omit
  scores$cell_line_name_ccle <- unlist(lapply(scores$CCLE_ID, function(x) strsplit(x,"_")[[1]][2]))
  
  scores <- left_join(scores, anno, by = c("cell_line_name_ccle"="StrippedCellLineName")) %>%
    mutate(`COSMIC ID` = (COSMICID)) %>%
    mutate(EMT_score = emt_score)
  scores <- left_join(scores, (matrix_resp%>%mutate(gdsc_resp = "yes"))[,c("COSMIC ID","TCGA Desc","gdsc_resp")])
  
  # PROBLEM only 474/728 (65%) cell lines overlap....
  if(T){
    EMTscores <- scores[,c("COSMIC ID","TCGA Desc","EMT_score")]
  }
  #EMT_gs <- list()
  #cancertypes <- cancertypes[cancertypes %in% names(which(table(matrix_exp$`TCGA Desc`) > 5))] # cancer types with more than 5 samples
  #names(EMT_gs) <- cancertypes
  #saveRDS(EMT_gs, file = "metadata/EMT_gs_gsva.rds")
  write_csv(EMTscores, file = "metadata/EMTscores_secrier_065.csv")
}
##############################################






