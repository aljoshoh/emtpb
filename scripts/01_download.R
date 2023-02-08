setwd('emtpb/')
library(readr)
library(dplyr)
library(tidyr)
library(tibble)

source("R/functions.R")
dir.create("data/")
dir.create("metadata/")


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
tmp <- tmp[,c("Cosmic ID","IC50","TCGA Classification","Drug Name","Drug ID","Dataset Version")] %>% 
  mutate(`DRUG ID` = paste0(`Drug ID`,"-",`Dataset Version`)) %>%
  mutate(`COSMIC ID` = `Cosmic ID`) %>%
  mutate(`TCGA Desc` = `TCGA Classification`) %>%
  dplyr::select(c("COSMIC ID","IC50","TCGA Desc","DRUG ID")) %>%
  pivot_wider(names_from = `DRUG ID`, values_from = `IC50`, ) %>%
  mutate(`COSMIC ID` = as.character(`COSMIC ID`))
matrix_resp <- tmp


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

# add EMT score
matrix_exp <- full_join(matrix_mut_pancan[,c("COSMIC ID","TCGA Desc")], matrix_exp)

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


# save
saveRDS(EMT_gs, file = "metadata/EMT_gs.rds")
write_csv(EMTscores, file = "metadata/EMTscores.csv")
write_csv(matrix_exp, file = "metadata/matrix_exp.csv")
write_csv(matrix_mut_pancan, file = "metadata/matrix_mut_pancan.csv")
write_csv(matrix_resp, file = "metadata/matrix_resp.csv")
for(cancertype in cancertypes){
  write_csv(list_mut_ct[[cancertype]], file = paste0("metadata/matrix_mut_",cancertype,".csv"))
}


#
#
#
