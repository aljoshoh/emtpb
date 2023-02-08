#download data #################
emtpb_download <- function(url,
                           path,
                           unzp = F){
  if(!file.exists(path)){
    download.file(url = url, 
                  destfile = path)
  }
  if(unzp){
    unzip(zipfile = path, exdir = "data")
  }
}
##################################


#EMT scoring ##################################
EMTscore_full <- function(gex = data.frame, marker_genes = c("CDH1", "CDH2", "VIM", "FN1")){
  #gex is a data frame with genes as columns and samples (cell lines) as rows, 
  #the first two columns are "COSMIC ID" and "TCGA Desc", genes are columns 3:ncol(gex)
  markers <- gex[,c("COSMIC ID", "TCGA Desc", marker_genes)] #marker genes
  correlates <- data.frame(tested_gene=character(), CDH1_r=numeric(), CDH1_pval=numeric(), CDH2_r=numeric(), CDH2_pval=numeric(), 
                           VIM_r=numeric(), VIM_pval=numeric(), FN1_r=numeric(), FN1_pval=numeric())
  gex <- gex[,!colnames(gex) %in% marker_genes] #don't use marker genes for score calculation
  for(i in 3:ncol(gex)){
    cor_cdh1 <- cor.test(markers$CDH1, gex[,i])
    cor_cdh2 <- cor.test(markers$CDH2, gex[,i])
    cor_vim <- cor.test(markers$VIM, gex[,i])
    cor_fn1 <- cor.test(markers$FN1, gex[,i])
    correlates[i,] <- c(tested_gene=colnames(gex)[i],
                        CDH1_r=cor_cdh1$estimate, 
                        CDH1_pval=cor_cdh1$p.value,
                        CDH2_r=cor_cdh2$estimate, 
                        CDH2_pval=cor_cdh2$p.value,
                        VIM_r=cor_vim$estimate,
                        VIMR_pval=cor_vim$p.value,
                        FN1_r=cor_fn1$estimate,
                        FN1_pval=cor_fn1$p.value)
  }
  for(j in 2:ncol(correlates)){
    correlates[,j] <- as.numeric(correlates[,j])
  }
  CDH1_rs <- correlates$CDH1_r
  names(CDH1_rs) <- correlates$tested_gene
  CDH1_rs <- sort(CDH1_rs, decreasing = T)
  
  CDH2_rs <- correlates$CDH2_r
  names(CDH2_rs) <- correlates$tested_gene
  CDH2_rs <- sort(CDH2_rs, decreasing = T)  
  
  VIM_rs <- correlates$VIM_r
  names(VIM_rs) <- correlates$tested_gene
  VIM_rs <- sort(VIM_rs, decreasing = T)
  
  FN1_rs <- correlates$FN1_r
  names(FN1_rs) <- correlates$tested_gene
  FN1_rs <- sort(FN1_rs, decreasing = T)
  
  EMTgenes <- c(names(CDH1_rs)[1:25], names(CDH2_rs)[1:25], names(VIM_rs)[1:25], names(FN1_rs)[1:25])
  EMTgenes <- unique(EMTgenes)
  
  e_m <- character()
  for(i in 1:length(EMTgenes)){
    e_m[i] <- case_when(EMTgenes[i] %in% names(CDH1_rs)[1:25] ~ "E",
                        TRUE ~ "M")
  }
  names(EMTgenes) <- e_m
  epgex <- gex[,which(colnames(gex) %in% EMTgenes[names(EMTgenes) == "E"])] #gex from epithelial marker genes
  mesgex <- gex[,which(colnames(gex) %in% EMTgenes[names(EMTgenes) == "M"])] #gex from mes. marker genes
  EMTscore <- numeric(length = nrow(gex)) #create vector for EMT score
  epscore <- rowSums(epgex)/ncol(epgex) #calculate mean expression of epithelial genes
  messcore <- rowSums(mesgex)/ncol(mesgex) #calculate mean expression of mes. genes
  EMTscore <- epscore - messcore #EMT score = epithelial - mesenchymal
  #create data frame with cosmic id, tcga label and emt score
  EMTscore_df <- data.frame(COSMIC_ID = gex$`COSMIC ID`, TCGA = gex$`TCGA Desc`, EMT_score = EMTscore)
  return(list(EMT_genes = EMTgenes, EMT_score = EMTscore_df))
}
##################################


