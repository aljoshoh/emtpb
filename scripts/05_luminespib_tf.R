library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggbeeswarm)
path <- "/Users/alexander.ohnmacht/research/marisa/emtpb/"
resp <- read_csv(paste0(path,"metadata/matrix_resp.csv")) #%>% dplyr::select(-c("COSMIC ID","TCGA Desc"))
EMTscores <- read_csv(paste0(path,"metadata/EMTscores.csv")) #%>% dplyr::select(-c("COSMIC ID","TCGA Desc"))
rna <- read_csv(paste0(path,"metadata/matrix_exp.csv")) %>% dplyr::filter(`COSMIC ID` %in% (EMTscores$`COSMIC ID`[EMTscores$`TCGA Desc` == "SKCM"]))
EMT_gs <- readRDS(paste0(path,"metadata/EMT_gs.rds"))

#Pearson:
if(!file.exists("metadata/cor_gex_luminespib.rds")){
  gex_gdsc_t <- rna
  skcm_emt <- EMTscores[EMTscores$`TCGA Desc` == "SKCM",]
  skcm_emt <- left_join(skcm_emt, resp[,c("COSMIC ID","1559-GDSC2"),drop = F])
  cor_gex_luminespib <- data.frame(gene=character(), R=numeric(), pval=numeric())
  for(gene in 3:ncol(gex_gdsc_t)){
    curr_gene <- colnames(gex_gdsc_t)[gene]
    tmp <- left_join(skcm_emt, gex_gdsc_t[, c(1,gene)], by = c("COSMIC ID" = "COSMIC ID"))
    tmp_cors <- cor.test(unlist(tmp[,"1559-GDSC2"]), unlist(tmp[,5]))
    cor_gex_luminespib <- rbind(cor_gex_luminespib, data.frame(gene=curr_gene, R=tmp_cors$estimate, pval=tmp_cors$p.value))
  }
  saveRDS(cor_gex_luminespib, file = "metadata/cor_gex_luminespib.rds")
}else{
  cor_gex_luminespib <- readRDS("metadata/cor_gex_luminespib.rds")
}

cor_gex_luminespib$fdr <- p.adjust(cor_gex_luminespib$pval, method = "BH")

cor_gex_luminespib$in_EMT_gs <- ifelse(cor_gex_luminespib$gene %in% EMT_gs$SKCM, "yes", "no")
cor_gex_luminespib$in_EMT_gs[cor_gex_luminespib$gene == "CDH1"] <- "seed"
cor_gex_luminespib$in_EMT_gs[cor_gex_luminespib$gene == "CDH2"] <- "seed"
cor_gex_luminespib$in_EMT_gs[cor_gex_luminespib$gene == "VIM"] <- "seed"
cor_gex_luminespib$in_EMT_gs[cor_gex_luminespib$gene == "FN1"] <- "seed"

logical_pos <- cor_gex_luminespib$fdr < 0.1 & cor_gex_luminespib$R > 0
logical_neg <- cor_gex_luminespib$fdr < 0.1 & cor_gex_luminespib$R < 0
cat(cor_gex_luminespib$gene[logical_pos])
cat(cor_gex_luminespib$gene[logical_neg])
ggplot(cor_gex_luminespib, aes(x=R, y=-log10(pval))) + geom_point(size = 0.5) + 
  geom_text_repel(aes(label=gene, color = in_EMT_gs),data=cor_gex_luminespib[logical,], #fdr < 0.01 and R > 0.6
                  size = 3) + theme_classic() + ggtitle("Correlation of luminespib drug response and gene expression") + 
  xlab("Pearson's R") + scale_color_manual(values = c("black", "darkolivegreen3", "darkolivegreen4")) +
  theme(legend.position = "none")



# supplement ? boxplot for mitf
skcm_emt <- skcm_emt %>%
  mutate(drugresp_strat = case_when(`1559-GDSC2` < median(na.omit(`1559-GDSC2`)) ~ "sensitive", T ~ "resistant"))
skcm_emt <- full_join(skcm_emt, rna[,c("COSMIC ID","MITF")])
ggplot(skcm_emt, aes(x = drugresp_strat, y = MITF, color = drugresp_strat)) + geom_beeswarm() + 
  stat_summary(fun = mean,geom = "crossbar", color = "black", width = 0.5) +
  scale_color_manual(values = c("resistant" = "darkolivegreen4", "sensitive" = "darkolivegreen2")) + theme_classic() +
  theme(legend.position = "none") + xlab("drug response") + ylab("MITF expression level") +
  ggtitle("MITF expression and response to luminespib",
          subtitle = paste("Welch t-test p-value:", format(round(t.test(skcm_emt$MITF[skcm_emt$drugresp_strat == "sensitive"], skcm_emt$MITF[skcm_emt$drugresp_strat == "resistant"])$p.value, digits = 6), scientific = F)))


#### limma for enrichments ?
rna <- read_csv(paste0(path,"metadata/matrix_exp.csv")) %>% dplyr::filter(`COSMIC ID` %in% (EMTscores$`COSMIC ID`[EMTscores$`TCGA Desc` == "SKCM"]))
gex_gdsc_t <- rna
skcm_emt <- EMTscores[EMTscores$`TCGA Desc` == "SKCM",]
skcm_emt <- left_join(skcm_emt, resp[,c("COSMIC ID","1559-GDSC2"),drop = F])
tmp <- left_join(skcm_emt, gex_gdsc_t, by = c("COSMIC ID" = "COSMIC ID"))
tmp <- na.omit(tmp)
exp_skcm <- tmp[,-c(1:5)]
luminespib_resp_skcm <- tmp$`1559-GDSC2`

my_design <- model.matrix(~luminespib_resp_skcm)
library(limma)
my_fit <- lmFit(exp_skcm%>%as.matrix%>%t, my_design)
my_ebayes <- eBayes(my_fit)
my_results <- topTable(my_ebayes, n=nrow(my_ebayes))%>%rownames_to_column("gene")


limma_luminespib_pos <- my_results$gene[(my_results$adj.P.Val < 0.1)&(my_results$logFC > 0)]
limma_luminespib_neg <- my_results$gene[(my_results$adj.P.Val < 0.1)&(my_results$logFC < 0)]

library(enrichR)


terms <- c("ARCHS4_Tissues","ChEA_2022")
enriched_pos <- enrichr(limma_luminespib_pos, terms); enriched_pos <- do.call(rbind, enriched_pos); enriched_pos$sign <- "positive"
enriched_neg <- enrichr(limma_luminespib_neg, terms); enriched_neg <- do.call(rbind, enriched_neg); enriched_neg$sign <- "negative"
enriched <- rbind(enriched_pos, enriched_neg)


strsplit(enriched$Genes[enriched$Term == "MITF 21258399 ChIP-Seq MELANOMA Human"],";")

colors_genes <- c(limma_luminespib_neg, limma_luminespib_pos)
names(colors_genes) = colors_genes
colors_genes[rep(T, length(colors_genes))] <- "black"
colors_genes[names(colors_genes) %in% strsplit(enriched$Genes[enriched$Term == "MITF 21258399 ChIP-Seq MELANOMA Human" & (enriched$sign == "positive")],";")[[1]]] <- "deepskyblue"
colors_genes[names(colors_genes) %in% strsplit(enriched$Genes[enriched$Term == "FIBROBLAST" & (enriched$sign == "negative")],";")[[1]]] <- "darkorange"
colors_genes[names(colors_genes) %in% c("CDH1","CDH2")] <- "deepskyblue4"
saveRDS(colors_genes, file = "metadata/colors_genes_enrichr_luminespib_melanoma.rds")

