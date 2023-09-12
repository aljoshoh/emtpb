### NEEDS: dat_skcm
library(org.Hs.eg.db)
library(httr)
library(jsonlite)
library(dplyr)
apiKey <- "6bd204bb24314605f5584fb1fb080669"

## Fucntions
query_lincs_id <- function(ID, # the InChIKey
                           type = "InChIKey", # which identifier 
                           which = "rep_samples",
                           start = 1,
                           amount = 1000
){ # which query
  base <- "https://api.clue.io"
  query <- paste0('/api/',which,'?filter={
             "where":{"',type,'":"',ID,'"},
             "limit":',amount,',"offset":',start,'
            }')
  call <- GET(gsub("\n","",gsub(" ","",paste0(base,query))), 
              add_headers(user_key = apiKey))
  print(call)
  data <- fromJSON(httr::content(call, "text"))
  return(data)
}

drugs_info_lincs <- read_delim("emtpb/data/lincs/compoundinfo_beta.txt", delim = "\t")
cl_info_lincs <- read_delim("emtpb/data/lincs/cellinfo_beta.txt", delim = "\t")

# pert_id: 
drugs_info_lincs[(drugs_info_lincs$cmap_name == "NVP-AUY922"),] #luminespib
drugs_info_lincs[(drugs_info_lincs$cmap_name == "staurosporine"),"pert_id"] #staurosporine
drugs_info_lincs[(drugs_info_lincs$cmap_name == "CHIR-99021"),"pert_id"] #staurosporine
ddd <- drugs_info_lincs[(drugs_info_lincs$cmap_name == "tanespimycin"),"pert_id"] #staurosporine





# QUERY
ID <- "BRD-K65182930" #luminespib (1) -> no data
ID <- "BRD-K70549064" #staurosporine (1) -> no cells
ID <- "BRD-K17953061" #staurosporine (2) -> only mes cells
ID <- "BRD-K76894938" #gsk (1) -> no cells
ID <- "BRD-K16189898" #gsk (2) -> no cells
ID <- "BRD-K11528507" # geldanamycin -> no cells
ID <- "BRD-A19500257" # geldanamycin hits !
ID <- "BRD-K41859756" #luminespib (2) -> no differential between epi and mes, but good !
#ID <- unlist(ddd[2,"pert_id"])


sig <- query_lincs_id(ID, type = "pert_id", which = "sigs") #only non-empty signature for this
sig <- sig[sig$is_gold,]

# cell lines that are skin cancer
lincs_anno <- read_csv(paste0("emtpb/data/Model.csv")) %>% dplyr::rename(`COSMIC ID` = COSMICID)
unique(sig$cell_id)

# mesenchymal signatures
mes <- dat_skcm[dat_skcm$mak > mean(na.omit(dat_skcm$mak)),]
mes <- mes[!is.na(mes$`COSMIC ID`),]
lincs_anno_mes <- dplyr::left_join(mes, lincs_anno, by = "COSMIC ID")
mes <- unique(sig$cell_id)[unique(sig$cell_id) %in% stringr::str_replace_all(lincs_anno_mes$CellLineName, "[[:punct:]]", "")]
message(length(mes)," mesenchymal cells")
mes <- sig[sig$cell_id %in% mes,]

# epithelial signatures
epi <- dat_skcm[dat_skcm$mak < mean(na.omit(dat_skcm$mak)),]
epi <- epi[!is.na(epi$`COSMIC ID`),]
lincs_anno_epi <- left_join(epi, lincs_anno, by = "COSMIC ID")
epi <- unique(sig$cell_id)[unique(sig$cell_id) %in% stringr::str_replace_all(lincs_anno_epi$CellLineName, "[[:punct:]]", "")]
message(length(epi)," epithelial cells")
epi <- sig[sig$cell_id %in% epi,]

sig <- mes
l <- list(); for(z in 1:nrow(sig)){
  d1 <- suppressMessages(ensembldb::select(org.Hs.eg.db, keys = sig$dn100_bing[[z]], keytype = "ENTREZID", columns = "SYMBOL"))
  d1$type <- "down"
  d2 <- suppressMessages(ensembldb::select(org.Hs.eg.db, keys = sig$up100_bing[[z]], keytype = "ENTREZID", columns = "SYMBOL"))
  d2$type <- "up"
  d <- rbind(d1, d2)
  d$cells <- list(unique(sig$cell_id[[z]]))
  d$sig_id <- z
  l[[i]] <- d
  i <- i+1
}; l <- do.call(rbind, l); l$emt <- "mes"; tmp <- l

sig <- epi
l <- list(); for(z in 1:nrow(sig)){
  d1 <- suppressMessages(ensembldb::select(org.Hs.eg.db, keys = sig$dn100_bing[[z]], keytype = "ENTREZID", columns = "SYMBOL"))
  d1$type <- "down"
  d2 <- suppressMessages(ensembldb::select(org.Hs.eg.db, keys = sig$up100_bing[[z]], keytype = "ENTREZID", columns = "SYMBOL"))
  d2$type <- "up"
  d <- rbind(d1, d2)
  d$cells <- list(unique(sig$cell_id[[z]]))
  d$sig_id <- z
  l[[i]] <- d
  i <- i+1
}; l <- do.call(rbind, l); l$emt <- "epi"; l <- rbind(l,tmp)


ll <- l
ll <- ll %>% 
  dplyr::select(c(type, emt, SYMBOL)) %>% 
  distinct %>%
  group_by(SYMBOL) %>%
  filter((n_distinct(paste(emt,type)) == 1)) %>%
  ungroup #%>% View
ll$id <- factor(paste(ll$emt, ll$type))

ll_epi_up <- ll$SYMBOL[ll$emt =="epi" & (ll$type == "up")]
ll_epi_dn <- ll$SYMBOL[ll$emt =="epi" & (ll$type == "down")]
ll_mes_up <- ll$SYMBOL[ll$emt =="mes" & (ll$type == "up")]
ll_mes_dn <- ll$SYMBOL[ll$emt =="mes" & (ll$type == "down")]

# check differential signatures
#cat(paste(ll_epi_dn,collapse = "\n"))
terms <- c("ChEA_2022","GO_Biological_Process_2023")
lll <- list(); for(id in levels(ll$id)){
  tmp <- enrichr(ll[ll$id == id,"SYMBOL",drop = TRUE], terms)
  tmp <- do.call(rbind, tmp)
  tmp$id <- id
  lll[[id]] <- tmp
}; lll_df <- do.call(rbind, lll)
lll_df <- lll_df[unlist(lapply(lll_df$Term, function(x) !grepl( "Mouse", x, fixed = TRUE))),]

file_lincs_enrichment <- "emtpb/metadata/lincs/luminespib_enrichments_diff.rds"
if(!file.exists(file_lincs_enrichment)){
  saveRDS(object = lll_df, file = file_lincs_enrichment)
  }


# check full signatures
l$id <- factor(paste(l$emt, l$type))
terms <- c("ChEA_2022","GO_Biological_Process_2023")
lll <- list(); for(id in levels(l$id)){
  tmp <- enrichr(l[l$id == id,"SYMBOL",drop = TRUE], terms)
  tmp <- do.call(rbind, tmp)
  tmp$id <- id
  lll[[id]] <- tmp
}; lll_df <- do.call(rbind, lll)
lll_df <- lll_df[unlist(lapply(lll_df$Term, function(x) !grepl( "Mouse", x, fixed = TRUE))),]

file_lincs_enrichment <- "emtpb/metadata/lincs/luminespib_enrichments_all.rds"
if(!file.exists(file_lincs_enrichment)){
  saveRDS(object = lll_df, file = file_lincs_enrichment)
}


# test for overall signature
file_lincs_enrichment <- "emtpb/metadata/lincs/luminespib_enrichments_all.rds"
t<- readRDS(file_lincs_enrichment)
t$emt <- unlist(lapply(t$id, function(x) strsplit(x," ")[[1]][1]))
t$updn <- unlist(lapply(t$id, function(x) strsplit(x," ")[[1]][2]))
t <- t[unlist(lapply(t$Term, function(x) grepl("GO:",x,fixed = TRUE))),]
View(t)





# get chir-99021 signature of A375 cells, not usable since only mesenchymal
test <- read_delim("/Users/alexander.ohnmacht/research/marisa/emtpb/metadata/lincs/lincs_chir-99021/sig_Tue_Sep_12_05_05_17_2023_6545018.xls")
chir_sig_up <-   (test %>% 
                 arrange(Significance_pvalue) %>% 
                 filter(Value_LogDiffExp > 0) %>% 
                 dplyr::select(c(Name_GeneSymbol)) %>%
                 distinct %>%
                 head(n = 100))$Name_GeneSymbol
chir_sig_dn <-   (test %>% 
                    arrange(Significance_pvalue) %>% 
                    filter(Value_LogDiffExp < 0) %>% 
                    dplyr::select(c(Name_GeneSymbol)) %>%
                    distinct %>%
                    head(n = 100))$Name_GeneSymbol
terms <- c("ChEA_2022","GO_Biological_Process_2023")
lll <- list(); for(id in c(1,2)){
  tmp <- enrichr(list(chir_sig_up, chir_sig_dn)[[id]], terms)
  tmp <- do.call(rbind, tmp)
  tmp$id <- list("up","dn")[[id]]
  lll[[id]] <- tmp
}; lll_df <- do.call(rbind, lll)
lll_df <- lll_df[unlist(lapply(lll_df$Term, function(x) !grepl( "Mouse", x, fixed = TRUE))),]


ID <- "BRD-K41859756" #luminespib (2) -> no differential between epi and mes, but good !
sig <- query_lincs_id(ID, type = "pert_id", which = "sigs") #only non-empty signature for this
sig <- sig[sig$is_gold,]
# cell lines that are skin cancer
lincs_anno <- read_csv(paste0("emtpb/data/Model.csv")) %>% dplyr::rename(`COSMIC ID` = COSMICID)
# mesenchymal signatures
mes <- dat_skcm[dat_skcm$mak > mean(na.omit(dat_skcm$mak)),]
mes <- mes[!is.na(mes$`COSMIC ID`),]
lincs_anno_mes <- dplyr::left_join(mes, lincs_anno, by = "COSMIC ID")
mes <- unique(sig$cell_id)[unique(sig$cell_id) %in% stringr::str_replace_all(lincs_anno_mes$CellLineName, "[[:punct:]]", "")]
message(length(mes)," mesenchymal cells")
mes <- sig[sig$cell_id %in% mes,]

mes_up <- suppressMessages(ensembldb::select(org.Hs.eg.db, keys = mes$up100_bing[[1]], keytype = "ENTREZID", columns = "SYMBOL"))$SYMBOL
mes_dn <- suppressMessages(ensembldb::select(org.Hs.eg.db, keys = mes$dn100_bing[[1]], keytype = "ENTREZID", columns = "SYMBOL"))$SYMBOL

intersect_dn <- enrichr(intersect(mes_dn, chir_sig_dn),terms);intersect_dn <- do.call(rbind, intersect_dn)

file_lincs_enrichment <- "emtpb/metadata/lincs/chir_990221_luminespib_enrichments_common_dn.rds"
if(!file.exists(file_lincs_enrichment)){
  saveRDS(object = intersect_dn, file = file_lincs_enrichment)
}


# differential signatures ? no...
ll <- l
ll %>% 
  dplyr::select(c(type, emt, SYMBOL)) %>% 
  distinct %>%
  group_by(SYMBOL) %>%
  filter((n_distinct(paste(emt, type)) == 2) & (n_distinct(emt) == 2) & (n_distinct(type) == 2)) %>%
  ungroup %>% View

  


