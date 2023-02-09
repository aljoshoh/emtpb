library(ggplot2)
library(gridExtra)
library(ggrepel)
library(readr)
library(dplyr)
#`194-GDSC1`,matrix_resp$`1559-GDSC2`

EMTscores <- read_csv("metadata/EMTscores.csv")
matrix_resp <- read_csv("metadata/matrix_resp.csv")

tmp1 <- read_csv("data/PANCANCER_Genetic_features_Sat Feb  4 15_33_58 2023.csv")
tmp2 <- read_csv("data/PANCANCER_Genetic_features_Sat Feb  4 15_42_40 2023.csv")
tmp <- full_join(tmp1, tmp2)
cell_line_names <- tmp[,c("COSMIC ID","Cell Line Name")] %>% distinct()
m <- full_join(EMTscores, matrix_resp)[,c("COSMIC ID","TCGA Desc","1559-GDSC2","194-GDSC1","EMT_score")]
m <- full_join(m, cell_line_names)
m$color <- "-"
m$color[m$`Cell Line Name` %in% ("COLO-783")] <- "mismatched"
m$color[m$`Cell Line Name` %in% ("UACC-257")] <- "mismatched"
m$color[m$`Cell Line Name` %in% c("SK-MEL-5","A375","RPMI-7951","IGR-37","SK-MEL-24")] <- "ordered before"


###
grid.arrange(
  ggplot(data=m[m$`TCGA Desc` == "SKCM",], aes(x= m$EMT_score[m$`TCGA Desc` == "SKCM"], 
           y = m$`1559-GDSC2`[m$`TCGA Desc` == "SKCM"]))+
  geom_point(alpha = 0.3, aes(color = color))+
  #geom_vline(xintercept = -2.380233)+
  #geom_hline(yintercept = -2.452889)+
  ylab("log(ic50)")+
  xlab("mesenchymal<>epithelial")+
  theme_minimal()+
  ggtitle("GDSC2")+
  geom_text_repel(aes(label = m$`Cell Line Name`[m$`TCGA Desc` == "SKCM"]),size = 2),
  ggplot(data=m[m$`TCGA Desc` == "SKCM",], aes(x= m$EMT_score[m$`TCGA Desc` == "SKCM"], 
                                             y = m$`194-GDSC1`[m$`TCGA Desc` == "SKCM"]))+
  geom_point(alpha = 0.3, aes(color = color))+
  #geom_vline(xintercept = -2.380233)+
  #geom_hline(yintercept = -1.569377)+
  ylab("log(ic50)")+
  xlab("mesenchymal<>epithelial")+
  theme_minimal()+
  ggtitle("GDSC1")+
  geom_text_repel(aes(label = m$`Cell Line Name`[m$`TCGA Desc` == "SKCM"]),size = 2) 
)
#plot(m$EMT_score[m$`TCGA Desc` == "SKCM"], m$`194-GDSC1`[m$`TCGA Desc` == "SKCM"])
