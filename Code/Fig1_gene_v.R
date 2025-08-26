
library(dplyr)
library(GeneNet)
library(Rgraphviz)
library(MCL)
library(ggplot2)
library(linkET)

################################################################################
### 0.0 ####P1 UC1 Link Genes##
rm(list = ls())
path <- getwd()
setwd(path)

input <- "input"
output <- "output"
if (!dir.exists(output)) {
  dir.create(output,recursive = TRUE)
}
genes_of_interest <- c("PTGR1","FOSB","CXCL5")
##network
algorithm_genes <- read.csv(file = paste0(input,"/CoreHD-DCI_genes.csv"), header = TRUE)
id_name <- read.csv(file = paste0(input,"/id_genename.csv"), header = TRUE)
P1_layer <- read.table(file = paste0(input,"/P1_layer.txt"), header = FALSE)
UC1_layer <- read.table(file = paste0(input,"/UC1_layer.txt"), header = FALSE)
P1_Conlayer <- read.table(file = paste0(input,"/P1_Conlayer.txt"), header = FALSE)
UC1_Conlayer <- read.table(file = paste0(input,"/UC1_Conlayer.txt"), header = FALSE)

################################################################################
# 1.0 mantel figures
# data read
P1_DEA <- read.csv(file = paste0(input,"/DEG_P1.csv"),
                   header = T)
colnames(P1_DEA)[1] <- "SYMBOL"
UC1_DEA <- read.csv(file = paste0(input,"/DEG_UC1.csv"),
                    header = T)
colnames(UC1_DEA)[1] <- "SYMBOL"

SORT_P1_DEA <- P1_DEA
SORT_P1_DEA <- P1_DEA %>%
  arrange(desc(abs(logFC))) %>%
  filter(adj.P.Val < 0.05) %>%
  head(round(0.09 * nrow(SORT_P1_DEA)))

SORT_UC1_DEA <- UC1_DEA
SORT_UC1_DEA <- UC1_DEA %>%
  arrange(desc(abs(logFC))) %>%
  filter(adj.P.Val < 0.05) %>%
  head(round(0.09 * nrow(SORT_UC1_DEA)))

GeneSymOlap <- intersect(SORT_P1_DEA$SYMBOL,SORT_UC1_DEA$SYMBOL)
length(GeneSymOlap)

P1_expr <- read.table(file = paste0(input,"/expr_P1.txt"),
                      sep = "\t",
                      header = T,
                      row.names = 1)
UC1_expr <- read.table(file = paste0(input,"/expr_UC1.txt"),
                       sep = "\t",
                       header = T,
                       row.names = 1)

matrix_P1 <- P1_expr[GeneSymOlap, ]
matrix_UC1 <- UC1_expr[GeneSymOlap,]

GeneSymOlap_DF <- data.frame(GeneSymOlap)
group_P1 <- read.csv(file = paste0(input,"/info_P1.csv"))
group_UC1 <- read.csv(file = paste0(input,"/info_UC1.csv"))

P1_expr_affected <- matrix_P1[,group_P1$Group==1]
P1_expr_unaffected <- matrix_P1[,group_P1$Group==0]

UC1_expr_affected <- matrix_UC1[,group_UC1$Group==1]
UC1_expr_unaffected <- matrix_UC1[,group_UC1$Group==0]
##############
### common DE -  expression matrix for diseases. 
P1_expr_affected_t <- t(P1_expr_affected)
UC1_expr_affected_t <- t(UC1_expr_affected)

P1_expr_affected_t_KEY <- P1_expr_affected_t[,c("PTGR1","FOSB","CXCL5")]
UC1_expr_affected_t_KEY <- UC1_expr_affected_t[,c("PTGR1","FOSB","CXCL5")]

P1_expr_affected_t_sel <- as.data.frame(P1_expr_affected_t) %>% select(-c("PTGR1","FOSB","CXCL5"))
UC1_expr_affected_t_sel <- as.data.frame(UC1_expr_affected_t) %>% select(-c("PTGR1","FOSB","CXCL5"))

P1_expr_affected_t_sel1 <- P1_expr_affected_t_sel[,1:50]
UC1_expr_affected_t_sel1 <- UC1_expr_affected_t_sel[,1:50]

###

################################################################################
# 2.0 ### select top 30 genes/ 10 highest r for each key genes USING correlation - linkET --  pearson  ##### 

################################################################################
P1_corr <- correlate(P1_expr_affected_t)
UC1_corr <- correlate(UC1_expr_affected_t)

P1_corr_r <- P1_corr$r
P1_corr_p <- P1_corr$p
UC1_corr_r <- UC1_corr$r
UC1_corr_p <- UC1_corr$p

P1_corr_r_sel <- as.data.frame(P1_corr_r[,genes_of_interest])
P1_corr_p_sel <- as.data.frame(P1_corr_p[,genes_of_interest])
UC1_corr_r_sel <- as.data.frame(UC1_corr_r[,genes_of_interest])
UC1_corr_p_sel <- as.data.frame(UC1_corr_p[,genes_of_interest])

all_genes <- rownames(P1_corr_r)
other_genes <- all_genes[!(all_genes %in% genes_of_interest)]
P1_corr_r_sel_sum <- as.data.frame(list(Sum = rowSums(P1_corr_r_sel[other_genes,]), Name = other_genes))  ### calculate the sum of correlation. 
UC1_corr_r_sel_sum <- as.data.frame(list(Sum = rowSums(UC1_corr_r_sel[other_genes,]), Name = other_genes))  ### calculate the sum of correlation. 

Number_sel = 30 
top_30_genes_p1 <- P1_corr_r_sel_sum[order(-P1_corr_r_sel_sum$`Sum`),"Name"][1:Number_sel]
top_30_genes_uc1 <- UC1_corr_r_sel_sum[order(-UC1_corr_r_sel_sum$`Sum`),"Name"][1:Number_sel]

#P1_selected_corr_matrix <- P1_corr_r[c(genes_of_interest, top_30_genes), c(genes_of_interest, top_30_genes)]
P1_expr_affected_t_sel2 <- P1_expr_affected_t_sel[,top_30_genes_p1]
UC1_expr_affected_t_sel2 <- UC1_expr_affected_t_sel[,top_30_genes_uc1]

#top_30_genes_p1 P1_corr_r_sel P1_corr_p_sel
#top_30_genes_uc1 UC1_corr_r_sel UC1_corr_p_sel
corr_p1_framedata <- as.data.frame(list(
  spec = c(rep(genes_of_interest[1], 30), rep(genes_of_interest[2], 30), rep(genes_of_interest[3], 30)),
  env = c(top_30_genes_p1,top_30_genes_p1,top_30_genes_p1),
  r = c(P1_corr_r_sel[top_30_genes_p1,genes_of_interest[1]], 
        P1_corr_r_sel[top_30_genes_p1,genes_of_interest[2]], 
        P1_corr_r_sel[top_30_genes_p1,genes_of_interest[3]]),
  p = c(P1_corr_p_sel[top_30_genes_p1,genes_of_interest[1]], 
        P1_corr_p_sel[top_30_genes_p1,genes_of_interest[2]], 
        P1_corr_p_sel[top_30_genes_p1,genes_of_interest[3]])
  ))%>% 
  mutate(R = cut(r, breaks = c(-Inf, -0.3, 0.3, Inf),
                  labels = c("< -0.3", "-0.3 - 0.3",">= 0.3")),
         P = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
  # mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.6, Inf),
  #                 labels = c("< 0.2", "0.2 - 0.6", ">= 0.6")),
  #        pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
  #                 labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

corr_uc1_framedata <- as.data.frame(list(
  spec = c(rep(genes_of_interest[1], 30), rep(genes_of_interest[2], 30), rep(genes_of_interest[3], 30)),
  env = c(top_30_genes_uc1,top_30_genes_uc1,top_30_genes_uc1),
  r = c(UC1_corr_r_sel[top_30_genes_uc1,genes_of_interest[1]], 
        UC1_corr_r_sel[top_30_genes_uc1,genes_of_interest[2]], 
        UC1_corr_r_sel[top_30_genes_uc1,genes_of_interest[3]]),
  p = c(UC1_corr_p_sel[top_30_genes_uc1,genes_of_interest[1]], 
        UC1_corr_p_sel[top_30_genes_uc1,genes_of_interest[2]], 
        UC1_corr_p_sel[top_30_genes_uc1,genes_of_interest[3]])
))%>% 
  mutate(R = cut(r, breaks = c(-Inf, -0.3, 0.3, Inf),
                  labels = c("< -0.3", "-0.3 - 0.3",">= 0.3")),
         P = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
# mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.6, Inf),
#                 labels = c("< 0.2", "0.2 - 0.6", ">= 0.6")),
#        pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
#                 labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))


################################################################################
# 3.0 ### visualization / generation of figure. 
# 
#set_corrplot_style(scale = ggplot2::scale_fill_viridis_c())
# set_default_style()

#######
# (qcorrplot1, qcorrplot2, ncol = 2)
#library(gridExtra)
#grid.arrange(qcorrplotP1, qcorrplotUC1, ncol = 2)

################################################################################
# 3.1 ### visualization / generation of figure. version 2

#set_corrplot_style()
qcorrplotP1 <- qcorrplot(correlate(P1_expr_affected_t_sel2), 
                         type = "upper",
                         diag = TRUE,
) +
  geom_square()+
  geom_mark(size=2,
            only_mark=T,
            sig_level=c(0.05,0.01,0.001),
            sig_thres=0.05)+
  geom_couple(aes(colour = R, size = P),data = corr_p1_framedata, curvature = 0.1,
              node.colour = c("blue", "blue"),
              node.fill = c("grey", "grey"),
              node.size = c(3.5, 3),
              label.size = 3,
              label.colour = "black"
  )+
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"),
                       limits = c(-1, 1),
                       breaks = seq(-1,1,0.5))+
  geom_mark(size=2.5,
            only_mark=T,
            sig_level=c(0.05,0.01,0.001),
            sig_thres=0.05)+
  scale_size_manual(values = c(1.,0.75,0.5)) +
  scale_colour_manual(values = c(">= 0.3" = "#2c7bb6","-0.3 - 0.3" ="#CCCCCC99","< -0.3" = "#fca3a3"))
ggsave("qcorrplotP1.pdf", plot = qcorrplotP1, width = 8,
       height = 6)
print(qcorrplotP1)

#######
qcorrplotUC1 <- qcorrplot(correlate(UC1_expr_affected_t_sel2), 
                          type = "lower",
                          diag = TRUE,
) +
  geom_square()+
  geom_mark(size=2.5,
            only_mark=T,
            sig_level=c(0.05,0.01,0.001),
            sig_thres=0.05)+
  geom_couple(aes(colour = R, size = P),data = corr_uc1_framedata, curvature = 0.1,
              node.colour = c("blue", "blue"),
              node.fill = c("grey", "grey"),
              node.size = c(3.5, 2.5),
              label.size = 3,
              label.colour = "black"
  )+
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"),
                       limits = c(-1, 1),
                       breaks = seq(-1,1,0.5))+
  geom_mark(size=2.5,
            only_mark=T,
            sig_level=c(0.05,0.01,0.001),
            sig_thres=0.05)+
  scale_size_manual(values = c(1.,0.75,0.5)) +
  scale_colour_manual(values = c(">= 0.3" = "#2c7bb6","-0.3 - 0.3" ="#CCCCCC99","< -0.3" = "#fca3a3"))

ggsave("qcorrplotUC1.pdf", plot = qcorrplotUC1,width = 8,
       height = 6)
print(qcorrplotUC1)

#"#D95F02""#1B9E77""#CCCCCC99"
################################################################################































################################################################################































################################################################################















































