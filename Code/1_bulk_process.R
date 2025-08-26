
# Step 1: Data Acquisition and Preprocessing ------------------------------
# Data Download & Clinical Info Extraction
# Gene Symbol Mapping & Filtering
# Log Transformation & Quality Check

# Output: Exported cleaned expression matrix 
# and filtered clinical information for downstream analysis.



##Periodontitis: GSE16134 GPL570
##UC1: GSE87466 GPL13158
##UC2: GSE59071 GPL6244
##UC3: GSE48958 GPL6244


library(idmap1)
library(devtools)
library(GEOmirror)
library(GEOquery)
library(AnnoProbe)
library(sva)
library(ggplot2)
library(ggforce)
library(dplyr)
library(httr)
library(jsonlite)

rm(list = ls())
path <- getwd()
setwd(path)

output <- "1_Processed"
if (!dir.exists(output)) {
  dir.create(output)
}

getOption("timeout")
options(timeout = 100000)

##Downloading
##P GSE16134
eset_P1_raw <- geoChina('GSE16134')
Anno_Info_GPL570 <- idmap('GPL570',type = 'soft')

##Clinical Annotation
##P1
eset_P1 <- eset_P1_raw[[1]]

phenoDat_P1 <- pData(eset_P1)
dim(phenoDat_P1)
colnames(phenoDat_P1)
clinical_data3 <- phenoDat_P1[,c("geo_accession","source_name_ch1","gingival tissue phenotype:ch1")] 

rownames(clinical_data3) <- NULL
colnames(clinical_data3) <- c("Samples","Source","Characteristics")
head(clinical_data3)

clinical_data3 <- clinical_data3 %>%
  mutate(Group = ifelse(Characteristics == "Healthy (no BoP + PPD of<5 mm + CAL of<3 2 mm", 0, 1))

write.csv(clinical_data3,
          file = paste0(output,"/Info_P1.csv"),
          row.names = F)





##Gene Symbol Mapping
data_P1 <- exprs(eset_P1)
data_P1 <- as.data.frame(data_P1)

##P1
same_P1 <- match(rownames(data_P1),Anno_Info_GPL570$ID)
data_P1$Gene_Symbol <- Anno_Info_GPL570[same_P1,c("symbol")]
data_P1$Gene_Symbol <- sub("///.*", "", data_P1$Gene_Symbol)
data_P1$Gene_Symbol <- trimws(data_P1$Gene_Symbol)  

url <- "https://rest.genenames.org/fetch/status/Approved"
response <- GET(url, accept("application/json"))

hgnc_data <- fromJSON(content(response, "text"))
hgnc_symbols <- hgnc_data$response$docs$symbol
head(hgnc_symbols)  

removed_genes <- setdiff(data_P1$Gene_Symbol, hgnc_symbols)
write.table(removed_genes, 
            file =paste0(output,"/Removed_symbols_P1.txt") , 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


data_P1 <- data_P1[data_P1$Gene_Symbol %in% hgnc_symbols, ]

## Handling Redundant Probes
data_P1$Mean_Expression <- rowMeans(data_P1[, !names(data_P1) %in% "Gene_Symbol"], na.rm = TRUE)
data_P1 <- data_P1 %>%
  group_by(Gene_Symbol) %>%
  slice_max(order_by = Mean_Expression, n = 1) %>%
  ungroup()
data_P1$Mean_Expression <- NULL
data_P1 <- as.data.frame(data_P1)
rownames(data_P1) <- data_P1$Gene_Symbol
data_P1 <- data_P1[,!colnames(data_P1)%in%c("Gene_Symbol")]


## Log Transformation
ex <- data_P1  
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))  # 计算分位数


LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)


if (LogC) { 
  ex[which(ex <= 0)] <- NaN 
  data_P1 <- log2(ex) 
  print("log2 transform is finished")
} else {
  print("log2 transform is not needed")
} 

sum(data_P1 == "")  
sum(data_P1 == 0)    
boxplot(data_P1, main = "Boxplot after cleaning", las = 2)

write.table(data.frame(gene_symbol=rownames(data_P1),data_P1),
            file = paste0(output,"/expr_P1.txt"),  
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)





##UC1
eset_UC1_raw <- geoChina('GSE87466')


Anno_Info_GPL13158 <- read.table("GPL13158-5065.txt",sep = "\t",quote = "",header = T,fill = T)
Anno_Info_GPL13158 <- Anno_Info_GPL13158[,c("ID","Gene.Symbol")]
colnames(Anno_Info_GPL13158) <- c("ID","symbol")


eset_UC1 <- eset_UC1_raw[[1]]

phenoDat_UC1 <- pData(eset_UC1)
dim(phenoDat_UC1)
colnames(phenoDat_UC1)
clinical_data1 <- phenoDat_UC1[,c("geo_accession","age:ch1","disease:ch1",
                                  "diseaseextent:ch1","tissue:ch1")] 

rownames(clinical_data1) <- NULL
colnames(clinical_data1) <- c("Samples","Age","Diagnosis",
                              "Disease extent","Tissue")
head(clinical_data1)
table(clinical_data1$Diagnosis)
clinical_data1 <- clinical_data1 %>%
  mutate(Group = ifelse(Diagnosis == "Normal", 0,1))

table(clinical_data1$Group)
write.csv(clinical_data1,
          file = paste0(output,"/Info_UC1.csv"),
          row.names = F)



data_UC1 <- exprs(eset_UC1)
data_UC1 <- as.data.frame(data_UC1)


boxplot(data_UC1,las=2)

same_UC1 <- match(rownames(data_UC1),Anno_Info_GPL13158$ID)
data_UC1$Gene_Symbol <- Anno_Info_GPL13158[same_UC1,c("symbol")]
data_UC1$Gene_Symbol <- sub("///.*", "", data_UC1$Gene_Symbol)
data_UC1$Gene_Symbol <- trimws(data_UC1$Gene_Symbol)  


url <- "https://rest.genenames.org/fetch/status/Approved"
response <- GET(url, accept("application/json"))

hgnc_data <- fromJSON(content(response, "text"))
hgnc_symbols <- hgnc_data$response$docs$symbol
head(hgnc_symbols)  

removed_genes <- setdiff(data_UC1$Gene_Symbol, hgnc_symbols)
write.table(removed_genes, 
            file =paste0(output,"/Removed_symbols_UC1.txt") , 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


data_UC1 <- data_UC1[data_UC1$Gene_Symbol %in% hgnc_symbols, ]


data_UC1$Mean_Expression <- rowMeans(data_UC1[, !names(data_UC1) %in% "Gene_Symbol"], na.rm = TRUE)
data_UC1 <- data_UC1 %>%
  group_by(Gene_Symbol) %>%
  slice_max(order_by = Mean_Expression, n = 1) %>%
  ungroup()
data_UC1$Mean_Expression <- NULL
data_UC1 <- as.data.frame(data_UC1)
rownames(data_UC1) <- data_UC1$Gene_Symbol
data_UC1 <- data_UC1[,!colnames(data_UC1)%in%c("Gene_Symbol")]


ex <- data_UC1 
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))  # 计算分位数


LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)


if (LogC) { 
  ex[which(ex <= 0)] <- NaN  # 将小于等于 0 的值设为 NaN
  data_UC1 <- log2(ex) 
  print("log2 transform is finished")
} else {
  print("log2 transform is not needed")
} 



sum(data_UC1 == "")  
sum(data_UC1 == 0)   

boxplot(data_UC1, main = "Boxplot after cleaning", las = 2)



write.table(data.frame(gene_symbol=rownames(data_UC1),data_UC1),
            file = paste0(output,"/expr_UC1.txt"),  
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)





##UC2

eset_UC2_raw <- geoChina('GSE59071')
Anno_Info_GPL6244 <- idmap('GPL6244',type = 'soft')


eset_UC2 <- eset_UC2_raw[[1]]

phenoDat_UC2 <- pData(eset_UC2)
dim(phenoDat_UC2)


colnames(phenoDat_UC2)
clinical_data2 <- phenoDat_UC2[,c("geo_accession","disease activity:ch1",
                                  "disease:ch1","tissue:ch1")] 

rownames(clinical_data2) <- NULL

colnames(clinical_data2) <- c("Samples","Disease activity","Diagnosis","Tissue")

head(clinical_data2)

clinical_data2 <- clinical_data2 %>%
  mutate(Group = ifelse(Diagnosis == "control", 0, 
                        ifelse(Diagnosis == "ulcerative colitis",1,2)))
table(clinical_data2$Diagnosis)
table(clinical_data2$Group)

clinical_data2$Label <- c(rep("CD Active",8),rep("Control",11),
                          rep("UC Active",74),rep("UC Inactive",23))
write.csv(clinical_data2,
          file = paste0(output,"/Info_UC2.csv"),
          row.names = F)


data_UC2 <- exprs(eset_UC2)
data_UC2 <- as.data.frame(data_UC2)


boxplot(data_UC2,las=2)

same_UC2 <- match(rownames(data_UC2),Anno_Info_GPL6244$ID)
data_UC2$Gene_Symbol <- Anno_Info_GPL6244[same_UC2,c("symbol")]
data_UC2$Gene_Symbol <- sub("///.*", "", data_UC2$Gene_Symbol)
data_UC2$Gene_Symbol <- trimws(data_UC2$Gene_Symbol) 


url <- "https://rest.genenames.org/fetch/status/Approved"
response <- GET(url, accept("application/json"))

hgnc_data <- fromJSON(content(response, "text"))
hgnc_symbols <- hgnc_data$response$docs$symbol
head(hgnc_symbols)  

removed_genes <- setdiff(data_UC2$Gene_Symbol, hgnc_symbols)
write.table(removed_genes, 
            file =paste0(output,"/Removed_symbols_UC2.txt") , 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


data_UC2 <- data_UC2[data_UC2$Gene_Symbol %in% hgnc_symbols, ]


data_UC2$Mean_Expression <- rowMeans(data_UC2[, !names(data_UC2) %in% "Gene_Symbol"], na.rm = TRUE)
data_UC2 <- data_UC2 %>%
  group_by(Gene_Symbol) %>%
  slice_max(order_by = Mean_Expression, n = 1, with_ties = FALSE) %>%
  ungroup()
data_UC2$Mean_Expression <- NULL
data_UC2 <- as.data.frame(data_UC2)
rownames(data_UC2) <- data_UC2$Gene_Symbol
data_UC2 <- data_UC2[,!colnames(data_UC2)%in%c("Gene_Symbol")]


ex <- data_UC2  
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))  # 计算分位数


LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)


if (LogC) { 
  ex[which(ex <= 0)] <- NaN  
  data_UC2 <- log2(ex) 
  print("log2 transform is finished")
} else {
  print("log2 transform is not needed")
} 



sum(data_UC2 == "") 
sum(data_UC2 == 0)   
boxplot(data_UC2, main = "Boxplot after cleaning", las = 2)


write.table(data.frame(gene_symbol=rownames(data_UC2),data_UC2),
            file = paste0(output,"/expr_UC2.txt"),  
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)


##UC3

eset_UC3_raw <- geoChina('GSE48958')


Anno_Info_GPL6244 <- idmap('GPL6244',type = 'soft')




eset_UC3 <- eset_UC3_raw[[1]]

phenoDat_UC3 <- pData(eset_UC3)
dim(phenoDat_UC3)


colnames(phenoDat_UC3)
clinical_data3 <- phenoDat_UC3[,c("geo_accession","disease activity:ch1",
                                  "disease:ch1","tissue:ch1" )] 

rownames(clinical_data3) <- NULL

colnames(clinical_data3) <- c(
  "Samples",
  "Disease Activity",
  "Diagnosis",
  "Tissue"
)

head(clinical_data3)
table(clinical_data3$`Disease Activity`)

clinical_data3 <- clinical_data3 %>%
  mutate(Group = ifelse(Diagnosis == "control", 0, 1)) 

table(clinical_data3$Group)


write.csv(clinical_data3,
          file = paste0(output,"/Info_UC3.csv"),
          row.names = F)


data_UC3 <- exprs(eset_UC3)
data_UC3 <- as.data.frame(data_UC3)


boxplot(data_UC3,las=2)
same_UC3 <- match(rownames(data_UC3),Anno_Info_GPL6244$ID)
data_UC3$Gene_Symbol <- Anno_Info_GPL6244[same_UC3,c("symbol")]
data_UC3$Gene_Symbol <- sub("///.*", "", data_UC3$Gene_Symbol)
data_UC3$Gene_Symbol <- trimws(data_UC3$Gene_Symbol)  #去除两端的空格

url <- "https://rest.genenames.org/fetch/status/Approved"
response <- GET(url, accept("application/json"))

hgnc_data <- fromJSON(content(response, "text"))
hgnc_symbols <- hgnc_data$response$docs$symbol
head(hgnc_symbols)  

removed_genes <- setdiff(data_UC3$Gene_Symbol, hgnc_symbols)
write.table(removed_genes, 
            file =paste0(output,"/Removed_symbols_UC3.txt") , 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


data_UC3 <- data_UC3[data_UC3$Gene_Symbol %in% hgnc_symbols, ]


data_UC3$Mean_Expression <- rowMeans(data_UC3[, !names(data_UC3) %in% "Gene_Symbol"], na.rm = TRUE)
data_UC3 <- data_UC3 %>%
  group_by(Gene_Symbol) %>%
  slice_max(order_by = Mean_Expression, n = 1, with_ties = FALSE) %>%
  ungroup()

data_UC3$Mean_Expression <- NULL
data_UC3 <- as.data.frame(data_UC3)
rownames(data_UC3) <- data_UC3$Gene_Symbol
data_UC3 <- data_UC3[,!colnames(data_UC3)%in%c("Gene_Symbol")]


ex <- data_UC3 
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))  # 计算分位数


LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)


if (LogC) { 
  ex[which(ex <= 0)] <- NaN  
  data_UC3 <- log2(ex) 
  print("log2 transform is finished")
} else {
  print("log2 transform is not needed")
} 



sum(data_UC3 == "") 
sum(data_UC3 == 0)    
boxplot(data_UC3, main = "Boxplot after cleaning", las = 2)


write.table(data.frame(gene_symbol=rownames(data_UC3),data_UC3),
            file = paste0(output,"/expr_UC3.txt"),  
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)




# Step 2: Differential Expression Analysis for Single Disease  --------

# Output:
# DEG result table with statistics and regulation status
# Volcano plot showing significantly up- and down-regulated genes

library(ggplot2)
library(limma)

rm(list = ls())
path <- getwd()
setwd(path)

input1 <- "1_Processed"
output <- "2_DE"
if (!dir.exists(output)) {
  dir.create(output)
}

perform_DE_analysis <- function(expr_file, info_file, group_names, contrast_str, output_prefix) {
  expr <- read.delim2(expr_file, header = TRUE, row.names = 1)
  group <- read.csv(info_file, header = TRUE)
  
  design <- model.matrix(~0 + factor(group$Group))
  colnames(design) <- group_names
  rownames(design) <- colnames(expr)
  
  contrast.matrix <- makeContrasts(contrasts = contrast_str, levels = design)
  
  expr_numeric <- as.data.frame(lapply(expr, as.numeric))
  rownames(expr_numeric) <- rownames(expr)
  colnames(expr_numeric) <- colnames(expr)
  
  fit <- lmFit(expr_numeric, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  DEG <- topTable(fit2, coef = 1, n = Inf)
  DEG$regulate <- ifelse(DEG$adj.P.Val > 0.05, "unchanged",
                         ifelse(DEG$logFC > 1, "up-regulated",
                                ifelse(DEG$logFC < -1, "down-regulated", "unchanged")))
  
  write.csv(DEG, file = paste0(output, "/DEG_", output_prefix, ".csv"))
  
  ## Volcano plot
  p <- ggplot(DEG, aes(x = logFC, y = -log10(adj.P.Val))) + 
    geom_point(alpha = 0.6, size = 3.5, aes(color = regulate)) + 
    ylab("-log10(adj.P.Val)") + 
    scale_color_manual(values = c("#0072B5", "grey", "#BC3C29")) + 
    geom_vline(xintercept = c(-1, 1), lty = 4, col = "black", lwd = 0.8) + 
    geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.8) + 
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14))
  
  ggsave(filename = paste0(output, "/Volcano_", output_prefix, ".png"),
         plot = p, width = 6, height = 6, dpi = 300)
  
  return(table(DEG$regulate))
}


##Data Preparation(from step1)
# Periodontitis
perform_DE_analysis(paste0(input1, "/expr_P1.txt"),
                    paste0(input1, "/Info_P1.csv"),
                    c("Control", "P"),
                    "P - Control",
                    "P1")

# UC1
perform_DE_analysis(paste0(input1, "/expr_UC1.txt"),
                    paste0(input1, "/Info_UC1.csv"),
                    c("Control", "UC"),
                    "UC - Control",
                    "UC1")

# UC2 (3 groups)
perform_DE_analysis(paste0(input1, "/expr_UC2.txt"),
                    paste0(input1, "/Info_UC2.csv"),
                    c("Control", "UC", "CD"),
                    "UC - Control",
                    "UC2")

# UC3
perform_DE_analysis(paste0(input1, "/expr_UC3.txt"),
                    paste0(input1, "/Info_UC3.csv"),
                    c("Control", "UC"),
                    "UC - Control",
                    "UC3")



# Step 3: Multiplex Network Construction----------------------------------------

# Building gene co-expression networks for Periodontitis and UC1,
# aligning network nodes
# Output: aligned edge lists and node mappings.

rm(list = ls())
path <- getwd()
setwd(path)


input1 <- "1_Processed"
input2 <- "2_DE"
output <- "4_Multiplex_Net/P1UC1"

if (!dir.exists(output)) {
  dir.create(output,recursive = TRUE)
}


# Expression Matrices and Group (from step1)
P1_expr <- read.table(file = paste0(input1,"/expr_P1.txt"),
                      sep = "\t",
                      header = T,
                      row.names = 1)
UC1_expr <- read.table(file = paste0(input1,"/expr_UC1.txt"),
                       sep = "\t",
                       header = T,
                       row.names = 1)
group_P1 <- read.csv(file = paste0(input1,"/info_P1.csv"))
group_UC1 <- read.csv(file = paste0(input1,"/info_UC1.csv"))


# DEA Results(from step2)
P1_DEA <- read.csv(file = paste0(input2,"/DEG_P1.csv"),
                   header = T)
colnames(P1_DEA)[1] <- "SYMBOL"
UC1_DEA <- read.csv(file = paste0(input2,"/DEG_UC1.csv"),
                    header = T)
colnames(UC1_DEA)[1] <- "SYMBOL"



# Obtain potential nodes
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

matrix_P1 <- P1_expr[GeneSymOlap, ]
matrix_UC1 <- UC1_expr[GeneSymOlap,]



#Construct Gene Co-expression Networks (via GeneNet)
GeneSymOlap_DF <- data.frame(GeneSymOlap)

P1_expr_affected <- matrix_P1[,group_P1$Group==1]
P1_expr_unaffected <- matrix_P1[,group_P1$Group==0]

UC1_expr_affected <- matrix_UC1[,group_UC1$Group==1]
UC1_expr_unaffected <- matrix_UC1[,group_UC1$Group==0]


## define network construction function     
#### BASED on GeneNet package. 
network_con <- function(temp_expr, ggmcutoff){
  pcor.dyn = ggm.estimate.pcor(t(temp_expr), method = "dynamic")
  arth.edges = network.test.edges(pcor.dyn,direct=FALSE)  ## 471*470/2 edges # 
  # pairwise corr based # 
  dim(arth.edges)
  #extract.network(network.all, method.ggm=c("prob", "qval","number"),
  #                cutoff.ggm=0.8, method.dir=c("prob","qval","number", "all"),
  #                cutoff.dir=0.8, verbose=TRUE)
  arth.net = extract.network(arth.edges, method.ggm="prob", cutoff.ggm=ggmcutoff)
  #arth.net
  el <- cbind(a=arth.net$node1, b=arth.net$node2, c=arth.net$prob)
  node.labels = as.character(row.names(temp_expr))  
  gr = network.make.graph(arth.net, node.labels, drop.singles=TRUE)  # drop single nodes
  
  return(list(out_graph = gr,out_edge = el, out_nodes = node.labels))
}



p_thres <- 0.95


P1_net <- network_con(P1_expr_affected,p_thres)
P1_net$out_graph
UC1_net <- network_con(UC1_expr_affected,p_thres)
UC1_net$out_graph

P1_netCon <- network_con(P1_expr_unaffected,p_thres)
P1_netCon$out_graph
UC1_netCon <- network_con(UC1_expr_unaffected,p_thres)
UC1_netCon$out_graph




#### network construction linking numbers  
P1_net
UC1_net
length(GeneSymOlap)


layerP1 <- P1_net$out_edge[,1:2]
P1_net$out_nodes # name and ids

layerUC1 <- UC1_net$out_edge[,1:2]

UC1_net$out_nodes # name and ids 


###
layerP1Con <- P1_netCon$out_edge[,1:2]
P1_netCon$out_nodes # name and ids

layerUC1Con <- UC1_netCon$out_edge[,1:2]

UC1_netCon$out_nodes # name and ids 


### Mapping ids ###    
layerP1_mapped <- layerP1 

P1_net_node_DF <- data.frame(P1_net$out_nodes)
UC1_net_node_DF <- data.frame(UC1_net$out_nodes)

### use pd as standard mappings 
P1_net_node_mappings <- P1_net_node_DF
layerP1_mapped <- layerP1
###

## ID UC1 <- gene <- ID P1
layerUC1_mapped <- layerUC1
layerUC1Con_mapped <- layerUC1Con


for (i in 1:nrow(layerUC1_mapped)) {
  for (j in 1:ncol(layerUC1_mapped)) {

    num_temp <- layerUC1_mapped[i, j]
    key_temp <- UC1_net$out_nodes[num_temp]
    
    # P1_net$out_nodes==key_temp
    # which(key_temp %in% P1_net$out_nodes)
    index_temp <- which(P1_net$out_nodes == key_temp)
    layerUC1_mapped[i, j] <- index_temp
  }
}



### 
for (i in 1:nrow(layerUC1Con_mapped)) {
  for (j in 1:ncol(layerUC1Con_mapped)) {
    
    num_temp <- layerUC1Con_mapped[i, j]
    key_temp <- UC1_netCon$out_nodes[num_temp]
    
    # P1_net$out_nodes==key_temp
    # which(key_temp %in% P1_net$out_nodes)
    index_temp <- which(P1_net$out_nodes == key_temp)
    layerUC1Con_mapped[i, j] <- index_temp
  }
}



### mapped ids ### for multinet ### 
layerP1_mapped
layerUC1_mapped
#layerUC12_mapped

write.table(layerP1_mapped, 
            file = paste0(output,"/P1_layer.txt"),
            row.names = FALSE,
            col.names = FALSE)
write.table(layerUC1_mapped, 
            file = paste0(output,"/UC1_layer.txt"),
            row.names = FALSE,
            col.names = FALSE)

write.table(layerP1Con, 
            file = paste0(output,"/P1_Conlayer.txt"),
            row.names = FALSE,
            col.names = FALSE)
write.table(layerUC1Con_mapped, 
            file = paste0(output,"/UC1_Conlayer.txt"),
            row.names = FALSE,
            col.names = FALSE)

write.csv(P1_net_node_DF, 
          file = paste0(output,"/id_genename.csv"))

