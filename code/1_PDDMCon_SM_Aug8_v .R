
#############################
### 1.0 Data input and preprocessing 
library('GEOquery')
library('sva')
library("AnnotationDbi")
library("hgu133plus2.db")    ##for Human
library('BatchQC')
library("graph")
library('GeneNet')
library("Rgraphviz")
library('MCL')

###############################
### 2.0 matrix and label extraction ### 
# PD data 
gsetPD <- getGEO("GSE16134", destdir = ".")   ### MicroArray ###
labelPD <- gsetPD$GSE16134_series_matrix.txt.gz$characteristics_ch1 # labels #

# T2DM data
gsetDM1 <- getGEO("GSE50397", destdir = ".")  ###  microarray ###GPL6244	
labelDM1 <- gsetDM1$GSE50397_series_matrix.txt.gz$characteristics_ch1.5  #labels #

##### labels pd
labelPD_N <- c() #120 patients. 
for (i_pd in 1:length(labelPD)) {
  if(substr(toString(labelPD[i_pd]),28,35)=='Diseased')
  {labelPD_N <- c(labelPD_N,1)} # 241
  else
  {labelPD_N <- c(labelPD_N,0)} # 69
}

##### labels dm
labelDM_1 <- c()
for (i_pd in 1:length(labelDM1)) {
  if(labelDM1[i_pd]!='')
  {
    if(as.numeric((substr(toString(labelDM1[i_pd]),8,100000)))>=5.7)  ### threshold 
    {labelDM_1 <- c(labelDM_1,1)} # 38
    else
    {labelDM_1 <- c(labelDM_1,0)} # 39
  }
  else
  {
    labelDM_1 <- c(labelDM_1,2)  # 12 missings 
  }
}

###############################
### 3.0 Expression extraction. 
# data matrix PD
gsetPD_one <- gsetPD[[1]]
exprSetPD <- exprs(gsetPD_one)
dim(exprSetPD)
colnames(exprSetPD)
row.names(exprSetPD)
#var(exprSetPD)
# platform GPL570

# data matrix DM
exprSetDM1 <- exprs(gsetDM1[[1]])
dim(exprSetDM1)
colnames(exprSetDM1)
row.names(exprSetDM1)
# platform GPL6244

## DE analysis based on GEO2R #### 
PD_table <- read.csv("GSE16134.top.table.csv")
DM_table1 <- read.csv("GSE50397.top.table.csv")
PD_DEG_genes <- PD_table[PD_table$P.Value<0.05&abs(PD_table$logFC)>0.2,]
dim(PD_DEG_genes)

library(ggplot2)
DEG <- PD_table
DEG$regulate <- ifelse(DEG$P.Value > 0.05, "unchanged",
                       ifelse(DEG$logFC > 0.2, "up-regulated",
                              ifelse(DEG$logFC < -0.2, "down-regulated", "unchanged")))
table(DEG$regulate)
pdf("volcano1.pdf")
ggplot(DEG,aes(x=logFC,y=-log10(P.Value)))+ 
  geom_point(alpha=0.6,size=3.5,aes(color=regulate))+ 
  ylab("-log10(P.Value)")+ #y轴的说明
  scale_color_manual(values = c("blue", "grey", "red"))+ 
  geom_vline(xintercept = c(-0.2,0.2),lty=4,col ="black",lwd=0.8)+ 
  geom_hline(yintercept=-log10(0.05),lty=4,col = "black",lwd=0.8)+ 
  theme_bw()  #火山图绘制
dev.off()

#volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=23,
#            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)) ,hl.col = 'cyan3')

PD_table_listfull <- PD_table$Gene.symbol
DM_DEG_genes <- DM_table1[DM_table1$P.Value<0.05&abs(DM_table1$logFC)>0.2,]
#DM_DEG_genes <- DM_table2[DM_table2$adj.P.Val<0.05&abs(DM_table2$logFC)>1,]
dim(DM_DEG_genes)

DEG <- DM_table1
DEG$regulate <- ifelse(DEG$P.Value > 0.05, "unchanged",
                       ifelse(DEG$logFC > 0.2, "up-regulated",
                              ifelse(DEG$logFC < -0.2, "down-regulated", "unchanged")))
table(DEG$regulate)
#max(-log(DEG$adj.P.Val))
pdf("volcano2.pdf")
ggplot(DEG,aes(x=logFC,y=-log10(P.Value)))+ 
  geom_point(alpha=0.6,size=3.5,aes(color=regulate))+
  ylim(-0.2, 4)+
  ylab("-log10(P.Value)")+ #
  scale_color_manual(values = c("blue", "grey", "red"))+ 
  geom_vline(xintercept = c(-0.2,0.2),lty=4,col ="black",lwd=0.8)+ 
  geom_hline(yintercept=-log10(0.05),lty=4,col = "black",lwd=0.8)+ 
  theme_bw()  #
dev.off()
#dm_list3 <- intersect(DM_table1_list3,DM_table2_list3)
#dm_list3

GeneSymOlap <- intersect(PD_DEG_genes$Gene.symbol, DM_DEG_genes$Gene.symbol)   ##### intersection of gene list #### 
length(GeneSymOlap)
#GeneSymOlap <- intersect(PD_table_list_p1, DM_table1_list3)  ### network went to large ### 
#GeneSymOlapfull <- intersect(PD_table_listfull, DM_table1_list3)  #### might need to optimize in future. #### // overlap. 

GeneSymOlap <- GeneSymOlap[(GeneSymOlap) %in% "" == FALSE] ####
length(GeneSymOlap)
# GeneSymOlapfull <- GeneSymOlapfull[(GeneSymOlapfull) %in% "" == FALSE] 
# my_list[names(my_list) %in% "number" == FALSE] 

write.csv(GeneSymOlap, file = 'PDDM_commonDEgenes.csv',
          quote=FALSE,col.names = FALSE,row.names=FALSE,
          fileEncoding = "UTF-8")

###############################
### 4.0 dimension reduction based on DE genes 
# PD_table, most of pd genes are DE.
# DM

### gene id and expression extraction 
# have ID overlappings ### expression with high prob ID# 
PD_table_id_list <- PD_table[unlist(lapply(PD_table$Gene.symbol, function(x) x %in% GeneSymOlap)),]$ID
PD_table_Genesymbol_list <- PD_table[unlist(lapply(PD_table$Gene.symbol, function(x) x %in% GeneSymOlap)),]$Gene.symbol
DM_table_id_list <- DM_table1[unlist(lapply(DM_table1$Gene.symbol, function(x) x %in% GeneSymOlap)),]$ID
DM_table_Genesymbol_list <- DM_table1[unlist(lapply(DM_table1$Gene.symbol, function(x) x %in% GeneSymOlap)),]$Gene.symbol

#length(intersect(DM_table_Genesymbol_list,DM_table_Genesymbol_list))
#### remove repeatation #### ##### 

#write.csv(PD_table_Genesymbol_list, '1.csv')
#write.csv(DM_table_Genesymbol_list, '2.csv')
## expression extraction ## 
## overlap gene - cal the mean # for PD. ### 

PD_expr_extr <- exprSetPD[PD_table_id_list,]
dim(PD_expr_extr)

# PD_expr_extr_Single 
# GeneSymOlap # for overlapping probeIDs select highest expressed probeID. #############################

select_HighProbe <- function(expr_matrix, id_mapping, genelist){  ### selected high probe to remove redundancies in gene IDs. ### 
  PD_expr_extr_Single_values <- c()
  PD_table <- id_mapping
  exprSetPD <- expr_matrix
  for (GeneCur in genelist) {
    print(GeneCur)
    PD_table_id_list_tmp <- PD_table[unlist(lapply(PD_table$Gene.symbol, function(x) x %in% GeneCur)),]$ID
    #print(PD_table_id_list_tmp)
    PD_expr_extr_tmp <- exprSetPD[PD_table_id_list_tmp,]
    
    #PD_expr_extr_tmp[which.max(rowMeans(PD_expr_extr_tmp)),]
    #PD_expr_extr_Single[PD_table_id_list_tmp,] <- PD_expr_extr_tmp[which.max(rowMeans(PD_expr_extr_tmp)),]
    #PD_expr_extr_Single[nrow(PD_expr_extr_Single) + 1,] <- as.numeric(PD_expr_extr_tmp[which.max(rowMeans(PD_expr_extr_tmp)),])
    
    if(length(PD_expr_extr_tmp)!=dim(exprSetPD)[2])
    {
      PD_expr_extr_Single_values <- rbind(PD_expr_extr_Single_values,as.numeric(PD_expr_extr_tmp[which.max(rowMeans(PD_expr_extr_tmp)),]))
    }
    else
    {
      PD_expr_extr_Single_values <- rbind(PD_expr_extr_Single_values,as.numeric(PD_expr_extr_tmp)) 
    }
  }
  
  return(PD_expr_extr_Single_values)
}
select_HighProbeDM <- function(expr_matrix, id_mapping, genelist){  ### selected high probe to remove redundancies in gene IDs. ### 
  PD_expr_extr_Single_values <- c()
  PD_table <- id_mapping
  exprSetPD <- expr_matrix
  for (GeneCur in genelist) {
    print(GeneCur)
    PD_table_id_list_tmp <- PD_table[unlist(lapply(PD_table$Gene.symbol, function(x) x %in% GeneCur)),]$ID
    #print(PD_table_id_list_tmp)
    PD_expr_extr_tmp <- exprSetPD[as.character(PD_table_id_list_tmp),]
    #print(PD_expr_extr_tmp)
    #PD_expr_extr_tmp[which.max(rowMeans(PD_expr_extr_tmp)),]
    #PD_expr_extr_Single[PD_table_id_list_tmp,] <- PD_expr_extr_tmp[which.max(rowMeans(PD_expr_extr_tmp)),]
    #PD_expr_extr_Single[nrow(PD_expr_extr_Single) + 1,] <- as.numeric(PD_expr_extr_tmp[which.max(rowMeans(PD_expr_extr_tmp)),])
    
    if(length(PD_expr_extr_tmp)!=dim(exprSetPD)[2])
    {
      PD_expr_extr_Single_values <- rbind(PD_expr_extr_Single_values,as.numeric(PD_expr_extr_tmp[which.max(rowMeans(PD_expr_extr_tmp)),]))
    }
    else
    {
      PD_expr_extr_Single_values <- rbind(PD_expr_extr_Single_values,as.numeric(PD_expr_extr_tmp)) 
    }
  }
  
  return(PD_expr_extr_Single_values)
}

PD_expr_extr_Single_values <- select_HighProbe(exprSetPD,PD_table,GeneSymOlap)
PD_expr_extr_Single <- data.frame(PD_expr_extr_Single_values, row.names = GeneSymOlap)
colnames(PD_expr_extr_Single)=colnames(exprSetPD)
dim(PD_expr_extr_Single)
write.csv(PD_expr_extr_Single, file = 'PD_expr_extr_Single.csv',
          quote=FALSE,
          fileEncoding = "UTF-8")

DM_expr_extr_Single_values <- select_HighProbeDM(exprSetDM1,DM_table1,GeneSymOlap)
DM_expr_extr_Single <- data.frame(DM_expr_extr_Single_values, row.names = GeneSymOlap)
colnames(DM_expr_extr_Single)=colnames(exprSetDM1)
dim(DM_expr_extr_Single)
write.csv(DM_expr_extr_Single, file = 'DM_expr_extr_Single.csv',
          quote=FALSE,
          fileEncoding = "UTF-8")

###############################
### 5.0 ID convertion and network construction ### PD, DM, Con ### 
# PD_expr_extr_Single 
# DM_expr_extr1
# DM_expr_extr2
GeneSymOlap_DF <- data.frame(GeneSymOlap)
write.csv(labelPD_N, file = 'label_pd.csv')
PD_expr_affected <- PD_expr_extr_Single[,labelPD_N==1]
PD_expr_unaffected <- PD_expr_extr_Single[,labelPD_N==0]

#DM_expr_extr1
write.csv(labelDM_1, file = 'labelDM_1.csv')
DM_expr_1 <- DM_expr_extr_Single[,labelDM_1==1]
DM_expr_con1 <- DM_expr_extr_Single[,!labelDM_1==0]

## define network construction function ##     #### BASED on GeneNet package. 
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
  node.labels = as.character(row.names(temp_expr))  ### node labels #### 
  gr = network.make.graph(arth.net, node.labels, drop.singles=TRUE)  # drop single nodes
  
  return(list(out_graph = gr,out_edge = el, out_nodes = node.labels))
}

p_thres <- 0.8 
PD_net <- network_con(PD_expr_affected,p_thres)
PD_net$out_graph
DM_net1 <- network_con(DM_expr_1,p_thres)
DM_net1$out_graph

PD_netCon <- network_con(PD_expr_unaffected,p_thres)
PD_netCon$out_graph
DM_netCon1 <- network_con(DM_expr_con1,p_thres)
DM_netCon1$out_graph

## define network visualization function ## 
net_visualization <- function(graph_temp){
  graph_temp <- graph_temp
  globalAttrs = list()
  globalAttrs$edge = list(color = "black", lty = "solid", lwd = 1, arrowsize=1)
  globalAttrs$node = list(fillcolor = gray(.95), shape = "ellipse", fixedsize = FALSE)
  #Set attributes of some particular nodes:
  nodeAttrs = list()
  #nodeAttrs$fillcolor = c('570' = "red", "81" = "red") # highlight hub nodes
  nodeAttrs$fillcolor = c('11000' = "red", "81" = "red") # highlight hub nodes
  #Set edge attributes:
  edi = edge.info(graph_temp) # edge directions and correlations
  edgeAttrs = list()
  edgeAttrs$dir = edi$dir # set edge directions
  cutoff = quantile(abs(edi$weight), c(0.2, 0.8)) # thresholds for line width / coloring
  edgeAttrs$lty = ifelse(edi$weight < 0, "dotted", "solid") # negative correlation
  #edgeAttrs$color = ifelse( abs(edi$weight <= cutoff[1]), "grey", "black") # lower 20% quantile
  edgeAttrs$color = ifelse( abs(edi$weight <= cutoff[1]), "blue", "red") # lower 20% quantile
  edgeAttrs$lwd = ifelse(abs(edi$weight >= cutoff[2]), 2, 1) # upper 20% quantile
  plot(graph_temp, attrs = globalAttrs, nodeAttrs = nodeAttrs, edgeAttrs = edgeAttrs, "fdp")
  
}

net_visualization(PD_net$out_graph)
net_visualization(DM_net1$out_graph)

net_visualization(PD_netCon$out_graph)
net_visualization(DM_netCon1$out_graph)

#### network construction linking numbers #### 
PD_net
DM_net1
#DM_net2

###############################
### 6.0 api for multinet and cytoscape
layerPD <- PD_net$out_edge[,1:2]
PD_net$out_nodes # name and ids

layerDM1 <- DM_net1$out_edge[,1:2]
#layerDM2 <- DM_net2$out_edge[,1:2]
DM_net1$out_nodes # name and ids 
#DM_net2$out_nodes

###
layerPDCon <- PD_netCon$out_edge[,1:2]
PD_netCon$out_nodes # name and ids

layerDMCon1 <- DM_netCon1$out_edge[,1:2]
#layerDMCon2 <- DM_netCon2$out_edge[,1:2]
DM_netCon1$out_nodes # name and ids 
#DM_netCon2$out_nodes

### Mapping ids ###    
layerPD_mapped <- layerPD 
#layerDM1_mapped
#layerDM2_mapped
PD_net_node_DF <- data.frame(PD_net$out_nodes)
DM_net1_node_DF <- data.frame(DM_net1$out_nodes)
#PD_netCon_node_DF <- data.frame(PD_netCon$out_nodes)
#DM_netCon1_node_DF <- data.frame(DM_netCon1$out_nodes)

### use pd as standard mappings 
PD_net_node_mappings <- PD_net_node_DF
layerPD_mapped <- layerPD
###

## ID DM <- gene <- ID PD
layerDM1_mapped <- layerDM1
#layerDM2_mapped <- layerDM2
layerDMCon1_mapped <- layerDMCon1
#layerDMCon2_mapped <- layerDMCon2

for (i in 1:nrow(layerDM1_mapped)) {
  for (j in 1:ncol(layerDM1_mapped)) {
    num_temp <- layerDM1_mapped[i, j]
    key_temp <- DM_net1$out_nodes[num_temp]
    
    # PD_net$out_nodes==key_temp
    # which(key_temp %in% PD_net$out_nodes)
    index_temp <- which(PD_net$out_nodes == key_temp)
    layerDM1_mapped[i, j] <- index_temp
  }
}

### 
for (i in 1:nrow(layerDMCon1_mapped)) {
  for (j in 1:ncol(layerDMCon1_mapped)) {
    #
    num_temp <- layerDMCon1_mapped[i, j]
    key_temp <- DM_netCon1$out_nodes[num_temp]
    
    # PD_net$out_nodes==key_temp
    # which(key_temp %in% PD_net$out_nodes)
    index_temp <- which(PD_net$out_nodes == key_temp)
    layerDMCon1_mapped[i, j] <- index_temp
  }
}

### mapped ids ### for multinet ### 
layerPD_mapped
layerDM1_mapped
#layerDM2_mapped

write.table(layerPD_mapped, file = 'pd_layer.txt',row.names = FALSE,
            col.names = FALSE)
write.table(layerDM1_mapped, file = 'dm1_layer.txt',row.names = FALSE,
            col.names = FALSE)
#write.table(layerDM2_mapped, file = 'dm2_layer.txt',row.names = FALSE, col.names = FALSE)

write.table(layerPDCon, file = 'pd_Conlayer.txt',row.names = FALSE,
            col.names = FALSE)
write.table(layerDMCon1_mapped, file = 'dm1_Conlayer.txt',row.names = FALSE,
            col.names = FALSE)
#write.table(layerDMCon2_mapped, file = 'dm2_Conlayer.txt',row.names = FALSE, col.names = FALSE)

write.csv(PD_net_node_DF, file = 'id_genename.csv')

### mapped ids ### for Cytoscape ### ####
# mapped genes. ###

layerPD_mappedGenes <- apply(layerPD_mapped, c(1,2), function(x) {PD_net_node_DF[x,]})
layerDM1_mappedGenes <- apply(layerDM1_mapped, c(1,2), function(x) {PD_net_node_DF[x,]})
#layerDM2_mappedGenes <- apply(layerDM2_mapped, c(1,2), function(x) {PD_net_node_DF[x,]})

layerPDCon_mappedGenes <- apply(layerPDCon, c(1,2), function(x) {PD_net_node_DF[x,]})
layerDM1Con_mappedGenes <- apply(layerDMCon1_mapped, c(1,2), function(x) {PD_net_node_DF[x,]})
#layerDM2Con_mappedGenes <- apply(layerDMCon2_mapped, c(1,2), function(x) {PD_net_node_DF[x,]})

write.table(layerPD_mappedGenes, file = 'pd_layer_genes.txt',row.names = FALSE,
            col.names = TRUE)
write.table(layerDM1_mappedGenes, file = 'dm1_layer_genes.txt',row.names = FALSE,
            col.names = TRUE)

write.table(layerPDCon_mappedGenes, file = 'pd_Conlayer_genes.txt',row.names = FALSE,
            col.names = TRUE)
write.table(layerDM1Con_mappedGenes, file = 'dm1_Conlayer_genes.txt',row.names = FALSE,
            col.names = TRUE)

mat  ### sparse matrix. 
