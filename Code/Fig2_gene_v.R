
library(dplyr)
library(GeneNet)
library(Rgraphviz)
library(MCL)
library(ggplot2)
library(linkET)
library(igraph)

library(tidyr)
library(reshape2)
library(ggrepel)
library(pheatmap)
library(ggpubr)

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggsci)
library(stringr)
library(extrafont)
library(DOSE)

library(data.table)
library('enrichplot')

################################################################################
### 0.0 ####Env network 
rm(list = ls())
path <- getwd()
setwd(path)
input <- "input"
output <- "output"
if (!dir.exists(output)) {
  dir.create(output,recursive = TRUE)
}

genes_of_interest <- c("PTGR1","FOSB","CXCL5")
algorithm_genes <- read.csv(file = paste0(input,"/CoreHD-DCI_genes.csv"), header = TRUE)
vital_genes <- algorithm_genes$x

id_name <- read.csv(file = paste0(input,"/id_genename.csv"), header = TRUE)
P1_layer <- read.table(file = paste0(input,"/P1_layer.txt"), header = FALSE)
UC1_layer <- read.table(file = paste0(input,"/UC1_layer.txt"), header = FALSE)
P1_Conlayer <- read.table(file = paste0(input,"/P1_Conlayer.txt"), header = FALSE)
UC1_Conlayer <- read.table(file = paste0(input,"/UC1_Conlayer.txt"), header = FALSE)

## Link genes
## network processing
process_network <- function(network_data, gene_id) {
  connected_nodes <- network_data[network_data$V1 == gene_id | network_data$V2 == gene_id, ]
  connected_genes <- unique(c(connected_nodes$V1, connected_nodes$V2))
  connected_genes <- connected_genes[connected_genes != gene_id]
  connected_gene_names <- id_name$P1_net.out_nodes[id_name$X %in% connected_genes]
  return(connected_gene_names)
}

# network plot function 
create_network_plot <- function(vital_gene, connected_gene_names, network_name) {
  if (length(connected_gene_names) == 0) {
    message(paste("Gene", vital_gene, "has no connections in", network_name, "- Skipping plot."))
    return()
  }
  
  edges <- cbind(vital_gene, connected_gene_names)
  graph <- graph_from_edgelist(as.matrix(edges), directed = FALSE)
  
  png(paste0(output, "/",vital_gene, "_", network_name, "_network.png"), width = 800, height = 600)
  plot(graph, 
       vertex.size = 30, 
       vertex.label.cex = 1.5, 
       vertex.label.color = "black", 
       vertex.color = c("#E64B35", rep("#4DBBD5", length(connected_gene_names))),
       main = paste(vital_gene, "Network in", network_name))
  dev.off()
}

# gene list saver
save_gene_list <- function(layer_list, file_name) {
  unique_genes <- unique(unlist(layer_list))
  write.csv(unique_genes, file = paste0(output,"/",file_name, ".csv"), row.names = FALSE)
}

# function for processing layers 
process_layer <- function(layer_data, layer_name) {
  gene_list_layer <- list()
  
  for (gene in vital_genes) {
    
    gene_id <- id_name$X[id_name$P1_net.out_nodes == gene]
    if (length(gene_id) == 0) {
      message(paste("Gene", gene, "not found in id_name. Skipping."))
      next
    }
    
    connected_gene_names <- process_network(layer_data, gene_id)
    gene_list_layer[[gene]] <- connected_gene_names
    
    write.csv(connected_gene_names, file = paste0(output, "/", gene, "_", layer_name, ".csv"), row.names = FALSE)
    
    create_network_plot(gene, connected_gene_names, layer_name)
  }
  return(gene_list_layer)
}

# save layers
P1_layer_genes <- process_layer(P1_layer, "P1_layer")
UC1_layer_genes <- process_layer(UC1_layer, "UC1_layer")

# save gene list
save_gene_list(P1_layer_genes, "P1_layer_gene_list")
save_gene_list(UC1_layer_genes, "UC1_layer_gene_list")

# network visualization
plot_vital_genes_layer <- function(vital_genes, gene_list_layer, layer_name) {
  edges <- c()
  all_genes <- unique(vital_genes)
  
  for (gene in vital_genes) {
    connected_genes <- gene_list_layer[[gene]]
    if (length(connected_genes) > 0) {
      gene_edges <- cbind(rep(gene, length(connected_genes)), connected_genes)
      edges <- rbind(edges, gene_edges)
      all_genes <- unique(c(all_genes, connected_genes))
    }
  }
  
  if (nrow(edges) > 0) {
    graph <- graph_from_edgelist(as.matrix(edges), directed = FALSE)
    vertex_colors <- ifelse(V(graph)$name %in% vital_genes, "#E64B35","#4DBBD5")
    
    layout <- layout_with_fr(graph)
    
    png(paste0(output,"/", layer_name, "_network.png"), width = 1000, height = 800)
    plot(graph, 
         layout = layout,              
         vertex.size = 15,             
         vertex.label.cex = 0.8,       
         vertex.color = vertex_colors, 
         vertex.label.color = "black", 
         edge.color = "gray", 
         edge.width = 1.5,             
         main = paste("Vital Genes Network -", layer_name))
    dev.off()
  } else {
    message(paste("No connections found in layer", layer_name, ". Skipping plot."))
  }
}

plot_vital_genes_layer(vital_genes, P1_layer_genes, "P1_layer")
plot_vital_genes_layer(vital_genes, UC1_layer_genes, "UC1_layer")

# all_gene_list
save_all_gene_list <- function(P1_layer_genes, UC1_layer_genes) {
  P1_genes <- unique(unlist(P1_layer_genes))
  UC1_genes <- unique(unlist(UC1_layer_genes))
  
  # fine the longest gene list
  max_length <- max(length(P1_genes), length(UC1_genes))
  
  # padding with NA
  P1_genes <- c(P1_genes, rep(NA, max_length - length(P1_genes)))
  UC1_genes <- c(UC1_genes, rep(NA, max_length - length(UC1_genes)))

  # create dataframe
  all_gene_list <- data.frame(
    P1_layer = P1_genes,
    UC1_layer = UC1_genes
  )
  # save all_gene_list as csv file
  write.csv(all_gene_list, 
            file = paste0(output,"/all_gene_list.csv"), row.names = FALSE)
}
# save gene list. 
save_all_gene_list(P1_layer_genes, UC1_layer_genes)

################################################################################
### 1.0 expression 

# P1_layer_genes
# UC1_layer_genes
UC1_expr <- read.table(file = paste0(input,"/expr_UC1.txt"),
                       sep = "\t",
                       header = T,
                       row.names = 1)
group_UC1 <- read.csv(file = paste0(input,"/info_UC1.csv"))

PTGR1_Link <- UC1_layer_genes$PTGR1
CXCL5_Link <- UC1_layer_genes$CXCL5
FOSB_Link <- UC1_layer_genes$FOSB

PTGR1_Link_expr <- UC1_expr[PTGR1_Link,]
CXCL5_Link_expr <- UC1_expr[CXCL5_Link,]
FOSB_Link_expr <- UC1_expr[FOSB_Link,]

PTGR1_expr <- UC1_expr["PTGR1",]
CXCL5_expr <- UC1_expr["CXCL5",]
FOSB_expr <- UC1_expr["FOSB",]

DEG_UC <- read.csv(file = paste0(input,"/DEG_UC1.csv"))

P1_expr <- read.table(file = paste0(input,"/expr_P1.txt"),
                      sep = "\t",
                      header = T,
                      row.names = 1)
group_P1 <- read.csv(file = paste0(input,"/info_P1.csv"))

P1PTGR1_Link <- P1_layer_genes$PTGR1
P1CXCL5_Link <- P1_layer_genes$CXCL5
P1FOSB_Link <- P1_layer_genes$FOSB

P1PTGR1_Link_expr <- P1_expr[P1PTGR1_Link,]
P1CXCL5_Link_expr <- P1_expr[P1CXCL5_Link,]
P1FOSB_Link_expr <- P1_expr[P1FOSB_Link,]

P1PTGR1_expr <- P1_expr["PTGR1",]
P1CXCL5_expr <- P1_expr["CXCL5",]
P1FOSB_expr <- P1_expr["FOSB",]

DEG_P1 <- read.csv(file = paste0(input,"/DEG_P1.csv"))

################################################################################
### 2.0 violin plot ### 

##表达量可视化
plot_boxplot <- function(data, group_info, title, output, widthB, heightB) {
  
  data_name <- deparse(substitute(data))
  data$gene <- rownames(data)
  
  melted_data <- melt(data, id.vars = "gene", variable.name = "sample", value.name = "expression")
  
  group_info_filtered <- group_info[group_info$Group %in% c(0, 1), ]  # 0 = Normal, 1 = Disease
  melted_data$Group <- group_info_filtered$Group[match(melted_data$sample, group_info_filtered$Samples)]
  
  melted_data <- melted_data[!is.na(melted_data$Group), ]
  melted_data$Group <- factor(melted_data$Group, levels = c(0, 1), labels = c("Normal", title))
  
  p <- ggplot(melted_data, aes(x = Group, y = expression, fill = Group)) +
    geom_violin(trim = FALSE, size = 0.8, alpha = 0.75, adjust = 1, scale = "area")+  ### estimate by density. 
    #geom_boxplot(trim = FALSE, size = 0.8, alpha = 0.75)+
    geom_jitter(position = position_jitter(width = 0.2), size = 0.5, alpha = 0.6) +
    geom_boxplot(outlier.shape = NA,width = 0.2,notch = FALSE, fill = "white") +
    #geom_boxplot(outlier.shape = NA) +  # 移除箱线图的离群点
    #geom_jitter(width = 0.2, shape = 16, color = "black", alpha = 0.5) +  # 标出数据点
    
    facet_wrap(~ gene, scales = "free_y") +
    theme_bw() +
    labs(title = title, x = "Group", y = "Expression") +
    #theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.y = unit(0, "lines"),
      #plot.title = element_blank(),
      strip.text = element_text(size = 12, margin = margin(0, 0, 2, 0)),
      strip.background = element_blank(),
      legend.position = "bottom",  
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(color = "black", size = 12),
      axis.ticks.y = element_line(color = "black", size = 0.25),
      axis.ticks.length.y = unit(0.05, "cm"),
      axis.title.y = element_text(color = "black", size = 12),
      axis.title.x = element_blank(),
      plot.margin = margin(1, 1, 1, 1)
    ) +
    #scale_fill_manual(values = c("#00A087FF", "#800080")) 
    # "#53a8e1", "#53e8e1", "red"
    scale_fill_manual(values = c("#53e8e1", "red")) 
  p <- p + stat_compare_means(aes(Group = Group),
                              method = "t.test",
                              label = "p.signif",
                              label.x = 1.4)
  print(data_name)
  ggsave(filename = paste0(output,"/", data_name, ".pdf"), plot = p,width = widthB, height = heightB)
  print(p)
}

################################################################################
### 2.1 boxplot in group

UC_expression <- UC1_expr[c("PTGR1","FOSB","CXCL5"),]
Periodontitis_expression <- P1_expr[c("PTGR1","FOSB","CXCL5"),]

plot_boxplot(UC_expression, group_UC1, "UC", output,5,3)           
plot_boxplot(Periodontitis_expression, group_P1, "Periodontitis", output,5,3)

######
################################################################################
### 3.0 volcano
#volcano()
#EnhancedVolcano()

uc1 <- ggplot(DEG_UC,aes(x=logFC,y=-log10(adj.P.Val)))+ 
  geom_point(alpha=0.6,size=3.5,aes(color=regulate))+ 
  ylab("-log10(adj.P.Val)")+ 
  scale_color_manual(values = c("#0072B5", "grey", "#BC3C29"))+ 
  geom_vline(xintercept = c(-1,1),lty=4,col ="black",lwd=0.8)+ 
  geom_hline(yintercept=-log10(0.05),lty=4,col = "black",lwd=0.8)+ 
  theme_bw() +
  theme(
    #text = element_text(family = "Arial", size = 14), 
    axis.title = element_text(size = 16),  
    legend.title = element_text(size = 16), 
    legend.text = element_text(size = 14)   
  ) 
print(uc1)
ggsave(filename = paste0(output,"/Volcano_UC1.pdf"), 
       plot = uc1, width = 6, height = 6)

p1 <- ggplot(DEG_P1,aes(x=logFC,y=-log10(adj.P.Val)))+ 
  geom_point(alpha=0.6,size=3.5,aes(color=regulate))+ 
  ylab("-log10(adj.P.Val)")+ 
  scale_color_manual(values = c("#0072B5", "grey", "#BC3C29"))+ 
  geom_vline(xintercept = c(-1,1),lty=4,col ="black",lwd=0.8)+ 
  geom_hline(yintercept=-log10(0.05),lty=4,col = "black",lwd=0.8)+ 
  theme_bw() +
  theme(
    #text = element_text(family = "Arial", size = 14), 
    axis.title = element_text(size = 16),  
    legend.title = element_text(size = 16), 
    legend.text = element_text(size = 14)   
  ) 
print(p1)
ggsave(filename = paste0(output,"/Volcano_p1.pdf"),
       plot = p1, width = 6, height = 6)

############~~~~~ version2 
uc1 <- ggplot(DEG_UC,aes(x=logFC,y=-log10(adj.P.Val)))+ 
  geom_point(alpha=0.6,size=3.5,aes(color=regulate))+ 
  ylab("-log10(adj.P.Val)")+ 
  #geom_text(data = top_genes, aes(label = gens),size = 3, vjust = -1)+
  #scale_color_manual(values = c("#0072B5", "grey", "#BC3C29"))+ 
  scale_color_manual(values = c("#BC3C29", "#3b6fd9", "#BC3C29"))+ 
  #scale_color_gradient(low = "blue", high = "red") + 
  geom_vline(xintercept = c(-1,1),lty=4,col ="black",lwd=0.8)+ 
  geom_hline(yintercept=-log10(0.05),lty=4,col = "black",lwd=0.8)+ 
  theme_minimal() +
  theme(
    #text = element_text(family = "Arial", size = 14), 
    axis.title = element_text(size = 16),  
    legend.title = element_text(size = 16), 
    legend.text = element_text(size = 14),
    panel.border = element_blank(),  # 去掉边框
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank(),  # 去掉次网格线
    axis.line = element_line(color = "black")  # 添加x和y轴
  ) 
print(uc1)
##"#ffa74f",
genes_of_interestUC <- subset(DEG_UC, DEG_UC$X =="CXCL5"|DEG_UC$X =="FOSB"|DEG_UC$X =="PTGR1") 
uc1 <- uc1 + geom_text(data = genes_of_interestUC, aes(label = X), color = "black", alpha=1, size = 5, , vjust = 1.5, hjust = - 0.0001)+
       geom_point(data = genes_of_interestUC, alpha = 0.7, size = 3.5, color = "black", shape = 19)
print(uc1)
## 
ggsave(filename = paste0(output,"/Volcano_UC1_2.pdf"), 
       plot = uc1, width = 6, height = 6)

p1 <- ggplot(DEG_P1,aes(x=logFC,y=-log10(adj.P.Val)))+ 
  geom_point(alpha=0.6,size=3.5,aes(color=regulate))+ 
  ylab("-log10(adj.P.Val)")+ 
  #scale_color_manual(values = c("#0072B5", "grey", "#BC3C29"))+ 
  scale_color_manual(values = c("#BC3C29", "#3b6fd9", "#BC3C29"))+ 
  geom_vline(xintercept = c(-1,1),lty=4,col ="black",lwd=0.8)+ 
  geom_hline(yintercept=-log10(0.05),lty=4,col = "black",lwd=0.8)+ 
  theme_bw() +
  theme(
    #text = element_text(family = "Arial", size = 14), 
    axis.title = element_text(size = 16),  
    legend.title = element_text(size = 16), 
    legend.text = element_text(size = 14),
    panel.border = element_blank(),  # 去掉边框
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank(),  # 去掉次网格线
    axis.line = element_line(color = "black")  # 添加x和y轴
  ) 
#print(p1)
##
genes_of_interestP <- subset(DEG_P1, DEG_P1$X =="CXCL5"|DEG_P1$X =="FOSB"|DEG_P1$X =="PTGR1") 
p1 <- p1 + geom_text(data = genes_of_interestP, aes(label = X), alpha=1, size = 5, , vjust = -1, hjust = - 0.0001)+
      geom_point(data = genes_of_interestP, alpha = 0.7, size = 3.5, color = "black", shape = 19)
print(p1)

ggsave(filename = paste0(output,"/Volcano_p1_2.pdf"),
       plot = p1, width = 6, height = 6)

#labs(title = "UC", x = bquote(~Log[2]~"(FC)"), y = bquote("-"~Log[10]~"(adj.P)")) +
########################~~~~~ version3  current in use
uc1 <- ggplot(DEG_UC,aes(x=logFC,y=-log10(adj.P.Val)))+ 
  geom_point(alpha=0.8,size=5,aes(color=-log10(adj.P.Val)))+ 
  labs(title = "UC", x = bquote(~Log[2]~"(FC)"), y = bquote("-"~Log[10]~"(adj.P)"))+ 
  #geom_text(data = top_genes, aes(label = gens),size = 3, vjust = -1)+
  #scale_color_manual(values = c("#0072B5", "grey", "#BC3C29"))+ 
  #scale_color_manual(values = c("#BC3C29", "#3b6fd9", "#BC3C29"))+ 
  #scale_color_gradient2(low = "blue", high = "red", midpoint = 15) + 
  #scale_color_gradientn(colors = c("#53a8e1", "red"), values = c(0, 0.65, 1)) + 
  scale_color_gradientn(colors = c("#53a8e1", "#53e8e1", "red"), values = c(0, 0.23, 1)) + 
  #geom_vline(xintercept = c(-1,1),lty=4,col ="black",lwd=0.8)+ 
  geom_hline(yintercept=-log10(0.05),lty=4,col = "black",lwd=0.8)+ 
  theme_minimal() +
  theme(
    #text = element_text(family = "Times New Roman", size = 14), 
    axis.title = element_text(size = 16),  
    legend.title = element_text(size = 16), 
    legend.text = element_text(size = 14),
    panel.border = element_blank(),  # 去掉边框
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank(),  # 去掉次网格线
    axis.line = element_line(color = "black"),  # 添加x和y轴
    plot.title = element_text(hjust = 0.5, size =  16 )
    ) 
print(uc1)
##"#ffa74f",
genes_of_interestUC <- subset(DEG_UC, DEG_UC$X =="CXCL5"|DEG_UC$X =="FOSB"|DEG_UC$X =="PTGR1") 
uc1 <- uc1 + geom_text(data = genes_of_interestUC, aes(label = X), color = "black", alpha=1, size = 5, , vjust = 1.5, hjust = - 0.0001)+
  geom_point(data = genes_of_interestUC, alpha = 0.7, size = 3.5, color = "black", shape = 19)+ggtitle("UC")
print(uc1)
## 
ggsave(filename = paste0(output,"/Volcano_UC1_3.pdf"), 
       plot = uc1, width = 8, height = 6)

p1 <- ggplot(DEG_P1,aes(x=logFC,y=-log10(adj.P.Val)))+ 
  geom_point(alpha=0.6,size=3.5,aes(color=-log10(adj.P.Val)))+ 
  labs(title = "Periodontitis", x = bquote(~Log[2]~"(FC)"), y = bquote("-"~Log[10]~"(adj.P)"))+ 
  #scale_color_manual(values = c("#0072B5", "grey", "#BC3C29"))+ 
  #scale_color_manual(values = c("#BC3C29", "#3b6fd9", "#BC3C29"))+ 
  scale_color_gradientn(colors = c("#53a8e1", "#53e8e1", "red"), values = c(0, 0.23, 1)) + 
  #geom_vline(xintercept = c(-1,1),lty=4,col ="black",lwd=0.8)+ 
  geom_hline(yintercept=-log10(0.05),lty=4,col = "black",lwd=0.8)+ 
  theme_bw() +
  theme(
    #text = element_text(family = "Arial", size = 14), 
    axis.title = element_text(size = 16),  
    legend.title = element_text(size = 16), 
    legend.text = element_text(size = 14),
    panel.border = element_blank(),  # 去掉边框
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank(),  # 去掉次网格线
    axis.line = element_line(color = "black"),  # 添加x和y轴
    plot.title = element_text(hjust = 0.5, size =  16)
  ) 
print(p1)
##
genes_of_interestP <- subset(DEG_P1, DEG_P1$X =="CXCL5"|DEG_P1$X =="FOSB"|DEG_P1$X =="PTGR1") 
p1 <- p1 + geom_text(data = genes_of_interestP, aes(label = X), alpha=1, size = 5, , vjust = -1, hjust = - 0.0001)+
  geom_point(data = genes_of_interestP, alpha = 0.7, size = 3.5, color = "black", shape = 19)+ggtitle("Periodontitis")
print(p1)

ggsave(filename = paste0(output,"/Volcano_p1_3.pdf"),
       plot = p1, width = 8, height = 6)

################################################################################
### 4.0 heatmap of key-link genes in UC. 

#pheatmap()  ConsensusClusterPlus survfit
#heatmap()
# UC1_expr
# P1_expr

###### p1
group_info_P1 <- group_P1
group_info_P1$Group <- ifelse(group_P1$Group == 1, "Periodontitis", "Control") 
row.names(group_info_P1) = group_info_P1$Samples
sort_index_P1 <- order(group_info_P1[["Group"]])
group_info_P1 <- group_info_P1[sort_index_P1, , drop = FALSE]
group_order_P1 <- group_info_P1$Samples
# expr_matrix_P1 <- matrix_P1[, group_order_P1]

##### UC1
group_info_UC1 <- group_UC1
group_info_UC1$Group <- ifelse(group_UC1$Group == 1, "UC", 
                               ifelse(group_UC1$Group == 0, "Control","Others")) 
row.names(group_info_UC1) = group_info_UC1$Samples

# UC1_layer_genes_selected = c("PTGR1","CXCL5","FOSB", UC1_layer_genes$PTGR1, UC1_layer_genes$CXCL5, UC1_layer_genes$FOSB)
GeneSymOlap <- union(UC1_layer_genes$PTGR1,UC1_layer_genes$CXCL5)
GeneSymOlap <- union(GeneSymOlap,UC1_layer_genes$FOSB)
UC1_layer_genes_selected = c("PTGR1","CXCL5","FOSB", GeneSymOlap)

# UC1_layer_genes_selected = c(UC1_layer_genes$PTGR1, UC1_layer_genes$CXCL5, UC1_layer_genes$FOSB)
UC1_layer_genes_selected_z = scale(UC1_layer_genes_selected)

UC1_expr_sel = UC1_expr[UC1_layer_genes_selected, ]
#UC1_expr_sel = UC1_expr[genes_of_interest, ]

normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

# adding more col annotation. 
#UC1_layer_genes_selected_z <- apply(UC1_expr_sel, 1, normalize)
# pheatmap(UC1_layer_genes_selected_z)
ann_colors = list(Group = c(Control = "#53a8e1",UC =  "#ca324c"))
#pheatmap(UC1_expr_sel, scale = 'row')
UC_p2 <- pheatmap(UC1_expr_sel, 
         scale = 'row', ## 归一化
         color = colorRampPalette(c("#0072B5", "white", "#BC3C29"))(100), #颜色
         #border_color = "black", #线的颜色
         fontsize = 14, #文字大小
         #fontfamily= "Arial",
         show_rownames = F,
         show_colnames = F,              
         cluster_cols = T, # 横向、纵向聚类
         cluster_rows = T,
         annotation_col = group_info_UC1[6],
         annotation_colors = ann_colors)
print(UC_p2)
ggsave(filename = paste0(output,"/UC1KeyLink_Pheatmap_in_UC.pdf"), 
       plot = UC_p2, width = 23, height = 15)

### optimizing 
library('RColorBrewer')
combined_annotations <- data.frame(group_info_UC1$Disease.extent, group_info_UC1[6])
colnames(combined_annotations) <- c("Subtype"," ")
UC_p2 <- pheatmap(UC1_expr_sel, 
                  scale = 'row', ## 归一化
                  cellwidth = 10, cellheight = 13,
                  #color = colorRampPalette(brewer.pal(11, "Spectral"), bias = 1)(100),
                  #colorRampPalette(c("#0072B5", "white", "white","#BC3C29"), bias = 1)(100), #颜色  # adjusting gradient of colors
                  colorRampPalette(c("#0072B5", "white", "#BC3C29"), bias = 1)(100), #颜色
                  border_color = "white", #线的颜色
                  fontsize = 16, #文字大小
                  #fontfamily= "Arial",
                  show_rownames = TRUE,
                  show_colnames = F,              
                  cluster_cols = T, # 横向、纵向聚类
                  cluster_rows = T,
                  annotation_col = group_info_UC1[6],#combined_annotations,
                  annotation_colors = ann_colors)
print(UC_p2)
ggsave(filename = paste0(output,"/UC1KeyLink_Pheatmap_in_UC1.pdf"), 
       plot = UC_p2, width = 23, height = 15)

################################################################################
### 5.0 function  
# UC1_genelist ----- 
# UC1_genelist <- unique(unlist(UC1_layer_genes))

UC1_genelist <- UC1_layer_genes_selected

perform_enrichment_analysis <- function(gene_list, ouput, prefix) {
  Genes <- bitr(gene_list,  
                fromType = "SYMBOL", 
                toType = c("ENTREZID"), 
                OrgDb = org.Hs.eg.db)
  
  if (is.null(Genes) || nrow(Genes) == 0) {
    stop(paste0("Gene mapping failed for ", prefix))
  }
  
  GO_enrichment <- enrichGO(gene = Genes$ENTREZID,
                            OrgDb = org.Hs.eg.db,
                            keyType = "ENTREZID",
                            ont = "all",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)
  
  if (!is.null(GO_enrichment) && nrow(as.data.frame(GO_enrichment)) > 0) {
    GO_result <- as.data.frame(GO_enrichment)
    write.table(GO_result, 
                file = paste0(ouput, "/", prefix, "_GO_enrichment.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    p1 <- dotplot(GO_enrichment, title = paste0(prefix, " GO Enrichment"), showCategory = 10) + 
      scale_color_gradientn(colors = c("#3C5488", "#E64B35"), name = "adj.p") +
      facet_grid(ONTOLOGY ~ ., scale = 'free')
    ggsave(filename = paste0(ouput, "/", prefix, "_GO_enrichment.png"), plot = p1, width = 12, height = 16, dpi = 300)
  } else {
    message("No significant GO enrichment results for ", prefix)
  }
  
  KEGG_enrichment <- enrichKEGG(gene = Genes$ENTREZID,
                                organism = "hsa", 
                                keyType = "kegg", 
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05)
  
  if (!is.null(KEGG_enrichment) && nrow(as.data.frame(KEGG_enrichment)) > 0) {
    KEGG_data <- DOSE::setReadable(KEGG_enrichment, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    KEGG_result <- as.data.frame(KEGG_data)
    write.table(KEGG_result, 
                file = paste0(ouput, "/", prefix, "_KEGG_enrichment.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    p2 <- dotplot(KEGG_enrichment, title = paste0(prefix, " KEGG Enrichment"), showCategory = 10) + 
      scale_color_gradientn(colors = c("#3C5488", "#E64B35"), name = "adj.p")
    ggsave(filename = paste0(ouput, "/", prefix, "_KEGG_enrichment.png"), plot = p2, width = 10, height = 8, dpi = 300)
  } else {
    message("No significant KEGG enrichment results for ", prefix)
  }
  DO_enrichment <- enrichDO(Genes$ENTREZID, 
                            pvalueCutoff = 0.05, 
                            qvalueCutoff = 0.05)
  
  if (!is.null(DO_enrichment) && nrow(as.data.frame(DO_enrichment)) > 0) {
    DO_result <- as.data.frame(DO_enrichment)
    # Convert gene IDs in DO results to gene symbols
    DO_result$geneID <- sapply(DO_result$geneID, function(gene_ids) {
      entrez_ids <- unlist(strsplit(gene_ids, "/"))
      gene_symbols <- mapIds(org.Hs.eg.db, keys = entrez_ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
      paste(gene_symbols, collapse = "/")
    })
    write.table(DO_result, 
                file = paste0(ouput, "/", prefix, "_DO_enrichment.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    p3 <- dotplot(DO_enrichment, title = paste0(prefix, " DO Enrichment"), showCategory = 10) + 
      scale_color_gradientn(colors = c("#3C5488", "#E64B35"), name = "adj.p")
    ggsave(filename = paste0(ouput, "/", prefix, "_DO_enrichment.png"), plot = p3, width = 10, height = 8, dpi = 300)
  } else {
    message("No significant DO enrichment results for ", prefix)
  }
  
  message("Enrichment analysis completed for ", prefix)
}

##### ##### ##### #####  optimized bubble plot
perform_enrichment_analysisOpt <- function(gene_list, ouput, prefix) {
  Genes <- bitr(gene_list,  
                fromType = "SYMBOL", 
                toType = c("ENTREZID"), 
                OrgDb = org.Hs.eg.db)
  
  if (is.null(Genes) || nrow(Genes) == 0) {
    stop(paste0("Gene mapping failed for ", prefix))
  }
  
  GO_enrichment <- enrichGO(gene = Genes$ENTREZID,
                            OrgDb = org.Hs.eg.db,
                            keyType = "ENTREZID",
                            ont = "all",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)
  
  if (!is.null(GO_enrichment) && nrow(as.data.frame(GO_enrichment)) > 0) {
    GO_result <- as.data.frame(GO_enrichment)
    write.table(GO_result, 
                file = paste0(ouput, "/", prefix, "_GO_enrichment.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    p1 <- dotplot(GO_enrichment, title = paste0(prefix, " GO Enrichment"), showCategory = 10) + 
      scale_color_gradientn(colors = c("#3C5488", "#E64B35"), name = "adj.p") +
      facet_grid(ONTOLOGY ~ ., scale = 'free')
    ggsave(filename = paste0(ouput, "/", prefix, "_GO_enrichment.png"), plot = p1, width = 12, height = 16, dpi = 300)
  } else {
    message("No significant GO enrichment results for ", prefix)
  }
  
  KEGG_enrichment <- enrichKEGG(gene = Genes$ENTREZID,
                                organism = "hsa", 
                                keyType = "kegg", 
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05)
  
  if (!is.null(KEGG_enrichment) && nrow(as.data.frame(KEGG_enrichment)) > 0) {
    KEGG_data <- DOSE::setReadable(KEGG_enrichment, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    KEGG_result <- as.data.frame(KEGG_data)
    write.table(KEGG_result, 
                file = paste0(ouput, "/", prefix, "_KEGG_enrichment.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    p2 <- dotplot(KEGG_enrichment, title = paste0(prefix, " KEGG Enrichment"), showCategory = 10) + 
      scale_color_gradientn(colors = c("#3C5488", "#E64B35"), name = "adj.p")
    ggsave(filename = paste0(ouput, "/", prefix, "_KEGG_enrichment.png"), plot = p2, width = 10, height = 8, dpi = 300)
  } else {
    message("No significant KEGG enrichment results for ", prefix)
  }
  DO_enrichment <- enrichDO(Genes$ENTREZID, 
                            pvalueCutoff = 0.05, 
                            qvalueCutoff = 0.05)
  
  if (!is.null(DO_enrichment) && nrow(as.data.frame(DO_enrichment)) > 0) {
    DO_result <- as.data.frame(DO_enrichment)
    # Convert gene IDs in DO results to gene symbols
    DO_result$geneID <- sapply(DO_result$geneID, function(gene_ids) {
      entrez_ids <- unlist(strsplit(gene_ids, "/"))
      gene_symbols <- mapIds(org.Hs.eg.db, keys = entrez_ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
      paste(gene_symbols, collapse = "/")
    })
    write.table(DO_result, 
                file = paste0(ouput, "/", prefix, "_DO_enrichment.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    p3 <- dotplot(DO_enrichment, title = paste0(prefix, " DO Enrichment"), showCategory = 10) + 
      scale_color_gradientn(colors = c("#3C5488", "#E64B35"), name = "adj.p")
    ggsave(filename = paste0(ouput, "/", prefix, "_DO_enrichment.png"), plot = p3, width = 10, height = 8, dpi = 300)
  } else {
    message("No significant DO enrichment results for ", prefix)
  }
  
  message("Enrichment analysis completed for ", prefix)
}

##### ##### ##### #####  optimized bar plot
######  output 
options(timeout = 120)
perform_enrichment_analysisOpt(gene_list=UC1_genelist,
                            ouput = output, 
                            prefix = "Key & link genes")

perform_enrichment_analysis(gene_list=UC1_genelist,
                            ouput = output, 
                            prefix = "Key & link genes")

################################################################################
### 5.0 function  v2
# UC1_genelist ----- 
# UC1_genelist <- unique(unlist(UC1_layer_genes))

gene_listUC=UC1_genelist
GenesUC <- bitr(gene_listUC,  
              fromType = "SYMBOL", 
              toType = c("ENTREZID"), 
              OrgDb = org.Hs.eg.db)

DEG_UC2 = copy(DEG_UC)
names(DEG_UC2)[1] <- "SYMBOL"
rownames(DEG_UC2) = DEG_UC2$SYMBOL
geneListUCFC <- DEG_UC2[gene_listUC,"logFC"]
names(geneListUCFC) <- gene_listUC

GO_enrichment <- enrichGO(gene = GenesUC$SYMBOL,
                          OrgDb = org.Hs.eg.db,
                          ont = "all",  ## BP, CC, MF
                          keyType       = 'SYMBOL',
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)  ## 0.05 p
summary(GO_enrichment)
barplot(GO_enrichment, showCategory=50)

GO_enrichmentBP <- enrichGO(gene = GenesUC$SYMBOL,
                          OrgDb = org.Hs.eg.db,
                          ont = "BP",  ## BP, CC, MF
                          keyType       = 'SYMBOL',
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)  ## 0.05 p
GO_enrichmentCC <- enrichGO(gene = GenesUC$SYMBOL,
                            OrgDb = org.Hs.eg.db,
                            ont = "CC",  ## BP, CC, MF
                            keyType       = 'SYMBOL',
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)  ## 0.05 p
GO_enrichmentMF <- enrichGO(gene = GenesUC$SYMBOL,
                            OrgDb = org.Hs.eg.db,
                            ont = "MF",  ## BP, CC, MF
                            keyType       = 'SYMBOL',
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)  ## 0.05 p
dotplot(GO_enrichmentBP, showCategory=20)
dotplot(GO_enrichmentCC, showCategory=20)
dotplot(GO_enrichmentMF, showCategory=20)

#### KEGG must using entrezid. #### #### #### #### #### #### #### #### #### #### 
KEGG_enrichment <- enrichKEGG(gene = GenesUC$ENTREZID,#GenesUC$SYMBOL,  ####
                              organism = "hsa", 
                              keyType = "kegg", 
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.05)
summary(KEGG_enrichment)   ### no significant terms. 
barplot(KEGG_enrichment, showCategory=20)

pdf(file="UC基因概念网络2.pdf",width = 5,height = 6)
cnetplot(GO_enrichment, foldChange=geneListUCFC,
         circular = TRUE,
         colorEdge = TRUE,
         color_category='firebrick', 
         color_gene='steelblue',
         layout="circle",
         edge.curved = 0.3
)
dev.off()

pdf(file="UC基因概念网络KEGG.pdf",width = 5,height = 6)
cnetplot(KEGG_enrichment, foldChange=geneListUCFC,
         circular = TRUE,
         colorEdge = TRUE,
         color_category='firebrick', 
         color_gene='steelblue',
         layout="circle",
         edge.curved = 0.3
)
dev.off()

barplot(GO_enrichment, showCategory=10)
barplot(KEGG_enrichment, showCategory=10)


p <- dotplot(GO_enrichment, showCategory=10)
p <- p + 
  scale_x_continuous(limits = c(0.05, 0.13))+ 
  theme(
  text = element_text(size = 12),
  
  strip.text = element_text(size = 20),
  axis.title = element_text(size = 12),  # 轴标题
  axis.text = element_text(size = 12), # 刻度标签
  aspect.ratio = 2  # 分面高度为宽度的60%
)
print(p)
ggsave("UCGO_dotplot.pdf", plot = p,width = 7,
       height = 5)

p2 <- dotplot(KEGG_enrichment, showCategory=10)
p2 <- p2 + 
  scale_x_continuous(limits = c(0.09, 0.17))+ 
  theme(
  text = element_text(size = 12),
  
  strip.text = element_text(size = 20),
  axis.title = element_text(size = 12),  # 轴标题
  axis.text = element_text(size = 12), # 刻度标签
  aspect.ratio = 1.8  # 分面高度为宽度的60%
)
print(p2)
ggsave("UCKEGG_dotplot.pdf", plot = p2,width = 7,
       height = 4)



###############################