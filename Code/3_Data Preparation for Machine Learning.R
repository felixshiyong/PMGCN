


# Identification and Visualization of Common DEGs between P and UC1 --------
library(pheatmap)
library(ggplot2)

rm(list = ls())

path <- getwd()
setwd(path)

input1 <- "1_Processed"
input2 <- "2_DE"
output <- "3_CommonDE/P1UC1"
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}


P1 <- read.delim2(paste0(input1, "/expr_P1.txt"), header = TRUE, row.names = 1)
Expr_P1 <- as.data.frame(lapply(P1, as.numeric))
rownames(Expr_P1) <- rownames(P1)

group_P1 <- read.csv(paste0(input1, "/Info_P1.csv"), header = TRUE)
DEG_P1 <- read.csv(file = paste0(input2,"/DEG_P1.csv"),row.names = 1)
P1_DE_genes <- DEG_P1[DEG_P1$adj.P.Val <0.05&abs(DEG_P1$logFC)> 1,]

UC1 <- read.delim2(paste0(input1, "/expr_UC1.txt"), header = TRUE, row.names = 1)
Expr_UC1 <- as.data.frame(lapply(UC1, as.numeric))
rownames(Expr_UC1) <- rownames(UC1)

group_UC1 <- read.csv(paste0(input1, "/Info_UC1.csv"), header = TRUE)
DEG_UC1 <- read.csv(file = paste0(input2,"/DEG_UC1.csv"),row.names = 1)
UC1_DE_genes <- DEG_UC1[DEG_UC1$adj.P.Val <0.05&abs(DEG_UC1$logFC)> 1,]

##Common DEGs
GeneSymOlap <- intersect(rownames(P1_DE_genes),rownames(UC1_DE_genes))
length(GeneSymOlap)

write.csv(GeneSymOlap,
          file = paste0(output,"/P1UC1_DEGs.csv"),
          quote=FALSE,
          row.names=FALSE,
          fileEncoding = "UTF-8")


#Com in P1
matrix_P1 <- Expr_P1[GeneSymOlap,]
write.table(data.frame(Symbol = rownames(matrix_P1),matrix_P1),
            file = paste0(output,"/P1UC1_DEGs_in_P1_matrix.txt"),
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

com_deg_in_P1 <- DEG_P1[GeneSymOlap,]
write.csv(com_deg_in_P1,
          file = paste0(output,"/P1UC1_deg_in_P1.csv"),
          quote = F,
          fileEncoding = "UTF-8")

#Com in UC1
matrix_UC1 <- Expr_UC1[GeneSymOlap,]
write.table(data.frame(Symbol = rownames(matrix_UC1),matrix_UC1),
            file = paste0(output,"/P1UC1_DEGs_in_UC1_matrix.txt"),
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

com_deg_in_UC1 <- DEG_UC1[GeneSymOlap,]
write.csv(com_deg_in_UC1,
          file = paste0(output,"/P1UC1_deg_in_UC1.csv"),
          quote = F,
          fileEncoding = "UTF-8")


###Heatmap plot
#P1
group_info_P1 <- group_P1
group_info_P1$Group <- ifelse(group_P1$Group == 1, "P", "Ctrl") 
row.names(group_info_P1) = group_info_P1$Samples

sort_index_P1 <- order(group_info_P1[["Group"]])
group_info_P1 <- group_info_P1[sort_index_P1, , drop = FALSE]

group_order_P1 <- group_info_P1$Samples
expr_matrix_P1 <- matrix_P1[, group_order_P1]


ann_colors = list(Group = c(Ctrl = "#20854E",P =  "#E18727FF"))

p1 <- pheatmap(expr_matrix_P1,
               color = colorRampPalette(c("#0072B5", "white", "#BC3C29"))(100), #颜色
               scale = "row", #归一化的方式
               border_color = NA, #线的颜色
               fontsize = 14, #文字大小
               fontfamily= "Arial",
               show_rownames = F,
               show_colnames = F,              
               cluster_cols = F, # 横向、纵向聚类
               cluster_rows = T,
               annotation_col = group_info_P1[4],
               annotation_colors = ann_colors)   

ggsave(filename = paste0(output,"/P1UC1_Pheatmap_in_P1.png"), 
       plot = p1, width = 12, height = 6, dpi = 300)

#UC1
group_info_UC1 <- group_UC1
group_info_UC1$Group <- ifelse(group_UC1$Group == 1, "UC", 
                               ifelse(group_UC1$Group == 0, "Ctrl","Others")) 
row.names(group_info_UC1) = group_info_UC1$Samples


select_UC1 <- group_info_UC1[group_info_UC1[["Group"]]%in%c("UC","Ctrl"),]
sort_index_UC1 <- order(select_UC1[["Group"]])
select_UC1 <- select_UC1[sort_index_UC1,,drop = F]
samples <- select_UC1[,1]

sub_matrix_UC1 <- matrix_UC1[,samples]
ann_colors = list(Group = c(Ctrl = "#20854E",UC =  "#E18727FF"))

p2 <- pheatmap(sub_matrix_UC1,
               color = colorRampPalette(c("#0072B5", "white", "#BC3C29"))(100), #颜色
               scale = "row", #归一化的方式
               border_color = NA, #线的颜色
               fontsize = 14, #文字大小
               fontfamily= "Arial",
               show_rownames = F,
               show_colnames = F,              
               cluster_cols = F, # 横向、纵向聚类
               cluster_rows = T,
               annotation_col = group_info_UC1[6],
               annotation_colors = ann_colors)   

ggsave(filename = paste0(output,"/P1UC1_Pheatmap_in_UC1.png"), 
       plot = p2, width = 12, height = 6, dpi = 300)




# Integrated Biomarker Selection from Network and DEG Analyses   ------------

rm(list = ls())
path <- getwd()
setwd(path)


input1 <- "1_Processed"
input2 <- "2_DE"
input3 <- "3_CommonDE/P1UC1"


output <- "7_Feature/P1UC1"
if (!dir.exists(output)) {
  dir.create(output,recursive = TRUE)
}



UC1_expr <- read.table(file = paste0(input1,"/expr_UC1.txt"),
                       sep = "\t",
                       header = T,
                       row.names = 1)
UC2_expr <- read.table(file = paste0(input1,"/expr_UC2.txt"),
                       sep = "\t",
                       header = T,
                       row.names = 1)
UC3_expr <- read.table(file = paste0(input1,"/expr_UC3.txt"),
                       sep = "\t",
                       header = T,
                       row.names = 1)



lines <- readLines("Network_results/filtered_genes_results_P1UC1.csv")
algorithm_genes <- list()
current_algorithm <- NULL

for (line in lines) {
  entries <- unlist(strsplit(line, ","))
  if (nchar(entries[1]) > 0) {
    current_algorithm <- entries[1]  # 记录当前算法名称
    algorithm_genes[[current_algorithm]] <- entries[-1]  # 初始化基因列表
  } else {
    algorithm_genes[[current_algorithm]] <- c(algorithm_genes[[current_algorithm]], entries)
  }
}


algorithm_genes <- lapply(algorithm_genes, function(genes) {
  genes[nchar(genes) > 0]
})


all_genes <- unique(unlist(algorithm_genes))


gene_matrix <- matrix(0, nrow = length(algorithm_genes), ncol = length(all_genes),
                      dimnames = list(names(algorithm_genes), all_genes))

for (algorithm in names(algorithm_genes)) {
  genes <- algorithm_genes[[algorithm]]
  gene_matrix[algorithm, genes] <- 1
}


# Extract key genes from network-based methods(percolation)
target_algorithm <- "CoreHD-DCI"
target_algorithm_gene <- lapply(target_algorithm, function(alg) algorithm_genes[[alg]])
target_algorithm_gene_list <- target_algorithm_gene[[1]]

ComKey_UC1_expr <- UC1_expr[target_algorithm_gene_list,]
write.table(data.frame(gene_symbol=rownames(ComKey_UC1_expr),ComKey_UC1_expr),
            file = paste0(output,"/ComKey_UC1_expr.txt"),
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)
ComKey_UC2_expr <- UC2_expr[target_algorithm_gene_list,]
write.table(data.frame(gene_symbol=rownames(ComKey_UC2_expr),ComKey_UC2_expr),
            file = paste0(output,"/ComKey_UC2_expr.txt"),
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

ComKey_UC3_expr <- UC3_expr[target_algorithm_gene_list,]
write.table(data.frame(gene_symbol=rownames(ComKey_UC3_expr),ComKey_UC3_expr),
            file = paste0(output,"/ComKey_UC3_expr.txt"),
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)



# Common DEGs
com_deg_in_UC1 <- read.csv(paste0(input3,"/P1UC1_deg_in_UC1.csv"),
                           header = T)


com_deg_genes <- com_deg_in_UC1[com_deg_in_UC1$adj.P.Val < 0.05, ]
com_deg_genes_sort <- com_deg_genes[order(abs(com_deg_genes$logFC), decreasing = TRUE), ]
com_deg_genes_top <- com_deg_genes_sort$X 


valid_genes <- com_deg_genes_top[com_deg_genes_top %in% rownames(UC2_expr) & com_deg_genes_top %in% rownames(UC3_expr)]


if (length(valid_genes) < 3) {
  additional_genes <- setdiff(com_deg_genes_top, valid_genes)  # UC1 中未被选中的基因
  additional_valid_genes <- additional_genes[additional_genes %in% rownames(UC2_expr) & additional_genes %in% rownames(UC3_expr)]
  valid_genes <- c(valid_genes, head(additional_valid_genes, 3 - length(valid_genes)))
}


com_degregulate_genes <- head(valid_genes, 3)


#UC1
com_degregulate_UC1_expr <- UC1_expr[com_degregulate_genes,]
write.table(data.frame(gene_symbol=rownames(com_degregulate_UC1_expr),com_degregulate_UC1_expr),
            file = paste0(output,"/CommonDEG_UC1_expr.txt"),
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

#UC2
com_degregulate_UC2_expr <- UC2_expr[com_degregulate_genes,]
write.table(data.frame(gene_symbol=rownames(com_degregulate_UC2_expr),com_degregulate_UC2_expr),
            file = paste0(output,"/CommonDEG_UC2_expr.txt"),
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

#UC3
com_degregulate_UC3_expr <- UC3_expr[com_degregulate_genes,]
write.table(data.frame(gene_symbol=rownames(com_degregulate_UC3_expr),com_degregulate_UC3_expr),
            file = paste0(output,"/CommonDEG_UC3_expr.txt"),
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)


##UC1 DEG
deg_in_UC1 <- read.csv(file = paste0(input2,"/DEG_UC1.csv"),header = T)
deg_genes <- deg_in_UC1[deg_in_UC1$adj.P.Val < 0.05,]
deg_genes_sort <- deg_genes[order(abs(deg_genes$logFC),decreasing = TRUE),]
degregulate_genes <- deg_genes_sort$X



valid_genes <- degregulate_genes[degregulate_genes %in% rownames(UC2_expr) & degregulate_genes %in% rownames(UC3_expr)]


if (length(valid_genes) < 3) {
  additional_genes <- setdiff(degregulate_genes, valid_genes)  
  additional_valid_genes <- additional_genes[additional_genes %in% rownames(UC2_expr) & additional_genes %in% rownames(UC3_expr)]
  valid_genes <- c(valid_genes, head(additional_valid_genes, 3 - length(valid_genes)))
}


degregulate_genes <- head(valid_genes, 3)


#UC1
degregulate_UC1_expr <- UC1_expr[degregulate_genes,]
write.table(data.frame(gene_symbol=rownames(degregulate_UC1_expr),degregulate_UC1_expr),
            file = paste0(output,"/deg_UC1_expr.txt"),
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)
#UC2
degregulate_UC2_expr <- UC2_expr[degregulate_genes,]
write.table(data.frame(gene_symbol=rownames(degregulate_UC2_expr),degregulate_UC2_expr),
            file = paste0(output,"/deg_UC2_expr.txt"),
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)
#UC3
degregulate_UC3_expr <- UC3_expr[degregulate_genes,]
write.table(data.frame(gene_symbol=rownames(degregulate_UC3_expr),degregulate_UC3_expr),
            file = paste0(output,"/deg_UC3_expr.txt"),
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)



