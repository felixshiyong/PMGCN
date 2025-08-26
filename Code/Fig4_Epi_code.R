library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(ggsci)
library(patchwork)
library(ggrepel)
library(ggpubr)
library(reshape2)
library(rstatix)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(AnnotationDbi)

#1.0 读取数据 ------------------------------------------------------------------
rm(list = ls())
gc()
path <- getwd()
setwd(path)

input13 <- "13_sc_cluster/Epi/UC"
output <- "15_Epi/expr/UC"

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

ss <- readRDS(file = paste0(input13, "/ss_anno3_Epi.rds"))
table(ss$cell_type)
table(ss$histo3)
table(ss$label)
# ss$cell_type <- ss$cell_type %>%
#   recode("LND" = "Mature Colonocytes")
ss$label <- recode(ss$label,
                   `healthy` = "Healthy",
                   `inflamed UC` = "Inflamed UC",
                   `non inflamed UC` = "Non inflamed UC")
table(ss$label)

# 修改 label 列的因子顺序
ss$label <- factor(ss$label, levels = c("Healthy", "Non inflamed UC", "Inflamed UC"))
# 修改 label 列的因子顺序
ss$histo3 <- factor(ss$histo3, levels = c("Healthy (N=14)", "Non inflamed UC (N=13)", "Inflamed UC (N=23)"))
table(ss$label)
table(ss$histo3)
#ss <- RunUMAP(ss, dims = 1:20,reduction = "harmony", min.dist = 0.3, spread =2)
ss <- RunUMAP(ss, dims = 1:20,reduction = "harmony", min.dist = 1, spread =7)

table(ss$cell_type)
ss$cell_type2 <- recode(ss$cell_type,
                    'undifferented cell' = "Undifferented cell"
)
table(ss$cell_type2)

# 2. 整体可视化 ----------------------------------------------------------------
type_colors <- c(
  "TA" = '#c8c7e1',
  "Early Colonocytes" = '#c4daec',
  "Intermediate Colonocytes" = '#3ca0cf',
  "Mature Colonocytes" = '#405993',
  "LND" = '#9a70a8',
  
  "BEST4/OTOP2" = '#d69971',   
  "Enteroendocrine" = '#64a776', 
  "Goblet" = '#df5734',
  
  
  "Tuft" = '#d25774',
  "Undifferented cell" = '#e6e2a3',
  "Stem" = '#696a6c'
)
p_umap <- DimPlot(ss, 
                  #pt.size = 1,
                  reduction = "umap", 
                  group.by = "cell_type2", 
                  label = F, 
                  #seed = 23,
                  repel = TRUE) +
  scale_color_manual(values = type_colors)+
  ggtitle("Annotation")+  # 指定颜色
  theme(
    panel.border = element_blank(),  # 去除边框
    panel.grid = element_blank(),    # 去除网格
    axis.text = element_blank(),     # 去除坐标轴文本
    axis.ticks = element_blank(),     # 去除坐标轴刻度
    aspect.ratio = 1,
    # 标题字体设置（主标题）
    plot.title = element_blank(),
    # 图例文字设置
    legend.text = element_text(
      size = 12
    ),
    # 图例标题设置（若需保留）
    legend.title = element_text(
      size = 12,          # 图例标题字号（如"cell_type"）
      face = "bold"
    )
  )   
p_umap

ggsave(filename = paste0(output,"/","Fig4.a_UMAP_Cell_Type1.pdf"), 
       plot = p_umap, 
       width = 8, 
       height = 6)

#### 按疾病状态可视化 UMAP
disease_colors <- c(
  "Healthy" = "#57b1ab",   
  "Non inflamed UC" = "#ffa74f", 
  "Inflamed UC" = "#9a72c7"
)

p_state <- DimPlot(ss, 
                   reduction = "umap", 
                   group.by = "label", 
                   label = F, 
                   repel = TRUE) +
  scale_color_manual(values = disease_colors)+  # 指定颜色
  theme(
    panel.border = element_blank(),  # 去除边框
    panel.grid = element_blank(),    # 去除网格
    axis.text = element_blank(),     # 去除坐标轴文本
    axis.ticks = element_blank()     # 去除坐标轴刻度
  ) +
  ggtitle("Histology")
p_state
ggsave(filename = paste0(output,"/UMAP_Disease_State.png"),
       plot = p_state,
       width = 8,
       height = 6)
# ggsave(filename = paste0(output,"/UMAP_Disease_State.pdf"), 
#        plot = p_state, 
#        width = 8, 
#        height = 6)

# 2.1 bubble plot ----------------------------------------------------------------
##注释结果
Epi_cell_markers <- list(
  "all" = c("EPCAM"),
  "Stem" = c("LGR5","OLFM4"),
  "TA" = c("MKI67", "TOP2A", "PCNA"),
  "Early Colonocytes" = c("CA2", "SLC26A2", "FABP1"),
  "Intermediate Colonocytes" = c("B3GNT7","ABR","ADH1C","STEAP3","ATP5G1","PCNP"),
  "Mature Colonocytes" = c("AQP8","GUCA2A","CA4","CEACAM1"),
  "LND" = c( "LCN2","NOS2" ,"DUOX2" ),
  "BEST4/OTOP2" = c("BEST4", "OTOP2", "CA7"),
  "Goblet" = c("CLCA1", "SPDEF", "FCGBP", "ZG16", "MUC2","TFF3"),
  "Enteroendocrine" = c("CHGA", "CHGB", "NEUROD1"),
  "Tuft" = c("POU2F3", "LRMP", "TRPM5")
)

all_marker <- unique(unlist(Epi_cell_markers))
p1 <- DotPlot(ss,
              features = all_marker,
              group.by = "cell_type",
              scale = F)+
  #cluster.idents = TRUE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_flip()
p1

table(ss$cell_type)
##按顺序画点图
##不展示某一类细胞
##ss_sub <- subset(ss, subset = cell_type != "undifferented cell")

dot_data <- DotPlot(ss,
                    features = unique(unlist(Epi_cell_markers)),
                    group.by = "cell_type",
                    cluster.idents = T,
                    scale = F)$data
dot_data
# 获取 y 轴（细胞类型）的顺序
cell_type_order <- levels(factor(dot_data$id))

# 按细胞类型顺序整理 marker_genes
ordered_markers <- unlist(Epi_cell_markers[cell_type_order])
# 补充没有的基因名称到完整序列中
all_markers <- unique(ordered_markers)
all_markers <- unique(c(all_markers,"EPCAM"))

# 重新绘制点图
p_filtered_marker <- DotPlot(ss,
                             features = all_markers,
                             group.by = "cell_type",
                             cluster.idents = T,
                             scale = F) +
  scale_color_gradient(low = "#FEE0D2", high = "#BC3C29FF") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(margin = margin(t = 0, r = 2, b = 0, l = 0)),  # 调整 y 轴标签间距
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加黑色边框
    panel.background = element_blank()  # 背景设为透明
  ) +
  #guides(size = guide_legend(title = "Dot Size")) + # 添加图例
  guides(size = guide_legend(title = "Percent Expressed", order = 1),
         color = guide_colorbar(title = "Average Expression", order = 2)) + # 添加图例
  scale_size(range = c(0, 10)) +
  coord_flip()

p_filtered_marker

ggsave(paste0(output,'/Fig4.b_dot_anno_marker.pdf'), 
       plot = p_filtered_marker, height = 12, width =7)


# 2.2 proportion plot ----------------------------------------------------------

# 计算每个疾病状态中每个细胞类型的比例
meta <- ss@meta.data 

# 计算每个疾病状态中每个细胞类型的比例
cell_proportions <- meta %>%
  group_by(label, cell_type) %>%
  summarise(count = n()) %>%
  group_by(label) %>%  # 按疾病状态重新分组
  mutate(proportion = count / sum(count)) %>%  # 计算每个疾病状态内细胞类型的比例
  arrange(label, cell_type)  # 排序

# 绘制比例条形图
p_proportions <- ggplot(cell_proportions, aes(x = label, y = proportion, fill = cell_type)) + 
  geom_bar(stat = "identity", position = "stack") +  # 堆积条形图
  # 注释掉 geom_text 以去除百分比标签
  # geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)),  # 显示百分比标签
  #           position = position_stack(vjust = 0.5), 
  #           size = 3) +  # 设置标签位置和大小
  labs(x = "Disease State", y = "Proportion", fill = "Cell Type") +  # 更新标签
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # 调整横坐标标签
        axis.text.y = element_text(margin = margin(t = 0, r = 2, b = 0, l = 0)),
        panel.border = element_blank(),  # 去除面板边框
        axis.line = element_line(color = "black", size = 0.5),  # 只保留外部的 x 和 y 轴线
        axis.ticks.y = element_line(color = "black", size = 0.5),  # 设置 y 轴刻度线
        axis.ticks.x = element_blank(),  # 去除 x 轴刻度线
        panel.grid = element_blank(),  # 去除网格
        axis.ticks.length = unit(0.2, "cm")) +  # 设置刻度线长度
  ggtitle("") +
  scale_fill_manual(values = type_colors)  +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05))
p_proportions
# 保存比例图
ggsave(filename = paste0(output, "/Fig4.c_Cell_Proportions_Stacked.pdf"),
       plot = p_proportions, width = 5, height = 8, dpi = 300)

# 3. PTGR1整体表达 ---------------------------------------------------------------
#vital gene 表达

# scale_fill_manual(values = c("#BC3C29FF", "#0072B5FF", "#E18727FF"))
#color = c("#B09C85FF", "#3C5488", "#E64B35"),
vital_genes <- c("PTGR1")
ss$anno_group <- paste(ss$cell_type, ss$group, sep = "_")
for (gene in vital_genes) {
  dotplot_data <- DotPlot(ss, 
                          features = gene, 
                          group.by = "anno_group",
                          scale = F)$data %>%
    separate(id, into = c("Cluster", "Group"), sep = "_")
  p_gene <- ggplot(dotplot_data, aes(x = Group, y = Cluster)) +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
    scale_color_gradient(low = "#FEE0D2", high = "#BC3C29FF", na.value = "grey") +
    scale_size(range = c(1, 10)) +
    labs(x = "Group", y = "Cluster", 
         color = "Average Expression",  # 修改颜色图例名称
         size = "Percent Expressed") +  # 修改大小图例名称
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.text.y = element_text(margin = margin(t = 0, r = 2, b = 0, l = 0)),
      #legend.position = "bottom",
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.background = element_blank()
      
    )
  dotplot_data
  ggsave(filename = paste0(output, "/", gene, "_Dot_group.png"), 
         plot = p_gene, width = 5, height = 6, dpi = 300)
}

###健康、非炎性、炎性
table(ss$label)
table(ss$cell_type)
ss$anno_label <- paste(ss$cell_type, ss$label, sep = "_")
# filtered_metadata <- ss@meta.data %>%
#   filter(label %in% c("healthy", "inflamed UC"))
# filtered_ss <- subset(ss, cells = rownames(filtered_metadata))
table(ss$anno_label)

for (gene in vital_genes) {
  dotplot_data <- DotPlot(ss, 
                          features = gene, 
                          group.by = "anno_label",
                          scale = F)$data %>%
    separate(id, into = c("Cluster", "Group"), sep = "_")
  
  # 修改 Group 列的因子顺序
  dotplot_data$Group <- factor(dotplot_data$Group,
                               levels = c("Healthy", "Non inflamed UC", "Inflamed UC"))
  
  p_gene <- ggplot(dotplot_data, aes(x = Group, y = Cluster)) +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
    scale_color_gradient(low = "#FEE0D2", high = "#BC3C29FF", na.value = "grey") +
    
    # 使用scale_color_gradient2，设置低值为蓝色，0值为白色，高值为红色
    #scale_color_gradient2(low = "#0072B5FF", mid = "white", high = "#BC3C29FF", midpoint = 0, na.value = "grey") +
    
    scale_size(range = c(1, 10)) +
    labs(x = "Group", y = "Cluster", 
         color = "Average Expression", 
         size = "Percent Expressed") +  
    ggtitle(gene) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(margin = margin(t = 0, r = 2, b = 0, l = 0)),
      #legend.position = "bottom",
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.background = element_blank()
      
    )+
    coord_flip()
  p_gene
  ggsave(filename = paste0(output, "/", gene, "_Dot_label.png"), 
         plot = p_gene, width = 8, height = 4, dpi = 300)
}

##### UMAP
p_vital <- FeaturePlot(ss, 
                       features = vital_genes,
                       reduction = "umap", 
                       #ncol = 3,
                       #split.by = "histo3",
                       cols = c("lightgrey", "#BC3C29FF")) 
p_vital
ggsave(filename = paste0(output,"/","UMAP_vital.png"), 
       plot = p_vital, 
       width = 4, 
       height = 4)
###
p_expr_group <- FeaturePlot(ss, 
                            features = vital_genes, 
                            split.by = "label",
                            cols = c("lightgrey", "#BC3C29FF"))
p_expr_group
ggsave(filename = paste0(output, "/", "UMAP_vital_group.png"), 
       plot = p_expr_group, width = 12, height = 4, dpi = 300)

# VlnPlot 
#c("#BC3C29FF", "#0072B5FF", "#E18727FF"))
disease_colors <- c(
  "Healthy" = "#57b1ab",   
  "Non inflamed UC" = "#ffa74f", 
  "Inflamed UC" = "#9a72c7"
)

p_expr_vln <- VlnPlot(ss, 
                      features = vital_genes, 
                      group.by = "cell_type2", 
                      split.by = "label",
                      cols = disease_colors,
                      log = TRUE,          # 对数变换
                      pt.size = 0,         # 不显示散点
                      layer = "data",   
                      combine = TRUE) 
p_expr_vln
ggsave(filename = paste0(output, "/", "vital_Vln_Expr.png"), 
       plot = p_expr_vln, width = 12, height = 4, dpi = 300)

p_expr_vln2 <- VlnPlot(ss, 
                      features = vital_genes, 
                      group.by = "cell_type2", 
                      #split.by = "label",
                      cols = type_colors,
                      log = TRUE,          # 对数变换
                      pt.size = 0,         # 不显示散点
                      layer = "data",   
                      combine = TRUE) 
p_expr_vln2
ggsave(filename = paste0(output, "/", "Fig4.d_vital_Vln.pdf"), 
       plot = p_expr_vln2, width = 8, height = 6, dpi = 300)


# 4. 细胞数量占比变化 ----------------------------------------------------------
ss$anno <- ss$cell_type
Idents(ss) <- "anno"
table(ss$anno)     # 注释信息

ss$NewCellType <- "Other"  
cell_types <- unique(ss$anno)
# 遍历每个细胞类型并标记 PTGR1 表达等于零的细胞
for (cell_type in cell_types) {
  
  cells_subset <- subset(ss, idents = cell_type)
  RNA_data <- GetAssayData(cells_subset, assay = "RNA", layer = "data")
  # 将 PTGR1 基因表达<= 0 的细胞标记为 "PTGR1- <细胞类型>"
  ss$NewCellType[Cells(cells_subset)[RNA_data["PTGR1", ] > 0]] <- paste0("PTGR1+ ", cell_type)
}

table(ss$NewCellType)
ptgr1_table <- as.data.frame(table(ss$NewCellType))
write.csv(ptgr1_table, file =paste0(output,"/PTGR1-_CellType_Table.csv") , row.names = FALSE)

#This function generate figures with selected cell type.
prop.plot <- function(input,title,color,hide.ns,method,legendLocation,stat.test,ct,strip.size=8,ncol=7,plot.margin=margin(0,0,0,2),title.y,vjust=-0.3,format=F,mult=c(0, 0.15)){
  input <- cbind(input[,ct],input$`sampleInfo[rownames(prop_ct), 2]`)
  tmp.plot <- melt(input,value.name="Prop",variable.name="cluster");
  colnames(tmp.plot)[1] <- "Groups";dim(tmp.plot)
  
  m <- aggregate(tmp.plot$Prop,by=list(tmp.plot$cluster),max)
  rownames(m) <- m$Group.1
  n <- length(unique(tmp.plot$Groups))
  
  stat.test <- stat.test[stat.test$cluster %in% ct,]
  stat.test$cluster <- factor(stat.test$cluster,levels = ct)
  
  for (i in levels(tmp.plot$cluster)) {
    if(n==2){
      stat.test$y.position[stat.test$cluster %in% i] <- c(m[i,2]+(m[i,2]/10))
    }else{
      stat.test$y.position[stat.test$cluster %in% i] <- c(m[i,2]+(m[i,2]/10)*3,m[i,2]+(m[i,2]/10),m[i,2]+(m[i,2]/10)*6)
    }
  }
  
  stat.test$y.position <- as.numeric(stat.test$y.position)
  
  if(isTRUE(format)){print("Formatted");stat.test$p.adj <- formatC(stat.test$p.adj,format="e",digits=1)}
  p1 <- ggplot(tmp.plot, aes(x=Groups, y=Prop,col=Groups)) + 
    
    # #散点图
    geom_jitter(shape=16,size=1)+facet_wrap(~cluster,scales="free",ncol = ncol)+
    
    # #小提琴图
    # geom_violin(trim=FALSE, size=0.8, alpha=0.7) +  
    # facet_wrap(~cluster, scales="free", ncol=3) +  
    # geom_jitter(position=position_jitter(width=0.2), size=0.5, alpha=0.6)+  # 添加样本点
    # 
    
    theme_bw()+labs(x=" ",y=title.y,title = title)+
    theme(
      panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
      panel.spacing.y = unit(0, "lines"),
      #plot.title = element_text(size=8,hjust = 0.5,face = "bold",margin = margin(0,0,0,0)),
      plot.title = element_blank(),
      strip.text = element_text(size=strip.size,margin = margin(0,0,2,0)),
      strip.background = element_blank(),
      legend.position = legendLocation,
      legend.text = element_text(colour = "black",size = 10),
      legend.title = element_blank(),legend.key.size = unit(0.5, "lines"),
      legend.margin = margin(0,0,0,0),legend.box.margin = margin(0,0,0,-5),legend.spacing.x = unit(0.1,"cm"),
      axis.text.x = element_blank(),axis.ticks.x = element_blank(),
      axis.text.y = element_text(color = "black",size = 5,margin = margin(r = 0)),
      axis.ticks.y = element_line(color = "black",size = 0.25),
      axis.ticks.length.y = unit(0.05,"cm"),
      axis.title.y = element_text(color = "black",size = 10),
      axis.title.x = element_blank(),plot.margin = plot.margin)+ 
    scale_color_manual(values=color)+
    scale_y_continuous(expand = expansion(mult = mult))+
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", color="black", width=0.1) +
    stat_summary(fun=mean, geom="point", color="black")
  
  p1 <- p1 + stat_pvalue_manual(stat.test, 
                                label = "p.adj",
                                hide.ns = T,
                                tip.length = 0,
                                size = 1.5,
                                color = "red",
                                vjust = vjust)
  return(p1)
}
p1

# ~4.1 PTGR1+细胞占比 -------------------------------------------------------------
#This function generate proportion of each cell type for each patient.
prop.fun <- function(bottom,ct.epi){
  if(bottom=="all") {
    meta <- ss@meta.data
    sampleInfo <- unique(data.frame(sample=meta$PatientID,histo3=meta$histo3))
    rownames(sampleInfo) <-sampleInfo$sample
    prop_ct <- table(meta$PatientID,meta$NewCellType)
    prop_ct<- as.data.frame.array(prop_ct/rowSums(prop_ct))
    prop_ct <- prop_ct[,ct.epi]
    input <- cbind(prop_ct,sampleInfo[rownames(prop_ct),2])
  }
  return(input)
}

# ss$histo3 <- factor(ss$histo3,levels = c("Control","Inactive CD","Active CD"))
# sample <- unique(ss@meta.data[,c("PatientID","histo3")]);table(sample$histo3)
# levels(ss$histo3) <- paste0(levels(ss$histo3)," (N=",table(sample$histo3)[levels(ss$histo3)],")")
ss$PatientID <- paste(ss$orig.ident, ss$batch, ss$label, sep = "_")
table(ss$PatientID)

input <- prop.fun(bottom = "all",ct.epi = unique(ss$NewCellType))
tmp.plot <- melt(input,value.name="Prop",variable.name="cluster")
colnames(tmp.plot)[1] <- "Groups"
dim(tmp.plot)

# 按细胞类型（cluster）和疾病状态（group）进行分组，并进行Wilcoxon检验
stat.test <- tmp.plot %>%
  group_by(cluster) %>%  # 按细胞类型分组
  do({
    # 当前细胞类型的数据
    group_data <- .
    
    # 获取不同的疾病状态组合
    disease_combinations <- combn(unique(group_data$Groups), 2, simplify = FALSE)
    
    # 对每个疾病状态组合执行 Wilcoxon 检验
    results <- lapply(disease_combinations, function(groups) {
      group1 <- groups[1]
      group2 <- groups[2]
      
      # 筛选出当前组合的两种疾病状态的数据
      test_data <- group_data %>% filter(Groups %in% c(group1, group2))
      
      # 执行 Wilcoxon 检验
      test_result <- wilcox.test(Prop ~ Groups, data = test_data)
      
      # 返回检验结果，包括细胞类型、疾病状态和p值
      data.frame(
        Group1 = group1,
        Group2 = group2,
        p_value = test_result$p.value,
        cluster = unique(group_data$cluster)
      )
    })
    
    # 合并所有的检验结果
    do.call(rbind, results)
  }) %>%
  ungroup()
stat.test
# FDR adjustment
# stat.test <- stat.test %>%
#   mutate(group = paste0(Group1, Group2)) %>%
#   mutate(p.adj = p.adjust(p_value, method = "fdr")) %>%
#   group_by(group) %>%
#   adjust_pvalue(method = "fdr") %>%
#   add_significance()

# FDR adjustment
#不要adj.p
stat.test <- stat.test %>%
  mutate(group = paste0(Group1, Group2)) %>%
  #mutate(p.adj = p.adjust(p_value, method = "fdr")) %>%
  mutate(p.adj = p_value) %>%
  group_by(group) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()

write.csv(stat.test, 
          file = paste0(output,"/Test_results_PTGR1-.csv"), 
          row.names = FALSE)
# 
all_cell_types <- unique(ss$NewCellType)
# 从 all_cell_types 中移除 "Other"
ct.imm <- setdiff(all_cell_types, "Other")
length(ct.imm)

# 修改 stat.test 格式，使其包含正确的列名和数据
stat.test <- stat.test %>%
  dplyr::rename(group1 = Group1, group2 = Group2, p.adj = p.adj)  # 重命名列
head(stat.test)
# "Healthy" = "#57b1ab",   
# "Non inflamed UC" = "#ffa74f", 
# "Inflamed UC" = "#9a72c7"
p_prop <- prop.plot(input, title = " ", 
                    color = c("#57b1ab", "#ffa74f", "#9a72c7"),
                    hide.ns = F, 
                    method = "fdr", 
                    legendLocation = "bottom", 
                    stat.test = stat.test, 
                    ct = ct.imm, 
                    strip.size = 10, 
                    ncol = 4, 
                    title.y = "Proportion", 
                    format = T)

p_prop
ggsave(filename = paste0(output,"/Proportion_PTGR1+.pdf"),
       plot = p_prop, width = 10, height = 5, dpi = 300)
ggsave(filename = paste0(output,"/Proportion_PTGR1+.png"),
       plot = p_prop, width = 10, height = 5, dpi = 300)
# ~4.2 all 细胞占比 -------------------------------------------------------------------
table(ss$anno)
prop.fun <- function(bottom, ct.epi) {
  
  if (bottom == "all") {
    meta <- ss@meta.data
    
    sampleInfo <- unique(data.frame(sample = meta$PatientID, histo3 = meta$histo3))
    
    rownames(sampleInfo) <- sampleInfo$sample
    
    # 计算每个患者的 vital+CellType 的比例
    prop_ct <- table(meta$PatientID, meta$anno)  # 计算患者ID和anno的频数表
    prop_ct <- as.data.frame.array(prop_ct / rowSums(prop_ct))  # 计算比例（除以每个患者的总数）
    
    # 选择 anno 列中指定的细胞类型（ct.epi 参数控制）
    prop_ct <- prop_ct[, ct.epi]
    
    # 将比例数据与样本信息结合，返回一个数据框（包含 anno 比例和 histo3 信息）
    input <- cbind(prop_ct, sampleInfo[rownames(prop_ct), 2])
  }
  
  return(input)
}

input <- prop.fun(bottom = "all",ct.epi = unique(ss$anno))
tmp.plot <- melt(input,value.name="Prop",variable.name="cluster")
colnames(tmp.plot)[1] <- "Groups"
dim(tmp.plot)

# 按细胞类型（cluster）和疾病状态（group）进行分组，并进行Wilcoxon检验
stat.test <- tmp.plot %>%
  group_by(cluster) %>%  # 按细胞类型分组
  do({
    # 当前细胞类型的数据
    group_data <- .
    
    # 获取不同的疾病状态组合
    disease_combinations <- combn(unique(group_data$Groups), 2, simplify = FALSE)
    
    # 对每个疾病状态组合执行 Wilcoxon 检验
    results <- lapply(disease_combinations, function(groups) {
      group1 <- groups[1]
      group2 <- groups[2]
      
      # 筛选出当前组合的两种疾病状态的数据
      test_data <- group_data %>% filter(Groups %in% c(group1, group2))
      
      # 执行 Wilcoxon 检验
      test_result <- wilcox.test(Prop ~ Groups, data = test_data, exact = FALSE)
      
      # 返回检验结果，包括细胞类型、疾病状态和p值
      data.frame(
        Group1 = group1,
        Group2 = group2,
        p_value = test_result$p.value,
        cluster = unique(group_data$cluster)
      )
    })
    
    # 合并所有的检验结果
    do.call(rbind, results)
  }) %>%
  ungroup()

stat.test
# FDR adjustment
#不要adj.p
stat.test <- stat.test %>%
  mutate(group = paste0(Group1, Group2)) %>%
  #mutate(p.adj = p.adjust(p_value, method = "fdr")) %>%
  mutate(p.adj = p_value) %>%
  group_by(group) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()

write.csv(stat.test, 
          file = paste0(output,"/Test_results_all.csv"), 
          row.names = FALSE)
# 
all_cell_types <- unique(ss$anno)
# 从 all_cell_types 中移除 "Other"
#ct.imm <- setdiff(all_cell_types, "Other")
ct.imm <- all_cell_types
length(ct.imm)
# 修改 stat.test 格式，使其包含正确的列名和数据
stat.test <- stat.test %>%
  dplyr::rename(group1 = Group1, group2 = Group2, p.adj = p.adj)  # 重命名列
head(stat.test)

p_prop <- prop.plot(input, 
                    title = " ", 
                    color = c("#B09C85FF", "#00A087FF", "#800080"),
                    hide.ns = T, 
                    method = "fdr", 
                    legendLocation = "bottom", 
                    stat.test = stat.test, 
                    ct = ct.imm, 
                    strip.size = 6, 
                    ncol = 4, 
                    title.y = "Proportion", 
                    format = TRUE)
p_prop
ggsave(filename = paste0(output,"/Proportion_all.png"),
       plot = p_prop, width = 6, height = 5, dpi = 300)

# 5. PTGR1表达量组间差异 ------------------------------------------------------------
##绘图 参考原文格式
##不删为0的 不标差异线
#cell_type2
#label
# for (gene in vital_genes) {
#   data_gene <- FetchData(ss, vars = c("anno", gene, "label"))
#   # # 过滤掉基因表达为0的细胞
#   # data_gene <- data_gene %>%
#   #   filter(!!sym(gene) > 0)
#   p <- ggplot(data_gene, 
#               aes(x = label, y = !!sym(gene), fill = label)) +
#     
#     geom_violin(trim = FALSE, size = 0.8, alpha = 0.7) +
#     geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 0.6, 
#                aes(color = label)) +
#     geom_boxplot(
#       position = position_dodge(width = 0.75),  # 调整箱子位置
#       width = 0.2,  # 控制箱体宽度
#       fill = "white",  # 白色箱体
#       color = "black",  # 箱体边框
#       alpha = 0.6,  # 设置透明度
#       outlier.shape = NA  # 移除箱外点
#     ) +
#     facet_wrap(~anno, scales = "free", ncol = 4) +
#     theme_bw() +
#     labs(
#       x = "Group",
#       y = paste(gene,"Expression Level" )
#       #title = paste("Expression of", gene)
#     ) +
#     theme(
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.spacing.y = unit(0, "lines"),
#       #plot.title = element_blank(),
#       strip.text = element_text(size = 10, margin = margin(0, 0, 2, 0)),
#       strip.background = element_blank(),
#       legend.position = "bottom",  # 图例
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       axis.text.y = element_text(color = "black", size = 10),
#       axis.ticks.y = element_line(color = "black", size = 0.25),
#       axis.ticks.length.y = unit(0.05, "cm"),
#       axis.title.y = element_text(color = "black", size = 10),
#       axis.title.x = element_blank(),
#       plot.margin = margin(1, 1, 1, 1)
#     ) +
#     #scale_fill_manual(values = c("#B09C85FF", "#00A087FF", "#800080")) +
#     scale_fill_manual(values = c("Healthy" = "#57b1ab", 
#                                  "Inflamed UC" = "#9a72c7", 
#                                  "Non inflamed UC" ="#ffa74f" ),name=NULL)+    
#     scale_color_manual(values = c("Healthy" = "#57b1ab", 
#                                   "Inflamed UC" = "#9a72c7", 
#                                   "Non inflamed UC" ="#ffa74f" ),name=NULL)
#   
#   #使用 `geom_signif` 添加星号和横线
#   p <- p + geom_signif(
#     comparisons = list(
#       c("Healthy", "Non inflamed UC"),
#       c("Healthy", "Inflamed UC"),
#       c("Non inflamed UC", "Inflamed UC")
#     ),
#     map_signif_level = TRUE,
#     #annotations = stat_test$p.adj.signif,  # 显示显著性符号（***, **, *）
#     y_position = c(3.5, 4, 4.5),
#     color = "red",
#     size = 0.2,
#     textsize = 3
#   )
#   ggsave(filename = paste0(output, "/Fig4.e_Expr2_", gene, ".pdf"),
#          plot = p, width = 9, height = 6, dpi = 300)
#   
#   cat("Saved violin plot with significance for gene:", gene, "\n")
# }

################
for (gene in vital_genes) {
  data_gene <- FetchData(ss, vars = c("cell_type2", gene, "label"))
  
  # 按分面计算y轴最大值
  max_values <- data_gene %>%
    group_by(cell_type2) %>%
    summarise(y_max = max(!!sym(gene), na.rm = TRUE)) 
  
  # 动态计算统计检验结果
  stat_test <- compare_means(
    formula = as.formula(paste(gene, "~ label")),
    data = data_gene,
    method = "wilcox.test",
    group.by = "cell_type2"
  ) %>% 
    left_join(max_values, by = "cell_type2") %>%  # 合并分面最大值
    mutate(
      # 按分面设置基准高度
      y_base = y_max * 1.,
      # 不同比较组的高度偏移
      y.position = case_when(
        group1 == "Healthy" & group2 == "Non inflamed UC" ~ y_base,
        group1 == "Healthy" & group2 == "Inflamed UC" ~ y_base * 0.85,
        group1 == "Non inflamed UC" & group2 == "Inflamed UC" ~ y_base * 0.7
      )
    )
  
  # 筛选显著结果
  sig_test <- stat_test %>% 
    filter(p < 0.05) %>%
    group_by(cell_type2) %>%
    arrange(desc(y.position)) 
  
  # 主绘图代码
  p <- ggplot(data_gene, aes(x = label, y = !!sym(gene)),fill = label, color = label) +
    geom_violin(aes(fill = label), trim = FALSE, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.3), 
               size = 0.3, 
               alpha = 0.6,
               aes(color = label)) + 
    facet_wrap(~cell_type2, scales = "free_y", nrow = 3) +  # 仅y轴自由缩放
    
    # 精准添加显著性标记
    ggsignif::geom_signif(
      data = sig_test,
      aes(xmin = group1, xmax = group2,
          annotations = case_when(
            p < 0.001 ~ "***",
            p < 0.01 ~ "**",
            p < 0.05 ~ "*"
          ),
          y_position = y.position),
      manual = TRUE,
      vjust = -0.2,      # 垂直位置微调
      tip_length = 0.02, # 横线末端长度
      color = 'red',
      textsize = 3,
      size = 0.5
    ) +    
    theme_bw() +
    labs(
      x = "Group",
      y = paste(gene, "Expression")
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.y = unit(0, "lines"),
      strip.text = element_text(size = 10, margin = margin(0, 0, 2, 0)),
      strip.background = element_blank(),
      legend.position = "bottom",  # 图例
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(color = "black", size = 10),
      axis.ticks.y = element_line(color = "black", size = 0.25),
      axis.ticks.length.y = unit(0.05, "cm"),
      axis.title.y = element_text(color = "black", size = 10),
      axis.title.x = element_blank(),
      plot.margin = margin(1, 1, 1, 1)
    ) +
    scale_fill_manual(values = c("Healthy" = "#57b1ab", 
                                 "Non inflamed UC" = "#ffa74f", 
                                 "Inflamed UC" = "#9a72c7"),name = NULL) +
    scale_color_manual(values = c("Healthy" = "#57b1ab", 
                                  "Non inflamed UC" = "#ffa74f", 
                                  "Inflamed UC" = "#9a72c7"),name = NULL)
  ggsave(
    filename = paste0(output, "/Fig4.Expr1_", gene, ".pdf"),
    plot = p, 
    width = 9,  # 适当增加宽度
    height = 6, # 动态高度
    dpi = 300
  )
}

# "Healthy" = "#57b1ab",   
# "Non inflamed UC" = "#ffa74f", 
# "Inflamed UC" = "#9a72c7"
# 6. 拟时序分析 ----------------------------------------------------------------
# 6.1 data reading -------------------------------------------------------------
library(monocle3)
library(cowplot)

rm(list = ls())
gc()
path <- getwd()
setwd(path)

input13 <- "13_sc_cluster/Epi/UC"
output <- "15_Epi/Monocle"

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

ss <- readRDS(file = paste0(input13, "/ss_anno3_Epi.rds"))
table(ss$cell_type)
table(ss$histo3)

# 修改 label 列的因子顺序
ss$label <- factor(ss$label, levels = c("healthy", "non inflamed UC", "inflamed UC"))
# 修改 label 列的因子顺序
ss$histo3 <- factor(ss$histo3, levels = c("Healthy (N=14)", "Non inflamed UC (N=13)", "Inflamed UC (N=23)"))

table(ss$histo3)
table(ss$label)

type_colors <- c(
  "TA" = '#c8c7e1',
  "Early Colonocytes" = '#c4daec',
  "Intermediate Colonocytes" = '#3ca0cf',
  "Mature Colonocytes" = '#405993',
  "LND" = '#9a70a8',
  
  "BEST4/OTOP2" = '#d69971',   
  "Enteroendocrine" = '#64a776', 
  "Goblet" = '#df5734',
  
  
  "Tuft" = '#d25774',
  "undifferented cell" = '#e6e2a3',
  "Stem" = '#696a6c'
)
# 6.1. 拟时序分析 -----------------------------------------------------------------
set.seed(42)
# 去批次，同上一个任务
ss <- RunUMAP(ss, dims = 1:20,reduction = "harmony", min.dist = 1, spread =7)

##创建CDS对象
data <- GetAssayData(ss, assay = 'RNA',layer = 'counts')
cell_metadata <- ss@meta.data 
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

##预处理：标准化和降维
cds <- preprocess_cds(cds,num_dim = 50)
monocle3::plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds)
#cds <- reduce_dimension(cds,preprocess_method = "PCA") #'PCA' or 'LSI'
plot_cells(cds)

colnames(colData(cds))
##以之前的seurat分群来添加颜色，和原来的分群对比
p1 <- plot_cells(cds,reduction_method = "UMAP",color_cells_by = "seurat_clusters")+
  ggtitle('cds.umap')
  #+scale_color_manual(values = type_colors)  # 指定颜色
p1
p1.1 <- plot_cells(cds,reduction_method = "UMAP",color_cells_by = "batch")+
  ggtitle('cds.umap')
#+scale_color_manual(values = type_colors)  # 指定颜色
p1.1

# 去批次 harmony. 
#
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(ss,reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds,reduction_method = "UMAP",color_cells_by = "seurat_clusters")+
  ggtitle('int.umap')
p2
p3 = p1|p2
#p3是monocle3重新降维的结果，p4是seurat降维的结果
p3
ggsave(paste0(output,"/umap_reduction_compare.png"),
       plot = p3,
       width = 10,
       height = 5)
##检查批次效应

p4 <- plot_cells(cds, 
                 color_cells_by="batch", 
                 label_cell_groups=FALSE)+
  ggtitle('before')
p4

###### already done by harmony. 
#align_cds() remove batch effects using mutual nearest neighbor alignment
# cds <- align_cds(cds, num_dim = 100, alignment_group = "batch")
# cds <- reduce_dimension(cds)
# p5 <- plot_cells(cds, 
#                  color_cells_by="batch", 
#                  label_cell_groups=FALSE)+
#   ggtitle('after')
# p5
# 
# p = p4|p5
# p
# ggsave(paste0(output,"/umap_remove_batch_effects.png"),
#        plot = p,
#        width = 10,
#        height = 5)
# colnames(colData(cds))

#cds <- reduce_dimension(cds, reduction_method="tSNE")
#plot_cells(cds, reduction_method="tSNE", color_cells_by="cell_type")
plot_cells(cds, reduction_method="UMAP", color_cells_by="cell_type")

#可视化指定基因
vital_gene <- c("PTGR1")
p6 <- plot_cells(cds,
                 genes = vital_gene,
                 #reduction_method="tSNE",
                 label_cell_groups = F,
                 group_cells_by="group",
                 show_trajectory_graph = F)+
  facet_wrap(~ group, ncol = 2)  # 按组分面，设置列数
p6
ggsave(paste0(output,"/umap_vital_expr.png"),
       plot = p6,
       width = 5,
       height = 5)

# ~6.2注释 ---------------------------------------------------------------------
##用之前的注释结果
p7 <- plot_cells(cds,reduction_method="UMAP", color_cells_by="cell_type")+
  scale_color_manual(values = type_colors)  # 应用自定义颜色
p7
ggsave(paste0(output,"/Fig4.a1_umap_seurat_anno.pdf"),
       plot = p7,
       width = 5,
       height = 5)

#Part 1. all cells Trajectory ------------------------------------------------------------
output_folder <- "15_Epi/Monocle/Trajectory_expression"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}
# ~ 6.2.1 拟合轨迹 --------------------------------------------------------------
cds2 <- cluster_cells(cds)
plot_cells(cds2,color_cells_by = "partition",label_cell_groups = F)

cds2 <- learn_graph(cds2)
p9 <- plot_cells(cds2,
                 color_cells_by = "partition",
                 label_groups_by_cluster = F,
                 label_leaves = T,
                 label_branch_points = T)
p9

ggsave(paste0(output,"/Trajectory.png"),
       plot = p9,
       width = 6,
       height = 6)

# ~ 6.3 确定起点 --------------------------------------------------------------
#手动选择起点
#cds2 <- order_cells(cds2)

# 获取“Stem”细胞作为根节点
root_cells <- rownames(pData(cds2))[pData(cds2)$cell_type == "Stem"]
cds2 <- order_cells(cds2, root_cells = root_cells)

p10 <- plot_cells(cds2,
                  color_cells_by = "pseudotime",
                  label_groups_by_cluster = T,
                  #trajectory_graph_segment_size = 2,
                  #labels_per_group=1,
                  label_cell_groups = F,
                  label_roots = F,
                  label_leaves = F,
                  label_branch_points = F)
p10
ggsave(paste0(output,"/Fig4.f_Trajectory_Pseudotime.pdf"),
       plot = p10,
       width = 8,
       height = 6)

# ~ 6.4 拟时序差异分析 -----------------------------------------------------------
##差异表达分析
# Track_genes <- graph_test(cds2,neighbor_graph = "principal_graph",cores = 6)   ### take time. 
# Track_genesWhole <- Track_genes
# Track_genes <- Track_genes %>%
#   filter(q_value < 1e-3)
# write.csv(Track_genes,
#           paste0(output,"/Trajectory_DEgenes.csv"),
#           row.names = F)

Track_genes_restored <- read.csv(
  file = paste0(output, "/Trajectory_DEgenes.csv"),
  header = TRUE,
  stringsAsFactors = FALSE   # 防止字符自动转为因子
) %>%
  # 可选：转换为tibble格式（若原始数据为tibble）
  # tibble::as_tibble()  
  # 恢复列类型（根据原始Track_genes结构调整）
  mutate(
    q_value = as.numeric(q_value),          # 确保数值列类型
    morans_I = as.numeric(morans_I),        # 若存在其他统计量列
    gene_short_name = as.character(gene_short_name)  # 确保基因名为字符型
  )
Track_genes <- Track_genes_restored

##标记核糖体和线粒体基因
ribosomal_genes <- grepl("^RPL|^RPS", rownames(Track_genes))
mitochondrial_genes <- grepl("^MT-", rownames(Track_genes))
# 合并，找到所有需要移除的基因
genes_to_remove <- ribosomal_genes | mitochondrial_genes
# 从 Track_genes 中移除这些基因
Track_genes_filtered <- Track_genes[!genes_to_remove, ]
any(grepl("^RPL|^RPS|^MT-", rownames(Track_genes_filtered), ignore.case = TRUE))

# 保存过滤后的结果
# write.csv(Track_genes_filtered,
#           paste0(output, "/Trajectory_DEgenes_filtered.csv"),
#           row.names = FALSE)

##挑选前n展示
# Track_genes_sig <- Track_genes_filtered %>%
#   top_n(n=50,morans_I) %>%
#   pull(gene_short_name) %>%
#   as.character()

##挑选方式二
Track_genes_sig <- Track_genes_filtered %>%
  filter(morans_I > 0.05) %>%     ####0.25  # 0.05 PTGR1
  pull(gene_short_name) %>%   
  as.character()              
#Track_genes_sig <- c(Track_genes_sig, "PTGR1")   ##### add PTGR1... 
##Moran's I > 0.25 基因在空间上有一定的正自相关性
##相邻的细胞中它们的表达模式较为相似。
##更高的 Moran's I 值表示基因在细胞群体中具有更强的聚集性
##表达相似的细胞更可能在空间上靠近

# ~ 6.5 结果可视化 -------------------------------------------------------------
###a.热图
library(dplyr)
library(ClusterGVis)
#https://blog.csdn.net/weixin_45822007/article/details/128566997
library(ComplexHeatmap)
library(factoextra) 
library(clusterProfiler)
library(NbClust)
library(Biobase)
library(TCseq)
library(Mfuzz)
library(circlize)

# 检查 Track_genes_sig 中哪些基因存在于 cds2 的行名中
valid_genes_in_cds <- Track_genes_sig %in% rownames(cds2)
# 过滤掉无效基因
Track_genes_sig_filtered <- Track_genes_sig[valid_genes_in_cds]
# 输出无效基因列表（可选）
if (sum(!valid_genes_in_cds) > 0) {
  invalid_genes <- Track_genes_sig[!valid_genes_in_cds]
  message("以下基因在 cds2 中不存在，将被过滤：\n", paste(invalid_genes, collapse = "\n"))
}

#列是按拟时序排列的细胞，行是筛选的基因
mat <- pre_pseudotime_matrix(cds_obj = cds2,
                             gene_list = Track_genes_sig_filtered)

# 检查基因符号是否有效
valid_genes <- rownames(mat) %in% keys(org.Hs.eg.db, keytype = "SYMBOL")
invalid_genes <- rownames(mat)[!valid_genes]

# 删去无效的
mat <- mat[valid_genes, ]
mat <- na.omit(mat) 
# check
#head(mat[1:5,1:5])
#dim(mat)

getClusters <- function(exp) {
  exp <- as.matrix(exp)
  factoextra::fviz_nbclust(exp, FUN = kmeans, method = "wss")
}
getClusters(exp = mat)
#px

#### 写入mat
write.csv(mat,paste0(output,"/cds2_mat.csv"),row.names = T)
#### 读取mat

mat <- as.data.frame(mat)
ck <- clusterData(obj = mat,
                  cluster.method = "kmeans",
                  cluster.num = 4) # 4 # 5

#saveRDS(ck, file = "cluster_results.rds")
loaded_ck <- readRDS("cluster_results.rds")  ### saving time 
ck <- loaded_ck

### 输出cluster4的基因列表 #### 
# write.csv(cluster_genes[[as.character(4)]] ,
#           paste0(output, "/genes_cluster4.csv"),
#           row.names = FALSE)

set.seed(123)
# 将PTGR1，手动添加mark基因
mark_genes <- c("PTGR1")   #### OR adding other genes? 
# add line annotation
pdf(paste0(output,'/Fig4.g1_heatmap_moran0.05.pdf'),height = 10,width = 8,onefile = F)
visCluster(object = ck,
           plot.type = "both",
           add.sampleanno = F,
           column_names_rot = 45,
           markGenes.side = "left",
           markGenes = mark_genes)
dev.off()

# ~ 6.6 结果可视化2 ------------------------------------------------------------
##富集分析
clusters <- unique(ck$wide.res$cluster)  # 获取所有 cluster
cluster_genes <- split(ck$wide.res$gene, ck$wide.res$cluster)  # 按 cluster 分组基因

# 创建一个列表存储富集分析结果
enrichment_results <- list()
# 对每个 cluster 做富集分析
for (cluster_id in clusters) {
  genes <- cluster_genes[[as.character(cluster_id)]]  # 当前 cluster 的基因
  
  # 富集分析
  ego <- enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,  
    keyType = "SYMBOL",
    ont = "all",  
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  # 保存结果
  enrichment_results[[as.character(cluster_id)]] <- data.frame(ego)
}

# 汇总富集结果
all_enrichments <- bind_rows(
  lapply(names(enrichment_results), function(cluster_id) {
    df <- enrichment_results[[cluster_id]]
    if (nrow(df) > 0) {
      df <- df %>% 
        dplyr::mutate(Cluster = cluster_id) %>% 
        dplyr::slice_head(n = 5)  # 选择每个 cluster 前 5 个显著结果
    }
    return(df)
  })
)

write.csv(all_enrichments,
          paste0(output, "/Trajectory_DEgenes_Enrichment.csv"),
          row.names = FALSE)

# enrich for clusters
enrich <- enrichCluster(object = ck,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        pvalueCutoff = 0.05,
                        topn = 5,
                        seed = 41)
head(enrich,3)

# 可视化富集结果与热图结合
pdf(paste0(output,'/Fig4.g2_moran0.05.pdf'),height = 9,width = 9,onefile = F)
visCluster(object = ck,
           plot.type = "both",
           add.sampleanno = F,
           markGenes = mark_genes,
           column_names_rot = 45,
           markGenes.side = "left",
           line.side = "left",
           annoTerm.data = enrich,
           go.col = rep(ggsci::pal_d3()(4),each = 5),
           go.size = "pval")
dev.off()

#---------------------------------------------------------------------------------------------------------------------------------------------------------



##挑选前n展示
# Track_genes_top <- Track_genes_filtered %>%
#   top_n(n=10,morans_I) %>%
#   pull(gene_short_name) %>%
#   as.character()
# output_folder <- paste0(output,"/Trajectory_expression")
# if (!dir.exists(output)) {
#   dir.create(output_folder, recursive = TRUE)
# }
# for (gene in Track_genes_top) {
#   p <- plot_cells(cds2,
#                   genes = gene,
#                   show_trajectory_graph = TRUE, 
#                   label_cell_groups = FALSE,
#                   label_leaves = FALSE,
#                   label_branch_points = FALSE) +
#     ggtitle(paste("Expression of", gene))
#   ggsave(paste0(output_folder, "/",gene,".png"), plot = p, width = 6, height = 5)
# }
# 
# p_all <- plot_cells(cds2,
#                     genes = Track_genes_top,
#                     show_trajectory_graph = TRUE,
#                     label_cell_groups = TRUE,
#                     label_leaves = TRUE,
#                     label_branch_points = TRUE) +
#   ggtitle("Highlighted Trajectory Genes")
# 
# p_all 
# ggsave(paste0(output, "/Highlighted Trajectory Genes.png"), 
#        plot = p_all, width = 10, height = 10)

##c.基因表达趋势图
library(ggplot2)
library(tidyr)
library(dplyr)

stem_markers <- c("LGR5","OLFM4","SOX9")
MCN_markers <- c("PGC","MUC6","AQP5")
gene2 <- c("DACH1","EGR1","FOSB","JUNB","ETS1","HES1","PITX2",
           "SPDEF","NR4A1","NME2","ATF3","ATF4","CREB3L1","CREB3L2")
vital_gene <- c("PTGR1")

genes <- c(stem_markers, MCN_markers,gene2, vital_gene)

#genes <- Track_genes_top
# 提取伪时序信息
pseudotime <- as.numeric(pseudotime(cds2))

# 提取显著基因的表达数据
#gene_expr_data <- as.data.frame(t(exprs(cds2)[Track_genes_top, ]))
gene_expr_data <- as.data.frame(t(exprs(cds2)[genes, ]))
gene_expr_data$pseudotime <- pseudotime

# 将数据转换为长格式
gene_expr_data_long <- gene_expr_data %>%
  gather(key = "gene", value = "expression", -pseudotime)

# 绘制趋势图
ggplot(gene_expr_data_long, aes(x = pseudotime, y = expression, color = gene)) +
  geom_line(size = 1) +
  theme_minimal() +
  labs(title = "Gene Expression Trends Over Pseudotime", x = "Pseudotime", y = "Expression") +
  theme(legend.position = "right")

ggsave(paste0(output, "/Gene Expression Trends Over Pseudotime.png"), 
       width = 8, height = 6)


##d.genes_in_pseudotime
for (gene in genes) {
  cds_subset <- cds2[rowData(cds2)$gene_short_name == gene, ]
  p <- monocle3::plot_genes_in_pseudotime( 
    cds_subset,
    min_expr = NULL,          # 可根据需要设置最小表达值过滤
    cell_size = 0.75,         # 点的大小
    panel_order = NULL,       # 面板顺序，默认按基因顺序
    color_cells_by = "pseudotime", # 根据伪时间着色
    trend_formula = "~ splines::ns(pseudotime, df=3)", # 趋势线公式
    label_by_short_name = TRUE, # 是否使用基因短名称作为标签
    vertical_jitter = NULL,   # 垂直抖动，默认无
    horizontal_jitter = NULL  # 水平抖动，默认无
  )
  
  ggsave(paste0(output_folder, "/",gene,".png"), 
         plot = p,
         width = 4, height = 2)
  
}

##看不同分化路径的表达结果
unique(pData(cds2)$cell_type)
for (gene in genes) {
  
  # 筛选包含目标基因的行
  cds_subset <- cds2[rowData(cds2)$gene_short_name == gene, ]
  cell_types <- pData(cds_subset)$cell_type
  
  # 获取每个路径对应的细胞
  path_LND <- cds_subset[, cell_types == "LND"]
  path_Early_Colonocytes <- cds_subset[, cell_types == "Early Colonocytes"]
  path_Goblet <- cds_subset[, cell_types == "Goblet"]
  
  # 定义一个列表，包含不同路径的子集
  path_list <- list("LND" = path_LND, "Early Colonocytes" = path_Early_Colonocytes, "Goblet" = path_Goblet)
  
  # 循环遍历每条路径，生成该路径下基因表达的趋势图
  for (path_name in names(path_list)) {
    
    # 提取当前路径的细胞子集
    cds_path <- path_list[[path_name]]
    
    # 生成基因在当前路径上的表达趋势图
    p <- monocle3::plot_genes_in_pseudotime( 
      cds_path,
      min_expr = NULL,                # 可根据需要设置最小表达值过滤
      cell_size = 0.75,               # 点的大小
      panel_order = NULL,             # 面板顺序，默认按基因顺序
      color_cells_by = "pseudotime",  # 根据伪时间着色
      trend_formula = "~ splines::ns(pseudotime, df=3)",  # 趋势线公式
      label_by_short_name = TRUE,     # 使用基因短名称作为标签
      vertical_jitter = NULL,         # 垂直抖动
      horizontal_jitter = NULL        # 水平抖动
    )
    
    ggsave(paste0(output_folder, "/", gene, "_", path_name, ".png"), 
           plot = p,
           width = 4, height = 2)
    
    #cat("Saved plot for gene:", gene, "in path:", path_name, "\n")
  }
}

# 7. 特定轨迹 ------------------------------------------------------------
# ~ 7.1 拟合轨迹 --------------------------------------------------------------
cds3 <- cluster_cells(cds)
plot_cells(cds3,color_cells_by = "partition",label_cell_groups = F)

cds3 <- learn_graph(cds3)

p9 <- plot_cells(cds3,
                 color_cells_by = "partition",
                 label_groups_by_cluster = F,
                 label_leaves = T,
                 label_branch_points = T)
p9

# ggsave(paste0(output,"/Trajectory.png"),
#        plot = p9,
#        width = 6,
#        height = 6)

# ~ 7.2 确定起点 --------------------------------------------------------------
#手动选择起点
#cds3 <- order_cells(cds3)

# 获取“Stem”细胞作为根节点
root_cells <- rownames(pData(cds3))[pData(cds3)$cell_type == "Stem"]
cds3 <- order_cells(cds3, root_cells = root_cells)

p10 <- plot_cells(cds3,
                  color_cells_by = "pseudotime",
                  label_groups_by_cluster = T,
                  label_leaves = F,
                  label_branch_points = F)
p10
# ggsave(paste0(output,"/Trajectory_Pseudotime.png"),
#        plot = p10,
#        width = 8,
#        height = 6)
