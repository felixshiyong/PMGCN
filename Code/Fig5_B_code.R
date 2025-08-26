
##提取 重新降维的B细胞
##P和UC都是三类：Cycling B  ,Memory B, Follicular B
##UC 只合并 Memory B

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(ggsci)
library(patchwork)
library(ggrepel)
library(ggpubr)
library(cowplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(reshape2)
library(rstatix)

##ClusterGVis画图需要程辑包‘monocle’
#在运行分析前，验证某函数来自 monocle3
#getAnywhere(preprocess_cds)

# 1. read data  ----------------------------------------------------------
rm(list = ls())
gc()
path <- getwd()
setwd(path)

input13_UC <- "13_sc_cluster/B/UC"
input13_P <- "13_sc_cluster/B/P"

output <- "14_B"
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

##UC
ss_UC <- readRDS(file = paste0(input13_UC, "/ss_anno_B.rds"))
ss_UC$label <- recode(ss_UC$label,
                      `healthy` = "Healthy",
                      `inflamed UC` = "Inflamed UC",
                      `non inflamed UC` = "Non inflamed UC")
# 修改 label 列的因子顺序
ss_UC$label <- factor(ss_UC$label, 
                      levels = c("Healthy", "Non inflamed UC", "Inflamed UC"))
# 修改 label 列的因子顺序
ss_UC$histo3 <- factor(ss_UC$histo3, levels = c("Healthy (N=14)", "Non inflamed UC (N=13)", "Inflamed UC (N=23)"))
ss_UC$cell_type <- recode(ss_UC$cell_type,
                          `IgA Memory B` = "Memory B",
                          `IgA/IgG Memory B` = "Memory B")
ss_UC$cell_type <- factor(ss_UC$cell_type, levels = c("Cycling B", "Follicular B", "Memory B"))

table(ss_UC$histo3)
table(ss_UC$label)
table(ss_UC$group)
table(ss_UC$cell_type)

##P
ss_P <- readRDS(file = paste0(input13_P, "/ss_anno_B.rds"))
table(ss_P$label)
ss_P$label <- recode(ss_P$label,  
                     `Mild Periodontitis` = "Mild P",
                     `Severe Periodontitis` = "Severe P")
table(ss_P$label)
table(ss_P$cell_type)
ss_P$cell_type <- recode(ss_P$cell_type,
                         #`IgA Memory B` = "Memory B",
                         `IgA/IgG Memory B` = "Memory B")
ss_P$cell_type <- factor(ss_P$cell_type, levels = c("Cycling B", "Follicular B", "Memory B"))

ss_P$histo3 <- factor(
  ss_P$label,
  levels = c("Healthy", "Mild P", "Severe P"), # 原始分组名称
  labels = c("Healthy (N=2)", "Mild P (N=1)", "Severe P (N=1)") # 新标签
)
table(ss_P$histo3)

#合并 "data" 层
data <- NULL
for (layer in c("data.1", "data.2", "data.3", "data.4")) {
  layer_data <- GetAssayData(ss_P[['RNA']], layer = layer)
  data <- cbind(data, layer_data)
}
# 将合并的数据存储为一个新的层
ss_P[['data']] <- CreateAssayObject(counts = data)

# 2.整体可视化 ----------------------------------------------------------------
# ~UC overall ---------------------------------------------------------------------
UC_colors <- c(
  "Healthy" = "#B09C85FF",
  "Non inflamed UC" = "#00A087FF",
  "Inflamed UC" = "#800080"
)

Type_colors_UC <- c(
  "Follicular B" = '#d69971',
  #"Memory B" = '#df5734',   
  "Cycling B" = '#61bada', 
  "Memory B" = '#9a70a8')
#"IgA/IgG Memory B" = '#64a776'

##umap
#ss_UC <- RunUMAP(ss_UC, dims = 1:20,reduction = "harmony", min.dist = 1, spread =10)
p_umap <- DimPlot(ss_UC, 
                  reduction = "umap", 
                  group.by = "cell_type", 
                  label = F, 
                  repel = TRUE) +
  scale_color_manual(values = Type_colors_UC)+  # 指定颜色
  theme(
    panel.border = element_blank(),  # 去除边框
    panel.grid = element_blank(),    # 去除网格
    axis.text = element_blank(),     # 去除坐标轴文本
    axis.ticks = element_blank(),     # 去除坐标轴刻度
    aspect.ratio = 1
  ) +
  ggtitle("B")  
p_umap
ggsave(filename = paste0(output,"/","FigS5.a_UMAP_UC_Type.pdf"), 
       plot = p_umap, 
       width = 8, 
       height = 6)

## annotation 
B_cell_markers <- list(
  "Follicular B" = c( "IGHD","TCL1A", "FCER2"),
  
  "Memory B" = c("CD27","CD69","CD44"),
  
  "Cycling B" = c("MKI67","TOP2A") 
  
  #"IgA Memory B" = c("IGHA1", "IGHA2"),
  #"IgA/IgG Memory B" = c("IGHA1", "IGHA2","IGHG1", "IGHG2", "IGHG3", "IGHG4")
)

all_marker <- unique(unlist(B_cell_markers))
table(ss_UC$cell_type)

##按顺序画点图
dot_data <- DotPlot(ss_UC,
                    features = unique(unlist(B_cell_markers)),
                    group.by = "cell_type",
                    cluster.idents = T,
                    scale = F)$data

# 获取 y 轴（细胞类型）的顺序
cell_type_order <- levels(factor(dot_data$id))
# 按细胞类型顺序整理 marker_genes
ordered_markers <- unlist(B_cell_markers[cell_type_order])

# 补充没有的基因名称到完整序列中
all_markers <- unique(ordered_markers)
#all_markers <- unique(c(all_markers,"EPCAM"))

# 重新绘制点图
p_filtered_marker <- DotPlot(ss_UC,
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
  guides(size = guide_legend(title = "Percent Expressed", order = 1),
         color = guide_colorbar(title = "Average Expression", order = 2)) + # 添加图例
  scale_size(range = c(0, 10)) +
  coord_flip()

p_filtered_marker
ggsave(paste0(output,'/FigS5.b_Dot_UC_marker.pdf'), 
       plot = p_filtered_marker, height = 6, width =4.5)

# 计算每个疾病状态中每个细胞类型的比例
meta <- ss_UC@meta.data 
# 计算每个疾病状态中每个细胞类型的比例
cell_proportions <- meta %>%
  group_by(label, cell_type) %>%
  summarise(count = n()) %>%
  group_by(label) %>%  # 按疾病状态重新分组
  mutate(proportion = count / sum(count)) %>%  # 计算每个疾病状态内细胞类型的比例
  arrange(label, cell_type)  # 排序

# 绘制比例条形图
p_proportions2 <- ggplot(cell_proportions, aes(x = label, y = proportion, fill = cell_type)) + 
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
  scale_fill_manual(values = Type_colors_UC)  +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05))

p_proportions2
# 保存比例图
ggsave(filename = paste0(output, "/FigS5.c_Cell_Proportions_Stacked_UC1.pdf"),
       plot = p_proportions2, width = 3.5, height = 6, dpi = 300)

vital_genes <- c("FOSB")
# ~P overall---------------------------------------------------------------------
table(ss_P$cell_type)
Type_colors_P <- c(
  "Follicular B" = '#d69971',
  "Cycling B" = '#61bada', 
  "Memory B" = '#9a70a8'
)
table(ss_P$group)
# P_colors <- c(
#   "Healthy" = "#B09C85FF",   
#   "Mild P" = "#00A087FF",
#   "Severe P" = "#800080"
# )
P_colors <- c(
  "Healthy" = "#20854EFF",   
  "Periodontitis" = "#E18727FF" 
)
##umap
#ss_P <- RunUMAP(ss_P, dims = 1:20,reduction = "harmony", min.dist = 1, spread =10)
p_umap <- DimPlot(ss_P, 
                  reduction = "umap", 
                  group.by = "cell_type", 
                  label = F, 
                  repel = TRUE) +
  scale_color_manual(values = Type_colors_P)+  # 指定颜色
  theme(
    panel.border = element_blank(),  # 去除边框
    panel.grid = element_blank(),    # 去除网格
    axis.text = element_blank(),     # 去除坐标轴文本
    axis.ticks = element_blank(),     # 去除坐标轴刻度
    aspect.ratio = 1
  ) +
  ggtitle("B")  
p_umap
ggsave(filename = paste0(output,"/","FigS5.f_UMAP_P_Type.pdf"), 
       plot = p_umap, 
       width = 8, 
       height = 6)

##注释结果
table(ss_P$cell_type)
# B_cell_markers <- list(
#   "Follicular B" = c( "IGHD","TCL1A", "FCER2"),
#   "Cycling B" = c("MKI67"), 
#   "IgA/IgG Memory B" = c("IGHA1", "IGHA2","IGHG1", "IGHG2", "IGHG3", "IGHG4")
# )
B_cell_markers <- list(
  #"IgG Plasma B" = c( "IGHG1", "IGHG2", "IGHG3", "IGHG4"),
  "Follicular B" = c( "IGHD","TCL1A", "FCER2"),
  
  "Memory B" = c("CD27","CD69","CD44"),
  
  "Cycling B" = c("MKI67","TOP2A")
  #"IgA/IgG Memory B" = c("IGHA1", "IGHA2","IGHG1", "IGHG2", "IGHG3", "IGHG4")
  #"Memory B" = c("CD27","CD69","CD44","MS4A1", "CD37")
)

all_marker <- unique(unlist(B_cell_markers))
table(ss_P$cell_type)
##按顺序画点图
##不展示某一类细胞
##ss_sub <- subset(ss, subset = cell_type != "undifferented cell")

dot_data <- DotPlot(ss_P,
                    features = unique(unlist(B_cell_markers)),
                    group.by = "cell_type",
                    cluster.idents = T,
                    scale = F)$data
# 获取 y 轴（细胞类型）的顺序
#cell_type_order <- levels(factor(dot_data$id))
# 按细胞类型顺序整理 marker_genes
#ordered_markers <- unlist(B_cell_markers[cell_type_order])

# 补充没有的基因名称到完整序列中
all_markers <- unique(ordered_markers)
#all_markers <- unique(c(all_markers,"IGHA1", "IGHA2","IGHG1", "IGHG2", "IGHG3", "IGHG4"))

# 重新绘制点图
p_filtered_marker <- DotPlot(ss_P,
                             features = all_markers,
                             group.by = "cell_type",
                             #cluster.idents = T,
                             scale = F) +
  scale_color_gradient(low = "#FEE0D2", high = "#BC3C29FF") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(margin = margin(t = 0, r = 2, b = 0, l = 0)),  # 调整 y 轴标签间距
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加黑色边框
    panel.background = element_blank()  # 背景设为透明
  ) +
  guides(size = guide_legend(title = "Percent Expressed", order = 1),
         color = guide_colorbar(title = "Average Expression", order = 2)) + # 添加图例
  scale_size(range = c(0, 10)) +
  coord_flip()

p_filtered_marker
ggsave(paste0(output,'/FigS5.g_dot_P_Type.pdf'), 
       plot = p_filtered_marker, height = 6, width =4.5)

# 计算每个疾病状态中每个细胞类型的比例
meta <- ss_P@meta.data 
# 计算每个疾病状态中每个细胞类型的比例
cell_proportions <- meta %>%
  group_by(group, cell_type) %>%
  summarise(count = n()) %>%
  group_by(group) %>%  # 按疾病状态重新分组
  mutate(proportion = count / sum(count)) %>%  # 计算每个疾病状态内细胞类型的比例
  arrange(group, cell_type)  # 排序

# 绘制比例条形图
p_proportions <- ggplot(cell_proportions, aes(x = group, y = proportion, fill = cell_type)) + 
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
  scale_fill_manual(values = Type_colors_P)  +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05))
p_proportions
# 保存比例图
ggsave(filename = paste0(output, "/FigS5.h_cellproportion_P_Type.pdf"),
       plot = p_proportions, width = 3, height = 6, dpi = 300)

vital_genes <- c("FOSB")
# ~UC proportion function ----------------------------------------------------------------
#vital+ 数量/表达量组间差异
###
#This function generate figures with selected cell type.
flag <- 0
prop.plot <- function(input,title,color,hide.ns,method,legendLocation,stat.test,ct,strip.size=8,ncol=7,
                      plot.margin=margin(0,0,0,2),title.y,vjust=-0.3,format=F,mult=c(0, 0.15)){
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
    theme_bw()+labs(x=" ",y=title.y,title = title)+
    theme(
      panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
      panel.spacing.y = unit(0, "lines"),
      #plot.title = element_text(size=8,hjust = 0.5,face = "bold",margin = margin(0,0,0,0)),
      plot.title = element_blank(),
      strip.text = element_text(size=strip.size,margin = margin(0,0,2,0)),
      strip.background = element_blank(),
      legend.position = legendLocation,
      legend.text = element_text(colour = "black",size = 6),
      legend.title = element_blank(),legend.key.size = unit(0.5, "lines"),
      legend.margin = margin(0,0,0,0),legend.box.margin = margin(0,0,0,-5),legend.spacing.x = unit(0.1,"cm"),
      axis.text.x = element_blank(),axis.ticks.x = element_blank(),
      axis.text.y = element_text(color = "black",size = 5,margin = margin(r = 0)),
      axis.ticks.y = element_line(color = "black",size = 0.25),
      axis.ticks.length.y = unit(0.05,"cm"),
      axis.title.y = element_text(color = "black",size = 8),
      axis.title.x = element_blank(),plot.margin = plot.margin)+ 
    scale_color_manual(values=color)+
    scale_y_continuous(expand = expansion(mult = mult))+
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", color="black", width=0.1) +
    stat_summary(fun=mean, geom="point", color="black")
  
  p1 <- p1 + stat_pvalue_manual(stat.test, 
                                label = "p.adj.signif",
                                hide.ns = T,
                                tip.length = 0,
                                size = 1.5,
                                color = "red",
                                vjust = vjust)
  return(p1)
}

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

# ~ UC FOSB+--------------------------------------------------------------------
ss <- ss_UC
ss$NewCellType <- "Other" 

Idents(ss) <- "cell_type"
cell_types <- c("Cycling B", "Follicular B","Memory B")
# 遍历每个细胞类型并标记 FOSB 表达等于零的细胞
for (cell_type in cell_types) {
  cells_subset <- subset(ss, idents = cell_type)
  RNA_data <- GetAssayData(cells_subset, assay = "SCT", layer = "data")
  # 将 FOSB 基因表达> 0 的细胞标记为 "FOSB+ <细胞类型>"
  ss$NewCellType[Cells(cells_subset)[RNA_data["FOSB", ] > 0]] <- paste0("FOSB+ ", cell_type)
}

table(ss$NewCellType)
FOSB_table <- as.data.frame(table(ss$NewCellType))
write.csv(FOSB_table, file =paste0(output,"/FOSB+_CellType_Table_UC.csv") , row.names = FALSE)

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
#不要adj.p
stat.test <- stat.test %>%
  mutate(group = paste0(Group1, Group2)) %>%
  #mutate(p.adj = p.adjust(p_value, method = "fdr")) %>%
  mutate(p.adj = p_value) %>%
  group_by(group) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()

write.csv(stat.test, 
          file = paste0(output,"/Test_results_FOSB+_UC.csv"), 
          row.names = FALSE)
# 
all_cell_types <- unique(ss$NewCellType)
# 从 all_cell_types 中移除 "Other"
ct.imm <- setdiff(all_cell_types, "Other")
length(ct.imm)
ct.imm <- c("FOSB+ Cycling B", "FOSB+ Follicular B", "FOSB+ Memory B" )
# 修改 stat.test 格式，使其包含正确的列名和数据
stat.test <- stat.test %>%
  dplyr::rename(group1 = Group1, group2 = Group2, p.adj = p.adj)  # 重命名列
head(stat.test)

p_prop1 <- prop.plot(input, title = " ", 
                     color = c("#57b1ab", "#ffa74f", "#9a72c7"),
                     hide.ns = T, 
                     method = "fdr", 
                     legendLocation = "bottom", 
                     stat.test = stat.test, 
                     ct = ct.imm, 
                     strip.size = 6, 
                     ncol = 3, 
                     title.y = "Proportion", 
                     format = TRUE)
p_prop1
ggsave(filename = paste0(output,"/FigS5.d_Proportion_FOSB_UC.pdf"),
       plot = p_prop1, width = 4, height = 2, dpi = 300)

## 3. 绘图 FOSB+的组间差异
#############################################################################
table(ss$cell_type)
vital_genes <- c("FOSB")

################
for (gene in vital_genes) {
  data_gene <- FetchData(ss, vars = c("cell_type", gene, "label"))
  
  # 按分面计算y轴最大值
  max_values <- data_gene %>%
    group_by(cell_type) %>%
    summarise(y_max = max(!!sym(gene), na.rm = TRUE)) 
  
  # 动态计算统计检验结果
  stat_test <- compare_means(
    formula = as.formula(paste(gene, "~ label")),
    data = data_gene,
    method = "wilcox.test",
    group.by = "cell_type"
  ) %>% 
    left_join(max_values, by = "cell_type") %>%  # 合并分面最大值
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
    group_by(cell_type) %>%
    arrange(desc(y.position)) 
  
  # 主绘图代码
  p <- ggplot(data_gene, aes(x = label, y = !!sym(gene)),fill = label, color = label) +
    geom_violin(aes(fill = label), trim = FALSE, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    # geom_point(position = position_jitter(width = 0.3), 
    #            size = 0.3, 
    #            alpha = 0.6,
    #            aes(color = label)) + 
    facet_wrap(~cell_type, scales = "free_y", nrow = 1) +  # 仅y轴自由缩放
    
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
    filename = paste0(output, "/FigS5.Expr1_", gene, ".pdf"),
    plot = p, 
    width = 5,  # 适当增加宽度
    height = 4, # 动态高度
    dpi = 300
  )
}

# ~ P FOSB+--------------------------------------------------------------------
###
ss <- ss_P
ss$NewCellType <- "Other" 

table(ss$label)
ss$PatientID <- paste(ss$orig.ident, ss$batch, ss$label, sep = "_")

table(ss$histo3)
Idents(ss) <- "cell_type"
cell_types <- unique(ss$cell_type)

# 遍历每个细胞类型并标记 FOSB 表达等于零的细胞
for (cell_type in cell_types) {
  cells_subset <- subset(ss, idents = cell_type)
  DefaultLayer(cells_subset[['RNA']])
  RNA_data <- GetAssayData(cells_subset, assay = "data")
  # 将 FOSB 基因表达> 0 的细胞标记为 "FOSB+ <细胞类型>"
  ss$NewCellType[Cells(cells_subset)[RNA_data["FOSB", ] > 0]] <- paste0("FOSB+ ", cell_type)
}

table(ss$NewCellType)

FOSB_table <- as.data.frame(table(ss$NewCellType))
write.csv(FOSB_table, file =paste0(output,"/FOSB+_CellType_Table_P.csv") , row.names = FALSE)

table(ss$PatientID)
table(meta$group)

##绘图 FOSB+的组间差异
table(ss$cell_type)
vital_genes <- c("FOSB")
table(ss$group)

for (gene in vital_genes) {
  data_gene <- FetchData(ss, vars = c("cell_type", gene, "group"))
  
  # 按分面计算y轴最大值
  max_values <- data_gene %>%
    group_by(cell_type) %>%
    summarise(y_max = max(!!sym(gene), na.rm = TRUE)) 
  
  # 动态计算统计检验结果
  stat_test <- compare_means(
    formula = as.formula(paste(gene, "~ group")),
    data = data_gene,
    method = "wilcox.test",
    group.by = "cell_type"
  ) %>% 
    left_join(max_values, by = "cell_type") %>%  # 合并分面最大值
    mutate(
      # 按分面设置基准高度
      y_base = y_max * 1.,
      # 不同比较组的高度偏移
      y.position = case_when(
        group1 == "Healthy" & group2 == "Periodontitis" ~ y_base
      )
    )
  
  # 筛选显著结果
  sig_test <- stat_test %>% 
    filter(p < 0.05) %>%
    group_by(cell_type) %>%
    arrange(desc(y.position)) 
  
  # 主绘图代码
  p <- ggplot(data_gene, aes(x = group, y = !!sym(gene)),fill = group, color = group) +
    geom_violin(aes(fill = group), trim = FALSE, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    # geom_point(position = position_jitter(width = 0.3), 
    #            size = 0.3, 
    #            alpha = 0.6) + 
    facet_wrap(~cell_type, scales = "free_y", nrow = 1) +  # 仅y轴自由缩放
    
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
                                 "Periodontitis" = "#9a72c7" ),name = NULL)+
    scale_color_manual(values = c("Healthy" = "#57b1ab", 
                                  "Periodontitis" = "#9a72c7" ),name = NULL)
  ggsave(
    filename = paste0(output, "/FigS5.Expr1_", gene, "Periodontitis.pdf"),
    plot = p, 
    width = 5,  # 适当增加宽度
    height = 4, # 动态高度
    dpi = 300
  )
}

########
##  4.0 KEGG Functions ###
######## in FigExt5B code





