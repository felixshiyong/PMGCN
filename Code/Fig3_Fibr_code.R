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

#1.0 读取数据 ------------------------------------------------------------------
rm(list = ls())
gc()
path <- getwd()
setwd(path)

input13 <- "13_sc_cluster/fibro"
output <- "16_Fibro"

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

ss_UC <- readRDS(file = paste0(input13, "/ss_UC_fibro3.rds"))
ss_P <- readRDS(file = paste0(input13, "/ss_P_fibro2.rds"))

table(ss_UC$label)

ss_UC$label <- recode(ss_UC$label,
                      `healthy` = "Healthy",
                      `inflamed UC` = "Inflamed UC",
                      `non inflamed UC` = "Non inflamed UC")
ss_UC$label <- factor(ss_UC$label, levels = c("Healthy", "Non inflamed UC", "Inflamed UC"))
# 修改 label 列的因子顺序
ss_UC$histo3 <- factor(ss_UC$histo3, levels = c("Healthy (N=14)", "Non inflamed UC (N=13)", "Inflamed UC (N=23)"))
table(ss_UC$label)
table(ss_UC$group)
ss_UC$group <- recode(ss_UC$group,  
                      `HC` = "Healthy",
                      `UC` = "UC")
table(ss_UC$group)
table(ss_UC$label)
table(ss_UC$histo3)

###
table(ss_P$label)
ss_P$label <- recode(ss_P$label,  
                     `Mild Periodontitis` = "Mild P",
                     `Severe Periodontitis` = "Severe P")
table(ss_P$label)
data <- NULL
for (layer in c("data.1", "data.2", "data.3", "data.4")) {
  layer_data <- GetAssayData(ss_P[['RNA']], layer = layer)
  data <- cbind(data, layer_data)
}
ss_P[['data']] <- CreateAssayObject(counts = data)
table(ss_P$group)
ss_P$group <- recode(ss_P$group,  
                     `HC` = "Healthy",
                     `P` = "Periodontitis")
table(ss_P$group)

#2.0 整体可视化 ----------------------------------------------------------------

UC_colors <- c(
  "Healthy" = "#57b1ab",   
  "Non inflamed UC" = "#ffa74f",
  "Inflamed UC" = "#9a72c7"
)

# P_colors <- c(
#   "Healthy" = "#B09C85FF",   
#   "Severe P" = "#800080", 
#   "Mild P" = "#00A087FF"
# )
# P_colors <- c(
#   "HC" = "#8ecfb0",   
#   "P" = "#f99655"
# )

P_colors <- c(
  "Healthy" = "#8ecfb0",   
  "Periodontitis" = "#f99655"
)
##umap
ss_UC <- RunUMAP(ss_UC, dims = 1:20,reduction = "harmony", min.dist = 1, spread =10)

p_umap1 <- DimPlot(ss_UC, 
                   pt.size = 0.5,
                   alpha = 0.8,
                   reduction = "umap", 
                   group.by = "label", 
                   label = F, 
                   repel = TRUE) +
  scale_color_manual(values = UC_colors)+  # 指定颜色
  theme(
    panel.border = element_blank(),  # 去除边框
    panel.grid = element_blank(),    # 去除网格
    axis.text = element_blank(),     # 去除坐标轴文本
    axis.ticks = element_blank(),     # 去除坐标轴刻度
    aspect.ratio = 1
  ) +
  ggtitle("UC Fibroblast")  
p_umap1
ggsave(filename = paste0(output,"/","Fig3.a_UMAP_UC_label.pdf"), 
       plot = p_umap1, 
       width = 8, 
       height = 6)

p_umap2 <- DimPlot(ss_P, 
                   pt.size = 0.5,
                   alpha = 0.8,
                   reduction = "umap", 
                   group.by = "group", 
                   label = F, 
                   repel = TRUE) +
  scale_color_manual(values = P_colors)+  # 指定颜色
  theme(
    panel.border = element_blank(),  # 去除边框
    panel.grid = element_blank(),    # 去除网格
    axis.text = element_blank(),     # 去除坐标轴文本
    axis.ticks = element_blank(),     # 去除坐标轴刻度
    aspect.ratio = 1
  ) +
  ggtitle("Periodontitis Fibroblast")  
p_umap2
ggsave(filename = paste0(output,"/","Fig3.c2_UMAP_P_label.pdf"), 
       plot = p_umap2, 
       width = 8, 
       height = 6)

#3.0 dot plot actived marker ---------------------------------------------------
##dot plot  actived marker 
actived_markers <- c("IL11","IL24","OSM","OSMR","WNT2B","WNT5B",
                     "CXCL1","CXCL5","CXCL8","CCL2","CCL7","TNFSF14",
                     "CLU","IL33","IL6","LOX","GREM1","MMP3","IL1B","WNT5A","IL6")

p_UC_marker <- DotPlot(ss_UC,
                       features = unique(actived_markers),
                       group.by = "group") +
  scale_color_gradient(low = "#FEE0D2", high = "#BC3C29FF") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(margin = margin(t = 0, r = 2, b = 0, l = 0)),  # 调整 y 轴标签间距
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加黑色边框
    panel.background = element_blank()  # 背景设为透明
  ) +
  scale_size(range = c(0, 7)) +
  coord_flip()

p_UC_marker
ggsave(filename = paste0(output,"/","Fig3.b_DOT_UC_actived.pdf"), 
       plot = p_UC_marker, 
       width = 4, 
       height = 6)

p_P_marker <- DotPlot(ss_P,
                      features = unique(actived_markers),
                      group.by = "group",
                      cluster.idents = T) +
  scale_color_gradient(low = "#FEE0D2", high = "#BC3C29FF") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(margin = margin(t = 0, r = 2, b = 0, l = 0)),  # 调整 y 轴标签间距
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加黑色边框
    panel.background = element_blank()  # 背景设为透明
  ) +
  scale_size(range = c(0, 7)) +
  coord_flip()
p_P_marker
ggsave(filename = paste0(output,"/","Fig3.d_DOT_P_actived.pdf"), 
       plot = p_P_marker, 
       width = 4, 
       height = 6)

#4.0  CXCL5+ 数量差异可视化 -----------------------------------------------

ss_UC_new <- ss_UC
ss_UC_new$NewCellType <- "Other" 

#expr_data <- FetchData(ss, vars = vital_genes, layer='counts')
##expr_data <- FetchData(ss, vars = vital_genes)
# RNA_data2 <- GetAssayData(ss_UC_new, assay = "RNA", layer = "scale.data")#####
RNA_data <- GetAssayData(ss_UC_new, assay = "RNA", layer = "data")
# RNA_data <- GetAssayData(ss_UC_new, assay = "RNA", layer = "data")

# 标记CXCL5+ Fibro (表达值 > 0)
ss_UC_new$NewCellType[Cells(ss_UC_new)[RNA_data["CXCL5", ] > 0]] <- "CXCL5+ Fibro"

#tempexpre2 = RNA_data2["CXCL5", ]
#tempexpre2 = as.data.frame(tempexpre2)
# 标记CXCL5- Fibro (表达值 == 0)
# ss_UC_new$NewCellType[Cells(ss_UC_new)[RNA_data["CXCL5", ] == 0]] <- "CXCL5- Fibro"
table(ss_UC_new$NewCellType)
ss_UC_new$PatientID <- paste(ss_UC_new$orig.ident, ss_UC_new$batch, ss_UC_new$label, sep = "_")

# 1. 提取CXCL5+ Fibro 细胞的数据并计算其比例
meta <- ss_UC_new@meta.data
meta$CXCL5_Fibro <- ifelse(meta$NewCellType == "CXCL5+ Fibro", 1, 0)  # 创建新的列，标记CXCL5+ Fibro细胞
table(meta$CXCL5_Fibro )

# 按照PatientID和histo3分组，计算每个组中CXCL5+ Fibro细胞的比例
prop_CXCL5_Fibro <- meta %>%
  group_by(PatientID, histo3) %>%
  summarise(CXCL5_Fibro_Prop = mean(CXCL5_Fibro))  # 计算每个组内CXCL5+ Fibro的比例

# # 2. 按histo3分组进行Wilcoxon检验，比较不同组之间CXCL5+ Fibro的比例
# stat.test <- prop_CXCL5_Fibro %>%
#   group_by(histo3) %>%
#   do({
#     group_data <- .
#     disease_combinations <- combn(unique(group_data$histo3), 2, simplify = FALSE)
# 
#     # 对每对疾病状态进行Wilcoxon检验
#     results <- lapply(disease_combinations, function(groups) {
#       group1 <- groups[1]
#       group2 <- groups[2]
# 
#       test_data <- group_data %>% filter(histo3 %in% c(group1, group2))
#       test_result <- wilcox.test(CXCL5_Fibro_Prop ~ histo3, data = test_data)
# 
#       data.frame(
#         Group1 = group1,
#         Group2 = group2,
#         p_value = test_result$p.value
#       )
#     })
# 
#     # 合并结果
#     do.call(rbind, results)
#   }) %>%
#   ungroup()
# # 
# # # # 3. FDR 校正
# stat.test <- stat.test %>%
#    mutate(group = paste0(Group1, Group2)) %>%
#    mutate(p.adj = p.adjust(p_value, method = "fdr"))

# 4. 绘制比例图，使用 `ggplot2` 并添加显著性标记
p_prop <- ggplot(prop_CXCL5_Fibro, aes(x = histo3, y = CXCL5_Fibro_Prop, col = histo3)) +
  geom_jitter(shape = 16, size = 2) +  # 添加散点
  #facet_wrap(~histo3, scales = "free", ncol = 3) +  # 按组分面
  theme_bw() +
  labs(x = "Disease State (histo3)", y = expression(paste("CXCL5"^"+"," Fb Ratio")), title = "CXCL5+ Fibro Cells in Different Histological States") +
  theme(
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    panel.spacing.y = unit(0, "lines"),
    #plot.title = element_text(size=8,hjust = 0.5,face = "bold",margin = margin(0,0,0,0)),
    plot.title = element_blank(),
    strip.text = element_text(size = 12,margin = margin(0,0,2,0)),
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(colour = "black",size = 12),
    legend.title = element_blank(),
    legend.key.size = unit(0.5, "lines"),
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(0,0,0,-5),
    legend.spacing.x = unit(0.1,"cm"),
    axis.text.x = element_blank(),axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black",size = 10,margin = margin(r = 0)),
    axis.ticks.y = element_line(color = "black",size = 0.25),
    axis.ticks.length.y = unit(0.05,"cm"),
    axis.title.y = element_text(color = "black",size = 12),
    axis.title.x = element_blank(),
    plot.margin = margin(1, 1, 1, 1)
  ) +
  scale_color_manual(values = c("#57b1ab", "#ffa74f", "#9a72c7")) +  # 设定颜色
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", color = "black", width = 0.1) +  # 错误条
  stat_summary(fun = mean, geom = "point", color = "black")  # 显示平均点

# 添加统计显著性标记
p_prop <- p_prop + geom_signif(
  comparisons = list(
    c("Healthy (N=14)", "Non inflamed UC (N=13)")
    #c("Healthy (N=14)", "Inflamed UC (N=23)"),
    #c("Non inflamed UC (N=13)", "Inflamed UC (N=23)")
  ),
  map_signif_level = TRUE,
  y_position = c(0.28),
  color = "red",textsize = 3
)+ geom_signif(
  comparisons = list(
    #c("Healthy (N=14)", "Non inflamed UC (N=13)"),
    c("Healthy (N=14)", "Inflamed UC (N=23)")
    #c("Non inflamed UC (N=13)", "Inflamed UC (N=23)")
  ),
  map_signif_level = TRUE,
  y_position = c(0.33),
  color = "red",textsize = 3
)
#+ stat_pvalue_manual(stat.test, label = "p.adj", hide.ns = TRUE, tip.length = 0, size = 1.5, color = "red")

p_prop
ggsave(filename = paste0(output,"/", "Fig3.e2_CXCL5+_Fibro_Proportion.pdf"), 
       plot = p_prop, width = 6, height =3.5, dpi = 300)

#5.0  CXCL5+ 表达差异可视化/表达量的组间差异 UC---------------------------------------------------
data_gene <- FetchData(ss_UC_new, vars = c("anno", "CXCL5", "label"))

#过滤掉CXCL5表达为0的细胞
# data_gene <- data_gene %>%
#   filter(CXCL5 > 0)
# UC_colors <- c(
#   "Healthy" = "#57b1ab",   
#   "Non inflamed UC" = "#ffa74f",
#   "Inflamed UC" = "#9a72c7"
# )
# 
# # P_colors <- c(
# #   "Healthy" = "#B09C85FF",   
# #   "Severe P" = "#800080", 
# #   "Mild P" = "#00A087FF"
# # )
# # P_colors <- c(
# #   "HC" = "#8ecfb0",   
# #   "P" = "#f99655"
# # )
# 
# P_colors <- c(
#   "Healthy" = "#8ecfb0",   
#   "Periodontitis" = "#f99655"
# )

p <- ggplot(data_gene, 
            aes(x = label, y = CXCL5, fill = label)) +
  geom_violin(trim = FALSE, size = 0.8, alpha = 0.7) +
  geom_point(position = position_jitter(width = 0.2), size = 2, alpha = 0.6, 
             aes(color = label)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),  # 调整箱体位置
    width = 0.2,  # 控制箱体宽度
    fill = "white",  # 白色箱体
    color = "black",  # 箱体边框
    alpha = 0.6,  # 设置透明度
    outlier.shape = NA  # 移除箱外点
  ) +
  facet_wrap(~anno, scales = "free", ncol = 1) +
  theme_bw() +
  labs(
    x = "Group",
    y = "CXCL5  Expression",
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.y = unit(0, "lines"),
    strip.text = element_blank(), #element_text(size = 20, margin = margin(0, 0, 2, 0)),
    strip.background = element_blank(),
    legend.position = "bottom",  # 图例
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 12),
    axis.ticks.y = element_line(color = "black", size = 0.25),
    axis.ticks.length.y = unit(0.05, "cm"),
    axis.title.y = element_text(color = "black", size = 12),
    axis.title.x = element_blank(),
    plot.margin = margin(1, 1, 1, 1)
  ) +
  scale_fill_manual(values = c("Healthy" = "#57b1ab", 
                               "Inflamed UC" = "#9a72c7", 
                               "Non inflamed UC" = "#ffa74f"),name = NULL)+
  scale_color_manual(values = c("Healthy" = "#57b1ab", 
                                "Inflamed UC" = "#9a72c7", 
                                "Non inflamed UC" = "#ffa74f"),name = NULL)

# 使用 `geom_signif` 添加星号和横线
p <- p + geom_signif(
  comparisons = list(
    c("Non inflamed UC", "Inflamed UC")
  ),
  map_signif_level = TRUE,
  y_position = c(4.5),
  color = "red",textsize = 3
)+ geom_signif(
  comparisons = list(
    c("Healthy", "Inflamed UC")
  ),
  map_signif_level = TRUE,
  y_position = c(4.9),
  color = "red",textsize = 3
)+
  # 控制y轴范围（根据数据范围调整 y_min 和 y_max）
  coord_cartesian(ylim = c(0, 5.2)) 
p
# 保存图片
ggsave(filename = paste0(output, "/Fig3.f_Expr_CXCL5.pdf"),
       plot = p, width = 5.5, height =3, dpi = 300)

#6.0  differential expression and functional analysis ---------------------------------------------------
#  CXCL5+/CXCL5-细胞之间的差异分析 -------------------------------------------------
# ~ UC --------------------------------------------------------------------
RNA_data_UC <- GetAssayData(ss_UC, assay = "RNA", layer = "data")
# 按CXCL5基因表达分组
ss_UC$CXCL5_status <- ifelse(RNA_data_UC["CXCL5", ] > 0, "Positive", "Negative")
table(ss_UC$CXCL5_status)

# 获取CXCL5基因的表达量数据
CXCL5_expr <- RNA_data_UC["CXCL5", ]
ss_UC$CXCL5_expr <- CXCL5_expr
table(ss_UC$CXCL5_expr)

# 绘制CXCL5表达量分布图
p <- ggplot(ss_UC@meta.data, aes(x = "", y = CXCL5_expr)) +
  geom_violin(fill = "lightblue", color = "black", alpha = 0.6) +
  geom_jitter(aes(color = CXCL5_expr), size = 1, width = 0.2, alpha = 0.6) +
  theme_minimal() +
  labs(title = "CXCL5 Expression Distribution in Fibroblasts", 
       x = "", y = "Expression Level") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

p
# 保存图像
ggsave(filename = paste0(output, "/CXCL5_UC.png"),
       plot = p,
       width = 8,
       height = 6,
       dpi = 300)

# 进行差异表达分析（对Positive与Negative组进行比较）  #######take time... ######
markers <- FindMarkers(ss_UC, 
                       ident.1 = "Positive", 
                       ident.2 = "Negative", 
                       group.by = "CXCL5_status", 
                       test.use = "MAST",
                       min.pct = 0.2)
markers
# 保存差异表达分析结果
write.csv(markers, file = paste0(output, "/CXCL5+_UC.csv"))

# 生成和保存CXCL5状态的分组统计表
status_table <- table(ss_UC$CXCL5_status)
write.csv(as.data.frame(status_table), file = paste0(output, "/Table_CXCL5_UC.csv"))

status_table
prop.table(status_table)

markers$Significance <- with(markers, 
                             ifelse(p_val_adj < 0.05 & avg_log2FC > 1, "Upregulated",
                                    ifelse(p_val_adj < 0.05 & avg_log2FC < -1, "Downregulated", "Not Significant")))


# 筛选显著的差异基因
significant_markers <- markers[abs(markers$avg_log2FC) > 1 & markers$p_val_adj < 0.05, ]
write.csv(significant_markers,
          file = paste0(output,"/DE_CXCL5+_UC.csv"))

#6.1 找出在差异基因中的相关因子  #############--------------------------------
# 从 GO terms 获取相关基因
library(org.Hs.eg.db)
secreted_genes <- unique(AnnotationDbi::select(org.Hs.eg.db, 
                                               keys = "GO:0005615", # 分泌到胞外空间的基因
                                               columns = "SYMBOL", 
                                               keytype = "GOALL")$SYMBOL)
inflammatory_genes <- unique(AnnotationDbi::select(org.Hs.eg.db, 
                                                   keys = "GO:0006954", # 炎性反应
                                                   columns = "SYMBOL", 
                                                   keytype = "GOALL")$SYMBOL)
cytokine_genes <- unique(AnnotationDbi::select(org.Hs.eg.db, 
                                               keys = "GO:0005125", # 细胞因子活性
                                               columns = "SYMBOL", 
                                               keytype = "GOALL")$SYMBOL)
chemokine_genes <- unique(AnnotationDbi::select(org.Hs.eg.db, 
                                                keys = "GO:0008009", # 趋化因子活性
                                                columns = "SYMBOL", 
                                                keytype = "GOALL")$SYMBOL)

#合并
factor_genes <- unique(c(secreted_genes, inflammatory_genes, 
                         cytokine_genes, chemokine_genes))

write.csv(factor_genes,
          file = paste0(output,"/factor_genes.csv"))

###
# 为每组基因添加类别标签
secreted_genes <- data.frame(Gene = secreted_genes, Factor = "Secreted")
inflammatory_genes <- data.frame(Gene = inflammatory_genes, Factor = "Inflammatory")
cytokine_genes <- data.frame(Gene = cytokine_genes, Factor = "Cytokine")
chemokine_genes <- data.frame(Gene = chemokine_genes, Factor = "Chemokine")

# 合并所有因子基因
factor_genes <- rbind(secreted_genes, inflammatory_genes, cytokine_genes, chemokine_genes)

# 保存到CSV文件
write.csv(factor_genes,
          file = paste0(output, "/factor_genes.csv"),
          row.names = FALSE)
###
highlight_genes <- rownames(significant_markers) %in% factor_genes$Gene
highlighted_markers <- significant_markers[highlight_genes, ]
sum(highlight_genes)

highlight_secreted_genes <- rownames(significant_markers) %in% secreted_genes$Gene
highlight_secreted_markers <- significant_markers[highlight_secreted_genes, ]
sum(highlight_secreted_genes)

inflammatory_genes_genes <- rownames(significant_markers) %in% inflammatory_genes$Gene
inflammatory_genes_markers <- significant_markers[inflammatory_genes_genes, ]
sum(inflammatory_genes_genes)

cytokine_genes <- rownames(significant_markers) %in% cytokine_genes$Gene
cytokine_genes_markers <- significant_markers[cytokine_genes, ]
sum(cytokine_genes)

chemokine_genes <- rownames(significant_markers) %in% chemokine_genes$Gene
chemokine_genes_markers <- significant_markers[chemokine_genes, ]
sum(chemokine_genes)

# 将因子标签添加到 highlighted_markers 数据框中
#highlighted_markers$Factor <- factor_genes$Factor[match(rownames(highlighted_markers), factor_genes$Gene)]
highlighted_markers$Factor <- factor_genes$Factor[match(rownames(highlighted_markers), factor_genes)]  ### problem. 
write.csv(highlighted_markers,
          file = paste0(output,"/Factor_CXCL5+_UC.csv"))

highlighted_markers
# secreted_genes
# 绘制火山图并标出这些因子
# scale_fill_manual(values = c('Up' = '#E64B35', 'Down' = '#3C5488'))
p1 <- ggplot(markers, 
             aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = Significance), alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Upregulated" = "#BC3C29FF", "Downregulated" = "#0072B5FF", "Not Significant" = "grey")) +
  theme(
    legend.position = "none",
    #axis.text.y = element_text(margin = margin(t = 0, r = 2, b = 0, l = 0)),  # 调整 y 轴标签间距
    axis.text.x = element_text(size = 12, color = "black"),  # X轴刻度数字
    axis.text.y = element_text(size = 12, color = "black"),  # Y轴刻度数字
    axis.title.x = element_text(size = 15, color = "black"),  # X轴刻度数字
    axis.title.y = element_text(size = 15, color = "black"),  # Y轴刻度数字
    # 分面标签（若有分面）
    plot.title = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加黑色边框
    panel.background = element_blank()  # 背景设为透明
  )+
  labs(title = "UC", x = bquote(~Log[2]~"(FC)"), y = bquote("-"~Log[10]~"(adj.P)")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  coord_cartesian(ylim = c(0, 300)) +
  coord_cartesian(xlim = c(-10, 18)) +
  #bquote("-"~Log[10]~"(adj.P)")
  # 添加分泌型因子、炎性因子、细胞因子、趋化因子的标签
  # geom_text_repel(data = highlighted_markers,
  #                 aes(label = rownames(highlighted_markers)),
  #                 size = 3,
  #                 box.padding = 0.5,
  #                 #max.overlaps = 20,  # 增大此值（如 100 或 Inf）
  #                 point.padding = 0.5,
  #                 max.overlaps = 10)   #### certain genes would be omitted.
  # geom_text_repel(data = highlight_secreted_markers,
  #               aes(label = rownames(highlight_secreted_markers)),
  #               size = 3,
  #               box.padding = 0.5,
  #               point.padding = 0.5,
  #               max.overlaps = 10)
  geom_text_repel(data = inflammatory_genes_markers,
                aes(label = rownames(inflammatory_genes_markers)),
                size = 3,
                box.padding = 0.5,
                point.padding = 0.5,
                max.overlaps = 100)
  # geom_text_repel(data = cytokine_genes_markers, 
  #               aes(label = rownames(cytokine_genes_markers)), 
  #               size = 3, 
  #               box.padding = 0.5, 
  #               point.padding = 0.5, 
  #               max.overlaps = 10)
  # geom_text_repel(data = chemokine_genes_markers, 
  #               aes(label = rownames(chemokine_genes_markers)), 
  #               size = 3, 
  #               box.padding = 0.5, 
  #               point.padding = 0.5, 
  #               max.overlaps = 10)
p1

ggsave(filename = paste0(output,"/","Fig3.h_Volcano_CXCL5+_UC.pdf"), 
       plot = p1, 
       width = 6.5, 
       height = 6.5)

#6.2 富集分析  #############--------------------------------
##
library(clusterProfiler)
gene_list <- rownames(significant_markers)

# 将基因符号转换为Entrez ID
Genes <- bitr(gene_list,  
              fromType = "SYMBOL", #输入数据的类型
              toType = c("ENTREZID"), #要转换的数据类型
              OrgDb = org.Hs.eg.db) #物种

# 进行GO富集分析
GO_enrichment <- enrichGO(gene = Genes$ENTREZID,
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID",
                          ont = "all",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,   #p
                          qvalueCutoff = 0.05, #fdr P
                          minGSSize = 10,
                          maxGSSize = 500,
                          readable = TRUE)

View(as.data.frame(GO_enrichment))
write.table(data.frame(ID = row.names(GO_enrichment@result), GO_enrichment@result),
            file = paste0(output, "/GO_CXCL5+_UC",  ".txt"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

p2 <- dotplot(GO_enrichment,
              title = "CXCL5+: Fibroblast",
              showCategory = 7,
              font.size = 22,
              split="ONTOLOGY") + 
  #aes(color = pvalue)+
  facet_grid(ONTOLOGY~., scale='free') + 
  scale_color_gradientn(colors = c("#3C5488", "#E64B35"),name = "adj.p")
p2
ggsave(filename = paste0(output,"/GO_CXCL5+_UC.png"),
       plot = p2,
       width = 14,
       height = 16,
       dpi = 300 )

# 读取本地文件 #################   KEGG    ############ 
kegg_gene <- read.delim("kegg_pathway_gene.txt", header = FALSE, col.names = c("pathway", "gene"))
kegg_name <- read.delim("kegg_pathway_name.txt", header = FALSE, col.names = c("pathway", "name"))
# 格式清理：移除冗余前缀 "path:"
kegg_gene$pathway <- gsub("path:", "", kegg_gene$pathway)
kegg_name$pathway <- gsub("path:", "", kegg_name$pathway)
# 构建本地数据库
TERM2GENE <- kegg_gene[, c("pathway", "gene")]  # 通路ID → 基因ID
TERM2GENE$gene <- gsub("hsa:", "", TERM2GENE$gene)
TERM2NAME <- kegg_name                          # 通路ID → 通路名称

KEGG_enrichment <- enricher(
  gene = Genes$ENTREZID,         # 输入基因需为 Entrez ID
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  TERM2GENE = TERM2GENE,         # 本地通路-基因映射表
  TERM2NAME = TERM2NAME          # 本地通路名称表
)

# KEGG_enrichment <- enrichKEGG(gene = Genes$ENTREZID,
#                               organism = "hsa",
#                               keyType = "kegg", #KEGG数据库
#                               pAdjustMethod = "BH",
#                               pvalueCutoff = 0.05,
#                               qvalueCutoff =0.05)

p3 <- dotplot(KEGG_enrichment,
              title = "CXCL5+: Fibroblast",
              showCategory = 7,
              font.size = 22) + 
  #aes(color = pvalue) +
  scale_color_gradientn(colors = c("#3C5488", "#E64B35"), name = "adj.p")
p3

ggsave(filename = paste0(output,"/KEGG_CXCL5+_UC1.png"),
       plot = p3,
       width = 10,
       height = 8,
       dpi = 300 )
KEGG_data <- DOSE::setReadable(KEGG_enrichment,
                               OrgDb = org.Hs.eg.db,
                               keyType = "ENTREZID") #将ENTREZID转换为genesymbol
View(as.data.frame(KEGG_data))
write.table(data.frame(ID=rownames(KEGG_data@result),KEGG_data@result),
            file = paste0(output, "/KEGG_CXCL5+_UC",  ".txt"),
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

####
##富集分析 区分上下调基因
#####~UC GO####
#筛选出上调、下调的
# 筛选上调和下调基因
upregulated_genes <- markers[markers$avg_log2FC > 1 & markers$p_val_adj < 0.05, ]
downregulated_genes <- markers[markers$avg_log2FC < -1 & markers$p_val_adj < 0.05, ]

# 保存上调和下调基因的结果
write.csv(upregulated_genes, file = paste0(output, "/Up_UC_Genes.csv"))
write.csv(downregulated_genes, file = paste0(output, "/Down_UC_Genes.csv"))

# 将上调基因的符号转换为Entrez ID
up_gene_list <- rownames(upregulated_genes)
up_genes <- bitr(up_gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# 进行GO富集分析
up_GO_enrichment <- enrichGO(gene = up_genes$ENTREZID,
                             OrgDb = org.Hs.eg.db,
                             keyType = "ENTREZID",
                             ont = "all",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,   # p值阈值
                             qvalueCutoff = 1,      # FDR阈值
                             minGSSize = 10,
                             maxGSSize = 500,
                             readable = TRUE)

# 查看GO分析结果并保存
View(as.data.frame(up_GO_enrichment))
write.table(data.frame(ID = row.names(up_GO_enrichment@result), up_GO_enrichment@result),
            file = paste0(output, "/GO_UC_Up.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# 将下调基因的符号转换为Entrez ID
down_gene_list <- rownames(downregulated_genes)
down_genes <- bitr(down_gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# 进行GO富集分析
down_GO_enrichment <- enrichGO(gene = down_genes$ENTREZID,
                               OrgDb = org.Hs.eg.db,
                               keyType = "ENTREZID",
                               ont = "all",
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 1,
                               minGSSize = 10,
                               maxGSSize = 500,
                               readable = TRUE)
# 查看GO分析结果并保存
View(as.data.frame(down_GO_enrichment))
write.table(data.frame(ID = row.names(down_GO_enrichment@result), down_GO_enrichment@result),
            file = paste0(output, "/GO_UC_Down.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# 合并上调和下调的GO分析结果
up_GO_result <- as.data.frame(up_GO_enrichment@result)
down_GO_result <- as.data.frame(down_GO_enrichment@result)

# 为上调和下调添加分组
up_GO_result$group <- "1"
down_GO_result$group <- "-1"

# 排序筛选前5条 BP, CC, MF 各自的前5条结果
up_GO_BP <- head(up_GO_result[up_GO_result$ONTOLOGY == "BP",], 5)
down_GO_BP <- head(down_GO_result[down_GO_result$ONTOLOGY == "BP",], 5)

up_GO_CC <- head(up_GO_result[up_GO_result$ONTOLOGY == "CC",], 5)
down_GO_CC <- head(down_GO_result[down_GO_result$ONTOLOGY == "CC",], 5)

up_GO_MF <- head(up_GO_result[up_GO_result$ONTOLOGY == "MF",], 5)
down_GO_MF <- head(down_GO_result[down_GO_result$ONTOLOGY == "MF",], 5)

# 合并所有的上调和下调的GO数据
dat <- rbind(up_GO_BP, down_GO_BP, up_GO_CC, down_GO_CC, up_GO_MF, down_GO_MF)

# 计算-log10(p.adjust)用于绘图
dat$p.adjust <- -log10(dat$p.adjust)

# 确保p.adjust是数值型
dat$group <- as.numeric(dat$group)

class(dat$p.adjust)
class(dat$group)
summary(dat$p.adjust)
# 将下调通路的p值乘以-1，使其显示为负值
dat$p.adjust <- dat$p.adjust * dat$group

# 排序
dat <- dat[order(dat$p.adjust, decreasing = F),]

# 添加threshold列，用于分组显示Up或Down
dat$threshold <- ifelse(dat$p.adjust > 0, 'Up', 'Down')
#把重复删除
dat <- dat[!duplicated(dat$Description), ]

dat$Description <- factor(dat$Description, levels = dat$Description)
#dat$Description <- factor(dat$Description, levels = unique(dat$Description[order(dat$Description)]))

# 设置GO分类颜色
color_map <- c('BP' ='#3b6fd9', 'CC' = '#f6cfd8', 'MF' = '#8d689d')
#color_map <- c('BP' ='#53a8e1', 'CC' = '#fca3a3', 'MF' = '#c8c7e1')
#ann_colors = list(Group = c(Control = "#53a8e1",UC =  "#ca324c"))
#color_map <- c('BP' ='#8ecfb0', 'CC' = '#ffa74f', 'MF' = '#c8c7e1')
dat$ONTOLOGY <- factor(dat$ONTOLOGY, levels = c("BP", "CC", "MF"))

# 绘制GO富集通路的条形图
p4 <- ggplot(data = dat, aes(x = Description, y = p.adjust, fill = ONTOLOGY)) +
  #geom_col() +  # 绘制条形图
  geom_col(width=0.6)+
  coord_flip() +  # 翻转坐标轴
  #geom_text(aes(label = GeneRatio), size = 6) +  # 显示每个通路的GeneRatio
  labs(x = NULL, y = bquote("-"~Log[10]~"(adj.P)"), title = "GO Enrichment") +
  scale_y_continuous(expand = expansion(add = c(0.1, 0.1)),  # 设置y轴扩展
                     limits = c(-12, 15), breaks = seq(-20, 20, 10),
                     labels = c(20,10,0,10,20))# 设置y轴标签
p4
p <- p4 + 
  # 上调通路标签
  geom_text(data = dat[dat$threshold == "Up",], 
            aes(x = Description, y =  -0.1, label = Description),
            hjust = 1, color = 'black', size = 5.5) + 
  # 下调通路标签
  geom_text(data = dat[dat$threshold == "Down",], 
            aes(x = Description, y = 0.1, label = Description),
            hjust = 0, color = 'black', size = 5.5) + 
  scale_x_discrete(labels = NULL) +  # 不显示x轴标签
  scale_fill_manual(values = color_map) +  # 添加BP, CC, MF颜色映射
  theme_bw() +  # 使用白色背景
  theme(legend.position = "bottom",  
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        axis.line.y = element_blank(),  # 删除纵坐标轴
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(),  # 删除网格线
        panel.grid.minor = element_blank(), 
        panel.border = element_blank()) +# 删除面板边框
  theme(plot.title = element_text(hjust = 0.5))  # 设置标题居中

p
ggsave(filename = paste0(output,"/Fig3.i_GO_UC_updown.pdf"),
       plot = p,
       width = 17,
       height = 8,
       dpi = 300 )

#####~UC KEGG####    ###########################################################
# 进行KEGG富集分析

up_KEGG_enrichment <- enricher(
  gene = up_genes$ENTREZID,         # 输入基因需为 Entrez ID
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  TERM2GENE = TERM2GENE,         # 本地通路-基因映射表
  TERM2NAME = TERM2NAME          # 本地通路名称表
)

# up_KEGG_enrichment <- enrichKEGG(gene = up_genes$ENTREZID,
#                                  organism = "hsa", 
#                                  keyType = "kegg", #KEGG数据库
#                                  pAdjustMethod = "BH",
#                                  pvalueCutoff = 0.05,
#                                  qvalueCutoff =1)
up_KEGG_enrichment <- DOSE::setReadable(up_KEGG_enrichment,
                                        OrgDb = org.Hs.eg.db,
                                        keyType = "ENTREZID") #将ENTREZID转换为genesymbol

# 查看KEGG分析结果并保存
View(as.data.frame(up_KEGG_enrichment))
write.table(data.frame(ID = row.names(up_KEGG_enrichment@result), up_KEGG_enrichment@result),
            file = paste0(output, "/KEGG_UC_Up.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# 进行KEGG富集分析
down_KEGG_enrichment <- enricher(
  gene = down_genes$ENTREZID,         # 输入基因需为 Entrez ID
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  TERM2GENE = TERM2GENE,         # 本地通路-基因映射表
  TERM2NAME = TERM2NAME          # 本地通路名称表
)

# down_KEGG_enrichment <- enrichKEGG(gene = down_genes$ENTREZID,
#                                    organism = "hsa", 
#                                    keyType = "kegg", #KEGG数据库
#                                    pAdjustMethod = "BH",
#                                    pvalueCutoff = 0.05,
#                                    qvalueCutoff =1)
down_KEGG_enrichment <- DOSE::setReadable(down_KEGG_enrichment,
                                          OrgDb = org.Hs.eg.db,
                                          keyType = "ENTREZID") #将ENTREZID转换为genesymbol

# 查看KEGG分析结果并保存
View(as.data.frame(down_KEGG_enrichment))
write.table(data.frame(ID = row.names(down_KEGG_enrichment@result), down_KEGG_enrichment@result),
            file = paste0(output, "/KEGG_UC_Down.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# 合并上调和下调的KEGG分析结果
up_KEGG_result <- as.data.frame(up_KEGG_enrichment@result)
down_KEGG_result <- as.data.frame(down_KEGG_enrichment@result)

# 为上调和下调添加分组
up_KEGG_result$group <- "1"
down_KEGG_result$group <- "-1"

# 排序筛选前n
up_KEGG_result <- up_KEGG_result[order(up_KEGG_result$p.adjust),]
down_KEGG_result <- down_KEGG_result[order(down_KEGG_result$p.adjust),]

up_KEGG <- head(up_KEGG_result, 10)
down_KEGG <- head(down_KEGG_result, 10)

# 合并上调和下调的KEGG数据
dat <- rbind(up_KEGG, down_KEGG)

# 计算-log10(p.adjust)用于绘图
dat$p.adjust <- -log10(dat$p.adjust)

# 确保p.adjust是数值型
dat$group <- as.numeric(dat$group)

class(dat$p.adjust)
class(dat$group)
summary(dat$p.adjust)
# 将下调通路的p值乘以-1，使其显示为负值
dat$p.adjust <- dat$p.adjust * dat$group

# 排序
dat <- dat[order(dat$p.adjust, decreasing = F),]

# 添加threshold列，用于分组显示Up或Down
dat$threshold <- ifelse(dat$p.adjust > 0, 'Up', 'Down')
dat$Description <- factor(dat$Description, levels = dat$Description)
# 绘制KEGG富集通路的条形图
p5 <- ggplot(data = dat, aes(x = Description, y = p.adjust, fill = threshold)) +
  geom_col(width = 0.6) +  # 绘制条形图
  coord_flip() +  # 翻转坐标轴
  #geom_text(aes(label = GeneRatio), size = 6) +  # 显示每个通路的GeneRatio
  scale_fill_manual(values = c('Up' = '#BC3C29FF', 'Down' = '#0072B5FF')) +  # 上调和下调的颜色 #3C5488", "#E64B35
  #"#BC3C29FF", "Downregulated" = "#0072B5FF"
  #scale_fill_manual(values = c('Up' = '#cea5c7', 'Down' = '#b7deea')) +  # 上调和下调的颜色 #3C5488", "#E64B35  
  labs(x = NULL, y = bquote("-"~Log[10]~"(adj.P)"), title = "KEGG Pathway Enrichment") +
  scale_y_continuous(expand = expansion(add = c(0.1, 0.1)),  # 设置y轴扩展
                     limits = c(-20, 30), breaks = seq(-20, 30, 10),
                     labels = c(20,10,0,10,20,30))  # 设置y轴标签
p5
p_kegg <- p5 +
  # 上调通路标签
  geom_text(data = dat[dat$threshold == "Up",], 
            aes(x = Description, y =  -0.1, label = Description),
            hjust = 1, color = 'black', size = 5.5) + 
  # 下调通路标签
  geom_text(data = dat[dat$threshold == "Down",], 
            aes(x = Description, y = 0.1, label = Description),
            hjust = 0, color = 'black', size = 5.5) + 
  scale_x_discrete(labels = NULL) +  # 不显示x轴标签
  #geom_text(x = 12, y = 10, label = "Up", size = 10, color = '#36638a') +  # 上调标记
  #geom_text(x = 2, y = -10, label = "Down", size = 10, color = '#7bcd7b') +  # 下调标记
  theme_bw() +  # 使用白色背景
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        axis.line.y = element_blank(),  # 删除纵坐标轴
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.position = "none",  # 不显示图例
        panel.grid.major = element_blank(),  # 删除网格线
        panel.grid.minor = element_blank(),
        panel.border = element_blank())+# 删除面板边框
  theme(plot.title = element_text(hjust = 0.5))  # 设置标题居中
p_kegg
ggsave(filename = paste0(output,"/Fig3.j_KEGG_UC_updown.pdf"),
       plot = p_kegg,
       width = 10,
       height = 5,
       dpi = 300 )

#7.0  periodontitis analysis Volcano plot --------------------------------------
####
# P_colors <- c(
#   "Healthy" = "#8ecfb0",   
#   "Periodontitis" = "#f99655"
# )

RNA_data_P <- GetAssayData(ss_P, assay = "data")
# 按CXCL5基因表达分组
ss_P$CXCL5_status <- ifelse(RNA_data_P["CXCL5", ] > 0, "Positive", "Negative")
table(ss_P$CXCL5_status)

# 获取CXCL5基因的表达量数据
CXCL5_expr <- RNA_data_P["CXCL5", ]
ss_P$CXCL5_expr <- CXCL5_expr
table(ss_P$CXCL5_expr)

# 绘制CXCL5表达量分布图
p <- ggplot(ss_P@meta.data, aes(x = "", y = CXCL5_expr)) +
  geom_violin(fill = "lightblue", color = "black", alpha = 0.6) +
  geom_jitter(aes(color = CXCL5_expr), size = 1, width = 0.2, alpha = 0.6) +
  theme_minimal() +
  labs(title = "CXCL5 Expression Distribution in Fibroblasts", 
       x = "", y = "Expression Level") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())
p
# 保存图像
ggsave(filename = paste0(output, "/CXCL5_P.png"),
       plot = p,
       width = 8,
       height = 6,
       dpi = 300)

# 进行差异表达分析（对Positive与Negative组进行比较）
ss_P <- JoinLayers(ss_P, assay = "RNA")
markersP <- FindMarkers(ss_P, 
                       ident.1 = "Positive", 
                       ident.2 = "Negative", 
                       group.by = "CXCL5_status", 
                       test.use = "MAST",
                       min.pct = 0.2)

# 保存差异表达分析结果
write.csv(markersP, file = paste0(output, "/CXCL5+_P.csv"))

# 生成和保存CXCL5状态的分组统计表
status_table <- table(ss_P$CXCL5_status)
write.csv(as.data.frame(status_table), file = paste0(output, "/Table_CXCL5_P.csv"))

status_table
prop.table(status_table)
markersP$Significance <- with(markersP, 
                             ifelse(p_val_adj < 0.05 & avg_log2FC > 1, "Upregulated",
                                    ifelse(p_val_adj < 0.05 & avg_log2FC < -1, "Downregulated", "Not Significant")))

# 筛选显著的差异基因
significant_markers <- markersP[abs(markersP$avg_log2FC) > 1 & markersP$p_val_adj < 0.05, ]
write.csv(significant_markers,
          file = paste0(output,"/DE_CXCL5+_P.csv"))

# 找出在差异基因中的相关因子
# 从 GO terms 获取相关基因
# library(org.Hs.eg.db)
# secreted_genes <- unique(AnnotationDbi::select(org.Hs.eg.db,
#                                                keys = "GO:0005615", # 分泌到胞外空间的基因
#                                                columns = "SYMBOL",
#                                                keytype = "GOALL")$SYMBOL)
# inflammatory_genes <- unique(AnnotationDbi::select(org.Hs.eg.db,
#                                                    keys = "GO:0006954", # 炎性反应
#                                                    columns = "SYMBOL",
#                                                    keytype = "GOALL")$SYMBOL)
# cytokine_genes <- unique(AnnotationDbi::select(org.Hs.eg.db,
#                                                keys = "GO:0005125", # 细胞因子活性
#                                                columns = "SYMBOL",
#                                                keytype = "GOALL")$SYMBOL)
# chemokine_genes <- unique(AnnotationDbi::select(org.Hs.eg.db,
#                                                 keys = "GO:0008009", # 趋化因子活性
#                                                 columns = "SYMBOL",
#                                                 keytype = "GOALL")$SYMBOL)
# 
# #合并
# factor_genes <- unique(c(secreted_genes, inflammatory_genes,
#                          cytokine_genes, chemokine_genes))
# # 
highlight_genes <- rownames(significant_markers) %in% factor_genes$Gene
highlighted_markers <- significant_markers[highlight_genes, ]
sum(highlight_genes)

inflammatory_genes_genes <- rownames(significant_markers) %in% inflammatory_genes$Gene
inflammatory_genes_markers <- significant_markers[inflammatory_genes_genes, ]
sum(inflammatory_genes_genes)

cytokine_genes <- rownames(significant_markers) %in% cytokine_genes$Gene
cytokine_genes_markers <- significant_markers[cytokine_genes, ]
sum(cytokine_genes)

chemokine_genes <- rownames(significant_markers) %in% chemokine_genes$Gene
chemokine_genes_markers <- significant_markers[chemokine_genes, ]
sum(chemokine_genes)

# 绘制火山图并标出这些因子
p1 <- ggplot(markersP, 
             aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = Significance), alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Upregulated" = "#BC3C29FF", "Downregulated" = "#0072B5FF", "Not Significant" = "grey")) +
  theme(
    legend.position = "none",
    #axis.text.y = element_text(margin = margin(t = 0, r = 2, b = 0, l = 0)),  # 调整 y 轴标签间距
    axis.text.x = element_text(size = 12, color = "black"),  # X轴刻度数字
    axis.text.y = element_text(size = 12, color = "black"),  # Y轴刻度数字
    axis.title.x = element_text(size = 15, color = "black"),  # X轴刻度数字
    axis.title.y = element_text(size = 15, color = "black"),  # Y轴刻度数字
    # 分面标签（若有分面）
    plot.title = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加黑色边框
    panel.background = element_blank()  # 背景设为透明
  )+
  labs(title = "UC", x = bquote(~Log[2]~"(FC)"), y = bquote("-"~Log[10]~"(adj.P)")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  coord_cartesian(ylim = c(0, 100)) +
  coord_cartesian(xlim = c(-7, 10)) +
  #coord_cartesian(xlim = c(-5, 10)) +
  # 添加分泌型因子、炎性因子、细胞因子、趋化因子的标签
  geom_text_repel(data = inflammatory_genes_markers, 
                  aes(label = rownames(inflammatory_genes_markers)), 
                  size = 3, 
                  box.padding = 0.5, 
                  point.padding = 0.5, 
                  max.overlaps = 20)
p1
ggsave(filename = paste0(output,"/","Fig3.k_Volcano_CXCL5+_P.pdf"), 
       plot = p1, 
       width = 6.5, 
       height = 6.5)

#8.0  periodontitis analysis prop ------------------------------------
ss_P_new <- ss_P
ss_P_new$NewCellType <- "Other" 
RNA_data <- GetAssayData(ss_P_new, assay = "data")
# 标记CXCL5+ Fibro (表达值 > 0)
ss_P_new$NewCellType[Cells(ss_P_new)[RNA_data["CXCL5", ] > 0]] <- "CXCL5+ Fibro"
table(ss_P_new$NewCellType)
ss_P_new$PatientID <- paste(ss_P_new$orig.ident, ss_P_new$batch, ss_P_new$label, sep = "_")

ss_P_new$histo3 <- ss_P_new$group
table(ss_P_new$histo3)

# 1. 提取CXCL5+ Fibro 细胞的数据并计算其比例
meta <- ss_P_new@meta.data
meta$CXCL5_Fibro <- ifelse(meta$NewCellType == "CXCL5+ Fibro", 1, 0)  # 创建新的列，标记CXCL5+ Fibro细胞
# 按照PatientID和histo3分组，计算每个组中CXCL5+ Fibro细胞的比例
prop_CXCL5_Fibro <- meta %>%
  group_by(PatientID, histo3) %>%
  summarise(CXCL5_Fibro_Prop = mean(CXCL5_Fibro))  # 计算每个组内CXCL5+ Fibro的比例

##CXCL5+ 表达量的组间差异
ss_P_new$anno <- ss_P_new$cell_type
data_gene <- FetchData(ss_P_new, vars = c("anno", "CXCL5", "histo3"))

# 过滤掉CXCL5表达为0的细胞
# data_gene <- data_gene %>%
#   filter(CXCL5 > 0)
p <- ggplot(data_gene, 
            aes(x = histo3, y = CXCL5, fill = histo3)) +
  geom_violin(trim = FALSE, size = 0.8, alpha = 0.7) +
  geom_boxplot(
    position = position_dodge(width = 0.75),  # 调整箱体位置
    width = 0.2,  # 控制箱体宽度
    fill = "white",  # 白色箱体
    color = "black",  # 箱体边框
    alpha = 0.6,  # 设置透明度
    outlier.shape = NA  # 移除箱外点
  ) +
  geom_point(position = position_jitter(width = 0.2), size = 2, alpha = 0.6, 
             aes(color = histo3)) +
  #facet_wrap(~anno, scales = "free", ncol = 5) +
  theme_bw() +
  labs(
    x = "Group",
    y = "CXCL5  Expression"
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.y = unit(0, "lines"),
    strip.text = element_text(size = 12, margin = margin(0, 0, 2, 0)),
    strip.background = element_blank(),
    legend.position = "bottom",  # 图例
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 12),
    axis.ticks.y = element_line(color = "black", size = 0.25),
    axis.ticks.length.y = unit(0.05, "cm"),
    axis.title.y = element_text(color = "black", size = 12),
    axis.title.x = element_blank(),
    plot.margin = margin(1, 1, 1, 1)
  ) +
  scale_fill_manual(values = c("Healthy" = "#20854EFF",   
                               "Periodontitis" = "#E18727FF"),name=NULL)+
  scale_color_manual(values = c("Healthy" = "#20854EFF",   
                                "Periodontitis" = "#E18727FF"),name=NULL)
# 使用 `geom_signif` 添加星号和横线
p <- p + geom_signif(
  comparisons = list(
    c("Healthy", "Periodontitis")
  ),
  map_signif_level = TRUE,
  #y_position = c( 3.5,4.5),
  y_position = c(2.9),
  color = "red",
  #size = 2,
  textsize = 3
)+
  # 控制y轴范围（根据数据范围调整 y_min 和 y_max）
  coord_cartesian(ylim = c(0, 3.1))
p
# 保存图片
ggsave(filename = paste0(output, "/Fig3.g_Expr_CXCL5_P_group.pdf"),
       plot = p, width = 4, height = 3, dpi = 300)

#9.0  periodontitis analysis functional ana ------------------------------------
library(clusterProfiler)
gene_list <- rownames(significant_markers)

# 将基因符号转换为Entrez ID
Genes <- bitr(gene_list,  
              fromType = "SYMBOL", #输入数据的类型
              toType = c("ENTREZID"), #要转换的数据类型
              OrgDb = org.Hs.eg.db) #物种

# 进行GO富集分析
GO_enrichment <- enrichGO(gene = Genes$ENTREZID,
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID",
                          ont = "all",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,   #p
                          qvalueCutoff = 0.05, #fdr P
                          minGSSize = 10,
                          maxGSSize = 500,
                          readable = TRUE)

View(as.data.frame(GO_enrichment))
write.table(data.frame(ID = row.names(GO_enrichment@result), GO_enrichment@result),
            file = paste0(output, "/GO_CXCL5+_P",  ".txt"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

p2 <- dotplot(GO_enrichment,
              title = "CXCL5+: Fibroblast",
              showCategory = 7,
              font.size = 22,
              split="ONTOLOGY") + 
  #aes(color = pvalue)+
  facet_grid(ONTOLOGY~., scale='free') + 
  scale_color_gradientn(colors = c("#3C5488", "#E64B35"),name = "adj.p")
p2
ggsave(filename = paste0(output,"/GO_CXCL5+_P.png"),
       plot = p2,
       width = 14,
       height = 16,
       dpi = 300 )

# KEGG富集分析  ########################################################################################
# up_KEGG_enrichment <- enricher(
#   gene = up_genes$ENTREZID,         # 输入基因需为 Entrez ID
#   pAdjustMethod = "BH",
#   pvalueCutoff = 0.05,
#   qvalueCutoff = 0.05,
#   TERM2GENE = TERM2GENE,         # 本地通路-基因映射表
#   TERM2NAME = TERM2NAME          # 本地通路名称表
# )
KEGG_enrichment <- enricher(gene = Genes$ENTREZID,
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.05,
                              qvalueCutoff =0.05,
                              TERM2GENE = TERM2GENE,         # 本地通路-基因映射表
                              TERM2NAME = TERM2NAME          # 本地通路名称表
                              )

p3 <- dotplot(KEGG_enrichment,
              title = "CXCL5+: Fibroblast",
              showCategory = 7,
              font.size = 22) + 
  #aes(color = pvalue) +
  scale_color_gradientn(colors = c("#3C5488", "#E64B35"), name = "adj.p")

p3
ggsave(filename = paste0(output,"/KEGG_CXCL5+_P.png"),
       plot = p3,
       width = 10,
       height = 8,
       dpi = 300 )
KEGG_data <- DOSE::setReadable(KEGG_enrichment,
                               OrgDb = org.Hs.eg.db,
                               keyType = "ENTREZID") #将ENTREZID转换为genesymbol
View(as.data.frame(KEGG_data))
write.table(data.frame(ID=rownames(KEGG_data@result),KEGG_data@result),
            file = paste0(output, "/KEGG_CXCL5+_P",  ".txt"),
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)
##富集分析 区分上下调基因
#####~P GO####  #####################
#筛选出上调、下调的
# 筛选上调和下调基因
markers <- markersP
upregulated_genes <- markers[markers$avg_log2FC > 1 & markers$p_val_adj < 0.05, ]
downregulated_genes <- markers[markers$avg_log2FC < -1 & markers$p_val_adj < 0.05, ]

# 保存上调和下调基因的结果
write.csv(upregulated_genes, file = paste0(output, "/Up_P_Genes.csv"))
write.csv(downregulated_genes, file = paste0(output, "/Down_P_Genes.csv"))

# 将上调基因的符号转换为Entrez ID
up_gene_list <- rownames(upregulated_genes)
up_genes <- bitr(up_gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# 进行GO富集分析
up_GO_enrichment <- enrichGO(gene = up_genes$ENTREZID,
                             OrgDb = org.Hs.eg.db,
                             keyType = "ENTREZID",
                             ont = "all",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,   # p值阈值
                             qvalueCutoff = 1,      # FDR阈值
                             minGSSize = 10,
                             maxGSSize = 500,
                             readable = TRUE)

# 查看GO分析结果并保存
View(as.data.frame(up_GO_enrichment))
write.table(data.frame(ID = row.names(up_GO_enrichment@result), up_GO_enrichment@result),
            file = paste0(output, "/GO_P_Up.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# 将下调基因的符号转换为Entrez ID
down_gene_list <- rownames(downregulated_genes)
down_genes <- bitr(down_gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# 进行GO富集分析
down_GO_enrichment <- enrichGO(gene = down_genes$ENTREZID,
                               OrgDb = org.Hs.eg.db,
                               keyType = "ENTREZID",
                               ont = "all",
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 1,
                               minGSSize = 10,
                               maxGSSize = 500,
                               readable = TRUE)

# 查看GO分析结果并保存
View(as.data.frame(down_GO_enrichment))
write.table(data.frame(ID = row.names(down_GO_enrichment@result), down_GO_enrichment@result),
            file = paste0(output, "/GO_P_Down.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# 合并上调和下调的GO分析结果
up_GO_result <- as.data.frame(up_GO_enrichment@result)
down_GO_result <- as.data.frame(down_GO_enrichment@result)

# 为上调和下调添加分组
up_GO_result$group <- "1"
down_GO_result$group <- "-1"

# 排序筛选前5条 BP, CC, MF 各自的前5条结果
up_GO_BP <- head(up_GO_result[up_GO_result$ONTOLOGY == "BP",], 5)
down_GO_BP <- head(down_GO_result[down_GO_result$ONTOLOGY == "BP",], 5)

up_GO_CC <- head(up_GO_result[up_GO_result$ONTOLOGY == "CC",], 5)
down_GO_CC <- head(down_GO_result[down_GO_result$ONTOLOGY == "CC",], 5)

up_GO_MF <- head(up_GO_result[up_GO_result$ONTOLOGY == "MF",], 5)
down_GO_MF <- head(down_GO_result[down_GO_result$ONTOLOGY == "MF",], 5)

# 合并所有的上调和下调的GO数据
dat <- rbind(up_GO_BP, down_GO_BP, up_GO_CC, down_GO_CC, up_GO_MF, down_GO_MF)
# 计算-log10(p.adjust)用于绘图
dat$p.adjust <- -log10(dat$p.adjust)
# 确保p.adjust是数值型
dat$group <- as.numeric(dat$group)

class(dat$p.adjust)
class(dat$group)
summary(dat$p.adjust)
# 将下调通路的p值乘以-1，使其显示为负值
dat$p.adjust <- dat$p.adjust * dat$group

# 排序
dat <- dat[order(dat$p.adjust, decreasing = F),]
# 添加threshold列，用于分组显示Up或Down
dat$threshold <- ifelse(dat$p.adjust > 0, 'Up', 'Down')
#dat$Description <- factor(dat$Description, levels = dat$Description)
dat$Description <- factor(dat$Description, levels = unique(dat$Description[order(dat$Description)]))

# 设置GO分类颜色
#color_map <- c('BP' ='#6d6fa0', 'CC' = '#8d689d', 'MF' = '#c8c7e1')
#color_map <- c('BP' ='#8ecfb0', 'CC' = '#ffa74f', 'MF' = '#c8c7e1')
color_map <- c('BP' ='#3b6fd9', 'CC' = '#f6cfd8', 'MF' = '#8d689d')
dat$ONTOLOGY <- factor(dat$ONTOLOGY, levels = c("BP", "CC", "MF"))

# 绘制GO富集通路的条形图
p4 <- ggplot(data = dat, aes(x = reorder(Description, p.adjust), y = p.adjust, fill = ONTOLOGY)) +
  geom_col(width = 0.6) +  # 绘制条形图
  coord_flip() +  # 翻转坐标轴
  #geom_text(aes(label = GeneRatio), size = 6) +  # 显示每个通路的GeneRatio
  labs(x = NULL, y = bquote("-"~Log[10]~"(adj.P)"), title = "GO Enrichment") +
  scale_y_continuous(expand = expansion(add = c(0.1, 0.1)),  # 设置y轴扩展
                     limits = c(-20, 25), breaks = seq(-20, 20, 10),
                     labels = c(20,10,0,10,20))# 设置y轴标签
p4
p <- p4 + 
  # 上调通路标签
  geom_text(data = dat[dat$threshold == "Up",], 
            aes(x = Description, y =  -0.1, label = Description),
            hjust = 1, color = 'black', size = 5.5) + 
  # 下调通路标签
  geom_text(data = dat[dat$threshold == "Down",], 
            aes(x = Description, y = 0.1, label = Description),
            hjust = 0, color = 'black', size = 5.5) + 
  scale_x_discrete(labels = NULL) +  # 不显示x轴标签
  scale_fill_manual(values = color_map) +  # 添加BP, CC, MF颜色映射
  theme_bw() +  # 使用白色背景
  theme(legend.position = "bottom",  
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        axis.line.y = element_blank(),  # 删除纵坐标轴
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(),  # 删除网格线
        panel.grid.minor = element_blank(), 
        panel.border = element_blank()) +# 删除面板边框
  theme(plot.title = element_text(hjust = 0.5))  # 设置标题居中

p
ggsave(filename = paste0(output,"/Fig3.l_GO_P_updown.pdf"),
       plot = p,
       width = 17,
       height = 8,
       dpi = 300 )

#####~P KEGG####
# 进行KEGG富集分析
up_KEGG_enrichment <- enricher(
  gene = up_genes$ENTREZID,         # 输入基因需为 Entrez ID
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  TERM2GENE = TERM2GENE,         # 本地通路-基因映射表
  TERM2NAME = TERM2NAME          # 本地通路名称表
)
# up_KEGG_enrichment <- enrichKEGG(gene = up_genes$ENTREZID,
#                                  organism = "hsa", 
#                                  keyType = "kegg", #KEGG数据库
#                                  pAdjustMethod = "BH",
#                                  pvalueCutoff = 0.05,
#                                  qvalueCutoff =1)
up_KEGG_enrichment <- DOSE::setReadable(up_KEGG_enrichment,
                                        OrgDb = org.Hs.eg.db,
                                        keyType = "ENTREZID") #将ENTREZID转换为genesymbol

# 查看KEGG分析结果并保存
View(as.data.frame(up_KEGG_enrichment))
write.table(data.frame(ID = row.names(up_KEGG_enrichment@result), up_KEGG_enrichment@result),
            file = paste0(output, "/KEGG_P_Up.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# 进行KEGG富集分析
down_KEGG_enrichment <- enricher(
  gene = down_genes$ENTREZID,         # 输入基因需为 Entrez ID
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 1,
  TERM2GENE = TERM2GENE,         # 本地通路-基因映射表
  TERM2NAME = TERM2NAME          # 本地通路名称表
)
# down_KEGG_enrichment <- enrichKEGG(gene = down_genes$ENTREZID,
#                                    organism = "hsa", 
#                                    keyType = "kegg", #KEGG数据库
#                                    pAdjustMethod = "BH",
#                                    pvalueCutoff = 0.05,
#                                    qvalueCutoff =1)
down_KEGG_enrichment <- DOSE::setReadable(down_KEGG_enrichment,
                                          OrgDb = org.Hs.eg.db,
                                          keyType = "ENTREZID") #将ENTREZID转换为genesymbol
# 查看KEGG分析结果并保存
View(as.data.frame(down_KEGG_enrichment))
write.table(data.frame(ID = row.names(down_KEGG_enrichment@result), down_KEGG_enrichment@result),
            file = paste0(output, "/KEGG_P_Down.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# 合并上调和下调的KEGG分析结果
up_KEGG_result <- as.data.frame(up_KEGG_enrichment@result)
down_KEGG_result <- as.data.frame(down_KEGG_enrichment@result)

# 为上调和下调添加分组
up_KEGG_result$group <- "1"
down_KEGG_result$group <- "-1"

# 排序筛选前n
up_KEGG_result <- up_KEGG_result[order(up_KEGG_result$p.adjust),]
down_KEGG_result <- down_KEGG_result[order(down_KEGG_result$p.adjust),]

up_KEGG <- head(up_KEGG_result, 10)
down_KEGG <- head(down_KEGG_result, 10)

# 合并上调和下调的KEGG数据
dat <- rbind(up_KEGG, down_KEGG)
# 计算-log10(p.adjust)用于绘图
dat$p.adjust <- -log10(dat$p.adjust)
# 确保p.adjust是数值型
dat$group <- as.numeric(dat$group)

class(dat$p.adjust)
class(dat$group)
summary(dat$p.adjust)
# 将下调通路的p值乘以-1，使其显示为负值
dat$p.adjust <- dat$p.adjust * dat$group
# 排序
dat <- dat[order(dat$p.adjust, decreasing = F),]

# 添加threshold列，用于分组显示Up或Down
dat$threshold <- ifelse(dat$p.adjust > 0, 'Up', 'Down')
dat$Description <- factor(dat$Description, levels = dat$Description)

# 绘制KEGG富集通路的条形图
p5 <- ggplot(data = dat, aes(x = Description, y = p.adjust, fill = threshold)) +
  geom_col(width = 0.6) +  # 绘制条形图
  coord_flip() +  # 翻转坐标轴
  #geom_text(aes(label = GeneRatio), size = 6) +  # 显示每个通路的GeneRatio
  scale_fill_manual(values = c('Up' = '#BC3C29FF', 'Down' = '#0072B5FF')) +  # 上调和下调的颜色
  labs(x = NULL, y = bquote("-"~Log[10]~"(adj.P)"), title = "KEGG Pathway Enrichment") +
  scale_y_continuous(expand = expansion(add = c(0.1, 0.1)),  # 设置y轴扩展
                     limits = c(-10, 10), breaks = seq(-10, 10, 5),
                     labels = c(10,5,0,5,10))  # 设置y轴标签
#"#BC3C29FF", "Downregulated" = "#0072B5FF"
p_kegg <- p5 +
  # 上调通路标签
  geom_text(data = dat[dat$threshold == "Up",], 
            aes(x = Description, y =  -0.1, label = Description),
            hjust = 1, color = 'black', size = 5.5) + 
  # 下调通路标签
  geom_text(data = dat[dat$threshold == "Down",], 
            aes(x = Description, y = 0.1, label = Description),
            hjust = 0, color = 'black', size = 5.5) + 
  scale_x_discrete(labels = NULL) +  # 不显示x轴标签
  #geom_text(x = 12, y = 10, label = "Up", size = 10, color = '#36638a') +  # 上调标记
  #geom_text(x = 2, y = -10, label = "Down", size = 10, color = '#7bcd7b') +  # 下调标记
  theme_bw() +  # 使用白色背景
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        axis.line.y = element_blank(),  # 删除纵坐标轴
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.position = "none",  # 不显示图例
        panel.grid.major = element_blank(),  # 删除网格线
        panel.grid.minor = element_blank(),
        panel.border = element_blank())+# 删除面板边框
  theme(plot.title = element_text(hjust = 0.5))  # 设置标题居中
p_kegg
ggsave(filename = paste0(output,"/Fig3.m_KEGG_P_updown.pdf"),
       plot = p_kegg,
       width = 13,
       height = 5,
       dpi = 300 )

