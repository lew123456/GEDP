# 设置工作目录（根据实际情况修改或删除）
setwd("F:/R")

# 加载所需包
library(limma)
library(sva)
library(Rtsne)
library(umap)
library(ggplot2)
library(ggforce)  # 用于绘制椭圆

# 读取数据 ---------------------------------------------------------------
# 读取第一个数据集（GSE56081组）
data_ss <- read.csv("GSE56081.csv", 
                    row.names = 1,  # 使用第一列作为行名
                    header = TRUE,  # 保留列标题
                    check.names = FALSE)  # 防止修改列名
# 读取第二个数据集（GSE23130组）
data_vv <- read.csv("GSE23130.csv",
                    row.names = 1,
                    header = TRUE,
                    check.names = FALSE)

# 数据预处理 -------------------------------------------------------------
# 获取共同基因列表
common_genes <- intersect(rownames(data_ss), rownames(data_vv))
# 提取共有基因数据
data_ss_common <- data_ss[common_genes, ]
data_vv_common <- data_vv[common_genes, ]
# 合并数据集（列合并）
combined_data <- cbind(data_ss_common, data_vv_common)

# 创建批次信息
batch <- factor(rep(c("GSE56081", "GSE23130"), 
                    times = c(ncol(data_ss_common), 
                              ncol(data_vv_common))))

# 定义颜色
colors <- c("#FF9999", "#66B2FF")  # 自定义颜色：红色和蓝色

# 绘制原始数据的彩色箱线图
png("boxplot_combined_data_colored.png", width = 1000, height = 600)  # 保存为 PNG 格式
boxplot(combined_data, col = rep(colors, times = c(ncol(data_ss_common), ncol(data_vv_common))), 
        main = "Boxplot of Combined Data", ylab = "Expression")
dev.off()

# 批次效应校正 -----------------------------------------------------------
# 使用ComBat进行校正
combat_corrected <- ComBat(
  dat = as.matrix(combined_data),
  batch = batch,
  mod = NULL,
  par.prior = TRUE,
  prior.plots = FALSE
)

# 绘制批次校正后的彩色箱线图
png("boxplot_data_colored.png", width = 1000, height = 600)  # 保存为 PNG 格式
boxplot(combat_corrected, col = rep(colors, times = c(ncol(data_ss_common), ncol(data_vv_common))), 
        main = "Boxplot of Batch-Corrected Data", ylab = "Expression")
dev.off()


# 将 combat_corrected 转换为数据框并保存为 CSV 文件
combat_corrected_df <- as.data.frame(combat_corrected)
write.csv(combat_corrected_df, file = "batch_corrected_data.csv", row.names = TRUE)

#========================================================================================================

# 绘图可视化
# PCA分析函数 ------------------------------------------------------------
perform_pca <- function(data_matrix, scale_data = TRUE) {
  # 转置矩阵使样本为行，基因为列
  pca_result <- prcomp(t(data_matrix), scale. = scale_data)
  return(pca_result)
}

# 执行PCA分析
pca_raw <- perform_pca(combined_data)      # 原始数据PCA

pca_corrected <- perform_pca(combat_corrected)  # 校正后数据PCA

# 绘图函数 ---------------------------------------------------------------
create_pca_plot <- function(pca_result, 
                            batch_info, 
                            title_text,
                            color_palette = c("#1F77B4", "#FF7F0E"),
                            point_size = 3.5,
                            ellipse_alpha = 0.15) {
  
  # 计算方差解释比例
  variance <- summary(pca_result)$importance[2, 1:2] * 100
  
  # 创建绘图数据框
  plot_data <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    Group = batch_info
  )
  
  # 生成坐标轴标签
  x_label <- sprintf("Dim1 (%.1f%%)", variance[1])
  y_label <- sprintf("Dim2 (%.1f%%)", variance[2])
  
  # 创建基础图形
  p <- ggplot(plot_data, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
    geom_point(size = point_size, alpha = 0.85) +
    geom_mark_ellipse(
      aes(fill = Group), 
      alpha = ellipse_alpha,
      show.legend = FALSE,
      linetype = "dashed",
      size = 0.6,
      expand = unit(0.65, "cm")
    ) +
    scale_color_manual(values = color_palette) +
    scale_fill_manual(values = color_palette) +
    scale_shape_manual(values = c(16, 17)) +  # 实心圆和三角
    labs(
      title = title_text,
      x = x_label,
      y = y_label
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 16,
        face = "bold",
        margin = margin(b = 15)
      ),
      axis.title = element_text(face = "bold", size = 13),
      axis.text = element_text(color = "gray30"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      aspect.ratio = 1  # 保持正方形画布
    ) +
    guides(
      color = guide_legend(
        nrow = 1,
        override.aes = list(size = 4, alpha = 1)
      ),
      shape = guide_legend(nrow = 1)
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dotted",
      color = "grey50",
      linewidth = 0.6
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dotted",
      color = "grey50",
      linewidth = 0.6
    )
  
  return(p)
}

# 生成图形 ---------------------------------------------------------------
# 设置统一配色方案
my_colors <- c("#3B9AB2", "#EBCC2A")  # 蓝金配色

# 原始数据PCA图
p_raw <- create_pca_plot(
  pca_result = pca_raw,
  batch_info = batch,
  title_text = "PCA Before Batch Correction",
  color_palette = my_colors
)

# 校正后数据PCA图
p_corrected <- create_pca_plot(
  pca_result = pca_corrected,
  batch_info = batch,
  title_text = "PCA After Batch Correction",
  color_palette = my_colors
)

# 显示图形
print(p_raw)
print(p_corrected)

# 保存图形 ---------------------------------------------------------------
# 设置保存参数
save_plot <- function(plot_obj, filename, width = 20, height = 18) {
  ggsave(
    filename = filename,
    plot = plot_obj,
    device = "png",
    units = "cm",
    width = width,
    height = height,
    dpi = 600,
    bg = "white"
  )
}

# 保存图片
save_plot(p_raw, "PCA_Before_Correction.png")
save_plot(p_corrected, "PCA_After_Correction.png")








