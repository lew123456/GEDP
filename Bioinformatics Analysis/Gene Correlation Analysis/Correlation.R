# 设置工作目录 (可选，如果数据文件在当前工作目录下则不需要)
setwd("F:/R")          # 设置工作目录路径，用于指定R的工作文件夹 

# 加载所需的R包
library(pheatmap)         # 用于绘制热图
library(ggcorrplot)       # 用于绘制相关性图
library(grid)             # 提供图形布局支持
library(corrplot)         # 用于绘制相关性图
library(RColorBrewer)     # 提供颜色方案
library(corrgram)         # 用于相关性分析
library(ggplot2)          # 用于数据可视化

# 读取CSV文件中的数据，并将第一列设置为行名（通常是基因名）
data <- read.csv("gene.csv", header = TRUE, row.names = 1)

#------------------------------------------------------------------------------------------

#绘制基因之间的相关性热图
# 转置数据，使得每一行代表一个样本，每一列代表一个基因
gene_data <- t(data)

# 计算基因表达数据的皮尔逊相关系数矩阵
cor_matrix <- cor(gene_data, method = "pearson")  # 使用皮尔逊相关方法

# 计算相关系数的显著性p值矩阵，置信水平为95%
pmtcars <- cor.mtest(gene_data, conf.level = 0.95)

# 自定义颜色方案：从紫色到白色再到橙色的渐变
col <- colorRampPalette(c("#008B8B", "white", "#FF9966"))(100)  # 生成100个颜色等级

# 打开PNG图形设备，准备保存图像
png("corrgram_plot.png", width = 2000, height = 3500, res = 300)

# 绘制相关性热图
corrplot(cor_matrix,
         col = col,               # 使用自定义的颜色方案
         tl.col = "black",        # 标签文字颜色为黑色
         method = 'ellipse',      # 使用椭圆表示相关性（正相关为圆形，负相关为扁平椭圆）
         order = 'hclust',        # 按层次聚类排序
         cl.pos = 'b',            # 颜色图例位置在底部
         cl.cex = 1.5,            # 颜色图例字体大小
         tl.cex = 1.5,            # 标签字体大小
         p.mat = pmtcars$p,       # 输入显著性p值矩阵
         sig.level = 0.01,        # 显著性阈值（p < 0.01的点会标记为显著）
         pch.cex = 5)             # 显著性标记的大小

# 关闭图形设备，保存图像
dev.off()

#------------------------------------------------------------------------------------------

#绘制一对基因之间的相关性图
genes <- c("MAPK1", "FOXO1", "ESR1", "SP1")
expr_data <- as.data.frame(t(data[genes, ]))

# 自定义主题
theme_custom <- theme_bw() +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 15),
    legend.position = "top",
    panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
    plot.background = element_rect(fill = "white", colour = NA)
  )

# 定义一个函数来计算相关性和绘制散点图
plot_gene_correlation <- function(gene1, gene2, data, output_file) {
  # 计算相关性和p值
  cor_result <- cor.test(data[[gene1]], data[[gene2]], method = "pearson")
  r_value <- round(cor_result$estimate, 2)
  p_value <- signif(cor_result$p.value, 2)
  
  # 绘制散点图
  plot <- ggplot(data, aes(x = .data[[gene1]], y = .data[[gene2]])) +
    geom_point(color = "#008B8B", size = 8, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "#FF9966", linewidth = 1, fill = "#ACE1AF", alpha = 0.3) +
    labs(
      title = paste(gene1, "&", gene2),
      x = paste(gene1, "Expression (log2)"),
      y = paste(gene2, "Expression (log2)")
    ) +
    annotate(
      "text",
      x = min(data[[gene1]]),
      y = max(data[[gene2]]),
      label = paste0("r = ", r_value, "\nP = ", p_value),
      size = 5,
      hjust = 0
    ) +
    theme_custom
  
  # 保存图像
  ggsave(output_file, plot, width = 6, height = 6, dpi = 300)
}

# 调用函数绘制每对基因的相关性图
plot_gene_correlation("MAPK1", "FOXO1", expr_data, "MAPK1_FOXO1_cor.png")
plot_gene_correlation("MAPK1", "ESR1", expr_data, "MAPK1_ESR1_cor.png")
plot_gene_correlation("SP1", "ESR1", expr_data, "SP1_ESR1_cor.png")
plot_gene_correlation("SP1", "FOXO1", expr_data, "SP1_FOXO1_cor.png")