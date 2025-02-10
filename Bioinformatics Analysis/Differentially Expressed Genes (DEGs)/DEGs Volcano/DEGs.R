# 设置工作目录 (可选，如果数据文件在当前工作目录下则不需要)
setwd("F:/R") 

# 加载必要的R包
library(limma)
library(ggplot2)

# 1. 数据读取
# 读取基因表达数据，假设文件名为 "Data_Processing.csv"，将第一行作为列名，将第一列作为行名。
#==========================================================================================================================================================
#header = TRUE (列名/表头):
#如果 CSV 文件第一行是列名,务必设置 header = TRUE,否则 R 会将列名当作数据。
#row.names = 1 (行名):
#row.names 可以是数字,表示第几列是行名；也可以是字符串,表示列名是什么。例如,如果行名在名为 "ID" 的列中,可以使用 row.names = "ID"。
#重要提示：当使用 row.names = 1 时,R 会将第一列作为行名,并且这一列将不再作为数据的一部分出现在数据框中。这意味着数据框的列数会比 CSV 文件中少一列。
#==========================================================================================================================================================
expression_data <- read.csv("Data_Processing.csv", header = TRUE, row.names = 1)

# 定义样本分组信息，假设前 5 个样本是 Group1，后 5 个样本是 Group2
#group <- factor(c(rep("Group1", 5), rep("Group2", 5)))
# 定义样本分组信息（按列名匹配的为一组，剩下的为一组）
group1_samples <- c("GSM1354764","GSM1354765","GSM1354766","GSM1354767","GSM1354768")
group <- factor(ifelse(colnames(expression_data) %in% group1_samples, "Group1", "Group2"))

# 检查分组信息长度是否与样本数一致
if (length(group) != ncol(expression_data)) {
  stop("错误: 分组信息的长度与样本数量不匹配！")
}

# 2. 差异表达分析
# 构建设计矩阵
design <- model.matrix(~ group)

# 拟合线性模型
fit <- lmFit(expression_data, design)

# 定义对比矩阵：比较 Group2 和 Group1
contrast.matrix <- makeContrasts(Group2vsGroup1 = groupGroup2, levels = design)

# 将对比应用于线性模型
fit2 <- contrasts.fit(fit, contrast.matrix)

# 计算经验贝叶斯统计量
fit2 <- eBayes(fit2)

# 获取差异表达分析结果，使用 BH 方法进行多重比较校正，并显示所有基因的结果
results <- topTable(fit2, coef = "Group2vsGroup1", adjust.method = "BH", number = Inf)

# 设置筛选标准：Fold Change 的绝对值大于 log2(2)=1，校正后的 p 值，即adj.P.Value小于 0.05
fold_change_cutoff <- 2
adj_p_cutoff <- 0.05

# 根据 Fold Change 和校正后的 p 值筛选显著差异表达基因 (DEGs)
degs <- results[
  (results$adj.P.Val < adj_p_cutoff) & (abs(results$logFC) > log2(fold_change_cutoff)),
]

# 将筛选出的 DEGs 保存到文件 "Gene_DEGs.csv"
write.csv(degs, file = "Gene_DEGs.csv", row.names = TRUE)

# 3. 火山图绘制
# 准备火山图数据
volcano_data <- data.frame(
  logFC = results$logFC,
  negLogPval = -log10(results$adj.P.Val),
  Gene = rownames(results),
  Regulation = ifelse(
    results$adj.P.Val >= adj_p_cutoff | abs(results$logFC) <= log2(fold_change_cutoff),
    "Not Significant",  # 非显著基因
    ifelse(
      results$logFC > log2(fold_change_cutoff),
      "Upregulated",  # 上调基因
      "Downregulated"  # 下调基因
    )
  )
)

# 绘制火山图
volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = negLogPval, color = Regulation)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c(
    "Not Significant" = "#b2b2b2",
    "Upregulated" = "#5cb07f",
    "Downregulated" = "#9793c6"
  )) +
  theme_classic() +
  labs(
    x = "log2(Fold Change)",
    y = "-log10(Adjusted p-value)",
    title = "Volcano Plot of Differential Expression",
    color = "Gene Regulation"
  ) +
  geom_hline(yintercept = -log10(adj_p_cutoff), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-log2(fold_change_cutoff), log2(fold_change_cutoff)), linetype = "dashed", color = "black") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 13),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.7),
    plot.margin = margin(1, 1, 1, 1, "cm") #增加图片边距，防止标签显示不完全
  ) +
  coord_cartesian(clip = "off") # 使用 coord_cartesian 函数并设置 clip = "off"，允许标签溢出坐标系

# 保存火山图
ggsave("Avolcano.png", plot = volcano_plot, width = 6, height = 7, dpi = 600)

# 4. 结果输出和保存
# 统计上调、下调和非显著基因的数量
num_upregulated_genes <- sum(volcano_data$Regulation == "Upregulated")
num_downregulated_genes <- sum(volcano_data$Regulation == "Downregulated")
num_not_significant_genes <- sum(volcano_data$Regulation == "Not Significant")

# 输出统计结果
cat("上调基因的数量为:", num_upregulated_genes, "\n")
cat("下调基因的数量为:", num_downregulated_genes, "\n")
cat("非显著基因的数量为:", num_not_significant_genes, "\n")

# 将上调和下调基因保存到文件
regulated_genes <- subset(volcano_data, Regulation %in% c("Upregulated", "Downregulated"))
write.csv(regulated_genes, file = "regulated_genes.csv", row.names = FALSE)
cat("上调基因和下调基因已保存到 regulated_genes.csv 文件中。\n")

# 将非显著基因保存到文件
not_significant_genes <- subset(volcano_data, Regulation == "Not Significant")
write.csv(not_significant_genes, file = "not_significant_genes.csv", row.names = FALSE) 
cat("非显著基因已保存到 not_significant_genes.csv 文件中。\n")