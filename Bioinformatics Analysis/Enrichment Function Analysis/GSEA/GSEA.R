# 设置工作目录（根据实际情况修改或删除）
setwd("F:/R")

# 加载包
library(msigdbr)          # 用于获取MSigDB中的基因集数据（如Hallmark基因集）
library(clusterProfiler)  # 用于执行GSEA分析
library(fgsea)
library(tidyverse)        # 用于数据处理和绘图
library(ggplot2)
library(gridExtra)
library(cowplot)

# 读取数据,从CSV文件中读取经过批次校正的基因表达数据，假设行名为基因名，列为样本
data <- read.csv("gene.csv", row.names = 1, header = TRUE)

# 定义分组,使用 setdiff 函数将不属于group1的样本归为group2
group1 <- c("GSM1354764","GSM1354765","GSM1354766","GSM1354767","GSM1354768")
group2 <- setdiff(colnames(data), group1)

# 计算基因排序指标
#----------------------------------------------------------------------------------------------------------
# 对每组样本计算基因表达的平均值（rowMeans）,以 group1 和 group2 的差异作为基因排序指标（gene_rank）,
# 随后依据该指标降序排列基因（sorted_genes）,为后续 GSEA 分析提供输入
#----------------------------------------------------------------------------------------------------------
avg_group1 <- rowMeans(data[, group1], na.rm = TRUE)
avg_group2 <- rowMeans(data[, group2], na.rm = TRUE)
gene_rank <- avg_group1 - avg_group2
sorted_genes <- sort(gene_rank, decreasing = TRUE)

# GSEA分析
#----------------------------------------------------------------------------------------------------------
# 使用 msigdbr 包下载人类（Homo sapiens）的Hallmark基因集
# 将基因集按名称拆分为列表（hallmark_list）
# 使用 fgsea 包执行快速GSEA分析
#----------------------------------------------------------------------------------------------------------
hallmark_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_gene_sets$gene_symbol, hallmark_gene_sets$gs_name)

set.seed(123)
fgsea_res <- fgsea(
  pathways = hallmark_list, # 基因集列表
  stats = sorted_genes,     # 基因排序指标
  minSize = 15,             # minSize 和 maxSize：限制基因集大小范围。
  maxSize = 500,
  eps = 0,                  # 控制数值稳定性
  nproc = 4                 # 并行计算的核心数
)

# 根据GSEA结果的p值排序，提取前8个最显著的通路名称。
top_pathways <- fgsea_res[order(pval), ]$pathway[1:8]  #可修改选择绘制排名前几的通路,如果想绘制前10的就修改为：pathway[1:10]

# 计算富集分数曲线函数
calc_es <- function(geneset, ranked_list) {
  indicator <- as.integer(names(ranked_list) %in% geneset)
  p_hit <- cumsum(abs(ranked_list) * indicator) / sum(abs(ranked_list) * indicator)
  p_miss <- cumsum(!indicator) / sum(!indicator)
  es <- p_hit - p_miss
  es[is.na(es)] <- 0
  return(es)
}

# 生成单个GSEA图的函数
generate_gsea_plot <- function(pathway_name) {
  pathway_genes <- hallmark_list[[pathway_name]]
  es <- calc_es(pathway_genes, sorted_genes)
  peak_idx <- which.max(abs(es))
  peak_es <- es[peak_idx]
  leading_edge <- 1:peak_idx
  
  # 构建绘图数据
  plot_df <- data.frame(
    rank = seq_along(sorted_genes),
    metric = sorted_genes,
    es = es,
    in_geneset = names(sorted_genes) %in% pathway_genes,
    in_leading_edge = seq_along(sorted_genes) %in% leading_edge
  )
  
  # 第一部分：ES曲线部分
  p1 <- ggplot(plot_df, aes(x = rank)) +
    geom_line(aes(y = es), color = "green", size = 1.0) +  # 加粗绿色曲线
    geom_vline(xintercept = peak_idx, linetype = "dashed", color = "red") +  # 红色竖线
    geom_hline(yintercept = 0, color = "grey30") +
    labs(y = "Enrichment Score", title = paste0(pathway_name, " ")) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          plot.title = element_text(size=10, hjust = 0.5, face = "bold"))
  
  # 第二部分：基因集位置
  p2 <- ggplot(plot_df) +
    geom_segment(aes(x = rank, xend = rank, y = 0, yend = 1, 
                     color = in_leading_edge), 
                 data = ~filter(.x, in_geneset)) +
    scale_color_manual(values = c("grey60", "darkgreen")) +
    scale_y_continuous(breaks = NULL) +
    labs(y = "Gene Set\nMembers", x = "") +
    theme_minimal() +
    theme(axis.line = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none")
  
  # 第三部分：排序指标分布
  p3 <- ggplot(plot_df, aes(x = rank)) +
    geom_area(aes(y = metric, fill = metric > 0)) +
    scale_fill_manual(values = c("blue", "red")) +
    geom_hline(yintercept = 0, color = "grey30") +
    labs(x = "Rank in Ordered Dataset", y = "Ranked List Metric") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # 组合图形
  combined_plot <- plot_grid(
    p1, p2, p3,
    ncol = 1,
    align = "v",
    rel_heights = c(2, 0.5, 1.5)
  )
  
  # 添加统计值
  pathway_stats <- fgsea_res[fgsea_res$pathway == pathway_name, ]
  final_plot <- ggdraw(combined_plot) + 
    draw_label(
      label = sprintf("ES = %.2f\np = %.2e\npadj = %.2e", 
                      pathway_stats$ES, 
                      pathway_stats$pval,
                      pathway_stats$padj),
      x = 0.85, y = 0.9,
      hjust = 1, vjust = 1,
      size = 8,
      color = "black"
    )
  return(final_plot)
}

# 生成前8通路的图形
plot_list <- lapply(top_pathways, generate_gsea_plot)

# 组合成每一行3个图形的布局
combined_plots <- plot_grid(plotlist = plot_list, ncol = 3)  #可修改ncol的值,控制每一行几个图形

# 保存图形
ggsave("GSEA_top8_pathways.png", combined_plots, 
       width = 16, height = 18, dpi = 300, bg = "white")

# 保存结果
fgsea_res %>%
  arrange(pval) %>%
  mutate(leadingEdge = map_chr(leadingEdge, ~ paste(.x, collapse = ","))) %>%
  write.csv("gsea_results.csv", row.names = FALSE)
































