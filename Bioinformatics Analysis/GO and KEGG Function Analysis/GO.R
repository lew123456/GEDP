# 设置工作目录 (可选，如果数据文件在当前工作目录下则不需要)
setwd("F:/R")          # 设置工作目录路径，用于指定R的工作文件夹

# 加载必要的R包
library(AnnotationDbi)    # 提供基因ID转换功能（底层依赖）
library(org.Hs.eg.db)     # 人类基因组注释数据库（提供基因ID映射关系）
library(clusterProfiler)  # 进行富集分析的核心包
library(dplyr)            # 数据处理（提供管道操作符和数据处理函数）
library(ggplot2)          # 绘图基础包
library(wordcloud)        # 生成词云图
library(stringr)          # 字符串处理（用于文本换行等操作）
library(enrichplot)       # 富集分析可视化扩展（支持多种富集结果绘图）
library(DOSE)             # 疾病本体分析（此代码中未直接使用）
library(pathview)         # KEGG通路可视化（此代码中未直接使用）
library(patchwork)        # 多图拼接（组合多个ggplot图形）
library(circlize)         # 绘制环形图（用于复杂可视化）
library(ComplexHeatmap)   # 绘制复杂热图
library(scales)           # 图形比例尺控制

# 1. 数据读取
#-------------------------------------------------------------------------------------------------------
# 读取想要进行富集分析的基因数据文件，必须包含SYMBOL列（基因符号）
# read.csv()：读取CSV格式文件
# file：文件路径
# header：是否包含表头（默认为TRUE）
# check.names：是否检查列名有效性（建议保持默认FALSE）
#-------------------------------------------------------------------------------------------------------
diff <- read.csv("gene.csv") 

# 2. 基因ID转换
#-------------------------------------------------------------------------------------------------------
# bitr()：生物学ID转换函数
# diff$SYMBOL：输入的基因符号向量
# fromType：输入ID类型（此处为基因符号SYMBOL）
# toType：目标ID类型（此处为Entrez基因ID）
# OrgDb：使用的注释数据库（人类基因org.Hs.eg.db）
# drop：是否删除无法转换的ID（默认TRUE，此处隐式使用）
#-------------------------------------------------------------------------------------------------------
gene.df <- bitr(diff$SYMBOL, 
                fromType = "SYMBOL", 
                toType = "ENTREZID", 
                OrgDb = org.Hs.eg.db)

gene <- gene.df$ENTREZID  # 提取转换后的Entrez ID向量

# 3. GO富集分析函数
#-------------------------------------------------------------------------------------------------------
# 定义函数：perform_go_enrichment
# 参数：
#   gene_list：基因列表（需为ENTREZID格式）
#   ontology：GO本体类型（BP/CC/MF）
#-------------------------------------------------------------------------------------------------------
perform_go_enrichment <- function(gene_list, ontology) {
  enrichGO(
    gene          = gene_list,       # 输入基因列表（必须为ENTREZID）
    OrgDb         = org.Hs.eg.db,    # 使用的注释数据库,通过这里指定物种,人类基因org.Hs.eg.db,小鼠用org.Mm.eg.db
    keyType       = "ENTREZID",      # 输入ID类型（与gene参数对应）
    ont           = ontology,        # GO本体类型：BP（生物过程）/CC（细胞组分）/MF（分子功能）
    pAdjustMethod = "BH",            # p值校正方法：可选"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
    pvalueCutoff  = 0.05,            # p值阈值（筛选显著项）
    qvalueCutoff  = 0.05,            # q值阈值（FDR校正后）
    readable      = TRUE             # 是否将ENTREZID转换为基因符号显示结果
  )
}

# 进行三个本体的GO富集分析
ego_BP <- perform_go_enrichment(gene, "BP")  # （生物过程）
ego_CC <- perform_go_enrichment(gene, "CC")  # （细胞组分）
ego_MF <- perform_go_enrichment(gene, "MF")  # （分子功能）


#4. 将结果保存到当前路径
ego_result_BP <- as.data.frame(ego_BP)
ego_result_CC <- as.data.frame(ego_CC)
ego_result_MF <- as.data.frame(ego_MF)
ego <- rbind(ego_result_BP,ego_result_CC,ego_result_MF)
write.csv(ego_result_BP,file = "ego_result_BP.csv",row.names = T)
write.csv(ego_result_CC,file = "ego_result_CC.csv",row.names = T)
write.csv(ego_result_MF,file = "ego_result_MF.csv",row.names = T)
write.csv(ego,file = "ego.csv",row.names = T)

# 5. 绘制GO富集分析结果图
#===================================================================================================================================
#' @param enrich_result 富集分析结果对象，enrichGO函数返回的对象
#' @param ontology 本体类型，字符串，可以是 "BP" (Biological Process), "CC" (Cellular Component), 或 "MF" (Molecular Function)
#' @param top_n 显示前多少个条目，数值，不能超过csv文件里面的富集通路个数，否则报错
#' @param output_filename 输出文件名，字符串，例如 "GO_BP_bubbleplot.png"
#' @return 无返回值，直接保存图片到指定路径
#===================================================================================================================================

# 自定义主题设置，用于绘制气泡图和柱状图
custom_theme <- theme_bw() +
  theme(
    plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),  # 图标题：14号加粗字体,水平居中(hjust=0.5)
    axis.title       = element_text(size = 13),                              # 坐标轴标题:(如"Gene Ratio")的字体大小设为13
    axis.text        = element_text(size = 12),                              # 坐标轴文本:(如基因名称)的字体大小设为12
    legend.title     = element_text(size = 12),                              # 图例标题:("Gene Count"等)使用12号字体
    legend.text      = element_text(size = 10),                              # 图例文本
    panel.grid.minor = element_blank(),                                      # 移除次要网格线
    panel.border     = element_rect(color = "black", linewidth = 0.7, fill = NA), # 设置边框:线宽0.7，不填充颜色
    plot.margin      = margin(t = 20, r = 20, b = 20, l = 100, unit = "pt")  # 设置图形边距
  )

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

generate_go_bubble <- function(enrich_result, ontology, top_n = 10, output_filename) {
  # 将富集结果转换为数据框
  data <- as.data.frame(enrich_result)
  
  # 计算GeneRatio，并对Description进行文本换行
  data <- data %>%
    mutate(
      GeneRatio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
      Description = str_wrap(Description, width = 35)  #对过长的描述文字自动换行，每行最多35字符，防止y轴标签重叠
    )
  
  # 截取前 top_n 个条目
  data <- data[1:top_n, ]
  
  # 绘制气泡图
  p <- ggplot(data, aes(x = GeneRatio, y = Description)) +
    geom_point(aes(size = Count, color = -log10(p.adjust))) +
    scale_color_gradient(low = "#9499c0", high = "#d6456c") +
    scale_size_continuous(range = c(4, 10)) +   #c(4,10)是描述气泡范围大小的,气泡最小4像素，最大10像素
    labs(
      title = paste("GO", ontology, sep = " "), # 标题中加入本体类型
      x = "Gene Ratio",
      size = "Gene Count",
      color = "-log10(p.adjust)"
    ) +
    ylab(NULL) +
    custom_theme +
    theme(axis.text.y = element_text(size = 10, hjust = 1, lineheight = 0.8))
  
  # 保存图片
  png(output_filename, width = 1300, height = 800, res = 150)  #保存图片的长、宽、分辨率
  print(p)
  dev.off()
}

# 绘制并保存气泡图
generate_go_bubble(ego_BP, "Biological Process (BP)", output_filename = "GO_BP_bubbleplot.png")
generate_go_bubble(ego_CC, "Cellular Component (CC)", output_filename = "GO_CC_bubbleplot.png")
generate_go_bubble(ego_MF, "Molecular Function (MF)", output_filename = "GO_MF_bubbleplot.png")

#--------------------------------------------------------------------------------------------------------

# 绘制柱状图
generate_go_bar <- function(enrich_result, ontology, top_n = 10, output_filename) {
  data <- as.data.frame(enrich_result) %>%
    arrange(p.adjust) %>%
    head(top_n) %>%
    mutate(
      Description = str_wrap(Description, width = 35),
      Description = factor(Description, levels = rev(Description))
    )
  
  p <- ggplot(data, aes(x = Count, y = Description, fill = -log10(p.adjust))) +
    geom_col(width = 0.8) +
    scale_fill_gradient(low = "#9499c0", high = "#d6456c") +
    labs(
      title = paste("GO", ontology, " "),
      x = "Gene Count",
      y = NULL,
      fill = "-log10(p.adjust)"
    ) +
    custom_theme +
    theme(axis.text.y = element_text(size = 10, lineheight = 0.8))
  
  png(output_filename, width = 1300, height = 800, res = 150)
  print(p)
  dev.off()
}

# 绘制并保存柱状图
generate_go_bar(ego_BP, "Biological Process (BP)", output_filename = "GO_BP_barplot.png")
generate_go_bar(ego_CC, "Cellular Component (CC)", output_filename = "GO_CC_barplot.png")
generate_go_bar(ego_MF, "Molecular Function (MF)", output_filename = "GO_MF_barplot.png")

#--------------------------------------------------------------------------------------------------------

# 绘制圈(弦)图
generate_go_circos <- function(enrich_result, output_filename, top_n = 8) {
  # 准备数据
  data <- as.data.frame(enrich_result)[1:top_n, ]
  genes <- strsplit(data$geneID, "/")
  links <- data.frame(
    Term = rep(data$Description, sapply(genes, length)),
    Gene = unlist(genes)
  )
  term_sizes <- setNames(sapply(genes, length), data$Description)
  
  # 创建颜色映射
  go_terms <- unique(links$Term)
  #grid.col <- setNames(hue_pal()(length(go_terms)), go_terms)  #颜色映射可选
  grid.col <- setNames(viridis_pal(option = "A")(length(go_terms)), go_terms)  # option可选A-D
  
  # 设置输出参数
  png(output_filename, width = 12000, height = 5500, res = 150) # **增加总宽度
  
  # **调整绘图边距（右边界大幅增加以容纳图例）
  par(mar = c(1, 1, 1, 100), xpd = TRUE) # xpd=TRUE允许在边距绘图
  
  # 初始化弦图
  circos.par(start.degree = 90, gap.degree = 4)
  chordDiagram(links, 
               grid.col = grid.col,
               transparency = 0.5,
               annotationTrack = "grid",
               preAllocateTracks = list(track.height = 0.1))
  
  # 基因标签轨道
  circos.trackPlotRegion(
    track.index = 1,
    panel.fun = function(x, y) {
      sector.name = get.cell.meta.data("sector.index")
      if (!sector.name %in% go_terms) {
        circos.text(CELL_META$xcenter, 
                    CELL_META$ylim[1], 
                    sector.name,
                    facing = "clockwise",
                    niceFacing = TRUE,
                    adj = c(0, 0.5),
                    cex = 6.0)   # 字体大小
      }
    }, 
    bg.border = NA
  )
  
  # 刻度轨道
  # 修正后的刻度轨道部分
  max_size <- max(term_sizes)
  circos.track(factors = links$Term, track.index = 2,
               ylim = c(0, max_size), track.height = 0.1,
               panel.fun = function(x, y) {
                 sector.index <- get.cell.meta.data("sector.index")
                 if (sector.index %in% names(term_sizes)) {  # 这里修改为%in%
                   current_size <- term_sizes[sector.index]
                   circos.axis(
                     major.at = c(0, current_size),
                     labels = c(0, current_size),
                     direction = "outside",
                     labels.cex = 5.5,    # 刻度大小
                     major.tick.length = 0.2
                   )
                 }
               }, bg.border = NA)
  
  # **计算图例位置（位于右侧边距中央）
  legend_x <- grconvertX(0.775, "npc")  # 水平位置（0-1范围）
  legend_y <- grconvertY(0.5, "npc")   # 垂直居中
  
  # 绘制图例
  legend(
    x = legend_x,
    y = legend_y,
    legend = go_terms,
    fill = grid.col,
    border = NA,
    bty = "n",
    cex = 5.5,   #图例字体大小
    xjust = 0,
    yjust = 0.5,
    x.intersp = 0.5,    # 图例项水平间距
    y.intersp = 2.0)    # 图例项垂直间距
  
  circos.clear()
  dev.off()
}
# 绘制并保存圈(弦)图
generate_go_circos(ego_BP, "GO_BP_circos.png")
generate_go_circos(ego_CC, "GO_CC_circos.png")
generate_go_circos(ego_MF, "GO_MF_circos.png")

#--------------------------------------------------------------------------------------------------------

# 绘制网络图
generate_go_cnet <- function(enrich_result, ontology, output_filename, top_n = 10) {
  # 提取前top_n个条目避免过于拥挤
  enrich_result <- setReadable(enrich_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  p <- cnetplot(enrich_result, showCategory = top_n, 
                colorEdge = TRUE, 
                node_label = "all",  # 显示所有节点标签
                cex_label_category = 1.2,
                cex_category = 1.5,
                color_category = "#d6456c",
                shadowtext = "category") +
    ggtitle(paste("GO", ontology, "Gene-Concept Network")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  png(output_filename, width = 2000, height = 1600, res = 150)
  print(p)
  dev.off()
}

# 绘制并保存网络图
generate_go_cnet(ego_BP, "Biological Process (BP)", output_filename = "GO_BP_cnetplot.png")
generate_go_cnet(ego_CC, "Cellular Component (CC)", output_filename = "GO_CC_cnetplot.png")
generate_go_cnet(ego_MF, "Molecular Function (MF)", output_filename = "GO_MF_cnetplot.png")

#--------------------------------------------------------------------------------------------------------

# 绘制词云图（Word cloud）
generate_go_wordcloud <- function(enrich_result, output_filename, top_n = 8) {
  data <- as.data.frame(enrich_result)[1:top_n, ]
  words <- strsplit(data$Description, " ")
  freq <- data$Count
  
  png(output_filename, width = 1000, height = 1000, res = 150)
  wordcloud(words = data$Description, 
            freq = freq,
            scale = c(3, 0.5),
            colors = brewer.pal(8, "Dark2"),
            random.order = FALSE)
  dev.off()
}

# 绘制并保存词云图
generate_go_wordcloud(ego_BP, "GO_BP_wordcloud.png")
generate_go_wordcloud(ego_CC, "GO_CC_wordcloud.png")
generate_go_wordcloud(ego_MF, "GO_MF_wordcloud.png")









