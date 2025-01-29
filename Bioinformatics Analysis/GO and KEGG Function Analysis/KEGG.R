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
library(pathview)         # KEGG通路可视化
library(patchwork)        # 多图拼接（组合多个ggplot图形）
library(circlize)         # 绘制环形图（用于复杂可视化）
library(ComplexHeatmap)   # 绘制复杂热图
library(scales)           # 图形比例尺控制
library(viridis)          # viridis 包，用于颜色

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
# OrgDb：使用的注释数据库（人类基因org.Hs.eg.db,小鼠用org.Mm.eg.db）
# drop：是否删除无法转换的ID（默认TRUE，此处隐式使用）
#-------------------------------------------------------------------------------------------------------
gene.df <- bitr(diff$SYMBOL, 
                fromType = "SYMBOL", 
                toType = "ENTREZID", 
                OrgDb = org.Hs.eg.db)    # 物种为人类

gene <- gene.df$ENTREZID  # 提取转换后的Entrez ID向量

# 3. KEGG富集分析函数
#-------------------------------------------------------------------------------------------------------
# 定义函数：perform_kegg_enrichment
# 参数：
#      genes：基因列表（需为ENTREZID格式）
#-------------------------------------------------------------------------------------------------------
perform_kegg_enrichment <- function(genes) {
  enrichKEGG(
    gene = genes,            # 输入基因列表
    organism = "hsa",        # 物种标识符,hsa人类,mmu小鼠,rno大鼠
    keyType = "ncbi-geneid", # 输入基因ID类型
    pvalueCutoff = 0.05,     # p值显著性阈值
    qvalueCutoff = 0.05      # q值（FDR校正后p值）阈值
  )
}

kegg <- perform_kegg_enrichment(gene)

#4. 将结果保存到当前路径
ego_result_kegg <- as.data.frame(kegg)
write.csv(ego_result_kegg,file = "ego_result_kegg.csv",row.names = T)

# 5. 绘制KEGG富集分析结果图
#===================================================================================================================================
#' @param enrich_result 富集分析结果对象，enrichKEGG函数返回的对象
#' @param ontology 本体类型，字符串
#' @param top_n 显示前多少个条目，数值，不能超过csv文件里面的富集通路个数，否则报错
#' @param output_filename 输出文件名，字符串，例如 "kegg_bubbleplot.png"
#' @return 无返回值，直接保存图片到指定路径
#===================================================================================================================================

# 绘制KEGG通路图和通路背景图,输出文件名格式为hsa+五位数字（KEGG官方ID）,hsa表示人类,其他物种需调整（如mmu为小鼠,rno为大鼠）
top_n <- 5  # 绘制前5个通路
for (i in 1:top_n) {
  pathway_id <- kegg@result$ID[i]
  pathway_name <- kegg@result$Description[i]
  clean_name <- gsub("[^[:alnum:]_]", "_", pathway_name)  # 更严格的清理
  pathview(gene.data = gene, pathway.id = pathway_id, species = "hsa")   # hsa表示人类
  old_file <- paste0("hsa", pathway_id, ".pathview.png")
  new_file <- paste0("Pathway_", i, "_", clean_name, ".png")
  if (file.exists(old_file)) file.rename(old_file, new_file)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 绘制气泡图
generate_kegg_bubble <- function(enrich_result, ontology, top_n = 10, output_filename) {
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
    geom_segment(aes(xend = 0, yend = Description), linetype = "dashed", color = "grey50") +
    scale_color_gradient(low = "#9869c9", high = "#b24175") +
    scale_size_continuous(range = c(4, 10)) +   #c(4,10)是描述气泡范围大小的,气泡最小4像素，最大10像素
    labs(
      title = paste("KEGG", ontology, sep = " "), # 标题中加入本体类型
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
generate_kegg_bubble(kegg, "Pathway", output_filename = "kegg_bubbleplot.png")

#--------------------------------------------------------------------------------------------------------

# 绘制柱状图
generate_kegg_bar <- function(enrich_result, ontology, top_n = 10, output_filename) {
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
      title = paste("KEGG", ontology, " "),
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
generate_kegg_bar(kegg, "Pathway", output_filename = "kegg_barplot.png")

#--------------------------------------------------------------------------------------------------------

# 绘制圈(弦)图
generate_kegg_circos <- function(enrich_result, output_filename, top_n = 6) {
  # 准备数据
  data <- as.data.frame(enrich_result)[1:top_n, ]
  genes_entrez <- strsplit(data$geneID, "/") # 先保留 Entrez ID
  
  # 创建 Entrez ID 到 SYMBOL 的映射
  entrez_to_symbol <- setNames(gene.df$SYMBOL, gene.df$ENTREZID)
  
  # 将 Entrez ID 转换为基因符号 (SYMBOL)
  genes_symbol_lists <- lapply(genes_entrez, function(entrez_ids) {
    sapply(entrez_ids, function(id) {
      if(id %in% names(entrez_to_symbol)) {
        return(entrez_to_symbol[id]) # 查找并返回基因符号
      } else {
        return(id) # 如果找不到对应的基因符号，则返回原始 Entrez ID (可以根据需求修改处理方式)
      }
    })
  })
  
  links <- data.frame(
    Term = rep(data$Description, sapply(genes_symbol_lists, length)),
    Gene = unlist(genes_symbol_lists) # 使用转换后的基因符号
  )
  term_sizes <- setNames(sapply(genes_entrez, length), data$Description) # term_sizes 仍然基于原始 Entrez ID 的数量
  
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
generate_kegg_circos(kegg, "kegg_circos.png")

#--------------------------------------------------------------------------------------------------------

# 绘制网络图
generate_kegg_cnet <- function(enrich_result, ontology, output_filename, top_n = 10) {
  # 提取前top_n个条目避免过于拥挤
  enrich_result <- setReadable(enrich_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  p <- cnetplot(enrich_result, showCategory = top_n, 
                colorEdge = TRUE, 
                node_label = "all",  # 显示所有节点标签
                cex_label_category = 1.2,
                cex_category = 1.5,
                color_category = "#d6456c",
                shadowtext = "category") +
    ggtitle(paste("KEGG", ontology, "Gene-Concept Network")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  png(output_filename, width = 2000, height = 1600, res = 150)
  print(p)
  dev.off()
}

# 绘制并保存网络图
generate_kegg_cnet(kegg, "Pathway", output_filename = "kegg_cnetplot.png")

#--------------------------------------------------------------------------------------------------------

# 绘制词云图（Word cloud）
generate_kegg_wordcloud <- function(enrich_result, output_filename, top_n = 6) {
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
generate_kegg_wordcloud(kegg, "kegg_wordcloud.png")








































